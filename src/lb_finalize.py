from __future__ import annotations

import json
import logging
import math
import time
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

from lb_common import bh_fdr, elapsed, lambda_gc_from_p, safe_neglog10, save_json
from lb_equivalence import (
    alias_aware_ld_clump,
    build_equivalence_tables,
    collapse_for_method,
    compute_or_load_equivalence,
)
from lb_meta import acat_pair, concordance_metrics, method_comparison, nearby_genes, parse_gff_genes, simes_pair
from lb_power import build_power_table

PRIMARY_TEST = "score"
SENSITIVITY_TEST = "lrt"
DIAGNOSTIC_TEST = "wald"


def build_test_table(
    inc: pd.DataFrame,
    sev: pd.DataFrame,
    test: str,
) -> pd.DataFrame:
    p_col = f"p_{test}"
    required = {"variant_index", "snp", "CHR", "pos", "MAF", "missing_rate", p_col}
    missing_inc = required.difference(inc.columns)
    missing_sev = required.difference(sev.columns)
    if missing_inc or missing_sev:
        raise ValueError(f"Missing columns for {test}: INC={sorted(missing_inc)}, SEV={sorted(missing_sev)}")

    base_cols = ["variant_index", "snp", "CHR", "pos", "MAF", "missing_rate"]
    inc_cols = base_cols + [c for c in ["beta", "se", "l_remle", p_col] if c in inc.columns]
    sev_cols = ["variant_index"] + [c for c in ["beta", "se", "l_remle", p_col] if c in sev.columns]
    left = inc[inc_cols].rename(
        columns={
            "beta": "beta_INC",
            "se": "se_INC",
            "l_remle": "l_remle_INC",
            p_col: "p_INC",
        }
    )
    right = sev[sev_cols].rename(
        columns={
            "beta": "beta_SEV",
            "se": "se_SEV",
            "l_remle": "l_remle_SEV",
            p_col: "p_SEV",
        }
    )
    out = left.merge(right, on="variant_index", how="inner", validate="one_to_one")
    if len(out) != len(inc):
        raise ValueError(f"Merged {len(out):,}/{len(inc):,} variants for {test}")
    out["p_ACAT"] = acat_pair(out["p_INC"].to_numpy(float), out["p_SEV"].to_numpy(float))
    out["p_SIMES"] = simes_pair(out["p_INC"].to_numpy(float), out["p_SEV"].to_numpy(float))
    for method in ["INC", "SEV", "ACAT", "SIMES"]:
        out[f"q_{method}"] = bh_fdr(out[f"p_{method}"].to_numpy(float))
    for trait in ["INC", "SEV"]:
        if f"beta_{trait}" in out.columns and f"se_{trait}" in out.columns:
            out[f"ci_low_{trait}"] = out[f"beta_{trait}"] - 1.96 * out[f"se_{trait}"]
            out[f"ci_high_{trait}"] = out[f"beta_{trait}"] + 1.96 * out[f"se_{trait}"]
    out["test"] = test
    return out


def build_single_trait_table(scan: pd.DataFrame, test: str, trait: str) -> pd.DataFrame:
    p_col = f"p_{test}"
    keep = [
        "variant_index", "snp", "CHR", "pos", "MAF", "missing_rate",
        *[c for c in ["beta", "se", "l_remle", p_col] if c in scan.columns],
    ]
    out = scan[keep].copy().rename(columns={p_col: "p"})
    out["q_BH"] = bh_fdr(out["p"].to_numpy(float))
    out["trait"] = trait
    out["test"] = test
    if "beta" in out and "se" in out:
        out["ci_low"] = out["beta"] - 1.96 * out["se"]
        out["ci_high"] = out["beta"] + 1.96 * out["se"]
    return out


def profile_test_table(table: pd.DataFrame, test: str, context: str = "loco") -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for method in ["INC", "SEV", "ACAT", "SIMES"]:
        p = table[f"p_{method}"].to_numpy(float)
        q = table[f"q_{method}"].to_numpy(float)
        rows.append(
            {
                "context": context,
                "test": test,
                "method": method,
                "n_variants": int(np.isfinite(p).sum()),
                "lambda_gc": lambda_gc_from_p(p),
                "minimum_p": float(np.nanmin(p)),
                "minimum_q": float(np.nanmin(q)),
                "BH_FDR_0p10_hits": int(np.nansum(q < 0.10)),
            }
        )
    return pd.DataFrame(rows)


def _jaccard_top(left: np.ndarray, right: np.ndarray, fraction: float) -> float:
    n = len(left)
    k = max(1, int(round(n * fraction)))
    a = set(np.argpartition(left, k - 1)[:k].tolist())
    b = set(np.argpartition(right, k - 1)[:k].tolist())
    return len(a & b) / len(a | b)


def compare_tests(tables: dict[str, pd.DataFrame]) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for left_name, right_name in [("score", "lrt"), ("score", "wald"), ("lrt", "wald")]:
        if left_name not in tables or right_name not in tables:
            continue
        left = tables[left_name]
        right = tables[right_name]
        for method in ["INC", "SEV", "ACAT", "SIMES"]:
            p_left = left[f"p_{method}"].to_numpy(float)
            p_right = right[f"p_{method}"].to_numpy(float)
            x = safe_neglog10(p_left)
            y = safe_neglog10(p_right)
            q_left = left[f"q_{method}"].to_numpy(float)
            q_right = right[f"q_{method}"].to_numpy(float)
            rows.append(
                {
                    "left_test": left_name,
                    "right_test": right_name,
                    "method": method,
                    "neglog10p_pearson": float(stats.pearsonr(x, y).statistic),
                    "neglog10p_spearman": float(stats.spearmanr(x, y).statistic),
                    "median_abs_neglog10p_difference": float(np.median(np.abs(x - y))),
                    "p99_abs_neglog10p_difference": float(np.quantile(np.abs(x - y), 0.99)),
                    "max_abs_neglog10p_difference": float(np.max(np.abs(x - y))),
                    "top_0p1pct_jaccard": _jaccard_top(p_left, p_right, 0.001),
                    "top_1pct_jaccard": _jaccard_top(p_left, p_right, 0.01),
                    "top_5pct_jaccard": _jaccard_top(p_left, p_right, 0.05),
                    "left_BH_FDR_0p10_hits": int(np.sum(q_left < 0.10)),
                    "right_BH_FDR_0p10_hits": int(np.sum(q_right < 0.10)),
                    "FDR_hit_overlap": int(np.sum((q_left < 0.10) & (q_right < 0.10))),
                }
            )
    return pd.DataFrame(rows)


def build_candidate_support(
    candidate_leads: pd.DataFrame,
    score_table: pd.DataFrame,
    lrt_table: pd.DataFrame,
    wald_table: pd.DataFrame,
    conditional_score: pd.DataFrame | None,
    positive_score: pd.DataFrame | None,
    global_classes: pd.DataFrame,
) -> pd.DataFrame:
    if candidate_leads.empty:
        return candidate_leads.copy()
    meta_cols = [
        "lead_rank", "lead_variant_index", "lead_snp", "CHR", "pos", "MAF",
        "candidate_criterion", "window_kb", "r2_threshold", "n_local_members",
        "global_alias_class_index", "test_class_index", "test_class_size",
        "candidate_source", "clump_method",
    ]
    meta_cols = [c for c in meta_cols if c in candidate_leads.columns]
    out = candidate_leads[meta_cols].copy()
    out = out.rename(columns={"lead_variant_index": "variant_index"})

    def add_table(
        source: pd.DataFrame,
        suffix: str,
        methods: list[str],
        effect_traits: list[str] | None = None,
        effect_suffix: str | None = None,
    ) -> None:
        nonlocal out
        keep = ["variant_index"]
        rename: dict[str, str] = {}
        for method in methods:
            for prefix in ["p", "q"]:
                col = f"{prefix}_{method}"
                if col in source.columns:
                    keep.append(col)
                    rename[col] = f"{col}_{suffix}"
        for trait in effect_traits or []:
            for col in [f"beta_{trait}", f"se_{trait}", f"ci_low_{trait}", f"ci_high_{trait}"]:
                if col in source.columns:
                    keep.append(col)
                    rename[col] = col if effect_suffix is None else f"{col}_{effect_suffix}"
        out = out.merge(source[keep].rename(columns=rename), on="variant_index", how="left")

    # GEMMA reports one beta/SE estimate alongside Wald/LRT/score p-values.
    # Store that estimate once, without implying that three distinct effects were fitted.
    add_table(
        score_table, "score", ["INC", "SEV", "ACAT", "SIMES"],
        effect_traits=["INC", "SEV"], effect_suffix=None,
    )
    add_table(lrt_table, "lrt", ["INC", "SEV", "ACAT", "SIMES"])
    add_table(wald_table, "wald_diagnostic", ["INC", "SEV", "ACAT", "SIMES"])
    if conditional_score is not None:
        add_table(
            conditional_score, "conditional_score", ["SEV", "ACAT", "SIMES"],
            effect_traits=["SEV"], effect_suffix="conditional",
        )
    if positive_score is not None:
        add_table(
            positive_score, "positive_only_score", ["SEV", "ACAT", "SIMES"],
            effect_traits=["SEV"], effect_suffix="positive_only",
        )

    aliases = global_classes[[
        "global_alias_class_index", "global_alias_class_id", "class_size",
        "n_chromosomes", "alias_mappings",
    ]].rename(columns={"class_size": "global_alias_class_size"})
    if "global_alias_class_index" in out.columns:
        out = out.merge(aliases, on="global_alias_class_index", how="left")
    return out.sort_values(["p_ACAT_score", "p_SEV_score", "CHR", "pos"], kind="mergesort").reset_index(drop=True)


def _plot_qq(p_values: np.ndarray, label: str, path: Path) -> None:
    p = np.sort(np.clip(np.asarray(p_values, dtype=float), 1e-300, 1.0))
    n = len(p)
    expected = -np.log10((np.arange(1, n + 1) - 0.5) / n)
    observed = -np.log10(p)
    fig, ax = plt.subplots(figsize=(5.5, 5.5))
    ax.scatter(expected, observed, s=5, alpha=0.45)
    limit = max(float(expected.max()), float(observed.max())) * 1.02
    ax.plot([0, limit], [0, limit], linestyle="--", linewidth=1)
    ax.set_xlim(0, limit)
    ax.set_ylim(0, limit)
    ax.set_xlabel("Expected −log10(p)")
    ax.set_ylabel("Observed −log10(p)")
    ax.set_title(label)
    fig.tight_layout()
    fig.savefig(path, dpi=300)
    plt.close(fig)


def _plot_test_diagnostic(comparison: pd.DataFrame, path: Path) -> None:
    subset = comparison.loc[comparison["method"].isin(["INC", "SEV", "ACAT"])]
    if subset.empty:
        return
    labels = [f"{a} vs {b}\n{m}" for a, b, m in subset[["left_test", "right_test", "method"]].itertuples(index=False)]
    values = subset["neglog10p_spearman"].to_numpy(float)
    fig, ax = plt.subplots(figsize=(max(8, len(labels) * 0.8), 4.8))
    ax.bar(np.arange(len(labels)), values)
    ax.set_xticks(np.arange(len(labels)))
    ax.set_xticklabels(labels, rotation=45, ha="right")
    ax.set_ylim(0, 1.05)
    ax.set_ylabel("Spearman correlation of −log10(p)")
    ax.set_title("GEMMA test-statistic sensitivity")
    fig.tight_layout()
    fig.savefig(path, dpi=300)
    plt.close(fig)


def write_workbook(path: Path, sheets: dict[str, pd.DataFrame]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with pd.ExcelWriter(path, engine="openpyxl") as writer:
        for name, table in sheets.items():
            safe_name = name[:31]
            # Excel's row limit is 1,048,576; the largest tables here are 567,758 rows.
            table.to_excel(writer, sheet_name=safe_name, index=False)


def finalize_results(
    project_root: str | Path,
    geno: np.memmap,
    metadata: pd.DataFrame,
    sample_indices: np.ndarray,
    scans: dict[tuple[str, str], pd.DataFrame],
    conditional_scans: dict[tuple[str, str], pd.DataFrame] | None,
    positive_scans: dict[tuple[str, str], pd.DataFrame] | None,
    gff_path: str | Path,
    logger: logging.Logger,
    candidate_top_n: int = 10,
    index_p: float = 1e-4,
    secondary_p: float = 1e-3,
    ld_window_kb: int = 1000,
    ld_r2: float = 0.20,
    power_sample_sizes: tuple[int, ...] = (100, 102),
) -> dict[str, Any]:
    root = Path(project_root)
    tables_dir = root / "06_tables"
    figures_dir = root / "07_figures"
    results_dir = root / "05_results"
    cache_dir = root / "04_intermediate" / "equivalence_cache"
    for d in [tables_dir, figures_dir, results_dir, cache_dir]:
        d.mkdir(parents=True, exist_ok=True)
    start = time.time()

    required = [("loco", "INC"), ("loco", "SEV")]
    for key in required:
        if key not in scans:
            raise ValueError(f"Missing required primary GEMMA scan: {key}")

    test_tables = {
        test: build_test_table(scans[("loco", "INC")], scans[("loco", "SEV")], test)
        for test in [PRIMARY_TEST, SENSITIVITY_TEST, DIAGNOSTIC_TEST]
    }
    score_table = test_tables[PRIMARY_TEST]
    lrt_table = test_tables[SENSITIVITY_TEST]
    wald_table = test_tables[DIAGNOSTIC_TEST]

    conditional_score = None
    if conditional_scans and ("loco", "SEV_COND") in conditional_scans:
        conditional_score = build_test_table(
            scans[("loco", "INC")], conditional_scans[("loco", "SEV_COND")], PRIMARY_TEST
        )
    positive_score = None
    if positive_scans and ("loco", "SEV_POS") in positive_scans:
        positive_score = build_test_table(
            scans[("loco", "INC")], positive_scans[("loco", "SEV_POS")], PRIMARY_TEST
        )

    profile_parts = [profile_test_table(table, test) for test, table in test_tables.items()]
    if conditional_score is not None:
        profile_parts.append(profile_test_table(conditional_score, PRIMARY_TEST, context="conditional_severity"))
    if positive_score is not None:
        profile_parts.append(profile_test_table(positive_score, PRIMARY_TEST, context="positive_only_severity"))
    profile = pd.concat(profile_parts, ignore_index=True)
    test_comparison = compare_tests(test_tables)
    correlations, enrichment = concordance_metrics(
        score_table["p_INC"].to_numpy(float), score_table["p_SEV"].to_numpy(float)
    )
    acat_simes = method_comparison(score_table, 0.10)

    eq = compute_or_load_equivalence(geno, metadata, sample_indices, cache_dir, logger)
    equivalence = build_equivalence_tables(metadata, score_table, eq)
    score_variant = equivalence["variant_table"]
    unique_acat = collapse_for_method(score_variant, "ACAT")
    unique_sev = collapse_for_method(score_variant, "SEV")

    acat_leads, acat_members = alias_aware_ld_clump(
        unique_acat, geno, sample_indices, "ACAT", index_p, secondary_p,
        ld_window_kb, ld_r2, 500, logger,
    )
    sev_leads, sev_members = alias_aware_ld_clump(
        unique_sev, geno, sample_indices, "SEV", index_p, secondary_p,
        ld_window_kb, ld_r2, 500, logger,
    )
    candidate_pool = pd.concat(
        [
            acat_leads.assign(candidate_source="score_ACAT"),
            sev_leads.assign(candidate_source="score_SEV"),
        ],
        ignore_index=True,
    )
    if candidate_pool.empty:
        candidate_pool = unique_acat.nsmallest(candidate_top_n, "p_ACAT").copy()
        candidate_pool["lead_rank"] = np.arange(1, len(candidate_pool) + 1)
        candidate_pool["lead_variant_index"] = candidate_pool["variant_index"]
        candidate_pool["lead_snp"] = candidate_pool["snp"]
        candidate_pool["candidate_source"] = "fallback_top_score_ACAT"
    candidate_pool = candidate_pool.sort_values(
        ["p_ACAT", "p_SEV", "CHR", "pos"], kind="mergesort"
    ).drop_duplicates("lead_variant_index", keep="first").head(candidate_top_n)
    candidates = build_candidate_support(
        candidate_pool, score_table, lrt_table, wald_table,
        conditional_score, positive_score, equivalence["global_alias_classes"],
    )

    genes = parse_gff_genes(gff_path)
    gene_input = candidates.rename(
        columns={"variant_index": "lead_variant_index", "lead_snp": "lead_snp"}
    )
    candidate_genes = nearby_genes(gene_input, genes, 50)

    unique_test_count = int(len(equivalence["test_classes"]))
    power = build_power_table(
        power_sample_sizes,
        [0.05, 0.10, 0.20, 0.30, 0.50],
        [0.80, 0.90],
        {
            "Bonferroni_SNP_rows": 0.05 / len(metadata),
            "Bonferroni_unique_LOCO_tests": 0.05 / unique_test_count,
            "Suggestive_p_1e-5": 1e-5,
        },
        n_covariates=1,
    )

    for test, table in test_tables.items():
        table.to_csv(tables_dir / f"gemma_loco_{test}_all_snps.csv.gz", index=False, compression="gzip")
        table.nsmallest(20000, "p_ACAT").to_csv(
            tables_dir / f"gemma_loco_{test}_top20000_by_ACAT.csv", index=False
        )
    if conditional_score is not None:
        conditional_score.to_csv(tables_dir / "conditional_severity_score_all_snps.csv.gz", index=False, compression="gzip")
    if positive_score is not None:
        positive_score.to_csv(tables_dir / "positive_only_severity_score_all_snps.csv.gz", index=False, compression="gzip")

    # Preserve full equivalence mappings as compressed CSVs; keep the review
    # workbook compact enough to open reliably in Excel.
    equivalence["test_classes"].to_csv(
        tables_dir / "Unique_Test_Classes_full.csv.gz", index=False, compression="gzip"
    )
    equivalence["global_alias_classes"].to_csv(
        tables_dir / "Global_Alias_Classes_full.csv.gz", index=False, compression="gzip"
    )
    equivalence["duplicate_mapping"].to_csv(
        tables_dir / "Duplicate_Marker_Mapping_full.csv.gz", index=False, compression="gzip"
    )

    compact_tables = {
        "Inference_Profile": profile,
        "Test_Sensitivity": test_comparison,
        "Score_Concordance": correlations,
        "Score_Enrichment": enrichment,
        "ACAT_vs_Simes": acat_simes,
        "Top_Ranked_Candidates": candidates,
        "Candidate_Genes": candidate_genes,
        "Power": power,
        "Equivalence_Summary": equivalence["summary"],
        "Top_Unique_Test_Classes": equivalence["test_classes"].nsmallest(20000, "minimum_p_ACAT"),
        "Largest_Alias_Classes": equivalence["global_alias_classes"].nlargest(5000, "class_size"),
        "ACAT_Clump_Leads": acat_leads,
        "ACAT_Clump_Members": acat_members,
        "SEV_Clump_Leads": sev_leads,
        "SEV_Clump_Members": sev_members,
    }
    write_workbook(results_dir / "LeafBlight_MultiEnvironment_Analysis_Results.xlsx", compact_tables)
    for name, table in compact_tables.items():
        table.to_csv(tables_dir / f"{name}.csv", index=False)

    _plot_qq(score_table["p_ACAT"].to_numpy(float), "Primary GEMMA LOCO score ACAT", figures_dir / "Figure_Score_ACAT_QQ.png")
    _plot_qq(score_table["p_SEV"].to_numpy(float), "Primary GEMMA LOCO score severity", figures_dir / "Figure_Score_SEV_QQ.png")
    _plot_test_diagnostic(test_comparison, figures_dir / "Figure_Test_Sensitivity.png")

    summary = {
        "primary_engine": "official GEMMA 0.98.5",
        "primary_context": "LOCO kinship",
        "primary_test": PRIMARY_TEST,
        "confirmatory_test": SENSITIVITY_TEST,
        "diagnostic_test": DIAGNOSTIC_TEST,
        "n_primary": int(len(sample_indices)),
        "n_variants": int(len(metadata)),
        "unique_loco_test_classes": unique_test_count,
        "score_acat_lambda_gc": float(lambda_gc_from_p(score_table["p_ACAT"])),
        "score_acat_minimum_p": float(score_table["p_ACAT"].min()),
        "score_acat_fdr_hits": int((score_table["q_ACAT"] < 0.10).sum()),
        "score_simes_fdr_hits": int((score_table["q_SIMES"] < 0.10).sum()),
        "score_inc_fdr_hits": int((score_table["q_INC"] < 0.10).sum()),
        "score_sev_fdr_hits": int((score_table["q_SEV"] < 0.10).sum()),
        "lrt_acat_fdr_hits": int((lrt_table["q_ACAT"] < 0.10).sum()),
        "wald_acat_fdr_hits_diagnostic": int((wald_table["q_ACAT"] < 0.10).sum()),
        "exploratory_candidate_regions": int(len(candidates)),
        "candidate_selection_rule": (
            f"alias-aware LD clumping of score-test ACAT/SEV signals at index p<={index_p:g}, "
            f"secondary p<={secondary_p:g}, +/-{ld_window_kb} kb, r2>={ld_r2:.2f}; "
            f"top {candidate_top_n} retained for summary"
        ),
        "interpretation_note": (
            "Only score-test and LRT results are inferential. Wald results are retained solely "
            "as a small-sample variance-component diagnostic. Top-ranked regions are exploratory, "
            "not genome-wide significant QTL."
        ),
        "elapsed": elapsed(start),
    }
    save_json(results_dir / "analysis_manifest.json", summary)
    lines = [
        "SORGHUM LEAF BLIGHT MULTI-ENVIRONMENT ANALYSIS",
        "=" * 52,
        "",
        f"Primary engine: {summary['primary_engine']}",
        f"Primary test: {summary['primary_test']} (confirmatory: {summary['confirmatory_test']}; diagnostic: {summary['diagnostic_test']})",
        f"Primary n: {summary['n_primary']}",
        f"Variants: {summary['n_variants']:,}",
        f"Unique LOCO test classes: {summary['unique_loco_test_classes']:,}",
        f"Score ACAT lambda_GC: {summary['score_acat_lambda_gc']:.4f}",
        f"Score ACAT minimum p: {summary['score_acat_minimum_p']:.6g}",
        f"Score ACAT BH-FDR 10% hits: {summary['score_acat_fdr_hits']}",
        f"Score Simes BH-FDR 10% hits: {summary['score_simes_fdr_hits']}",
        f"Score incidence BH-FDR 10% hits: {summary['score_inc_fdr_hits']}",
        f"Score severity BH-FDR 10% hits: {summary['score_sev_fdr_hits']}",
        f"LRT ACAT BH-FDR 10% hits: {summary['lrt_acat_fdr_hits']}",
        f"Wald ACAT BH-FDR 10% hits (diagnostic only): {summary['wald_acat_fdr_hits_diagnostic']}",
        "",
        "INTERPRETATION",
        "----------",
        summary["interpretation_note"],
    ]
    (results_dir / "ANALYSIS_SUMMARY.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
    logger.info("Final result consolidation complete (%s)", elapsed(start))
    return summary
