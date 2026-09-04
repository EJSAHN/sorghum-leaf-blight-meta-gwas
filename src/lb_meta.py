from __future__ import annotations

import gzip
import logging
import math
import re
from pathlib import Path
from typing import Any, Protocol

import numpy as np
import pandas as pd
from scipy import stats

from lb_common import bh_fdr, safe_neglog10
from lb_genotype import chromosome_to_int, load_variant_dosages


class AssociationResult(Protocol):
    beta: np.ndarray
    se: np.ndarray
    p: np.ndarray
    n: int
    sample_ids: list[str]


def acat_pair(p1: np.ndarray, p2: np.ndarray) -> np.ndarray:
    a = np.clip(np.asarray(p1, dtype=float), 1e-15, 1 - 1e-15)
    b = np.clip(np.asarray(p2, dtype=float), 1e-15, 1 - 1e-15)
    t = 0.5 * (np.tan((0.5 - a) * np.pi) + np.tan((0.5 - b) * np.pi))
    return np.clip(0.5 - np.arctan(t) / np.pi, 1e-300, 1.0)


def simes_pair(p1: np.ndarray, p2: np.ndarray) -> np.ndarray:
    a = np.asarray(p1, dtype=float)
    b = np.asarray(p2, dtype=float)
    lo = np.minimum(a, b)
    hi = np.maximum(a, b)
    return np.clip(np.minimum(2.0 * lo, hi), 1e-300, 1.0)


def build_meta_table(metadata: pd.DataFrame, inc: AssociationResult, sev: AssociationResult) -> pd.DataFrame:
    if inc.sample_ids != sev.sample_ids:
        raise ValueError("Incidence and severity association scans used different sample sets")
    out = metadata[["variant_index", "snp", "CHR", "pos", "MAF", "missing_rate"]].copy()
    out["beta_INC"] = inc.beta
    out["se_INC"] = inc.se
    out["p_INC"] = inc.p
    out["beta_SEV"] = sev.beta
    out["se_SEV"] = sev.se
    out["p_SEV"] = sev.p
    out["p_ACAT"] = acat_pair(inc.p, sev.p)
    out["q_ACAT"] = bh_fdr(out["p_ACAT"])
    out["p_SIMES"] = simes_pair(inc.p, sev.p)
    out["q_SIMES"] = bh_fdr(out["p_SIMES"])
    return out


def association_table(metadata: pd.DataFrame, result: AssociationResult, trait_label: str) -> pd.DataFrame:
    out = metadata[["variant_index", "snp", "CHR", "pos", "MAF", "missing_rate"]].copy()
    out["beta"] = result.beta
    out["se"] = result.se
    out["ci_low"] = result.beta - 1.96 * result.se
    out["ci_high"] = result.beta + 1.96 * result.se
    out["p"] = result.p
    out["q_BH"] = bh_fdr(result.p)
    out["n"] = result.n
    out["trait"] = trait_label
    return out


def fixed_distance_leads(meta: pd.DataFrame, p_col: str, q_col: str, fdr: float, bin_kb: int) -> pd.DataFrame:
    candidates = meta.loc[meta[q_col] < fdr].sort_values(p_col, kind="mergesort")
    if candidates.empty:
        return pd.DataFrame(columns=["lead_rank", "lead_variant_index", "lead_snp", "CHR", "pos", p_col, q_col])
    bin_size = int(bin_kb * 1000)
    seen: set[tuple[int, int]] = set()
    rows = []
    for idx, row in candidates.iterrows():
        key = (int(row["CHR"]), int(row["pos"] // bin_size))
        if key in seen:
            continue
        seen.add(key)
        item = row.to_dict()
        item["lead_variant_index"] = int(idx)
        item["lead_snp"] = row["snp"]
        rows.append(item)
    result = pd.DataFrame(rows).reset_index(drop=True)
    result.insert(0, "lead_rank", np.arange(1, len(result) + 1))
    return result


def _r2_against_lead(dosages: np.ndarray, lead_row: int) -> np.ndarray:
    centered = dosages - dosages.mean(axis=1, keepdims=True)
    lead = centered[lead_row]
    lead_ss = float(lead @ lead)
    if lead_ss <= 1e-12:
        return np.zeros(len(dosages), dtype=float)
    ss = np.sum(centered * centered, axis=1)
    numer = centered @ lead
    denom = ss * lead_ss
    return np.clip(np.divide(numer * numer, denom, out=np.zeros_like(numer), where=denom > 1e-12), 0, 1)


def _candidate_set(
    meta: pd.DataFrame,
    p_col: str,
    q_col: str | None,
    fdr: float | None,
    p_threshold: float | None,
    fallback_top_n: int | None,
) -> tuple[pd.DataFrame, str]:
    if q_col and fdr is not None:
        candidates = meta.loc[meta[q_col] < fdr].sort_values(p_col, kind="mergesort").copy()
        criterion = f"{q_col} < {fdr}"
    elif p_threshold is not None:
        candidates = meta.loc[meta[p_col] <= p_threshold].sort_values(p_col, kind="mergesort").copy()
        criterion = f"{p_col} <= {p_threshold}"
    else:
        candidates = pd.DataFrame()
        criterion = "none"
    if candidates.empty and fallback_top_n:
        candidates = meta.nsmallest(int(fallback_top_n), p_col).copy()
        criterion = f"fallback top {fallback_top_n} by {p_col}"
    return candidates, criterion


def ld_clump(
    meta: pd.DataFrame,
    geno: np.memmap,
    sample_indices: np.ndarray,
    p_col: str,
    q_col: str | None,
    fdr: float | None,
    window_kb: int,
    r2_threshold: float,
    logger: logging.Logger,
    p_threshold: float | None = None,
    fallback_top_n: int | None = None,
    secondary_p: float | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    candidates, criterion = _candidate_set(meta, p_col, q_col, fdr, p_threshold, fallback_top_n)
    lead_columns = [
        "lead_rank", "lead_variant_index", "lead_snp", "CHR", "pos", p_col,
        "MAF", "clump_member_count", "clump_window_kb", "clump_r2_threshold",
        "selection_criterion",
    ]
    if q_col:
        lead_columns.insert(6, q_col)
    member_columns = [
        "lead_snp", "member_snp", "member_variant_index", "CHR", "pos", p_col,
        "r2_to_lead", "distance_bp",
    ]
    if candidates.empty:
        logger.info("LD clumping %s: no candidates under %s", p_col, criterion)
        return pd.DataFrame(columns=lead_columns), pd.DataFrame(columns=member_columns)

    candidates = candidates.copy()
    candidate_indices = candidates.index.to_numpy(dtype=int)
    assigned = np.zeros(len(candidates), dtype=bool)
    positions = candidates["pos"].to_numpy(dtype=int)
    chroms = candidates["CHR"].to_numpy(dtype=int)
    window = int(window_kb * 1000)
    lead_rows: list[dict[str, Any]] = []
    member_rows: list[dict[str, Any]] = []
    logger.info(
        "LD clumping %s: %d candidates (%s), window=%dkb, r2>=%.2f",
        p_col, len(candidates), criterion, window_kb, r2_threshold,
    )

    for order_pos in range(len(candidates)):
        if assigned[order_pos]:
            continue
        lead_meta_index = int(candidate_indices[order_pos])
        lead = candidates.iloc[order_pos]
        local_mask = (~assigned) & (chroms == int(lead["CHR"])) & (np.abs(positions - int(lead["pos"])) <= window)
        if secondary_p is not None:
            local_mask &= candidates[p_col].to_numpy(dtype=float) <= float(secondary_p)
        local_positions = np.flatnonzero(local_mask)
        # The lead itself must always be included even if secondary_p is stricter than its p-value.
        if order_pos not in local_positions:
            local_positions = np.unique(np.append(local_positions, order_pos))
        local_meta_indices = candidate_indices[local_positions]
        dosage = load_variant_dosages(geno, local_meta_indices)[:, sample_indices]
        lead_local = int(np.flatnonzero(local_meta_indices == lead_meta_index)[0])
        r2 = _r2_against_lead(dosage, lead_local)
        clumped = local_positions[r2 >= r2_threshold]
        assigned[clumped] = True

        row = {
            "lead_rank": len(lead_rows) + 1,
            "lead_variant_index": lead_meta_index,
            "lead_snp": lead["snp"],
            "CHR": int(lead["CHR"]),
            "pos": int(lead["pos"]),
            p_col: float(lead[p_col]),
            "MAF": float(lead["MAF"]),
            "clump_member_count": int(len(clumped)),
            "clump_window_kb": int(window_kb),
            "clump_r2_threshold": float(r2_threshold),
            "selection_criterion": criterion,
        }
        if q_col:
            row[q_col] = float(lead[q_col])
        lead_rows.append(row)

        for local_pos, r2_value in zip(local_positions, r2):
            if r2_value < r2_threshold:
                continue
            member = candidates.iloc[int(local_pos)]
            member_rows.append(
                {
                    "lead_snp": lead["snp"],
                    "member_snp": member["snp"],
                    "member_variant_index": int(candidate_indices[int(local_pos)]),
                    "CHR": int(member["CHR"]),
                    "pos": int(member["pos"]),
                    p_col: float(member[p_col]),
                    "r2_to_lead": float(r2_value),
                    "distance_bp": int(member["pos"] - lead["pos"]),
                }
            )
    return pd.DataFrame(lead_rows, columns=lead_columns), pd.DataFrame(member_rows, columns=member_columns)


def add_trait_effects(
    leads: pd.DataFrame,
    results: dict[str, AssociationResult],
) -> pd.DataFrame:
    out = leads.copy()
    if "lead_variant_index" not in out.columns:
        out["lead_variant_index"] = pd.Series(dtype="int64")
    for label, result in results.items():
        for suffix in ["beta", "se", "ci_low", "ci_high", "p"]:
            out[f"{suffix}_{label}"] = pd.Series(dtype="float64")
        out[f"n_{label}"] = pd.Series(dtype="int64")
    if out.empty:
        return out
    idx = out["lead_variant_index"].to_numpy(dtype=int)
    for label, result in results.items():
        beta = result.beta[idx]
        se = result.se[idx]
        out[f"beta_{label}"] = beta
        out[f"se_{label}"] = se
        out[f"ci_low_{label}"] = beta - 1.96 * se
        out[f"ci_high_{label}"] = beta + 1.96 * se
        out[f"p_{label}"] = result.p[idx]
        out[f"n_{label}"] = result.n
    return out


def concordance_metrics(p_inc: np.ndarray, p_sev: np.ndarray) -> tuple[pd.DataFrame, pd.DataFrame]:
    x = safe_neglog10(p_inc)
    y = safe_neglog10(p_sev)
    pearson_r, pearson_p = stats.pearsonr(x, y)
    spearman_r, spearman_p = stats.spearmanr(x, y)
    correlations = pd.DataFrame(
        [
            {"metric": "Pearson", "coefficient": pearson_r, "p_value": pearson_p, "n_snps": len(x)},
            {"metric": "Spearman", "coefficient": spearman_r, "p_value": spearman_p, "n_snps": len(x)},
        ]
    )
    rows = []
    n = len(x)
    for fraction in [0.01, 0.05, 0.10]:
        k = max(1, int(round(n * fraction)))
        idx_inc = set(np.argpartition(p_inc, k - 1)[:k].tolist())
        idx_sev = set(np.argpartition(p_sev, k - 1)[:k].tolist())
        overlap = len(idx_inc & idx_sev)
        union = len(idx_inc | idx_sev)
        expected = (k * k) / n
        table = np.array([[overlap, k - overlap], [k - overlap, n - 2 * k + overlap]], dtype=int)
        odds_ratio, fisher_p = stats.fisher_exact(table, alternative="greater")
        rows.append(
            {
                "top_fraction": fraction,
                "top_n_each": k,
                "overlap": overlap,
                "jaccard": overlap / union if union else np.nan,
                "expected_overlap": expected,
                "fold_enrichment": overlap / expected if expected else np.nan,
                "odds_ratio": odds_ratio,
                "fisher_p": fisher_p,
            }
        )
    return correlations, pd.DataFrame(rows)


def method_comparison(meta: pd.DataFrame, fdr: float) -> pd.DataFrame:
    rank_r, rank_p = stats.spearmanr(meta["p_ACAT"], meta["p_SIMES"])
    acat_sig = set(meta.index[meta["q_ACAT"] < fdr])
    simes_sig = set(meta.index[meta["q_SIMES"] < fdr])
    rows = [
        {
            "comparison": "Genome-wide p-value rank correlation",
            "value": float(rank_r),
            "secondary_value": float(rank_p),
            "details": "Spearman correlation; secondary=p-value",
        },
        {"comparison": "FDR-significant SNP count ACAT", "value": len(acat_sig), "secondary_value": np.nan, "details": f"q < {fdr}"},
        {"comparison": "FDR-significant SNP count Simes", "value": len(simes_sig), "secondary_value": np.nan, "details": f"q < {fdr}"},
    ]
    for fraction in [0.01, 0.05, 0.10]:
        n = len(meta)
        k = max(1, int(round(n * fraction)))
        a = set(np.argpartition(meta["p_ACAT"].to_numpy(), k - 1)[:k].tolist())
        s = set(np.argpartition(meta["p_SIMES"].to_numpy(), k - 1)[:k].tolist())
        rows.append(
            {
                "comparison": f"Top {int(fraction*100)}% overlap",
                "value": len(a & s),
                "secondary_value": len(a & s) / len(a | s),
                "details": "secondary=Jaccard",
            }
        )
    return pd.DataFrame(rows)


def lead_window_support(
    leads: pd.DataFrame,
    meta: pd.DataFrame,
    window_kb: int,
    p_threshold: float,
) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    window = int(window_kb * 1000)
    for _, lead in leads.iterrows():
        region = meta.loc[
            (meta["CHR"].astype(int) == int(lead["CHR"]))
            & meta["pos"].between(int(lead["pos"]) - window, int(lead["pos"]) + window)
        ]
        min_inc = float(region["p_INC"].min()) if not region.empty else np.nan
        min_sev = float(region["p_SEV"].min()) if not region.empty else np.nan
        rows.append(
            {
                "lead_snp": lead["lead_snp"],
                "CHR": int(lead["CHR"]),
                "pos": int(lead["pos"]),
                "min_p_INC_window": min_inc,
                "min_p_SEV_window": min_sev,
                "INC_support": bool(min_inc < p_threshold),
                "SEV_support": bool(min_sev < p_threshold),
                "both_support": bool(min_inc < p_threshold and min_sev < p_threshold),
                "window_kb": int(window_kb),
                "p_threshold": float(p_threshold),
            }
        )
    return pd.DataFrame(rows)


def compare_scans(full: pd.DataFrame, comparator: pd.DataFrame, label: str) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for p_col in ["p_INC", "p_SEV", "p_ACAT", "p_SIMES"]:
        rho, p_value = stats.spearmanr(full[p_col], comparator[p_col])
        rows.append(
            {
                "comparison": label,
                "metric": p_col,
                "spearman_rho": rho,
                "p_value": p_value,
                "full_min_p": float(full[p_col].min()),
                "comparator_min_p": float(comparator[p_col].min()),
            }
        )
    return pd.DataFrame(rows)


def leave_one_location_out_stability(
    leads: pd.DataFrame,
    full_meta: pd.DataFrame,
    loo_meta: dict[str, pd.DataFrame],
    window_kb: int,
) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    window = int(window_kb * 1000)
    for _, lead in leads.iterrows():
        idx = int(lead["lead_variant_index"])
        for location, meta in loo_meta.items():
            region = meta.loc[
                (meta["CHR"].astype(int) == int(lead["CHR"]))
                & meta["pos"].between(int(lead["pos"]) - window, int(lead["pos"]) + window)
            ]
            rows.append(
                {
                    "lead_snp": lead["lead_snp"],
                    "CHR": int(lead["CHR"]),
                    "pos": int(lead["pos"]),
                    "omitted_location": location,
                    "full_p_ACAT": float(full_meta.loc[idx, "p_ACAT"]),
                    "loo_exact_p_INC": float(meta.loc[idx, "p_INC"]),
                    "loo_exact_p_SEV": float(meta.loc[idx, "p_SEV"]),
                    "loo_exact_p_ACAT": float(meta.loc[idx, "p_ACAT"]),
                    "loo_min_p_ACAT_window": float(region["p_ACAT"].min()) if not region.empty else np.nan,
                    "loo_min_p_INC_window": float(region["p_INC"].min()) if not region.empty else np.nan,
                    "loo_min_p_SEV_window": float(region["p_SEV"].min()) if not region.empty else np.nan,
                    "beta_INC_direction_same": bool(np.sign(meta.loc[idx, "beta_INC"]) == np.sign(full_meta.loc[idx, "beta_INC"])),
                    "beta_SEV_direction_same": bool(np.sign(meta.loc[idx, "beta_SEV"]) == np.sign(full_meta.loc[idx, "beta_SEV"])),
                }
            )
    return pd.DataFrame(rows)


def parse_gff_genes(gff_path: str | Path) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    opener = gzip.open if str(gff_path).lower().endswith(".gz") else open
    with opener(gff_path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2].lower() != "gene":
                continue
            chrom = chromosome_to_int(parts[0])
            if chrom is None:
                continue
            try:
                start, end = int(parts[3]), int(parts[4])
            except ValueError:
                continue
            attrs = {}
            for token in parts[8].split(";"):
                if "=" in token:
                    key, value = token.split("=", 1)
                    attrs[key.strip()] = value.strip()
            gene = attrs.get("Name") or attrs.get("ID") or attrs.get("gene_id") or ""
            gene = re.sub(r"^(gene:|gene-)", "", gene)
            rows.append({"CHR": chrom, "start": start, "end": end, "gene": gene, "attributes": parts[8]})
    return pd.DataFrame(rows)


def nearby_genes(leads: pd.DataFrame, genes: pd.DataFrame, flank_kb: int) -> pd.DataFrame:
    columns = ["lead_snp", "CHR", "lead_pos", "gene", "gene_start", "gene_end", "distance_bp", "within_flank"]
    if genes.empty or leads.empty:
        return pd.DataFrame(columns=columns)
    flank = int(flank_kb * 1000)
    rows: list[dict[str, Any]] = []
    for _, lead in leads.iterrows():
        chrom_genes = genes.loc[genes["CHR"].astype(int) == int(lead["CHR"])].copy()
        if chrom_genes.empty:
            continue
        pos = int(lead["pos"])
        distance = np.where(
            pos < chrom_genes["start"], chrom_genes["start"] - pos,
            np.where(pos > chrom_genes["end"], pos - chrom_genes["end"], 0),
        )
        chrom_genes["distance_bp"] = distance.astype(int)
        selected = chrom_genes.loc[chrom_genes["distance_bp"] <= flank].copy()
        if selected.empty:
            selected = chrom_genes.nsmallest(1, "distance_bp").copy()
        for _, gene in selected.sort_values("distance_bp").iterrows():
            rows.append(
                {
                    "lead_snp": lead["lead_snp"],
                    "CHR": int(lead["CHR"]),
                    "lead_pos": pos,
                    "gene": gene["gene"],
                    "gene_start": int(gene["start"]),
                    "gene_end": int(gene["end"]),
                    "distance_bp": int(gene["distance_bp"]),
                    "within_flank": bool(gene["distance_bp"] <= flank),
                }
            )
    return pd.DataFrame(rows, columns=columns)
