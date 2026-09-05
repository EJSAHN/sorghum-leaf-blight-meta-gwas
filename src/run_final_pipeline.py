from __future__ import annotations

import argparse
import json
import shutil
import sys
import time
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy import stats

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from lb_common import (
    ensure_directories,
    load_config,
    save_json,
    setup_logger,
    sha256_file,
    software_versions,
    git_commit,
)
from lb_finalize import build_test_table, finalize_results
from lb_genotype import build_or_load_genotype_cache
from lb_gemma import (
    load_scan_tables,
    prepare_kinship_text,
    run_gemma_scans,
    save_scan_tables,
    write_covariates,
)
from lb_phenotype import (
    build_accession_traits,
    build_canonical_phenotype,
    build_leave_one_location_out_traits,
    build_location_specific_traits,
    build_population_table,
)
from lb_plink import write_unique_plink_files
from lb_structure import build_or_load_population_structure, select_ld_pruned_structure_markers

VERSION = "2.0.1"


def require(path: Path, label: str) -> Path:
    if not path.exists():
        raise FileNotFoundError(f"Required {label} not found: {path}")
    return path


def read_traits(paths: dict[str, Path]) -> pd.DataFrame:
    path = paths["intermediate"] / "phenotype_all_accession_traits.csv"
    if not path.exists():
        raise FileNotFoundError(f"Prepared phenotype traits not found: {path}")
    return pd.read_csv(path)


def prepare_stage(paths: dict[str, Path], cfg: dict[str, Any], logger, force: bool) -> dict[str, Any]:
    inputs = paths["inputs"]
    phenotype = require(inputs / "Phenotype.xlsx", "corrected phenotype workbook")
    vcf = require(inputs / "700k.vcf", "VCF")
    gff = require(inputs / "Sbicolor_454_v3.1.1.gene.gff3.gz", "Sorghum GFF")
    del phenotype, gff

    expected_hash = str(cfg.get("expected_vcf_sha256", "")).upper().strip()
    observed_hash = sha256_file(vcf)
    if expected_hash and observed_hash != expected_hash:
        raise ValueError(f"VCF SHA-256 mismatch: expected {expected_hash}, observed {observed_hash}")

    canonical, phenotype_qc = build_canonical_phenotype(inputs, paths["intermediate"], logger)
    traits, phenotype_models = build_accession_traits(canonical, paths["intermediate"], logger)
    locations = [str(x) for x in cfg["locations"]]
    lolo_traits, lolo_diagnostics = build_leave_one_location_out_traits(
        canonical, locations, paths["intermediate"], logger
    )
    location_traits = build_location_specific_traits(canonical, locations, paths["intermediate"])
    traits = traits.merge(lolo_traits, on="ID_std", how="left").merge(location_traits, on="ID_std", how="left")
    traits_path = paths["intermediate"] / "phenotype_all_accession_traits.csv"
    traits.to_csv(traits_path, index=False)

    expected_rows = cfg.get("expected_canonical_rows")
    expected_ids = cfg.get("expected_phenotype_ids")
    if expected_rows is not None and len(canonical) != int(expected_rows):
        raise ValueError(f"Canonical rows={len(canonical)}, expected {expected_rows}")
    if expected_ids is not None and canonical["ID_std"].nunique() != int(expected_ids):
        raise ValueError(
            f"Phenotype IDs={canonical['ID_std'].nunique()}, expected {expected_ids}"
        )

    geno, metadata, matched_ids, genotype_manifest = build_or_load_genotype_cache(
        vcf,
        traits["ID_std"].astype(str).tolist(),
        paths["cache"],
        float(cfg["maf_min"]),
        float(cfg["missing_max"]),
        logger,
    )
    if cfg.get("expected_matched_ids") is not None and len(matched_ids) != int(cfg["expected_matched_ids"]):
        raise ValueError(f"Matched IDs={len(matched_ids)}, expected {cfg['expected_matched_ids']}")
    if cfg.get("expected_post_qc_snps") is not None and len(metadata) != int(cfg["expected_post_qc_snps"]):
        raise ValueError(f"Post-QC SNPs={len(metadata)}, expected {cfg['expected_post_qc_snps']}")

    population_table, primary_ids, full_ids = build_population_table(
        traits,
        matched_ids,
        int(cfg["primary_min_locations"]),
        int(cfg["expected_primary_n"]) if cfg.get("expected_primary_n") is not None else None,
        paths["intermediate"],
        logger,
    )

    all_indices = np.arange(len(full_ids), dtype=int)
    common_markers = select_ld_pruned_structure_markers(
        geno,
        metadata,
        all_indices,
        paths["structure"],
        int(cfg["structure_physical_thin_kb"]),
        int(cfg["structure_ld_window_kb"]),
        float(cfg["structure_ld_r2"]),
        int(cfg["structure_max_candidate_markers"]),
        int(cfg["structure_min_markers"]),
        logger,
        force=force,
    )
    n_pc_save = max(5, int(cfg.get("pc_sensitivity_count", 3)))
    primary_structure = build_or_load_population_structure(
        "PRIMARY_N100",
        geno,
        metadata,
        matched_ids,
        primary_ids,
        common_markers,
        paths["structure"],
        n_pc_save,
        int(cfg["variant_block_size"]),
        logger,
        force=force,
    )
    full_structure = build_or_load_population_structure(
        "FULL_N102",
        geno,
        metadata,
        matched_ids,
        full_ids,
        common_markers,
        paths["structure"],
        n_pc_save,
        int(cfg["variant_block_size"]),
        logger,
        force=force,
    )

    indexed = traits.set_index("ID_std")

    def phenotype_vector(ids: list[str], column: str) -> np.ndarray:
        values = pd.to_numeric(indexed.reindex(ids)[column], errors="coerce")
        if values.isna().any():
            missing = values.index[values.isna()].tolist()
            raise ValueError(f"Missing {column} values for: {missing[:10]}")
        return values.to_numpy(float)

    primary_phenotypes: dict[str, np.ndarray] = {
        "INC": phenotype_vector(primary_ids, "INC_ADJ_INT"),
        "SEV": phenotype_vector(primary_ids, "SEV_ADJ_INT"),
        "SEV_COND": phenotype_vector(primary_ids, "SEV_COND_ADJ_INT"),
        "SEV_POS": phenotype_vector(primary_ids, "SEV_POS_ADJ_INT"),
    }
    for location in locations:
        label = "".join(ch if ch.isalnum() else "_" for ch in location).strip("_").upper()
        primary_phenotypes[f"INC_LOO_{label}"] = phenotype_vector(primary_ids, f"INC_LOO_{label}_INT")
        primary_phenotypes[f"SEV_LOO_{label}"] = phenotype_vector(primary_ids, f"SEV_LOO_{label}_INT")

    primary_plink = write_unique_plink_files(
        geno,
        metadata,
        primary_structure.sample_indices,
        primary_ids,
        primary_phenotypes,
        paths["intermediate"] / "plink_primary_n100",
        common_markers["variant_index"].to_numpy(dtype=int),
        logger,
        force=force,
    )
    primary_kinship = prepare_kinship_text(
        primary_structure,
        paths["intermediate"] / "kinship_primary_n100",
        logger,
    )
    pc3_path = write_covariates(
        primary_structure,
        paths["intermediate"] / "kinship_primary_n100" / "covariates_pc3.txt",
        3,
    )

    full_phenotypes = {
        "INC": phenotype_vector(full_ids, "INC_ADJ_INT"),
        "SEV": phenotype_vector(full_ids, "SEV_ADJ_INT"),
    }
    full_plink = write_unique_plink_files(
        geno,
        metadata,
        full_structure.sample_indices,
        full_ids,
        full_phenotypes,
        paths["intermediate"] / "plink_full_n102",
        common_markers["variant_index"].to_numpy(dtype=int),
        logger,
        force=force,
    )
    full_kinship = prepare_kinship_text(
        full_structure,
        paths["intermediate"] / "kinship_full_n102",
        logger,
    )

    prep_manifest = {
        "pipeline_version": VERSION,
        "vcf_sha256": observed_hash,
        "canonical_rows": int(len(canonical)),
        "phenotype_ids": int(canonical["ID_std"].nunique()),
        "matched_ids": int(len(matched_ids)),
        "primary_n": int(len(primary_ids)),
        "full_n": int(len(full_ids)),
        "post_qc_snps": int(len(metadata)),
        "structure_markers": int(len(common_markers)),
        "trait_file": str(traits_path),
        "primary_mapping": str(primary_plink["variant_mapping"]),
        "full_mapping": str(full_plink["variant_mapping"]),
        "pc3_covariates": str(pc3_path),
        "phenotype_qc": phenotype_qc,
        "genotype_manifest": genotype_manifest,
        "software_versions": software_versions(),
        "source_git_commit": git_commit(Path(__file__).resolve().parents[1]),
        "config_sha256": sha256_file(paths["config"] / "final_config.json"),
    }
    save_json(paths["intermediate"] / "prepare_manifest.json", prep_manifest)
    logger.info("Preparation stage complete")
    # Explicitly release large memmaps.
    del primary_structure.spectral_genotypes
    del full_structure.spectral_genotypes
    del geno
    return prep_manifest


def _prefixes(plink_dir: Path, labels: list[str]) -> dict[str, Path]:
    return {label: plink_dir / f"leaf_blight_{label.lower()}" for label in labels}


def _chromosome_lists(plink_dir: Path) -> dict[int, Path]:
    return {chrom: plink_dir / "chromosome_snp_lists" / f"chr{chrom}_snps.txt" for chrom in range(1, 11)}


def _kinship_paths(kinship_dir: Path) -> dict[str, Path]:
    return {
        "global": kinship_dir / "global_kinship.cXX.txt",
        **{f"chr{chrom}": kinship_dir / f"loco_chr{chrom}.cXX.txt" for chrom in range(1, 11)},
    }


def gemma_stage(paths: dict[str, Path], cfg: dict[str, Any], gemma_wsl: str, logger, force: bool) -> None:
    if not gemma_wsl:
        raise ValueError("--gemma-wsl is required for the GEMMA stage")

    primary_plink_dir = paths["intermediate"] / "plink_primary_n100"
    primary_mapping = require(primary_plink_dir / "variant_id_mapping.csv.gz", "primary variant mapping")
    primary_labels = ["INC", "SEV", "SEV_COND", "SEV_POS"]
    primary_scans, primary_status = run_gemma_scans(
        gemma_wsl,
        _prefixes(primary_plink_dir, primary_labels),
        primary_mapping,
        _kinship_paths(paths["intermediate"] / "kinship_primary_n100"),
        _chromosome_lists(primary_plink_dir),
        paths["intermediate"] / "gemma_primary_n100",
        logger,
        lmm_mode=4,
        run_global=True,
        run_loco=True,
        force=force,
    )
    save_scan_tables(primary_scans, paths["results"] / "gemma_primary_n100")
    primary_status.to_csv(paths["tables"] / "gemma_primary_n100_status.csv", index=False)

    if bool(cfg.get("run_pc3_sensitivity", True)):
        pc3_scans, pc3_status = run_gemma_scans(
            gemma_wsl,
            _prefixes(primary_plink_dir, ["INC", "SEV"]),
            primary_mapping,
            _kinship_paths(paths["intermediate"] / "kinship_primary_n100"),
            _chromosome_lists(primary_plink_dir),
            paths["intermediate"] / "gemma_primary_n100_pc3",
            logger,
            lmm_mode=4,
            run_global=False,
            run_loco=True,
            covariates=paths["intermediate"] / "kinship_primary_n100" / "covariates_pc3.txt",
            force=force,
        )
        save_scan_tables(pc3_scans, paths["results"] / "gemma_primary_n100_pc3")
        pc3_status.to_csv(paths["tables"] / "gemma_primary_n100_pc3_status.csv", index=False)

    if bool(cfg.get("run_full_n102_sensitivity", True)):
        full_plink_dir = paths["intermediate"] / "plink_full_n102"
        full_mapping = require(full_plink_dir / "variant_id_mapping.csv.gz", "full-n102 variant mapping")
        full_scans, full_status = run_gemma_scans(
            gemma_wsl,
            _prefixes(full_plink_dir, ["INC", "SEV"]),
            full_mapping,
            _kinship_paths(paths["intermediate"] / "kinship_full_n102"),
            _chromosome_lists(full_plink_dir),
            paths["intermediate"] / "gemma_full_n102",
            logger,
            lmm_mode=4,
            run_global=False,
            run_loco=True,
            force=force,
        )
        save_scan_tables(full_scans, paths["results"] / "gemma_full_n102")
        full_status.to_csv(paths["tables"] / "gemma_full_n102_status.csv", index=False)

    if bool(cfg.get("run_leave_one_location_out_score", True)):
        lolo_labels = []
        for location in cfg["locations"]:
            label = "".join(ch if ch.isalnum() else "_" for ch in str(location)).strip("_").upper()
            lolo_labels.extend([f"INC_LOO_{label}", f"SEV_LOO_{label}"])
        lolo_scans, lolo_status = run_gemma_scans(
            gemma_wsl,
            _prefixes(primary_plink_dir, lolo_labels),
            primary_mapping,
            _kinship_paths(paths["intermediate"] / "kinship_primary_n100"),
            _chromosome_lists(primary_plink_dir),
            paths["intermediate"] / "gemma_lolo_score",
            logger,
            lmm_mode=3,
            run_global=False,
            run_loco=True,
            force=force,
        )
        save_scan_tables(lolo_scans, paths["results"] / "gemma_lolo_score")
        lolo_status.to_csv(paths["tables"] / "gemma_lolo_score_status.csv", index=False)

    logger.info("Official GEMMA stage complete")


def _reference_validation(summary: dict[str, Any], cfg: dict[str, Any]) -> pd.DataFrame:
    spec = cfg.get("reference_validation", {})
    if not spec or not bool(spec.get("enabled", False)):
        return pd.DataFrame(columns=["metric", "observed", "expected", "tolerance", "status", "kind"])
    rows: list[dict[str, Any]] = []
    failures: list[str] = []
    for metric, expected in spec.get("exact", {}).items():
        observed = summary.get(metric)
        ok = observed == expected
        rows.append({
            "metric": metric, "observed": observed, "expected": expected,
            "tolerance": "exact", "status": "PASS" if ok else "FAIL", "kind": "primary",
        })
        if not ok:
            failures.append(f"{metric}: observed={observed}, expected={expected}")
    for metric, definition in spec.get("numeric", {}).items():
        observed = float(summary.get(metric, float("nan")))
        expected = float(definition["value"])
        if "absolute_tolerance" in definition:
            tolerance = float(definition["absolute_tolerance"] )
            ok = np.isfinite(observed) and abs(observed - expected) <= tolerance
            tolerance_label = f"absolute <= {tolerance:g}"
        elif "log10_absolute_tolerance" in definition:
            tolerance = float(definition["log10_absolute_tolerance"] )
            ok = (
                np.isfinite(observed) and observed > 0 and expected > 0
                and abs(np.log10(observed) - np.log10(expected)) <= tolerance
            )
            tolerance_label = f"|delta log10| <= {tolerance:g}"
        else:
            raise ValueError(f"Reference validation metric {metric} lacks a tolerance")
        rows.append({
            "metric": metric, "observed": observed, "expected": expected,
            "tolerance": tolerance_label, "status": "PASS" if ok else "FAIL", "kind": "primary",
        })
        if not ok:
            failures.append(f"{metric}: observed={observed}, expected={expected}, {tolerance_label}")
    for metric, expected in spec.get("diagnostic_expected", {}).items():
        observed = summary.get(metric)
        ok = observed == expected
        rows.append({
            "metric": metric, "observed": observed, "expected": expected,
            "tolerance": "diagnostic exact", "status": "PASS" if ok else "WARN", "kind": "diagnostic",
        })
    result = pd.DataFrame(rows)
    result.attrs["primary_failures"] = failures
    result.attrs["strict_primary"] = bool(spec.get("strict_primary", False))
    return result


def _candidate_lolo_support(
    candidates: pd.DataFrame,
    primary_score: pd.DataFrame,
    lolo_scans: dict[tuple[str, str], pd.DataFrame],
    locations: list[str],
) -> pd.DataFrame:
    if candidates.empty or "variant_index" not in candidates.columns:
        return pd.DataFrame()
    candidate_ids = candidates[["variant_index", "lead_snp", "CHR", "pos"]].drop_duplicates("variant_index")
    primary = primary_score.set_index("variant_index")
    rows: list[dict[str, Any]] = []
    for location in locations:
        label = "".join(ch if ch.isalnum() else "_" for ch in str(location)).strip("_").upper()
        inc_key = ("loco", f"INC_LOO_{label}")
        sev_key = ("loco", f"SEV_LOO_{label}")
        if inc_key not in lolo_scans or sev_key not in lolo_scans:
            continue
        table = build_test_table(lolo_scans[inc_key], lolo_scans[sev_key], "score").set_index("variant_index")
        for row in candidate_ids.itertuples(index=False):
            variant_index = int(row.variant_index)
            if variant_index not in table.index or variant_index not in primary.index:
                continue
            current = table.loc[variant_index]
            base = primary.loc[variant_index]
            rows.append({
                "variant_index": variant_index,
                "lead_snp": row.lead_snp,
                "CHR": int(row.CHR),
                "pos": int(row.pos),
                "omitted_location": str(location),
                "p_INC_score": float(current["p_INC"]),
                "p_SEV_score": float(current["p_SEV"]),
                "p_ACAT_score": float(current["p_ACAT"]),
                "p_SIMES_score": float(current["p_SIMES"]),
                "beta_INC_score": float(current.get("beta_INC", np.nan)),
                "beta_SEV_score": float(current.get("beta_SEV", np.nan)),
                "primary_p_INC_score": float(base["p_INC"]),
                "primary_p_SEV_score": float(base["p_SEV"]),
                "primary_p_ACAT_score": float(base["p_ACAT"]),
                "primary_beta_INC_score": float(base.get("beta_INC", np.nan)),
                "primary_beta_SEV_score": float(base.get("beta_SEV", np.nan)),
                "INC_effect_direction_preserved": bool(
                    np.sign(current.get("beta_INC", np.nan)) == np.sign(base.get("beta_INC", np.nan))
                ),
                "SEV_effect_direction_preserved": bool(
                    np.sign(current.get("beta_SEV", np.nan)) == np.sign(base.get("beta_SEV", np.nan))
                ),
            })
    return pd.DataFrame(rows)


def _append_workbook_sheets(workbook: Path, sheets: dict[str, pd.DataFrame]) -> None:
    if not sheets or not workbook.exists():
        return
    with pd.ExcelWriter(workbook, engine="openpyxl", mode="a", if_sheet_exists="replace") as writer:
        for name, table in sheets.items():
            if table is not None and not table.empty:
                table.to_excel(writer, sheet_name=name[:31], index=False)


def _scan_comparison(primary: pd.DataFrame, comparator: pd.DataFrame, label: str) -> pd.DataFrame:
    required = {"variant_index"}
    if not required.issubset(primary.columns) or not required.issubset(comparator.columns):
        raise ValueError("Sensitivity tables must contain variant_index")
    columns = ["variant_index"]
    for method in ["INC", "SEV", "ACAT", "SIMES"]:
        columns.extend([f"p_{method}", f"q_{method}"])
    left = primary[columns].copy()
    right = comparator[columns].copy()
    merged = left.merge(right, on="variant_index", how="inner", suffixes=("_primary", "_comparator"), validate="one_to_one")
    if len(merged) != len(primary) or len(merged) != len(comparator):
        raise ValueError(
            f"Sensitivity comparison {label!r} matched {len(merged):,} variants; "
            f"primary={len(primary):,}, comparator={len(comparator):,}"
        )
    rows = []
    for method in ["INC", "SEV", "ACAT", "SIMES"]:
        x = merged[f"p_{method}_primary"].to_numpy(float)
        y = merged[f"p_{method}_comparator"].to_numpy(float)
        rows.append(
            {
                "comparison": label,
                "method": method,
                "n_matched_variants": int(len(merged)),
                "spearman_rho": float(stats.spearmanr(x, y).statistic),
                "primary_min_p": float(np.min(x)),
                "comparator_min_p": float(np.min(y)),
                "primary_fdr_hits": int((merged[f"q_{method}_primary"] < 0.10).sum()),
                "comparator_fdr_hits": int((merged[f"q_{method}_comparator"] < 0.10).sum()),
            }
        )
    return pd.DataFrame(rows)


def finalize_stage(paths: dict[str, Path], cfg: dict[str, Any], logger) -> dict[str, Any]:
    primary_scans = load_scan_tables(paths["results"] / "gemma_primary_n100")
    require(paths["cache"] / "cache_manifest.json", "genotype cache manifest")
    genotype_manifest = json.loads((paths["cache"] / "cache_manifest.json").read_text(encoding="utf-8-sig"))
    metadata = pd.read_csv(paths["cache"] / "variant_metadata.csv.gz")
    geno = np.memmap(
        paths["cache"] / "genotypes_int8.bin",
        dtype=np.int8,
        mode="r",
        shape=(int(genotype_manifest["n_variants"]), int(genotype_manifest["n_samples"])),
    )
    population = pd.read_csv(paths["intermediate"] / "population_definition.csv")
    primary_ids = population.loc[population["PRIMARY_ELIGIBLE"].astype(bool), "ID_std"].astype(str).tolist()
    all_ids = [str(x) for x in genotype_manifest["sample_ids"]]
    id_to_index = {sid: idx for idx, sid in enumerate(all_ids)}
    sample_indices = np.array([id_to_index[sid] for sid in primary_ids], dtype=int)

    summary = finalize_results(
        paths["root"],
        geno,
        metadata,
        sample_indices,
        primary_scans,
        primary_scans,
        primary_scans,
        paths["inputs"] / "Sbicolor_454_v3.1.1.gene.gff3.gz",
        logger,
        candidate_top_n=int(cfg["candidate_top_n"]),
        index_p=float(cfg["candidate_index_p"]),
        secondary_p=float(cfg["candidate_secondary_p"]),
        ld_window_kb=int(cfg["ld_clump_window_kb"]),
        ld_r2=float(cfg["ld_clump_r2"]),
        power_sample_sizes=(int(cfg["expected_primary_n"]), int(cfg["expected_matched_ids"])),
    )

    primary_score = build_test_table(primary_scans[("loco", "INC")], primary_scans[("loco", "SEV")], "score")
    sensitivity_rows = []
    pc3_dir = paths["results"] / "gemma_primary_n100_pc3"
    if pc3_dir.exists():
        pc3 = load_scan_tables(pc3_dir)
        if ("loco", "INC") in pc3 and ("loco", "SEV") in pc3:
            pc3_score = build_test_table(pc3[("loco", "INC")], pc3[("loco", "SEV")], "score")
            sensitivity_rows.append(_scan_comparison(primary_score, pc3_score, "K-only primary vs K+PC3"))
    full_dir = paths["results"] / "gemma_full_n102"
    if full_dir.exists():
        full = load_scan_tables(full_dir)
        if ("loco", "INC") in full and ("loco", "SEV") in full:
            full_score = build_test_table(full[("loco", "INC")], full[("loco", "SEV")], "score")
            sensitivity_rows.append(_scan_comparison(primary_score, full_score, "n=100 primary vs n=102"))
    if sensitivity_rows:
        sensitivity = pd.concat(sensitivity_rows, ignore_index=True)
        sensitivity.to_csv(paths["tables"] / "Model_Sensitivity.csv", index=False)

    lolo_dir = paths["results"] / "gemma_lolo_score"
    workbook_additions: dict[str, pd.DataFrame] = {}
    if lolo_dir.exists():
        lolo = load_scan_tables(lolo_dir)
        lolo_rows = []
        for location in cfg["locations"]:
            label = "".join(ch if ch.isalnum() else "_" for ch in str(location)).strip("_").upper()
            inc_key = ("loco", f"INC_LOO_{label}")
            sev_key = ("loco", f"SEV_LOO_{label}")
            if inc_key not in lolo or sev_key not in lolo:
                continue
            table = build_test_table(lolo[inc_key], lolo[sev_key], "score")
            comp = _scan_comparison(primary_score, table, f"leave out {location}")
            lolo_rows.append(comp)
        if lolo_rows:
            lolo_genomewide = pd.concat(lolo_rows, ignore_index=True)
            lolo_genomewide.to_csv(paths["tables"] / "Leave_One_Location_Out.csv", index=False)
            workbook_additions["LOO_Genomewide"] = lolo_genomewide
        candidate_path = paths["tables"] / "Top_Ranked_Candidates.csv"
        if candidate_path.exists():
            candidates = pd.read_csv(candidate_path)
            lolo_candidates = _candidate_lolo_support(candidates, primary_score, lolo, list(cfg["locations"]))
            if not lolo_candidates.empty:
                lolo_candidates.to_csv(
                    paths["tables"] / "Leave_One_Location_Out_Candidates.csv", index=False
                )
                workbook_additions["LOO_Candidates"] = lolo_candidates

    final_manifest = paths["results"] / "analysis_manifest.json"
    if final_manifest.exists():
        payload = json.loads(final_manifest.read_text(encoding="utf-8"))
        payload["pipeline_version"] = VERSION
        payload["sensitivity_tables_written"] = bool(sensitivity_rows)
        payload["lolo_table_written"] = (paths["tables"] / "Leave_One_Location_Out.csv").exists()
        payload["lolo_candidate_table_written"] = (
            paths["tables"] / "Leave_One_Location_Out_Candidates.csv"
        ).exists()
        validation = _reference_validation(payload, cfg)
        validation.to_csv(paths["tables"] / "Reference_Reproduction_Validation.csv", index=False)
        workbook_additions["Reference_Validation"] = validation
        primary_failures = list(validation.attrs.get("primary_failures", []))
        strict_reference = bool(validation.attrs.get("strict_primary", False))
        payload["reference_validation_passed"] = bool(
            validation.empty or (validation.loc[validation["kind"] == "primary", "status"] == "PASS").all()
        )
        payload["reference_validation_failures"] = primary_failures
        save_json(final_manifest, payload)
        summary = payload
        _append_workbook_sheets(
            paths["results"] / "LeafBlight_MultiEnvironment_Analysis_Results.xlsx",
            workbook_additions,
        )
        if primary_failures and strict_reference:
            raise RuntimeError(
                "Reference values did not reproduce within the configured tolerances. Validation output: "
                + "; ".join(primary_failures)
            )
    del geno
    logger.info("Finalization stage complete")
    return summary


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Sorghum leaf blight multi-environment genomic analysis pipeline"
    )
    parser.add_argument("--version", action="version", version=VERSION)
    parser.add_argument("--project-root", required=True)
    parser.add_argument(
        "--stage", choices=["prepare", "gemma", "finalize", "all"], default="all"
    )
    parser.add_argument("--gemma-wsl", default="")
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--config", default="")
    args = parser.parse_args()

    paths = ensure_directories(args.project_root)
    config_path = Path(args.config) if args.config else paths["config"] / "final_config.json"
    if not config_path.exists():
        bundled = Path(__file__).resolve().parents[1] / "config" / "final_config.json"
        if not bundled.exists():
            raise FileNotFoundError(f"Configuration not found: {config_path}")
        config_path.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(bundled, config_path)
    cfg = load_config(config_path)
    stamp = time.strftime("%Y%m%d_%H%M%S")
    logger = setup_logger(paths["logs"] / f"analysis_pipeline_{args.stage}_{stamp}.log")
    logger.info("Sorghum leaf blight pipeline v%s", VERSION)
    logger.info("Project root: %s", paths["root"])

    if args.stage in {"prepare", "all"}:
        prepare_stage(paths, cfg, logger, args.force)
    if args.stage in {"gemma", "all"}:
        gemma_stage(paths, cfg, args.gemma_wsl, logger, args.force)
    if args.stage in {"finalize", "all"}:
        summary = finalize_stage(paths, cfg, logger)
        print("FINAL_SUMMARY=" + str(paths["results"] / "ANALYSIS_SUMMARY.txt"))
        print("FINAL_WORKBOOK=" + str(paths["results"] / "LeafBlight_MultiEnvironment_Analysis_Results.xlsx"))
        print("FINAL_MANIFEST=" + str(paths["results"] / "analysis_manifest.json"))
        logger.info("Final summary: %s", summary)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception:
        import traceback

        traceback.print_exc()
        raise
