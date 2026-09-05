from __future__ import annotations

import hashlib
import logging
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from lb_common import bh_fdr, elapsed, load_json, save_json


@dataclass
class EquivalenceResult:
    orientation: np.ndarray
    global_inverse: np.ndarray
    global_counts: np.ndarray
    global_digests: np.ndarray
    chrom_inverse: np.ndarray
    chrom_counts: np.ndarray
    chrom_chromosomes: np.ndarray
    chrom_digests: np.ndarray


def _canonical_digest(row: np.ndarray) -> tuple[bytes, int]:
    """Return an allele-orientation-invariant digest and orientation flag.

    ``row`` is an int8 dosage vector in {-1, 0, 1, 2}. The alternative
    orientation is 2-dosage for called genotypes; missing values remain -1.
    Exact and allele-complemented genotype patterns therefore share a class.
    """
    x = np.asarray(row, dtype=np.int8)
    comp = x.copy()
    called = comp >= 0
    comp[called] = 2 - comp[called]
    raw = x.tobytes(order="C")
    alt = comp.tobytes(order="C")
    if alt < raw:
        return hashlib.blake2b(alt, digest_size=16).digest(), -1
    return hashlib.blake2b(raw, digest_size=16).digest(), 1


def _chrom_keys(chromosomes: np.ndarray, digests: np.ndarray) -> np.ndarray:
    dtype = np.dtype([("chrom", np.int16), ("digest", "S16")])
    out = np.empty(len(chromosomes), dtype=dtype)
    out["chrom"] = np.asarray(chromosomes, dtype=np.int16)
    out["digest"] = digests
    return out


def compute_or_load_equivalence(
    geno: np.memmap,
    metadata: pd.DataFrame,
    sample_indices: np.ndarray,
    cache_dir: str | Path,
    logger: logging.Logger,
    force: bool = False,
) -> EquivalenceResult:
    cache = Path(cache_dir)
    cache.mkdir(parents=True, exist_ok=True)
    npz_path = cache / "primary_n100_genotype_equivalence.npz"
    json_path = cache / "primary_n100_genotype_equivalence.json"
    expected = {
        "n_variants": int(len(metadata)),
        "sample_indices": [int(x) for x in np.asarray(sample_indices, dtype=int)],
        "algorithm": "BLAKE2b-128; exact-or-allele-complement dosage vector",
    }

    if not force and npz_path.exists() and json_path.exists():
        try:
            old = load_json(json_path)
            if all(old.get(k) == v for k, v in expected.items()):
                z = np.load(npz_path, allow_pickle=False)
                logger.info(
                    "Reusing equivalence cache: %s global alias classes; %s chromosome-specific test classes",
                    f"{len(z['global_counts']):,}", f"{len(z['chrom_counts']):,}",
                )
                return EquivalenceResult(
                    orientation=z["orientation"],
                    global_inverse=z["global_inverse"],
                    global_counts=z["global_counts"],
                    global_digests=z["global_digests"],
                    chrom_inverse=z["chrom_inverse"],
                    chrom_counts=z["chrom_counts"],
                    chrom_chromosomes=z["chrom_chromosomes"],
                    chrom_digests=z["chrom_digests"],
                )
        except Exception as exc:
            logger.warning("Equivalence cache could not be reused: %s", exc)

    start = time.time()
    n = len(metadata)
    digests = np.empty(n, dtype="S16")
    orientation = np.empty(n, dtype=np.int8)
    logger.info(
        "Computing exact/complement genotype signatures for %s variants over %d primary samples",
        f"{n:,}", len(sample_indices),
    )
    for i in range(n):
        digest, orient = _canonical_digest(np.asarray(geno[i, sample_indices], dtype=np.int8))
        digests[i] = digest
        orientation[i] = orient
        if i and i % 50000 == 0:
            logger.info("  signatures %s/%s (%s)", f"{i:,}", f"{n:,}", elapsed(start))

    global_digests, global_inverse, global_counts = np.unique(
        digests, return_inverse=True, return_counts=True
    )
    keys = _chrom_keys(metadata["CHR"].to_numpy(dtype=np.int16), digests)
    chrom_unique, chrom_inverse, chrom_counts = np.unique(
        keys, return_inverse=True, return_counts=True
    )
    chrom_chromosomes = np.asarray(chrom_unique["chrom"], dtype=np.int16)
    chrom_digests = np.asarray(chrom_unique["digest"], dtype="S16")

    np.savez_compressed(
        npz_path,
        orientation=orientation,
        global_inverse=global_inverse.astype(np.int32),
        global_counts=global_counts.astype(np.int32),
        global_digests=global_digests,
        chrom_inverse=chrom_inverse.astype(np.int32),
        chrom_counts=chrom_counts.astype(np.int32),
        chrom_chromosomes=chrom_chromosomes,
        chrom_digests=chrom_digests,
    )
    save_json(
        json_path,
        {
            **expected,
            "global_alias_class_count": int(len(global_counts)),
            "chromosome_specific_test_class_count": int(len(chrom_counts)),
            "variants_in_global_duplicate_classes": int(np.sum(global_counts[global_inverse] > 1)),
            "variants_in_chromosome_duplicate_classes": int(np.sum(chrom_counts[chrom_inverse] > 1)),
            "largest_global_class": int(global_counts.max()),
            "largest_chromosome_class": int(chrom_counts.max()),
            "elapsed": elapsed(start),
        },
    )
    logger.info(
        "Equivalence analysis complete: %s global alias classes; %s chromosome-specific test classes (%s)",
        f"{len(global_counts):,}", f"{len(chrom_counts):,}", elapsed(start),
    )
    return EquivalenceResult(
        orientation=orientation,
        global_inverse=global_inverse.astype(np.int32),
        global_counts=global_counts.astype(np.int32),
        global_digests=global_digests,
        chrom_inverse=chrom_inverse.astype(np.int32),
        chrom_counts=chrom_counts.astype(np.int32),
        chrom_chromosomes=chrom_chromosomes,
        chrom_digests=chrom_digests,
    )


def _hex(values: np.ndarray) -> list[str]:
    return [bytes(x).hex().upper() for x in values]


def _representatives(values: np.ndarray, inverse: np.ndarray, n_classes: int) -> np.ndarray:
    frame = pd.DataFrame(
        {
            "class_index": inverse,
            "p": np.asarray(values, dtype=float),
            "variant_index": np.arange(len(values), dtype=np.int64),
        }
    )
    frame["p_sort"] = frame["p"].where(np.isfinite(frame["p"]), np.inf)
    reps = (
        frame.sort_values(["p_sort", "variant_index"], kind="mergesort")
        .drop_duplicates("class_index", keep="first")
        .set_index("class_index")["variant_index"]
        .reindex(np.arange(n_classes))
        .to_numpy(dtype=np.int64)
    )
    return reps


def build_equivalence_tables(
    metadata: pd.DataFrame,
    primary_meta: pd.DataFrame,
    eq: EquivalenceResult,
) -> dict[str, pd.DataFrame]:
    metadata = metadata.reset_index(drop=True)
    scan = primary_meta.reset_index(drop=True)
    if not (len(metadata) == len(scan) == len(eq.global_inverse)):
        raise ValueError("Metadata, primary scan, and equivalence arrays have different lengths")
    methods = ["INC", "SEV", "ACAT", "SIMES"]
    for method in methods:
        if f"p_{method}" not in scan.columns:
            raise ValueError(f"Primary scan lacks p_{method}")

    n_chrom_classes = len(eq.chrom_counts)
    reps = {
        method: _representatives(scan[f"p_{method}"].to_numpy(), eq.chrom_inverse, n_chrom_classes)
        for method in methods
    }
    acat_reps = reps["ACAT"]
    rep_meta = metadata.iloc[acat_reps].reset_index(drop=True)
    chrom_table = pd.DataFrame(
        {
            "test_class_index": np.arange(n_chrom_classes, dtype=np.int32),
            "test_class_id": [
                f"CHR{int(ch):02d}_{digest}" for ch, digest in zip(eq.chrom_chromosomes, _hex(eq.chrom_digests))
            ],
            "CHR": eq.chrom_chromosomes.astype(int),
            "class_size": eq.chrom_counts.astype(int),
            "representative_variant_index_ACAT": acat_reps,
            "representative_snp_ACAT": rep_meta["snp"].astype(str).to_numpy(),
            "representative_pos_ACAT": rep_meta["pos"].astype(int).to_numpy(),
            "global_alias_class_index": eq.global_inverse[acat_reps].astype(int),
        }
    )
    chrom_table["global_alias_class_size"] = eq.global_counts[
        chrom_table["global_alias_class_index"].to_numpy(dtype=int)
    ]
    for method in methods:
        r = reps[method]
        p_min = scan[f"p_{method}"].to_numpy(dtype=float)[r]
        chrom_table[f"representative_variant_index_{method}"] = r
        chrom_table[f"minimum_p_{method}"] = p_min
        chrom_table[f"q_unique_{method}"] = bh_fdr(p_min)

    # Same genotype vector on the same chromosome is the same LOCO test and
    # should produce the same statistic. Record deviations as validation metrics.
    sort_idx = np.argsort(eq.chrom_inverse, kind="stable")
    sorted_cls = eq.chrom_inverse[sort_idx]
    starts = np.r_[0, np.flatnonzero(sorted_cls[1:] != sorted_cls[:-1]) + 1]
    for method in methods:
        vals = scan[f"p_{method}"].to_numpy(dtype=float)[sort_idx]
        safe_min = np.minimum.reduceat(np.where(np.isfinite(vals), vals, np.inf), starts)
        safe_max = np.maximum.reduceat(np.where(np.isfinite(vals), vals, -np.inf), starts)
        spread = safe_max - safe_min
        chrom_table[f"within_class_p_spread_{method}"] = spread

    n_global = len(eq.global_counts)
    global_rep_idx = _representatives(scan["p_ACAT"].to_numpy(), eq.global_inverse, n_global)
    global_rep = metadata.iloc[global_rep_idx].reset_index(drop=True)
    global_table = pd.DataFrame(
        {
            "global_alias_class_index": np.arange(n_global, dtype=np.int32),
            "global_alias_class_id": [f"G_{x}" for x in _hex(eq.global_digests)],
            "class_size": eq.global_counts.astype(int),
            "representative_variant_index": global_rep_idx,
            "representative_snp": global_rep["snp"].astype(str).to_numpy(),
            "representative_CHR": global_rep["CHR"].astype(int).to_numpy(),
            "representative_pos": global_rep["pos"].astype(int).to_numpy(),
        }
    )
    meta_map = metadata[["variant_index", "snp", "CHR", "pos"]].copy()
    meta_map["global_alias_class_index"] = eq.global_inverse
    grouped = meta_map.groupby("global_alias_class_index", sort=False)
    global_table = global_table.merge(
        grouped.agg(
            n_chromosomes=("CHR", "nunique"),
            minimum_position=("pos", "min"),
            maximum_position=("pos", "max"),
        ).reset_index(),
        on="global_alias_class_index", how="left",
    )
    # Build long mapping strings only for duplicate classes. Singleton classes
    # use their representative mapping; this avoids a very expensive groupby
    # apply over hundreds of thousands of one-row groups.
    global_table["alias_mappings"] = (
        global_table["representative_snp"].astype(str)
        + "@chr" + global_table["representative_CHR"].astype(str)
        + ":" + global_table["representative_pos"].astype(str)
    )
    duplicate_rows = meta_map.loc[eq.global_counts[eq.global_inverse] > 1]
    if not duplicate_rows.empty:
        duplicate_strings = duplicate_rows.groupby("global_alias_class_index", sort=False).apply(
            lambda x: "; ".join(
                f"{s}@chr{int(c)}:{int(pos)}"
                for s, c, pos in x[["snp", "CHR", "pos"]].itertuples(index=False, name=None)
            )
        )
        mapped = global_table["global_alias_class_index"].map(duplicate_strings)
        global_table.loc[mapped.notna(), "alias_mappings"] = mapped[mapped.notna()]

    variant_table = scan.copy()
    variant_table["global_alias_class_index"] = eq.global_inverse.astype(np.int32)
    variant_table["global_alias_class_size"] = eq.global_counts[eq.global_inverse].astype(np.int32)
    variant_table["test_class_index"] = eq.chrom_inverse.astype(np.int32)
    variant_table["test_class_size"] = eq.chrom_counts[eq.chrom_inverse].astype(np.int32)
    q_map = chrom_table.set_index("test_class_index")
    for method in methods:
        variant_table[f"q_unique_{method}"] = variant_table["test_class_index"].map(
            q_map[f"q_unique_{method}"]
        ).to_numpy(dtype=float)

    duplicates = metadata.loc[eq.global_counts[eq.global_inverse] > 1, ["variant_index", "snp", "CHR", "pos", "MAF"]].copy()
    duplicates["orientation_to_canonical"] = eq.orientation[eq.global_counts[eq.global_inverse] > 1]
    duplicates["global_alias_class_index"] = eq.global_inverse[eq.global_counts[eq.global_inverse] > 1]
    duplicates["global_alias_class_size"] = eq.global_counts[eq.global_inverse[eq.global_counts[eq.global_inverse] > 1]]
    duplicates["test_class_index"] = eq.chrom_inverse[eq.global_counts[eq.global_inverse] > 1]
    duplicates["test_class_size"] = eq.chrom_counts[eq.chrom_inverse[eq.global_counts[eq.global_inverse] > 1]]
    duplicates = duplicates.merge(
        global_table[["global_alias_class_index", "global_alias_class_id", "n_chromosomes"]],
        on="global_alias_class_index", how="left",
    ).merge(
        chrom_table[["test_class_index", "test_class_id"]],
        on="test_class_index", how="left",
    )

    summary = pd.DataFrame(
        {
            "metric": [
                "variant_rows",
                "global_exact_or_complement_alias_classes",
                "chromosome_specific_unique_test_classes",
                "variants_in_global_duplicate_classes",
                "variants_in_chromosome_duplicate_classes",
                "global_alias_classes_spanning_multiple_chromosomes",
                "largest_global_alias_class_size",
                "largest_chromosome_test_class_size",
                "max_within_test_class_ACAT_p_spread",
            ],
            "value": [
                len(metadata),
                n_global,
                n_chrom_classes,
                int(np.sum(eq.global_counts[eq.global_inverse] > 1)),
                int(np.sum(eq.chrom_counts[eq.chrom_inverse] > 1)),
                int((global_table["n_chromosomes"] > 1).sum()),
                int(eq.global_counts.max()),
                int(eq.chrom_counts.max()),
                float(chrom_table["within_class_p_spread_ACAT"].max()),
            ],
        }
    )
    return {
        "summary": summary,
        "test_classes": chrom_table,
        "global_alias_classes": global_table,
        "duplicate_mapping": duplicates,
        "variant_table": variant_table,
    }


def collapse_for_method(variant_table: pd.DataFrame, method: str) -> pd.DataFrame:
    """Return one representative per chromosome-specific test, then per global alias.

    The first collapse defines unique LOCO tests. The second prevents identical or
    allele-complemented patterns on distant chromosomes from being counted as
    independent biological candidate regions.
    """
    p_col = f"p_{method}"
    x = variant_table.copy()
    x["_p_sort"] = pd.to_numeric(x[p_col], errors="coerce").fillna(np.inf)
    x = x.sort_values(["_p_sort", "variant_index"], kind="mergesort")
    x = x.drop_duplicates("test_class_index", keep="first")
    x = x.sort_values(["_p_sort", "variant_index"], kind="mergesort")
    x = x.drop_duplicates("global_alias_class_index", keep="first")
    return x.drop(columns="_p_sort").reset_index(drop=True)


def _r2_to_lead(dosages: np.ndarray, lead_row: int) -> np.ndarray:
    x = np.asarray(dosages, dtype=float)
    x = x - x.mean(axis=1, keepdims=True)
    lead = x[lead_row]
    lead_ss = float(lead @ lead)
    ss = np.sum(x * x, axis=1)
    numer = x @ lead
    denom = ss * lead_ss
    return np.clip(
        np.divide(numer * numer, denom, out=np.zeros_like(numer), where=denom > 1e-12),
        0.0, 1.0,
    )


def alias_aware_ld_clump(
    unique_table: pd.DataFrame,
    geno: np.memmap,
    sample_indices: np.ndarray,
    method: str,
    index_p: float,
    secondary_p: float,
    window_kb: int,
    r2_threshold: float,
    fallback_top_n: int,
    logger: logging.Logger,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    p_col = f"p_{method}"
    index = unique_table.loc[pd.to_numeric(unique_table[p_col], errors="coerce") <= index_p].copy()
    criterion = f"{p_col}<={index_p:g}"
    if index.empty:
        index = unique_table.nsmallest(fallback_top_n, p_col).copy()
        criterion = f"top_{fallback_top_n}_global_alias_classes_by_{p_col}"
    secondary = unique_table.loc[pd.to_numeric(unique_table[p_col], errors="coerce") <= secondary_p].copy()
    if secondary.empty:
        secondary = index.copy()
    index = index.sort_values([p_col, "CHR", "pos"], kind="mergesort")
    secondary = secondary.sort_values([p_col, "CHR", "pos"], kind="mergesort")

    assigned: set[int] = set()
    lead_rows: list[dict[str, Any]] = []
    member_rows: list[dict[str, Any]] = []
    window = int(window_kb * 1000)
    rank = 0
    for _, candidate in index.iterrows():
        v = int(candidate["variant_index"])
        if v in assigned:
            continue
        rank += 1
        chrom = int(candidate["CHR"])
        pos = int(candidate["pos"])
        local = secondary.loc[
            (secondary["CHR"].astype(int) == chrom)
            & secondary["pos"].astype(int).between(pos - window, pos + window)
        ].copy()
        local = local.loc[~local["variant_index"].astype(int).isin(assigned)]
        if v not in set(local["variant_index"].astype(int)):
            local = pd.concat([candidate.to_frame().T, local], ignore_index=True)
        local = local.drop_duplicates("variant_index").reset_index(drop=True)
        idx = local["variant_index"].to_numpy(dtype=int)
        dosage = np.asarray(geno[idx][:, sample_indices], dtype=float)
        lead_row = int(np.flatnonzero(idx == v)[0])
        local["r2_to_lead"] = _r2_to_lead(dosage, lead_row)
        members = local.loc[local["r2_to_lead"] >= r2_threshold].copy()
        for _, member in members.iterrows():
            assigned.add(int(member["variant_index"]))
            row = member.to_dict()
            row.update(
                {
                    "lead_rank": rank,
                    "lead_variant_index": v,
                    "lead_snp": str(candidate["snp"]),
                    "clump_method": method,
                    "candidate_criterion": criterion,
                    "window_kb": int(window_kb),
                    "r2_threshold": float(r2_threshold),
                }
            )
            member_rows.append(row)
        lead = candidate.to_dict()
        lead.update(
            {
                "lead_rank": rank,
                "lead_variant_index": v,
                "lead_snp": str(candidate["snp"]),
                "clump_method": method,
                "candidate_criterion": criterion,
                "window_kb": int(window_kb),
                "r2_threshold": float(r2_threshold),
                "n_local_members": int(len(members)),
            }
        )
        lead_rows.append(lead)
    leads = pd.DataFrame(lead_rows)
    members = pd.DataFrame(member_rows)
    logger.info(
        "Alias-aware %s clumping: %d index classes -> %d candidate signals (%s; r2>=%.2f)",
        method, len(index), len(leads), criterion, r2_threshold,
    )
    return leads, members


def enrich_candidates(
    candidates: pd.DataFrame,
    scan_tables: dict[str, pd.DataFrame],
    global_classes: pd.DataFrame,
) -> pd.DataFrame:
    if candidates.empty:
        return candidates.copy()
    primary = scan_tables["primary"]
    cols = [
        "variant_index", "snp", "CHR", "pos", "MAF", "missing_rate",
        "beta_INC", "se_INC", "p_INC", "beta_SEV", "se_SEV", "p_SEV",
        "p_ACAT", "q_ACAT", "p_SIMES", "q_SIMES",
        "q_unique_INC", "q_unique_SEV", "q_unique_ACAT", "q_unique_SIMES",
        "global_alias_class_index", "global_alias_class_size", "test_class_index", "test_class_size",
    ]
    cols = [c for c in cols if c in primary.columns]
    # Clumping input carries a copy of most scan columns. Preserve only
    # clump/selection metadata before joining the primary scan,
    # otherwise pandas creates ambiguous _x/_y effect columns.
    overlap = [c for c in candidates.columns if c in primary.columns and c != "variant_index"]
    out = candidates.drop(columns=overlap, errors="ignore").merge(
        primary[cols], on="variant_index", how="left"
    )
    for trait in ["INC", "SEV"]:
        out[f"ci_low_{trait}"] = out[f"beta_{trait}"] - 1.96 * out[f"se_{trait}"]
        out[f"ci_high_{trait}"] = out[f"beta_{trait}"] + 1.96 * out[f"se_{trait}"]

    extras = {
        "conditional": ["p_SEV", "p_ACAT", "p_SIMES"],
        "positive": ["p_SEV", "p_ACAT", "p_SIMES"],
        "pc3": ["p_INC", "p_SEV", "p_ACAT", "p_SIMES"],
        "full102": ["p_INC", "p_SEV", "p_ACAT", "p_SIMES"],
        "raw": ["p_INC", "p_SEV", "p_ACAT", "p_SIMES"],
    }
    for label, wanted in extras.items():
        table = scan_tables.get(label)
        if table is None:
            continue
        keep = ["variant_index"] + [c for c in wanted if c in table.columns]
        subset = table[keep].rename(columns={c: f"{c}_{label}" for c in keep if c != "variant_index"})
        out = out.merge(subset, on="variant_index", how="left")

    g = global_classes[[
        "global_alias_class_index", "global_alias_class_id", "class_size", "n_chromosomes", "alias_mappings"
    ]].rename(columns={"class_size": "global_alias_count"})
    out = out.drop(columns="global_alias_class_size", errors="ignore").merge(
        g, on="global_alias_class_index", how="left"
    )
    out = out.sort_values(["p_ACAT", "p_SEV", "CHR", "pos"], kind="mergesort").reset_index(drop=True)
    out.insert(0, "candidate_rank", np.arange(1, len(out) + 1))
    return out
