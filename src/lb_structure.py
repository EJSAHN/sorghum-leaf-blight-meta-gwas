from __future__ import annotations

import logging
import math
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from lb_common import bytes_hash, dataframe_hash, elapsed, load_json, save_json, slugify
from lb_genotype import load_variant_dosages


@dataclass
class PopulationStructure:
    label: str
    sample_ids: list[str]
    sample_indices: np.ndarray
    pcs: pd.DataFrame
    global_kinship: np.ndarray
    global_eigenvalues: np.ndarray
    global_eigenvectors: np.ndarray
    loco_eigenvalues: np.ndarray
    loco_eigenvectors: np.ndarray
    marker_table: pd.DataFrame
    structure_hash: str
    spectral_genotypes: np.memmap


def _symmetrize(matrix: np.ndarray) -> np.ndarray:
    return (np.asarray(matrix, dtype=float) + np.asarray(matrix, dtype=float).T) / 2.0


def _eigh_desc(matrix: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    values, vectors = np.linalg.eigh(_symmetrize(matrix))
    order = np.argsort(values)[::-1]
    values = np.maximum(values[order], 0.0)
    vectors = vectors[:, order]
    return values, vectors


def _candidate_marker_table(metadata: pd.DataFrame, physical_thin_kb: int, max_candidates: int) -> pd.DataFrame:
    meta = metadata[["variant_index", "snp", "CHR", "pos", "MAF", "missing_rate"]].copy()
    meta = meta.sort_values(["CHR", "pos"], kind="mergesort")
    meta["thin_bin"] = (meta["pos"].astype(np.int64) // int(physical_thin_kb * 1000)).astype(np.int64)
    meta["maf_priority"] = -pd.to_numeric(meta["MAF"], errors="coerce").fillna(0.0)
    meta = meta.sort_values(
        ["CHR", "thin_bin", "missing_rate", "maf_priority", "pos"],
        kind="mergesort",
    )
    candidates = meta.groupby(["CHR", "thin_bin"], sort=True, as_index=False).head(1).copy()
    candidates = candidates.sort_values(["CHR", "pos"], kind="mergesort").reset_index(drop=True)
    if len(candidates) > max_candidates:
        # Deterministic genome-wide thinning if an unusually dense input exceeds the cap.
        keep = np.linspace(0, len(candidates) - 1, max_candidates, dtype=int)
        candidates = candidates.iloc[np.unique(keep)].reset_index(drop=True)
    return candidates.drop(columns=["maf_priority"])


def select_ld_pruned_structure_markers(
    geno: np.memmap,
    metadata: pd.DataFrame,
    sample_indices: np.ndarray,
    cache_dir: str | Path,
    physical_thin_kb: int,
    ld_window_kb: int,
    ld_r2: float,
    max_candidates: int,
    min_markers: int,
    logger: logging.Logger,
    force: bool = False,
) -> pd.DataFrame:
    cache = Path(cache_dir)
    cache.mkdir(parents=True, exist_ok=True)
    marker_path = cache / "common_ld_pruned_structure_markers.csv"
    manifest_path = cache / "common_ld_pruned_structure_markers.json"
    metadata_signature = dataframe_hash(metadata[["variant_index", "snp", "CHR", "pos", "MAF", "missing_rate"]])
    config = {
        "metadata_signature": metadata_signature,
        "sample_indices": [int(x) for x in sample_indices],
        "physical_thin_kb": int(physical_thin_kb),
        "ld_window_kb": int(ld_window_kb),
        "ld_r2": float(ld_r2),
        "max_candidates": int(max_candidates),
    }
    if not force and marker_path.exists() and manifest_path.exists():
        old = load_json(manifest_path)
        if all(old.get(k) == v for k, v in config.items()):
            markers = pd.read_csv(marker_path)
            if len(markers) >= min_markers:
                logger.info("Reusing %d LD-pruned structure markers", len(markers))
                return markers

    start_time = time.time()
    candidates = _candidate_marker_table(metadata, physical_thin_kb, max_candidates)
    logger.info(
        "Structure marker selection: %d physically thinned candidates (one per %d kb bin)",
        len(candidates), physical_thin_kb,
    )
    dosage = load_variant_dosages(geno, candidates["variant_index"].to_numpy(dtype=int))[:, sample_indices]
    centered = dosage - dosage.mean(axis=1, keepdims=True)
    ss = np.sum(centered * centered, axis=1)

    selected_local: list[int] = []
    selected_by_chr: dict[int, list[int]] = {chrom: [] for chrom in range(1, 11)}
    window_bp = int(ld_window_kb * 1000)
    chrom_values = candidates["CHR"].to_numpy(dtype=int)
    pos_values = candidates["pos"].to_numpy(dtype=int)

    for i in range(len(candidates)):
        if ss[i] <= 1e-12:
            continue
        chrom = int(chrom_values[i])
        nearby = [j for j in selected_by_chr[chrom] if pos_values[i] - pos_values[j] <= window_bp]
        if nearby:
            numer = centered[nearby] @ centered[i]
            denom = ss[nearby] * ss[i]
            r2_values = np.divide(numer * numer, denom, out=np.zeros_like(numer), where=denom > 1e-12)
            if np.any(r2_values >= ld_r2):
                continue
        selected_local.append(i)
        selected_by_chr[chrom].append(i)

    markers = candidates.iloc[selected_local].copy().reset_index(drop=True)
    if len(markers) < min_markers:
        raise ValueError(
            f"Only {len(markers)} LD-pruned structure markers remained; expected at least {min_markers}."
        )
    markers.insert(0, "structure_marker_rank", np.arange(1, len(markers) + 1))
    markers.to_csv(marker_path, index=False)
    marker_hash = bytes_hash(markers["variant_index"].to_numpy(dtype=np.int64).tobytes())
    save_json(
        manifest_path,
        {
            **config,
            "candidate_count": int(len(candidates)),
            "selected_count": int(len(markers)),
            "marker_hash": marker_hash,
            "elapsed": elapsed(start_time),
        },
    )
    logger.info("Selected %d LD-pruned structure markers in %s", len(markers), elapsed(start_time))
    return markers


def _standardized_marker_matrix(
    geno: np.memmap,
    marker_indices: np.ndarray,
    sample_indices: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    dosage = load_variant_dosages(geno, marker_indices)[:, sample_indices]
    mean = dosage.mean(axis=1)
    p = np.clip(mean / 2.0, 1e-6, 1 - 1e-6)
    scale = np.sqrt(2.0 * p * (1.0 - p))
    valid = np.isfinite(scale) & (scale > 1e-8) & (np.var(dosage, axis=1) > 1e-12)
    z = (dosage[valid] - mean[valid, None]) / scale[valid, None]
    return z.astype(np.float64), valid


def _normalize_kinship(matrix: np.ndarray) -> tuple[np.ndarray, float]:
    matrix = _symmetrize(matrix)
    scale = float(np.mean(np.diag(matrix)))
    if not np.isfinite(scale) or scale <= 1e-12:
        raise ValueError("Kinship matrix has a non-positive diagonal scale")
    return matrix / scale, scale


def _structure_cache_valid(manifest_path: Path, expected: dict[str, Any]) -> bool:
    if not manifest_path.exists():
        return False
    try:
        old = load_json(manifest_path)
    except Exception:
        return False
    return all(old.get(key) == value for key, value in expected.items())


def build_or_load_population_structure(
    label: str,
    geno: np.memmap,
    metadata: pd.DataFrame,
    all_sample_ids: list[str],
    population_ids: list[str],
    common_markers: pd.DataFrame,
    cache_dir: str | Path,
    n_pc_save: int,
    block_size: int,
    logger: logging.Logger,
    force: bool = False,
) -> PopulationStructure:
    cache = Path(cache_dir) / slugify(label)
    cache.mkdir(parents=True, exist_ok=True)
    id_to_index = {sid: i for i, sid in enumerate(all_sample_ids)}
    missing = [sid for sid in population_ids if sid not in id_to_index]
    if missing:
        raise ValueError(f"Population {label} contains IDs absent from genotype cache: {missing}")
    sample_indices = np.array([id_to_index[sid] for sid in population_ids], dtype=int)
    marker_indices = common_markers["variant_index"].to_numpy(dtype=int)
    marker_hash = bytes_hash(marker_indices.astype(np.int64).tobytes())
    structure_hash = bytes_hash(
        ("|".join(population_ids) + "|" + marker_hash + "|LOCO_KINSHIP_V02").encode("utf-8")
    )
    expected = {
        "label": label,
        "sample_ids": population_ids,
        "marker_hash": marker_hash,
        "n_variants": int(len(metadata)),
        "structure_hash": structure_hash,
    }
    manifest_path = cache / "structure_manifest.json"
    required = [
        cache / "global_kinship.npy",
        cache / "global_eigenvalues.npy",
        cache / "global_eigenvectors.npy",
        cache / "loco_eigenvalues.npy",
        cache / "loco_eigenvectors.npy",
        cache / "pcs.csv",
        cache / "markers_used.csv",
    ]

    if not force and _structure_cache_valid(manifest_path, expected) and all(p.exists() for p in required):
        logger.info("Reusing genomic structure for %s", label)
        global_k = np.load(cache / "global_kinship.npy")
        global_values = np.load(cache / "global_eigenvalues.npy")
        global_vectors = np.load(cache / "global_eigenvectors.npy")
        loco_values = np.load(cache / "loco_eigenvalues.npy")
        loco_vectors = np.load(cache / "loco_eigenvectors.npy")
        pcs = pd.read_csv(cache / "pcs.csv").set_index("ID_std")
        markers_used = pd.read_csv(cache / "markers_used.csv")
    else:
        start_time = time.time()
        z, valid = _standardized_marker_matrix(geno, marker_indices, sample_indices)
        markers_used = common_markers.loc[valid].copy().reset_index(drop=True)
        markers_used.to_csv(cache / "markers_used.csv", index=False)
        chroms = markers_used["CHR"].to_numpy(dtype=int)
        n_markers = len(markers_used)
        if n_markers < 500:
            raise ValueError(f"Only {n_markers} polymorphic structure markers in population {label}")

        sum_total = z.T @ z
        sum_by_chr = np.zeros((10, len(population_ids), len(population_ids)), dtype=np.float64)
        count_by_chr = np.zeros(10, dtype=int)
        for chrom in range(1, 11):
            mask = chroms == chrom
            count_by_chr[chrom - 1] = int(mask.sum())
            if mask.any():
                sum_by_chr[chrom - 1] = z[mask].T @ z[mask]

        global_k, global_scale = _normalize_kinship(sum_total / n_markers)
        global_values, global_vectors = _eigh_desc(global_k)
        pcs_array = global_vectors[:, : min(n_pc_save, len(population_ids) - 1)] * np.sqrt(
            np.maximum(global_values[: min(n_pc_save, len(population_ids) - 1)], 0.0)
        )[None, :]
        pcs = pd.DataFrame(
            pcs_array,
            index=population_ids,
            columns=[f"PC{i+1}" for i in range(pcs_array.shape[1])],
        )
        pcs.index.name = "ID_std"

        loco_values = np.zeros((10, len(population_ids)), dtype=np.float64)
        loco_vectors = np.zeros((10, len(population_ids), len(population_ids)), dtype=np.float64)
        loco_scales: list[float] = []
        loco_marker_counts: list[int] = []
        for chrom in range(1, 11):
            m_loco = n_markers - int(count_by_chr[chrom - 1])
            if m_loco <= 0:
                raise ValueError(f"No markers remain for LOCO kinship excluding chromosome {chrom}")
            k_loco, scale = _normalize_kinship((sum_total - sum_by_chr[chrom - 1]) / m_loco)
            values, vectors = _eigh_desc(k_loco)
            loco_values[chrom - 1] = values
            loco_vectors[chrom - 1] = vectors
            loco_scales.append(scale)
            loco_marker_counts.append(m_loco)

        np.save(cache / "global_kinship.npy", global_k)
        np.save(cache / "global_eigenvalues.npy", global_values)
        np.save(cache / "global_eigenvectors.npy", global_vectors)
        np.save(cache / "loco_eigenvalues.npy", loco_values)
        np.save(cache / "loco_eigenvectors.npy", loco_vectors)
        pcs.to_csv(cache / "pcs.csv")
        save_json(
            manifest_path,
            {
                **expected,
                "n_samples": len(population_ids),
                "common_marker_count": int(len(common_markers)),
                "population_polymorphic_marker_count": int(n_markers),
                "markers_per_chromosome": {str(c): int((chroms == c).sum()) for c in range(1, 11)},
                "loco_marker_counts": {str(c): int(loco_marker_counts[c - 1]) for c in range(1, 11)},
                "global_kinship_scale": global_scale,
                "loco_kinship_scales": {str(c): float(loco_scales[c - 1]) for c in range(1, 11)},
                "global_eigenvalues": global_values.tolist(),
                "elapsed": elapsed(start_time),
            },
        )
        logger.info(
            "Built LD-pruned genomic structure for %s: n=%d, markers=%d (%s)",
            label, len(population_ids), n_markers, elapsed(start_time),
        )

    spectral_path = cache / "loco_spectral_genotypes.float32.bin"
    spectral_manifest_path = cache / "loco_spectral_genotypes.json"
    expected_bytes = int(len(metadata) * len(population_ids) * np.dtype(np.float32).itemsize)
    spectral_expected = {
        "structure_hash": structure_hash,
        "sample_ids": population_ids,
        "n_variants": int(len(metadata)),
        "n_samples": int(len(population_ids)),
    }
    spectral_valid = (
        not force
        and spectral_path.exists()
        and spectral_path.stat().st_size == expected_bytes
        and _structure_cache_valid(spectral_manifest_path, spectral_expected)
    )
    if not spectral_valid:
        start_time = time.time()
        logger.info(
            "Building LOCO spectral genotype cache for %s (%d variants x %d samples; one-time)",
            label, len(metadata), len(population_ids),
        )
        spectral = np.memmap(
            spectral_path,
            dtype=np.float32,
            mode="w+",
            shape=(len(metadata), len(population_ids)),
        )
        for chrom in range(1, 11):
            indices = metadata.index[metadata["CHR"].astype(int) == chrom].to_numpy(dtype=int)
            u = loco_vectors[chrom - 1]
            for start in range(0, len(indices), block_size):
                idx = indices[start : start + block_size]
                dosage = load_variant_dosages(geno, idx)[:, sample_indices]
                spectral[idx, :] = (dosage @ u).astype(np.float32)
            spectral.flush()
            logger.info(
                "  spectral cache %s chromosome %d complete (%s)",
                label, chrom, elapsed(start_time),
            )
        save_json(spectral_manifest_path, {**spectral_expected, "elapsed": elapsed(start_time)})
        del spectral
    spectral = np.memmap(
        spectral_path,
        dtype=np.float32,
        mode="r",
        shape=(len(metadata), len(population_ids)),
    )

    return PopulationStructure(
        label=label,
        sample_ids=population_ids,
        sample_indices=sample_indices,
        pcs=pcs,
        global_kinship=global_k,
        global_eigenvalues=global_values,
        global_eigenvectors=global_vectors,
        loco_eigenvalues=loco_values,
        loco_eigenvectors=loco_vectors,
        marker_table=markers_used,
        structure_hash=structure_hash,
        spectral_genotypes=spectral,
    )
