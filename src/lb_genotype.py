from __future__ import annotations

import gzip
import logging
import os
import re
import time
from pathlib import Path
from typing import Iterator, Any

import numpy as np
import pandas as pd

from lb_common import normalize_id, save_json, load_json, sha256_file, elapsed


def open_text(path: str | Path):
    path_str = str(path)
    if path_str.lower().endswith((".gz", ".bgz")):
        return gzip.open(path_str, "rt", encoding="utf-8", errors="replace")
    return open(path_str, "rt", encoding="utf-8", errors="replace")


def chromosome_to_int(value: str) -> int | None:
    text = re.sub(r"^chr", "", str(value).strip(), flags=re.IGNORECASE).lstrip("0")
    try:
        chrom = int(text)
    except Exception:
        return None
    return chrom if 1 <= chrom <= 10 else None


def parse_gt_token(token: str) -> int:
    if not token or token == "." or "." in token:
        return -1
    alleles = re.split(r"[|/]", token)
    if len(alleles) != 2:
        return -1
    try:
        a0, a1 = int(alleles[0]), int(alleles[1])
    except ValueError:
        return -1
    if a0 not in (0, 1) or a1 not in (0, 1):
        return -1
    return a0 + a1


def read_vcf_header(path: str | Path) -> tuple[list[str], int]:
    with open_text(path) as handle:
        line_number = 0
        for line in handle:
            line_number += 1
            if line.startswith("#CHROM"):
                fields = line.rstrip("\r\n").split("\t")
                if len(fields) < 10:
                    raise ValueError("VCF header contains no sample columns")
                return fields[9:], line_number
    raise ValueError(f"No #CHROM header found in {path}")


def _cache_valid(cache_dir: Path, vcf_hash: str, matched_ids: list[str]) -> bool:
    manifest_path = cache_dir / "cache_manifest.json"
    required = [
        manifest_path,
        cache_dir / "genotypes_int8.bin",
        cache_dir / "variant_metadata.csv.gz",
        cache_dir / "sample_ids.txt",
    ]
    if not all(path.exists() for path in required):
        return False
    try:
        manifest = load_json(manifest_path)
    except Exception:
        return False
    return (
        manifest.get("vcf_sha256") == vcf_hash
        and manifest.get("sample_ids") == matched_ids
        and int(manifest.get("n_variants", 0)) > 0
    )


def build_or_load_genotype_cache(
    vcf_path: str | Path,
    phenotype_ids: list[str],
    cache_dir: str | Path,
    maf_min: float,
    missing_max: float,
    logger: logging.Logger,
) -> tuple[np.memmap, pd.DataFrame, list[str], dict[str, Any]]:
    cache = Path(cache_dir)
    cache.mkdir(parents=True, exist_ok=True)
    vcf_path = Path(vcf_path)

    hash_note = vcf_path.parent / "input_hashes.csv"
    vcf_hash = None
    if hash_note.exists():
        try:
            hashes = pd.read_csv(hash_note)
            row = hashes.loc[hashes["Name"].astype(str).str.lower() == vcf_path.name.lower()]
            if not row.empty:
                vcf_hash = str(row.iloc[0]["SHA256"]).upper()
        except Exception:
            vcf_hash = None
    if not vcf_hash:
        logger.info("Computing SHA-256 for VCF (one-time cache validation)...")
        vcf_hash = sha256_file(vcf_path)

    raw_samples, _ = read_vcf_header(vcf_path)
    normalized_to_raw: dict[str, str] = {}
    normalized_to_col: dict[str, int] = {}
    duplicate_norm: list[str] = []
    for idx, raw in enumerate(raw_samples):
        norm = normalize_id(raw)
        if norm in normalized_to_col:
            duplicate_norm.append(norm)
        else:
            normalized_to_col[norm] = idx
            normalized_to_raw[norm] = raw
    if duplicate_norm:
        raise ValueError(f"Duplicate normalized VCF sample IDs: {sorted(set(duplicate_norm))[:10]}")

    matched_ids = [pid for pid in phenotype_ids if pid in normalized_to_col]
    missing_ids = [pid for pid in phenotype_ids if pid not in normalized_to_col]
    if len(matched_ids) < 30:
        raise ValueError(f"Only {len(matched_ids)} phenotype IDs matched VCF samples")

    if _cache_valid(cache, vcf_hash, matched_ids):
        manifest = load_json(cache / "cache_manifest.json")
        metadata = pd.read_csv(cache / "variant_metadata.csv.gz")
        n_variants = int(manifest["n_variants"])
        n_samples = int(manifest["n_samples"])
        geno = np.memmap(cache / "genotypes_int8.bin", dtype=np.int8, mode="r", shape=(n_variants, n_samples))
        logger.info("Reusing genotype cache: %s variants x %s samples", n_variants, n_samples)
        return geno, metadata, matched_ids, manifest

    logger.info(
        "Building genotype cache from %s; matched %s/%s phenotype IDs",
        vcf_path,
        len(matched_ids),
        len(phenotype_ids),
    )
    logger.info("Phenotype IDs without VCF genotypes: %s", ", ".join(missing_ids) if missing_ids else "none")

    selected_vcf_columns = [9 + normalized_to_col[sid] for sid in matched_ids]
    metadata_rows: list[tuple[int, str, int, int, float, float]] = []
    binary_path = cache / "genotypes_int8.bin"
    start = time.time()
    total_variant_lines = 0
    kept = 0
    skip_non_autosome = 0
    skip_multiallelic = 0
    skip_maf = 0
    skip_missing = 0
    skip_no_calls = 0

    with open_text(vcf_path) as handle, open(binary_path, "wb") as binary:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            total_variant_lines += 1
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) < 10:
                continue
            chrom = chromosome_to_int(fields[0])
            if chrom is None:
                skip_non_autosome += 1
                continue
            alt = fields[4]
            if "," in alt:
                skip_multiallelic += 1
                continue
            try:
                pos = int(fields[1])
            except ValueError:
                continue
            fmt = fields[8].split(":")
            try:
                gt_index = fmt.index("GT")
            except ValueError:
                gt_index = 0

            dosages = np.empty(len(selected_vcf_columns), dtype=np.int8)
            for out_idx, field_idx in enumerate(selected_vcf_columns):
                sample_field = fields[field_idx]
                parts = sample_field.split(":")
                gt_token = parts[gt_index] if gt_index < len(parts) else "."
                dosages[out_idx] = parse_gt_token(gt_token)

            called = dosages >= 0
            n_called = int(called.sum())
            if n_called == 0:
                skip_no_calls += 1
                continue
            missing_rate = 1.0 - n_called / len(dosages)
            if missing_rate > missing_max:
                skip_missing += 1
                continue
            alt_freq = float(dosages[called].sum() / (2.0 * n_called))
            maf = min(alt_freq, 1.0 - alt_freq)
            if maf < maf_min:
                skip_maf += 1
                continue

            snp_id = fields[2]
            if snp_id in {"", "."}:
                snp_id = f"S{chrom:02d}_{pos}"
            binary.write(dosages.tobytes(order="C"))
            metadata_rows.append((kept, snp_id, chrom, pos, maf, missing_rate))
            kept += 1

            if total_variant_lines % 100000 == 0:
                logger.info(
                    "VCF scan: %s variant lines, %s retained (%s)",
                    f"{total_variant_lines:,}",
                    f"{kept:,}",
                    elapsed(start),
                )

    metadata = pd.DataFrame(
        metadata_rows,
        columns=["variant_index", "snp", "CHR", "pos", "MAF", "missing_rate"],
    )
    metadata.to_csv(cache / "variant_metadata.csv.gz", index=False, compression="gzip")
    (cache / "sample_ids.txt").write_text("\n".join(matched_ids) + "\n", encoding="utf-8")
    sample_map = pd.DataFrame(
        {
            "ID_std": matched_ids,
            "VCF_sample": [normalized_to_raw[x] for x in matched_ids],
            "VCF_sample_index": [normalized_to_col[x] for x in matched_ids],
        }
    )
    sample_map.to_csv(cache / "sample_map.csv", index=False)

    manifest = {
        "vcf_path": str(vcf_path),
        "vcf_sha256": vcf_hash,
        "vcf_size_bytes": int(vcf_path.stat().st_size),
        "raw_vcf_sample_count": len(raw_samples),
        "phenotype_id_count": len(phenotype_ids),
        "n_samples": len(matched_ids),
        "sample_ids": matched_ids,
        "missing_phenotype_ids": missing_ids,
        "n_variant_lines": total_variant_lines,
        "n_variants": kept,
        "maf_min": maf_min,
        "missing_max": missing_max,
        "skip_non_autosome": skip_non_autosome,
        "skip_multiallelic": skip_multiallelic,
        "skip_maf": skip_maf,
        "skip_missing": skip_missing,
        "skip_no_calls": skip_no_calls,
    }
    save_json(cache / "cache_manifest.json", manifest)

    if kept == 0:
        raise ValueError("No variants passed VCF quality control")
    expected_bytes = kept * len(matched_ids)
    actual_bytes = binary_path.stat().st_size
    if actual_bytes != expected_bytes:
        raise RuntimeError(f"Genotype cache size mismatch: expected {expected_bytes}, got {actual_bytes}")

    geno = np.memmap(binary_path, dtype=np.int8, mode="r", shape=(kept, len(matched_ids)))
    logger.info("Genotype cache complete: %s variants x %d samples in %s", f"{kept:,}", len(matched_ids), elapsed(start))
    return geno, metadata, matched_ids, manifest


def iter_genotype_blocks(geno: np.memmap, block_size: int) -> Iterator[tuple[int, int, np.ndarray]]:
    n_variants = geno.shape[0]
    for start in range(0, n_variants, block_size):
        end = min(n_variants, start + block_size)
        block = np.asarray(geno[start:end], dtype=np.float64)
        called = block >= 0
        counts = called.sum(axis=1)
        sums = np.where(called, block, 0.0).sum(axis=1)
        means = np.divide(sums, counts, out=np.zeros_like(sums), where=counts > 0)
        if np.any(~called):
            row_idx, col_idx = np.where(~called)
            block[row_idx, col_idx] = means[row_idx]
        yield start, end, block


def compute_pcs(
    geno: np.memmap,
    sample_ids: list[str],
    cache_dir: str | Path,
    n_components: int,
    block_size: int,
    logger: logging.Logger,
) -> pd.DataFrame:
    cache = Path(cache_dir)
    pcs_path = cache / "sample_pcs.csv"
    eig_path = cache / "pca_eigenvalues.csv"
    if pcs_path.exists() and eig_path.exists():
        pcs = pd.read_csv(pcs_path).set_index("ID_std")
        if len(pcs) == len(sample_ids) and list(pcs.index) == sample_ids and pcs.shape[1] >= n_components:
            logger.info("Reusing cached principal components")
            return pcs

    logger.info("Computing PCA from the matched, QC-filtered genotype matrix")
    n_samples = geno.shape[1]
    covariance = np.zeros((n_samples, n_samples), dtype=np.float64)
    total_centered_ss = 0.0
    start_time = time.time()
    for start, end, block in iter_genotype_blocks(geno, block_size):
        block -= block.mean(axis=1, keepdims=True)
        covariance += block.T @ block
        total_centered_ss += float(np.square(block).sum())
        if start and start % (block_size * 10) == 0:
            logger.info("PCA covariance: %s/%s variants (%s)", f"{end:,}", f"{geno.shape[0]:,}", elapsed(start_time))

    eigenvalues, eigenvectors = np.linalg.eigh(covariance)
    order = np.argsort(eigenvalues)[::-1]
    eigenvalues = np.maximum(eigenvalues[order], 0.0)
    eigenvectors = eigenvectors[:, order]
    keep = min(n_components, n_samples - 1)
    scores = eigenvectors[:, :keep] * np.sqrt(eigenvalues[:keep])[None, :]
    pcs = pd.DataFrame(scores, index=sample_ids, columns=[f"PC{i+1}" for i in range(keep)])
    pcs.index.name = "ID_std"
    pcs.to_csv(pcs_path)
    explained = eigenvalues / eigenvalues.sum() if eigenvalues.sum() > 0 else np.zeros_like(eigenvalues)
    pd.DataFrame(
        {
            "component": np.arange(1, len(eigenvalues) + 1),
            "eigenvalue": eigenvalues,
            "explained_fraction": explained,
        }
    ).to_csv(eig_path, index=False)
    logger.info("PCA complete in %s", elapsed(start_time))
    return pcs


def load_variant_dosages(geno: np.memmap, indices: np.ndarray | list[int]) -> np.ndarray:
    block = np.asarray(geno[np.asarray(indices, dtype=int)], dtype=np.float64)
    called = block >= 0
    counts = called.sum(axis=1)
    sums = np.where(called, block, 0.0).sum(axis=1)
    means = np.divide(sums, counts, out=np.zeros_like(sums), where=counts > 0)
    if np.any(~called):
        rows, cols = np.where(~called)
        block[rows, cols] = means[rows]
    return block
