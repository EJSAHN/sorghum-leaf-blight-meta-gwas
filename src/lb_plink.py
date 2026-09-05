from __future__ import annotations

import logging
import os
import shutil
import time
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from lb_common import elapsed, save_json

# PLINK 1 BED codes in SNP-major mode:
# 00 = homozygous allele 1; 01 = missing; 10 = heterozygous; 11 = homozygous allele 2.
_CODE = np.array([0b00, 0b10, 0b11, 0b01], dtype=np.uint8)


def _pack_variant(dosage: np.ndarray) -> bytes:
    x = np.asarray(dosage, dtype=np.int8)
    lookup_index = np.where(x < 0, 3, x).astype(np.int8)
    if np.any((lookup_index < 0) | (lookup_index > 3)):
        raise ValueError("Genotype dosages must be in {-1, 0, 1, 2}")
    codes = _CODE[lookup_index]
    n_bytes = (len(codes) + 3) // 4
    packed = np.zeros(n_bytes, dtype=np.uint8)
    for offset in range(4):
        vals = codes[offset::4]
        packed[: len(vals)] |= vals << (2 * offset)
    return packed.tobytes()


def _unpack_variant(block: bytes, n_samples: int) -> np.ndarray:
    raw = np.frombuffer(block, dtype=np.uint8)
    out = np.full(n_samples, -1, dtype=np.int8)
    reverse = {0b00: 0, 0b10: 1, 0b11: 2, 0b01: -1}
    for i in range(n_samples):
        code = int((raw[i // 4] >> (2 * (i % 4))) & 0b11)
        out[i] = reverse[code]
    return out


def hardlink_or_copy(source: str | Path, destination: str | Path) -> None:
    src = Path(source)
    dst = Path(destination)
    dst.parent.mkdir(parents=True, exist_ok=True)
    if dst.exists():
        dst.unlink()
    try:
        os.link(src, dst)
    except OSError:
        shutil.copy2(src, dst)


def write_unique_plink_files(
    geno: np.memmap,
    metadata: pd.DataFrame,
    sample_indices: np.ndarray,
    sample_ids: list[str],
    phenotypes: dict[str, np.ndarray],
    output_dir: str | Path,
    structure_variant_indices: np.ndarray,
    logger: logging.Logger,
    force: bool = False,
) -> dict[str, Path]:
    """Write reproducible PLINK BED/BIM/FAM files with guaranteed-unique marker IDs.

    The VCF cache contains dosages but not REF/ALT alleles. Dosage 0/1/2 is
    therefore encoded as A/A, A/G, G/G. This preserves association p-values;
    effect direction is relative to the synthetic A/G coding and is documented
    in the output manifest.
    """
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    n_variants = len(metadata)
    n_samples = len(sample_indices)
    if n_variants != geno.shape[0]:
        raise ValueError("Metadata and genotype cache have different variant counts")
    if n_samples != len(sample_ids):
        raise ValueError("sample_indices and sample_ids have different lengths")

    shared = out / "leaf_blight_primary"
    bed = shared.with_suffix(".bed")
    bim = shared.with_suffix(".bim")
    mapping_path = out / "variant_id_mapping.csv.gz"
    manifest_path = out / "plink_export_manifest.json"
    bytes_per_variant = (n_samples + 3) // 4
    expected_bed_size = 3 + n_variants * bytes_per_variant
    start = time.time()

    if force or not bed.exists() or bed.stat().st_size != expected_bed_size:
        logger.info(
            "Writing PLINK BED: %s variants x %d samples (%.1f MB expected)",
            f"{n_variants:,}", n_samples, expected_bed_size / 1e6,
        )
        tmp = bed.with_suffix(".bed.tmp")
        with open(tmp, "wb") as handle:
            handle.write(bytes([0x6C, 0x1B, 0x01]))
            for start_idx in range(0, n_variants, 20000):
                end_idx = min(start_idx + 20000, n_variants)
                block = np.asarray(geno[start_idx:end_idx][:, sample_indices], dtype=np.int8)
                for row in block:
                    handle.write(_pack_variant(row))
                if start_idx and start_idx % 100000 == 0:
                    logger.info(
                        "  PLINK BED %s/%s (%s)",
                        f"{start_idx:,}", f"{n_variants:,}", elapsed(start),
                    )
        os.replace(tmp, bed)
    else:
        logger.info("Reusing PLINK BED: %s", bed)

    map_table = metadata[["variant_index", "snp", "CHR", "pos", "MAF", "missing_rate"]].copy()
    map_table = map_table.reset_index(drop=True)
    if not np.array_equal(map_table["variant_index"].to_numpy(dtype=int), np.arange(n_variants)):
        raise ValueError("variant_index must be a contiguous zero-based sequence")
    map_table.insert(1, "gemma_id", [f"G{i + 1:09d}" for i in range(n_variants)])
    if map_table["gemma_id"].duplicated().any():
        raise RuntimeError("Generated GEMMA marker IDs are not unique")
    map_table.to_csv(mapping_path, index=False, compression="gzip")

    bim_table = pd.DataFrame(
        {
            "CHR": map_table["CHR"].astype(int),
            "gemma_id": map_table["gemma_id"],
            "cm": 0,
            "pos": map_table["pos"].astype(int),
            "A1": "A",
            "A2": "G",
        }
    )
    with open(bim, "w", encoding="ascii", newline="\n") as handle:
        bim_table.to_csv(handle, sep="\t", header=False, index=False, lineterminator="\n")

    paths: dict[str, Path] = {
        "shared_bed": bed,
        "shared_bim": bim,
        "variant_mapping": mapping_path,
    }
    for label, values in phenotypes.items():
        y = np.asarray(values, dtype=float)
        if len(y) != n_samples or not np.isfinite(y).all():
            raise ValueError(f"Phenotype {label} must contain {n_samples} finite values")
        prefix = out / f"leaf_blight_{label.lower()}"
        hardlink_or_copy(bed, prefix.with_suffix(".bed"))
        hardlink_or_copy(bim, prefix.with_suffix(".bim"))
        fam = pd.DataFrame(
            {
                "FID": sample_ids,
                "IID": sample_ids,
                "PAT": 0,
                "MAT": 0,
                "SEX": 0,
                "PHENO": y,
            }
        )
        with open(prefix.with_suffix(".fam"), "w", encoding="ascii", newline="\n") as handle:
            fam.to_csv(handle, sep="\t", header=False, index=False, float_format="%.12g", lineterminator="\n")
        paths[f"prefix_{label}"] = prefix

    structure_ids = map_table.iloc[np.asarray(structure_variant_indices, dtype=int)]["gemma_id"].astype(str)
    structure_list = out / "kinship_structure_snps.txt"
    with open(structure_list, "w", encoding="ascii", newline="\n") as handle:
        handle.write("\n".join(structure_ids.tolist()) + "\n")
    paths["structure_snp_list"] = structure_list

    chromosome_dir = out / "chromosome_snp_lists"
    chromosome_dir.mkdir(parents=True, exist_ok=True)
    chromosome_rows: list[dict[str, Any]] = []
    for chrom in range(1, 11):
        ids = map_table.loc[map_table["CHR"].astype(int) == chrom, "gemma_id"].astype(str).tolist()
        if not ids:
            raise ValueError(f"No variants found on chromosome {chrom}")
        path = chromosome_dir / f"chr{chrom}_snps.txt"
        with open(path, "w", encoding="ascii", newline="\n") as handle:
            handle.write("\n".join(ids) + "\n")
        if b"\r" in path.read_bytes():
            raise RuntimeError(f"CR byte found in Linux SNP list: {path}")
        paths[f"chromosome_{chrom}_snps"] = path
        chromosome_rows.append({"CHR": chrom, "n_snps": len(ids), "path": str(path)})
    pd.DataFrame(chromosome_rows).to_csv(chromosome_dir / "manifest.csv", index=False)

    # Deterministic read-back validation.
    check_rows = sorted(set([0, n_variants - 1] + np.linspace(0, n_variants - 1, 12, dtype=int).tolist()))
    with open(bed, "rb") as handle:
        if handle.read(3) != bytes([0x6C, 0x1B, 0x01]):
            raise ValueError("PLINK BED magic bytes are invalid")
        for idx in check_rows:
            handle.seek(3 + idx * bytes_per_variant)
            observed = _unpack_variant(handle.read(bytes_per_variant), n_samples)
            expected = np.asarray(geno[idx, sample_indices], dtype=np.int8)
            if not np.array_equal(observed, expected):
                raise ValueError(f"PLINK BED read-back mismatch at variant index {idx}")

    save_json(
        manifest_path,
        {
            "n_variants": n_variants,
            "n_samples": n_samples,
            "bytes_per_variant": bytes_per_variant,
            "bed_size_bytes": bed.stat().st_size,
            "synthetic_marker_id_format": "G#########",
            "allele1": "A",
            "allele2": "G",
            "dosage_encoding": "0=A/A, 1=A/G, 2=G/G, -1=missing",
            "readback_checked_variant_indices": check_rows,
            "phenotype_prefixes": {
                key: str(value) for key, value in paths.items() if key.startswith("prefix_")
            },
            "elapsed": elapsed(start),
        },
    )
    paths["manifest"] = manifest_path
    logger.info("PLINK export complete (%s)", elapsed(start))
    return paths
