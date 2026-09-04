from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
import pandas as pd

from lb_plink import _pack_variant, _unpack_variant, write_unique_plink_files


def _logger() -> logging.Logger:
    logger = logging.getLogger("test_plink")
    logger.handlers.clear()
    logger.addHandler(logging.NullHandler())
    return logger


def test_pack_unpack_round_trip() -> None:
    dosage = np.array([0, 1, 2, -1, 2, 0, 1], dtype=np.int8)
    packed = _pack_variant(dosage)
    observed = _unpack_variant(packed, len(dosage))
    np.testing.assert_array_equal(observed, dosage)


def test_unique_plink_export_and_readback(tmp_path: Path) -> None:
    # Include one marker on each sorghum chromosome because the production
    # exporter deliberately requires chromosomes 1-10 for LOCO lists.
    raw = np.array(
        [[(i + j) % 3 for j in range(5)] for i in range(10)],
        dtype=np.int8,
    )
    mmap_path = tmp_path / "geno.bin"
    mm = np.memmap(mmap_path, dtype=np.int8, mode="w+", shape=raw.shape)
    mm[:] = raw
    mm.flush()
    metadata = pd.DataFrame(
        {
            "variant_index": np.arange(10),
            "snp": ["dup", "dup", *[f"s{i}" for i in range(2, 10)]],
            "CHR": np.arange(1, 11),
            "pos": np.arange(1, 11) * 100,
            "MAF": [0.2] * 10,
            "missing_rate": [0.0] * 10,
        }
    )
    paths = write_unique_plink_files(
        mm,
        metadata,
        np.arange(5),
        [f"S{i}" for i in range(5)],
        {"INC": np.arange(5, dtype=float), "SEV": np.arange(5, dtype=float)[::-1]},
        tmp_path / "plink",
        np.array([0, 2]),
        _logger(),
        force=True,
    )
    mapping = pd.read_csv(paths["variant_mapping"])
    assert mapping["gemma_id"].is_unique
    assert mapping["gemma_id"].iloc[0] == "G000000001"
    assert mapping["gemma_id"].iloc[-1] == "G000000010"
    assert paths["shared_bed"].stat().st_size == 3 + 10 * 2
    for chrom in range(1, 11):
        content = paths[f"chromosome_{chrom}_snps"].read_bytes()
        assert b"\r" not in content
        assert content.endswith(b"\n")
    del mm
