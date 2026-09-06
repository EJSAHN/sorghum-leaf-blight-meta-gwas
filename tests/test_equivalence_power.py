from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
import pandas as pd

from lb_equivalence import compute_or_load_equivalence
from lb_power import build_power_table


def _logger() -> logging.Logger:
    logger = logging.getLogger("test_equivalence")
    logger.handlers.clear()
    logger.addHandler(logging.NullHandler())
    return logger


def test_exact_and_complement_equivalence(tmp_path: Path) -> None:
    raw = np.array(
        [
            [0, 0, 1, 2],
            [0, 0, 1, 2],
            [2, 2, 1, 0],
            [0, 1, 1, 2],
            [0, 0, 1, 2],
        ],
        dtype=np.int8,
    )
    path = tmp_path / "g.bin"
    mm = np.memmap(path, dtype=np.int8, mode="w+", shape=raw.shape)
    mm[:] = raw
    mm.flush()
    metadata = pd.DataFrame(
        {
            "variant_index": np.arange(5),
            "snp": [f"s{i}" for i in range(5)],
            "CHR": [1, 1, 1, 1, 2],
            "pos": [10, 20, 30, 40, 50],
            "MAF": [0.25] * 5,
            "missing_rate": [0.0] * 5,
        }
    )
    result = compute_or_load_equivalence(mm, metadata, np.arange(4), tmp_path / "eq", _logger(), force=True)
    assert result.global_inverse[0] == result.global_inverse[1]
    assert result.global_inverse[0] == result.global_inverse[2]
    assert result.global_inverse[0] == result.global_inverse[4]
    assert result.chrom_inverse[0] == result.chrom_inverse[1]
    assert result.chrom_inverse[0] == result.chrom_inverse[2]
    assert result.chrom_inverse[0] != result.chrom_inverse[4]
    assert result.global_inverse[3] != result.global_inverse[0]
    del mm


def test_detectable_effect_increases_for_lower_maf_and_higher_power() -> None:
    table = build_power_table(
        sample_sizes=[100],
        mafs=[0.05, 0.50],
        target_powers=[0.80, 0.90],
        alpha_definitions={"a": 0.05 / 500000},
        n_covariates=1,
    )
    pivot = table.set_index(["maf", "target_power"])["minimum_detectable_beta_residual_SD_per_allele"]
    assert pivot.loc[(0.05, 0.80)] > pivot.loc[(0.50, 0.80)]
    assert pivot.loc[(0.50, 0.90)] > pivot.loc[(0.50, 0.80)]
