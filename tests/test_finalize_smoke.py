from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
import pandas as pd

from lb_finalize import finalize_results


def test_compact_finalization_smoke(tmp_path: Path) -> None:
    rng = np.random.default_rng(12)
    n_variants, n_samples = 50, 12
    raw = rng.integers(0, 3, size=(n_variants, n_samples), dtype=np.int8)
    geno_path = tmp_path / "geno.bin"
    geno = np.memmap(geno_path, dtype=np.int8, mode="w+", shape=raw.shape)
    geno[:] = raw
    geno.flush()
    metadata = pd.DataFrame(
        {
            "variant_index": np.arange(n_variants),
            "snp": [f"S{(i % 10) + 1:02d}_{1000 + i}" for i in range(n_variants)],
            "CHR": (np.arange(n_variants) % 10) + 1,
            "pos": 1000 + np.arange(n_variants) * 1000,
            "MAF": np.clip(raw.mean(axis=1) / 2.0, 0.05, 0.5),
            "missing_rate": np.zeros(n_variants),
        }
    )

    def scan() -> pd.DataFrame:
        return metadata.assign(
            beta=rng.normal(size=n_variants),
            se=np.full(n_variants, 0.2),
            l_remle=np.ones(n_variants),
            p_score=np.linspace(1e-5, 0.9, n_variants),
            p_lrt=np.linspace(2e-5, 0.95, n_variants),
            p_wald=np.linspace(3e-5, 0.99, n_variants),
        )

    scans = {
        ("loco", "INC"): scan(),
        ("loco", "SEV"): scan(),
        ("loco", "SEV_COND"): scan(),
        ("loco", "SEV_POS"): scan(),
    }
    gff = tmp_path / "genes.gff3"
    gff.write_text(
        "##gff-version 3\n"
        "Chr01\ttest\tgene\t1\t50000\t.\t+\t.\tID=Sobic.001G000100.v3.2;Name=Sobic.001G000100\n",
        encoding="utf-8",
    )
    logger = logging.getLogger("finalize-smoke")
    logger.handlers.clear()
    logger.addHandler(logging.NullHandler())
    summary = finalize_results(
        tmp_path,
        geno,
        metadata,
        np.arange(n_samples),
        scans,
        scans,
        scans,
        gff,
        logger,
        candidate_top_n=5,
        index_p=1e-4,
        secondary_p=1e-3,
        ld_window_kb=1000,
        ld_r2=0.2,
        power_sample_sizes=(n_samples,),
    )
    assert summary["primary_test"] == "score"
    assert (tmp_path / "05_results" / "LeafBlight_MultiEnvironment_Analysis_Results.xlsx").exists()
    assert (tmp_path / "05_results" / "analysis_manifest.json").exists()
    del geno
