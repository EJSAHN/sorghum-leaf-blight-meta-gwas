from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from lb_gemma import map_assoc, read_assoc, windows_to_wsl


def test_windows_to_wsl() -> None:
    windows_path = "X:" + r"\data\study\file.txt"
    expected = "/mnt/" + "x/data/study/file.txt"
    assert windows_to_wsl(windows_path) == expected


def test_read_and_map_assoc(tmp_path: Path) -> None:
    assoc = tmp_path / "x.assoc.txt"
    assoc.write_text(
        "chr rs ps n_miss allele1 allele0 af beta se logl_H1 l_remle p_wald p_lrt p_score\n"
        "1 G000000001 100 0 A G 0.2 0.5 0.1 -10 1.2 1e-4 2e-4 3e-4\n"
        "2 G000000002 200 0 A G 0.3 -0.2 0.2 -11 0.8 0.1 0.2 0.3\n",
        encoding="ascii",
    )
    parsed = read_assoc(assoc)
    mapping = pd.DataFrame(
        {
            "variant_index": [0, 1],
            "gemma_id": ["G000000001", "G000000002"],
            "snp": ["s1", "s2"],
            "CHR": [1, 2],
            "pos": [100, 200],
            "MAF": [0.2, 0.3],
            "missing_rate": [0.0, 0.0],
        }
    )
    mapped = map_assoc(parsed, mapping)
    assert mapped["snp"].tolist() == ["s1", "s2"]
    np.testing.assert_allclose(mapped["p_score"], [3e-4, 0.3])
