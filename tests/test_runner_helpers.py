from __future__ import annotations

import numpy as np
import pandas as pd

from run_final_pipeline import _reference_audit, _scan_comparison


def _meta_table(offset: float = 0.0) -> pd.DataFrame:
    p_inc = np.array([0.01, 0.2, 0.7]) + offset
    p_sev = np.array([0.02, 0.3, 0.8]) + offset
    p_acat = np.minimum(1.0, (p_inc + p_sev) / 2)
    p_simes = np.minimum(1.0, 2 * np.minimum(p_inc, p_sev))
    out = pd.DataFrame(
        {
            "variant_index": [0, 1, 2],
            "p_INC": p_inc,
            "p_SEV": p_sev,
            "p_ACAT": p_acat,
            "p_SIMES": p_simes,
        }
    )
    for method in ["INC", "SEV", "ACAT", "SIMES"]:
        out[f"q_{method}"] = np.minimum(1.0, out[f"p_{method}"] * 3)
    return out


def test_scan_comparison_matches_by_variant_index_not_row_order() -> None:
    primary = _meta_table()
    comparator = _meta_table(0.001).iloc[::-1].reset_index(drop=True)
    result = _scan_comparison(primary, comparator, "x")
    assert set(result["n_matched_variants"]) == {3}
    assert np.all(result["spearman_rho"] > 0.99)


def test_reference_audit_passes_and_fails() -> None:
    summary = {"n_primary": 100, "score_acat_lambda_gc": 1.014}
    cfg = {
        "reference_audit": {
            "enabled": True,
            "strict_primary": True,
            "exact": {"n_primary": 100},
            "numeric": {
                "score_acat_lambda_gc": {"value": 1.014163, "absolute_tolerance": 0.02}
            },
        }
    }
    result = _reference_audit(summary, cfg)
    assert (result["status"] == "PASS").all()
    bad = dict(summary, n_primary=99)
    failed = _reference_audit(bad, cfg)
    assert failed.attrs["strict_primary"] is True
    assert failed.attrs["primary_failures"]
    assert "FAIL" in failed["status"].tolist()
