from __future__ import annotations

import numpy as np

from lb_common import bh_fdr, normalize_id, rank_int
from lb_meta import acat_pair, simes_pair


def test_normalize_id_aliases() -> None:
    assert normalize_id(" SC748-5 ") == "SC748_5"
    assert normalize_id("BTx6") == "BTX623"
    assert normalize_id(" n-42 ") == "N_42"


def test_rank_int_preserves_order_and_missing() -> None:
    values = np.array([3.0, 1.0, np.nan, 2.0, 2.0])
    transformed = rank_int(values)
    assert np.isnan(transformed[2])
    assert transformed[1] < transformed[3]
    assert transformed[3] == transformed[4]
    assert transformed[3] < transformed[0]
    assert abs(np.nanmean(transformed)) < 1e-12


def test_acat_and_simes_are_bounded_and_monotone() -> None:
    p1 = np.array([1e-8, 0.02, 0.8])
    p2 = np.array([0.5, 0.03, 0.9])
    acat = acat_pair(p1, p2)
    simes = simes_pair(p1, p2)
    assert np.all((acat > 0) & (acat <= 1))
    assert np.all((simes > 0) & (simes <= 1))
    assert acat[0] < acat[1] < acat[2]
    assert simes[0] < simes[1] < simes[2]


def test_bh_fdr_is_order_aware() -> None:
    p = np.array([0.001, 0.01, 0.2, np.nan])
    q = bh_fdr(p)
    assert q[0] <= q[1] <= q[2]
    assert q[3] == 1.0
