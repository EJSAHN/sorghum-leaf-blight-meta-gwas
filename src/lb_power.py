from __future__ import annotations

import math
from typing import Iterable

import numpy as np
import pandas as pd
from scipy import optimize, stats


def two_sided_power(beta_std: float, n: int, maf: float, alpha: float, n_covariates: int = 1) -> float:
    """Approximate power for an additive SNP effect on a standardized trait.

    The approximation uses a noncentral t test with genotype variance 2p(1-p).
    It is intended as a transparent detectable-effect analysis for a small GWAS,
    not as a substitute for simulation under the full mixed model.
    """
    df = max(int(n - n_covariates - 1), 2)
    geno_var = 2.0 * maf * (1.0 - maf)
    ncp = abs(float(beta_std)) * math.sqrt(n * geno_var)
    critical = stats.t.ppf(1.0 - alpha / 2.0, df=df)
    lower = float(stats.nct.cdf(-critical, df, ncp))
    upper = float(stats.nct.sf(critical, df, ncp))
    if not np.isfinite(lower):
        lower = 0.0
    if not np.isfinite(upper):
        # At large positive noncentrality, right-tail power tends to one.
        upper = 1.0 if ncp > critical else 0.0
    return float(np.clip(lower + upper, 0.0, 1.0))


def minimum_detectable_beta(
    n: int,
    maf: float,
    alpha: float,
    target_power: float,
    n_covariates: int = 1,
) -> float:
    if not (0 < alpha < 1 and 0 < target_power < 1 and 0 < maf <= 0.5):
        raise ValueError("Invalid alpha, target power, or MAF")
    fn = lambda beta: two_sided_power(beta, n, maf, alpha, n_covariates) - target_power
    upper = 0.25
    while fn(upper) < 0 and upper < 20:
        upper *= 2
    if fn(upper) < 0:
        return float("nan")
    return float(optimize.brentq(fn, 0.0, upper, xtol=1e-10, rtol=1e-10))


def build_power_table(
    sample_sizes: Iterable[int],
    mafs: Iterable[float],
    target_powers: Iterable[float],
    alpha_definitions: dict[str, float],
    n_covariates: int = 1,
) -> pd.DataFrame:
    rows = []
    for n in sample_sizes:
        for maf in mafs:
            for alpha_label, alpha in alpha_definitions.items():
                for power in target_powers:
                    beta = minimum_detectable_beta(n, maf, alpha, power, n_covariates)
                    marginal_r2 = beta * beta * 2.0 * maf * (1.0 - maf)
                    rows.append(
                        {
                            "n_samples": int(n),
                            "maf": float(maf),
                            "alpha_definition": alpha_label,
                            "alpha": float(alpha),
                            "target_power": float(power),
                            "minimum_detectable_standardized_per_allele_beta": beta,
                            "approximate_marginal_variance_explained": marginal_r2,
                            "approximate_marginal_variance_explained_percent": 100.0 * marginal_r2,
                            "residual_df_assumption": int(max(n - n_covariates - 1, 2)),
                            "method_note": "Noncentral-t approximation; standardized phenotype; additive genotype variance 2p(1-p)",
                        }
                    )
    return pd.DataFrame(rows)
