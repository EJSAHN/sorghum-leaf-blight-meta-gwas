"""Independent-observation detectable effects in residual-standard-deviation units."""
from __future__ import annotations

import math
from typing import Iterable
import numpy as np
import pandas as pd
from scipy import optimize, stats


def two_sided_power(beta_std: float, n: int, maf: float, alpha: float, n_covariates: int = 1) -> float:
    """Power for beta/sigma_e with assumed genotype variance 2f(1-f).

    n_covariates includes the intercept. Genomic covariance and uncertainty in
    estimated accession phenotypes are not included in this approximation.
    """
    if n_covariates < 1 or n <= n_covariates + 1 or not (0 < maf <= .5 and 0 < alpha < 1):
        raise ValueError('Require positive residual degrees of freedom, valid MAF and alpha')
    df = int(n - n_covariates - 1)
    critical = stats.t.ppf(1 - alpha / 2, df)
    delta = abs(float(beta_std)) * math.sqrt(n * 2 * maf * (1 - maf))
    value = stats.nct.sf(critical, df, -delta) + stats.nct.sf(critical, df, delta)
    if not np.isfinite(value):
        raise ArithmeticError('Non-finite noncentral-t probability')
    return float(np.clip(value, 0, 1))


def minimum_detectable_beta(n: int, maf: float, alpha: float, target_power: float,
                            n_covariates: int = 1) -> float:
    """Minimum absolute beta/sigma_e; this quantity is not phenotype R-squared."""
    if not (0 < target_power < 1):
        raise ValueError('Target power must be between zero and one')
    two_sided_power(0, n, maf, alpha, n_covariates)
    if target_power <= alpha:
        return 0.0
    fun = lambda beta: two_sided_power(beta, n, maf, alpha, n_covariates) - target_power
    high = .25
    while fun(high) < 0 and high < 32:
        high *= 2
    if fun(high) < 0:
        raise ArithmeticError('Unable to bracket detectable effect')
    return float(optimize.brentq(fun, 0, high, xtol=1e-10, rtol=1e-10))


def build_power_table(sample_sizes: Iterable[int], mafs: Iterable[float],
                      target_powers: Iterable[float], alpha_definitions: dict[str, float],
                      n_covariates: int = 1) -> pd.DataFrame:
    """Return residual-SD effect thresholds without inferring total variance explained."""
    rows = []
    sample_sizes, mafs, target_powers = list(sample_sizes), list(mafs), list(target_powers)
    for n in sample_sizes:
        for maf in mafs:
            for label, alpha in alpha_definitions.items():
                for power in target_powers:
                    beta = minimum_detectable_beta(n, maf, alpha, power, n_covariates)
                    rows.append(dict(
                        n_samples=int(n), maf=float(maf), alpha_definition=label,
                        alpha=float(alpha), target_power=float(power),
                        minimum_detectable_beta_residual_SD_per_allele=beta,
                        residual_df_assumption=int(n-n_covariates-1),
                        method_note='Independent observations; residual SD per allele; additive genotype '
                                    'variance 2f(1-f); no genomic covariance or phenotype-estimation uncertainty.'
                    ))
    return pd.DataFrame(rows)
