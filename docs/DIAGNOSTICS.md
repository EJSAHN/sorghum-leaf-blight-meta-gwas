# Model diagnostics

The primary association workflow and the auxiliary diagnostics have separate
entry points. The latter reads an existing analysis project; it does not alter
its input files or regenerate the main association scans.

```text
python src/run_sensitivity.py --project-root PROJECT --output-dir NEW_OUTPUT --stage all
```

The output must be a new directory outside `PROJECT`. Individual stages are
`ols`, `effects`, and `power`. Existing NumPy, SciPy and pandas installations are
used. The `effects` stage also requires GEMMA 0.98.5, available through WSL on
Windows, or supplied using `--gemma-executable` on native Linux.

## OLS

OLS diagnostics use the four accession phenotypes `INC_ADJ_RAW`, `INC_ADJ_INT`,
`SEV_ADJ_RAW`, and `SEV_ADJ_INT`, and 0, 1, 2, 3, or 5 genomic PCs. The intercept
and the same PC covariates are projected out of both phenotype and genotype.
The residual degrees of freedom are n minus the covariate rank minus one.
Missing dosages are mean-imputed across all 102 QC samples before selecting
the 100 primary accessions. All QC markers are scanned; FDR is recalculated
for each of the 20 trait/covariate combinations.

The routine uses the structure marker list and PC file in the current project.
It reconstructs the PCs from genotype rows and checks nested covariate spaces,
allowing eigenvector sign changes. No replacement PC basis is silently used.
The associated diagnostic dataset uses 8,566 structure markers.

`Transformation_Checks.csv` contains lambda GC, minimum p, and the number of
markers with BH-adjusted q < 0.10. This is a diagnostic table, not the primary
mixed-model discovery set.

## Candidate effects

Primary candidate coefficients are obtained with GEMMA `-lmm 3`, using the same
sample order, phenotype and chromosome-specific kinship as the primary scans.
The primary score p-values are checked against the existing candidate table.
These coefficients are compared only with saved `-lmm 3` coefficients from the
location-omission scans. Both therefore use covariance estimated under the null
by maximum likelihood. Files produced with `-lmm 4` retain alternative-model
REML coefficients and are not used as the primary side of this comparison.

`Candidate_Effect_Direction_Check.csv` and `Candidate_Direction_Summary.csv`
report matched-estimator comparisons. Effects are in inverse-normal phenotype
units per synthetic allele1 A. A/G are dosage labels, not original nucleotide
identities. Direction retention is not independent replication. Wald intervals
in the main candidate table are unadjusted, alternative-model intervals, not
confidence intervals obtained by inverting a score test.

## Detectable effects

`Power.csv` reports `minimum_detectable_beta_residual_SD_per_allele`.
The noncentral-t approximation uses beta/sigma_e, residual degrees of freedom
n minus the number of covariates minus one, and genotype variance 2f(1-f).
For the intercept-only approximation this is n - 2. Kinship, inbreeding-specific
genotype variance, and uncertainty in estimated phenotypes are not modeled.
No observed or total-phenotype explained-variance estimate is inferred.

## Optional independent reference tables

To compare numerical outputs against separately supplied supporting tables:

```text
python src/run_sensitivity.py --project-root PROJECT --output-dir NEW_OUTPUT --reference-dir TABLE_DIRECTORY
```

Expected filenames are `Accession_Phenotypes.csv`, `Population_Definition.csv`,
`Structure_PCs.csv`, `Structure_Markers.csv`, `Top_Ranked_Candidates.csv`,
`LOO_Candidates.csv`, `Power.csv`, and `Transformation_Checks.csv`. Only the
files relevant to the selected stage are used. These study tables are not
redistributed in the software source archive. Without `--reference-dir`,
inputs come from the local project and no independent OLS reference comparison
is claimed. The summary records which type of comparison was performed.
