# Analysis specification

## Data flow

```text
Plot-level leaf-blight observations
  -> missing-value and range validation
  -> field-stratum-adjusted accession estimates on the raw scale
  -> rank-based inverse-normal transformation
  -> genotype matching
  -> primary population represented in at least four field locations
  -> official GEMMA 0.98.5 LOCO association tests
```

The full genotype-matched population is retained as a prespecified sensitivity analysis.

## Phenotypes

- `INC_ADJ_INT`: adjusted incidence.
- `SEV_ADJ_INT`: adjusted standard severity.
- `SEV_COND_ADJ_INT`: severity adjusted for plot-level incidence.
- `SEV_POS_ADJ_INT`: severity estimated only from disease-positive plots.

The latter two phenotypes evaluate symptom burden after accounting for disease establishment.

## Association tests

Official GEMMA 0.98.5 is supplied with validated global and chromosome-specific kinship matrices.

1. LOCO score test: primary inference.
2. LOCO likelihood-ratio test: confirmatory sensitivity.
3. LOCO Wald test: diagnostic.

Wald results are retained to characterize small-sample variance-component sensitivity and are not used to define genome-wide discoveries.

## Dual-trait analysis

Incidence and severity p-values are combined with both ACAT and Simes. A combined p-value is not interpreted as proof that a marker affects both traits. Shared signal is evaluated separately using genome-wide concordance, enrichment, incidence-adjusted severity, and leave-one-location-out analyses.

## Multiple testing and region summaries

- BH-FDR is calculated across all 567,758 post-QC SNP rows.
- A sensitivity analysis groups exact and allele-complement patterns into chromosome-specific unique LOCO test classes.
- Region summaries begin with primary score-test rankings, group global genotype aliases, and apply genotype-based local LD clumping.
- A top-ranked region is labelled exploratory unless it satisfies the prespecified FDR criterion.
- Bayesian credible sets and formal fine-mapping are not performed.

## Power and uncertainty

Main candidate effects and standard errors are alternative-model REML estimates from GEMMA mode 4. The beta +/- 1.96 SE intervals are unadjusted Wald intervals. Matched-estimator location-omission comparisons use GEMMA mode 3 on both sides and are produced by the auxiliary diagnostic command. The noncentral-t detectable-effect calculation reports beta in residual SD per allele, under independent observations and genotype variance 2f(1-f); no total-phenotype explained variance is inferred.
