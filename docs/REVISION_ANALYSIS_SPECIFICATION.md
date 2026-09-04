# Revision analysis specification

## Study-specific data flow

```text
Original plot-level leaf-blight records
  -> explicit missing-value and correction audit
  -> field-stratum-adjusted accession estimates on the raw scale
  -> rank-based inverse-normal transformation
  -> genotype matching
  -> primary n=100 population (at least four field locations)
  -> official GEMMA 0.98.5 LOCO association tests
```

The all-matched n=102 population is a prespecified sensitivity analysis. The
workflow never substitutes the earlier first-record-per-accession phenotype.

## Phenotypes

- `INC_ADJ_INT`: adjusted incidence.
- `SEV_ADJ_INT`: adjusted standard severity.
- `SEV_COND_ADJ_INT`: severity adjusted for plot-level incidence.
- `SEV_POS_ADJ_INT`: severity estimated only from disease-positive plots.

The latter two are sensitivity analyses intended to separate post-infection
symptom burden from disease establishment.

## Inferential hierarchy

Official GEMMA 0.98.5 is supplied with audited global and chromosome-specific
kinship matrices.

1. LOCO score test: primary inference.
2. LOCO likelihood-ratio test: confirmatory sensitivity.
3. LOCO Wald test: diagnostic only.

The Wald test produced an unstable extreme tail in this small panel that was not
supported by score or likelihood-ratio tests. The software preserves the Wald
output for auditability but excludes it from discovery claims.

## Dual-trait analysis

Incidence and severity p-values are combined with both ACAT and Simes. ACAT is
not interpreted as proof that a locus affects both traits; shared architecture
is evaluated separately using genome-wide concordance, enrichment, conditional
severity, and leave-one-location-out analyses.

## Multiple testing and candidate regions

- BH-FDR is calculated across all 567,758 post-QC SNP rows.
- A sensitivity analysis collapses exact/allele-complement patterns into
  chromosome-specific unique LOCO test classes.
- Candidate regions are selected from primary score-test rankings, collapsed
  for global genotype aliases, and grouped by genotype-based local LD.
- A top-ranked region is labelled exploratory unless it passes the prespecified
  FDR criterion.
- No Bayesian credible sets or formal fine-mapping are performed.

## Power and uncertainty

Trait-specific beta estimates, standard errors, and 95% confidence intervals
are retained. A noncentral-t detectable-effect analysis quantifies the minimum
standardized per-allele effect detectable at selected allele frequencies and
multiplicity thresholds. This calculation is presented as an approximation,
not as a simulation of the full mixed model.
