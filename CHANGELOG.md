# Changelog

## 2.0.2

- Add raw/INT-by-PC OLS diagnostics and genotype-based PC reconstruction.
- Add matched-GEMMA-mode-3 primary/location-omission coefficient comparisons.
- Remove mixed-estimator direction comparisons from the main output.
- State detectable effects in residual SD per allele and remove incorrect
  total-variance-explained labels.
- Label candidate effects and Wald intervals by estimator and dosage coding.
- Verify the complete source inventory before publishing a tagged release.

## 2.0.1 — 2026-09-04

### Distribution and metadata

- Uses explicit project and input-root parameters throughout the source distribution.
- Removes obsolete comparison fields and repository-maintenance helpers.
- Uses neutral output names for summaries, workbooks, tables, figures, and logs.
- Adds Zenodo metadata and updates citation metadata for version 2.0.1.
- Adds tests that verify consistency among release metadata files.

### Analysis implementation

- Preserves the statistical models, thresholds, and inferential hierarchy of version 2.0.0.
- Retains official GEMMA 0.98.5 LOCO score tests as primary inference.
- Retains LOCO likelihood-ratio tests as confirmatory sensitivity analyses.
- Retains LOCO Wald tests as diagnostics only.

## 2.0.0 — 2026-09-04

### Analysis workflow

- Constructs accession-level phenotypes from plot-level field observations before transformation.
- Models field-stratum effects and distinguishes numeric zero from missing observations.
- Uses 100 accessions represented in at least four locations as the primary population and 102 genotype-matched accessions as a sensitivity population.
- Computes standard severity, incidence-adjusted severity, and disease-positive-only severity.
- Uses LD-pruned structure markers and chromosome-specific LOCO kinship matrices.
- Compares ACAT and Simes dual-trait p-value combinations.
- Groups genotype-equivalent markers before genotype-based LD clumping.
- Reports effects, standard errors, confidence intervals, detectable-effect calculations, and leave-one-location-out stability.

### Reproducibility

- Uses a deterministic configuration and recorded input checksums.
- Verifies sample order, PLINK BED encoding, kinship matrices, and principal output dimensions.
- Excludes raw data, derived study results, and binary caches from the source repository.
