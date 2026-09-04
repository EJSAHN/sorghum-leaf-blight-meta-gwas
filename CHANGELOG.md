# Changelog

## 2.0.0 — 2026-09-04

### Scientific corrections

- Supersedes the pre-revision exploratory workflow.
- Corrects phenotype construction so plot-level field structure is modeled
  before accession-level transformation.
- Replaces first-record retention with adjusted accession estimates.
- Defines the primary population as the 100 genotype-matched accessions observed
  in at least four locations; retains n=102 as sensitivity.
- Uses official GEMMA 0.98.5 LOCO score tests as primary inference and LRT as
  confirmatory sensitivity.
- Restricts Wald tests to diagnostics because their extreme tail was not
  supported by score/LRT tests in this small panel.
- Removes the prior 37-locus result and all credible-set/fine-mapping claims from
  the reproducible inferential workflow.

### New analyses

- Standard, incidence-adjusted, and disease-positive-only severity phenotypes.
- K+PC3, n=102, and leave-one-location-out sensitivities.
- ACAT and Simes dual-trait comparisons.
- Exact/allele-complement genotype-alias audit.
- Alias-aware genotype-based LD clumping.
- SNP effects, standard errors, and 95% confidence intervals.
- Formal detectable-effect analysis.
- Unique synthetic marker IDs and LF-only chromosome lists for official GEMMA.

### Reproducibility

- Versioned configuration and deterministic seed.
- Input hashes, sample-order manifests, BED read-back, kinship audits, and final
  machine-readable manifest.
- Raw data, results, caches, and manuscript files remain excluded from GitHub.
