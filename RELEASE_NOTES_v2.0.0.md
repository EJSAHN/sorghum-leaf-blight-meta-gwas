# Release v2.0.0 — revised multi-environment analysis

This release supersedes the exploratory workflow that accompanied the initial
submission. It rebuilds leaf-blight phenotypes from the original plot-level
records and implements the inferential hierarchy adopted during revision.

## Inferential hierarchy

- **Primary:** official GEMMA 0.98.5 LOCO score test.
- **Confirmatory sensitivity:** official GEMMA 0.98.5 LOCO likelihood-ratio test.
- **Diagnostic only:** official GEMMA 0.98.5 LOCO Wald test.

The Wald tail is retained for transparency but is not used to define reported
discoveries. The revised workflow does not restore the submitted 37-locus
result and does not claim credible sets or formal fine-mapping.

## Key corrections and additions

- Reconstructs 1,339 plot-level observations from five West African field sites.
- Models field-stratum effects before accession-level transformation.
- Uses 100 accessions represented at four or more locations as the primary
  analysis population and 102 genotype-matched accessions as sensitivity.
- Includes standard, incidence-adjusted, and disease-positive-only severity.
- Uses LD-pruned structure markers and chromosome-specific LOCO kinship.
- Runs official GEMMA score, LRT, and Wald tests with unique synthetic marker IDs.
- Compares ACAT with Simes.
- Audits exact and allele-complement genotype aliases.
- Uses genotype-based LD clumping for exploratory candidate-region summaries.
- Reports effects, standard errors, 95% confidence intervals, model sensitivity,
  leave-one-location-out stability, and detectable-effect calculations.

## Data and local execution

Raw VCF and phenotype files are not distributed in this repository. The
Windows helper `SETUP_AND_RUN_LOCAL.ps1` locates the audited local inputs,
verifies the VCF SHA-256, prepares a project on `D:`, installs dependencies,
downloads and verifies official GEMMA 0.98.5 in WSL, and runs the workflow.
