# Release v2.0.1

Version 2.0.1 provides a public distribution of the sorghum leaf blight multi-environment genomic analysis workflow.

## Analysis scope

- Reconstructs 1,339 plot-level observations from five West African field locations.
- Models field-stratum effects before accession-level transformation.
- Uses 100 accessions represented in at least four locations as the primary analysis population and 102 genotype-matched accessions as a sensitivity population.
- Includes standard, incidence-adjusted, and disease-positive-only severity phenotypes.
- Uses official GEMMA 0.98.5 LOCO score tests for primary inference.
- Uses GEMMA LOCO likelihood-ratio tests as confirmatory sensitivity analyses.
- Retains GEMMA LOCO Wald tests as diagnostics.
- Compares ACAT and Simes p-value combinations.
- Groups genotype-equivalent markers before genotype-based LD clumping.
- Reports effect estimates, standard errors, 95% confidence intervals, detectable-effect calculations, and leave-one-location-out stability.

## Distribution changes

- Uses generic project paths and neutral public-facing terminology.
- Uses generic project and input-root parameters in the Windows helper.
- Uses neutral result filenames and directory labels.
- Adds explicit Zenodo metadata and updated citation metadata.

Raw genotype and phenotype files are not included in the source archive.
