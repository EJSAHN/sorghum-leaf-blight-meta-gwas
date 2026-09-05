# Reproducibility

## Phenotype processing order

The workflow follows this sequence:

1. parse the five field workbooks;
2. preserve plot, replicate, block, site, country, and year metadata;
3. distinguish missing tokens from numeric zero;
4. apply documented corrections from the leaf-blight phenotype workbook;
5. fit field-stratum-adjusted accession effects on the raw scale; and
6. apply rank-based inverse-normal transformation to the accession effects.

## Severity branches

Three severity phenotypes are generated:

- `SEV_ADJ_INT`: standard adjusted severity;
- `SEV_COND_ADJ_INT`: severity adjusted for plot-level incidence; and
- `SEV_POS_ADJ_INT`: severity estimated only from plots with incidence greater than zero.

## Association hierarchy

- Score test: primary inference.
- Likelihood-ratio test: confirmatory sensitivity.
- Wald test: diagnostic.

The hierarchy is defined in `config/final_config.json` and recorded in the analysis manifest.

## Multiplicity and region summaries

BH-FDR is computed across all post-QC SNP rows. A separate sensitivity groups exact and allele-complement genotype patterns into chromosome-specific unique LOCO test classes. Exploratory region summaries are generated after global genotype-equivalence grouping and genotype-based LD clumping. A region is called genome-wide significant only when it satisfies the prespecified FDR criterion.

## Versioned preservation

Each software version is identified by a Git tag and preserved as a Zenodo software record. The source tag and version-specific DOI provide an immutable reference for the exact software distribution.
