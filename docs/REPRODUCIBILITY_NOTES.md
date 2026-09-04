# Reproducibility and interpretation notes

## Phenotype order of operations

The earlier exploratory code rank-transformed plot rows and then retained one
row per accession. Version 2.0.0 instead performs the following sequence:

1. parse the five original field workbooks;
2. preserve plot, replicate, block, site, country, and year metadata;
3. distinguish missing tokens from numeric zero;
4. apply documented corrections from the corrected leaf-blight workbook;
5. fit field-stratum-adjusted accession effects on the raw scale; and
6. apply rank-based inverse-normal transformation to accession effects.

This order is the basis of the revised analysis.

## Severity sensitivity branches

The workflow creates three severity phenotypes:

- `SEV_ADJ_INT`: standard adjusted severity;
- `SEV_COND_ADJ_INT`: severity adjusted for plot-level incidence; and
- `SEV_POS_ADJ_INT`: severity estimated only from plots with incidence > 0.

Their interpretation should follow the final confirmation of the field scoring
protocol. All three are computed so the numerical analysis does not depend on
that wording decision.

## Association tests

Official GEMMA 0.98.5 is run with supplied, audited global and
chromosome-specific kinship matrices.

- Score test: primary inference.
- Likelihood-ratio test: confirmatory sensitivity.
- Wald test: diagnostic only.

This hierarchy is prespecified in `config/final_config.json` and repeated in the
final result manifest.

## Multiplicity and candidate regions

BH-FDR is computed across all post-QC SNP rows. A second sensitivity collapses
exact/complement genotype patterns to chromosome-specific unique LOCO test
classes. Exploratory candidate regions are generated only after global alias
collapse and genotype-based LD clumping. They must not be described as
significant QTL unless they pass the prespecified FDR criterion.

## Permanent archive

After the revised code and manuscript are frozen, create a GitHub release and a
permanent DOI archive from the exact same tagged commit. Add the commit hash,
release tag, and DOI to the manuscript Data Availability statement.
