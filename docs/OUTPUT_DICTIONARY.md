# Output dictionary

## Main analysis

`05_results/LeafBlight_MultiEnvironment_Analysis_Results.xlsx` contains compact
results. Full-marker probability tables remain in the analysis project.

| Output table | Meaning |
|---|---|
| Inference_Profile | Number of tested SNP rows, lambda GC, minimum p/q and FDR counts for each test and trait/combination. |
| Test_Sensitivity | Contrasts of score, likelihood-ratio and Wald probability profiles. |
| Score_Concordance | Incidence/severity association-strength correlations. |
| Score_Enrichment | Overlap, enrichment and set statistics at specified rank fractions. |
| ACAT_vs_Simes | Comparison of the two p-value combination methods. |
| Top_Ranked_Candidates | Exploratory alias-aware LD-clumped regions with separately identified p/q fields. |
| Candidate_Genes | Nearby annotated genes; proximity is not evidence of causality. |
| Power | Noncentral-t detectable effects in residual SD per allele. |
| Model_Sensitivity | K-only/K+PC3 and primary/full-sample probability comparisons. |
| Leave_One_Location_Out | Genome-wide comparisons after individual location omissions. |
| Leave_One_Location_Out_Candidates | Candidate score probabilities and omission-mode-3 coefficients. No comparison to mode-4 primary coefficients is made. |

Main candidate beta/SE values come from GEMMA mode 4 alternative-model REML.
Their beta +/- 1.96 SE intervals are unadjusted Wald intervals, not score-test
confidence intervals. Effect metadata describe the estimator and synthetic
allele1 A coding; A/G do not identify the original nucleotide alleles.

## Auxiliary command

`python src/run_sensitivity.py --project-root PROJECT --output-dir NEW_OUTPUT`
produces the following in a separate output folder.

| Output | Meaning |
|---|---|
| Transformation_Checks.csv | Twenty OLS raw/INT-by-PC diagnostics using the project's primary PC basis. |
| PC_Reconstruction_Columns.csv | Correlation checks for stored versus genotype-reconstructed PCs. |
| PC_Reconstruction_Subspaces.csv | Distances and angles between nested PC spaces. |
| OLS_Reproduction.csv | Optional comparisons with independent expected numerical tables. |
| Candidate_Effect_Direction_Check.csv | Primary and omission coefficients both from GEMMA mode 3, paired with score probabilities. |
| Candidate_Direction_Summary.csv | Counts of evaluated and direction-preserved candidate/omission comparisons. |
| Power.csv | Independently recalculated residual-SD effect thresholds; no total variance explained is inferred. |
| Calculation_Summary.json | Requested calculations, completion status and whether independent reference tables were supplied. |

See `docs/DIAGNOSTICS.md` for required local files and optional reference inputs.
