# Output dictionary

## Primary result workbook

`05_results/LeafBlight_PEIR1_Final_Analysis_Results.xlsx` contains compact,
reviewable tables. Full all-SNP tables remain as compressed CSV files under
`06_tables`.

### Principal sheets and CSV files

- `Inference_Profile`: lambda GC, minimum p/q, and FDR-hit counts for score,
  likelihood-ratio, and Wald tests. Wald rows are diagnostic only.
- `Test_Sensitivity`: genome-wide agreement and extreme-tail disagreement among
  score, LRT, and Wald tests.
- `Score_Concordance`: Pearson and Spearman concordance of incidence and severity
  score-test association profiles.
- `Score_Enrichment`: overlap and enrichment among top-ranked incidence and
  severity SNPs.
- `ACAT_vs_Simes`: agreement between the two dual-trait combination methods.
- `Top_Ranked_Candidates`: exploratory alias-aware, LD-clumped score-test regions
  with trait-specific p/q values, effects, standard errors, and confidence
  intervals.
- `Candidate_Genes`: genes within 50 kb of each exploratory lead, or the nearest
  gene if none falls within that interval.
- `Power`: detectable standardized effects by n, MAF, power, and alpha.
- `Equivalence_Summary`: counts of SNP rows, chromosome-specific unique LOCO
  tests, and global exact/complement genotype aliases.
- `ACAT_Clump_Leads` / `SEV_Clump_Leads`: exploratory score-test clumping output.

### Additional sensitivity files

- `Final_Model_Sensitivity.csv`: K-only primary versus K+PC3 and n=100 versus
  n=102 score-test comparisons.
- `Final_Leave_One_Location_Out.csv`: genome-wide score-test comparisons after
  omitting each field location.
- `Unique_Test_Classes_full.csv.gz`: full chromosome-specific unique-test table.
- `Global_Alias_Classes_full.csv.gz`: full global genotype-alias table.
- `Duplicate_Marker_Mapping_full.csv.gz`: SNP-to-equivalence-class mapping.

## Machine-readable manifest

`05_results/final_analysis_manifest.json` records the inferential hierarchy,
analysis dimensions, principal result counts, candidate-selection rule, and
interpretation guardrail.
