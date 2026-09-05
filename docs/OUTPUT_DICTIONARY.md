# Output dictionary

## Result workbook

`05_results/LeafBlight_MultiEnvironment_Analysis_Results.xlsx` contains compact result tables. Full all-SNP tables remain as compressed CSV files under `06_tables`.

### Principal sheets and CSV files

- `Inference_Profile`: genomic inflation factor, minimum p/q, and FDR-hit counts for score, likelihood-ratio, and Wald tests. Wald rows are diagnostic.
- `Test_Sensitivity`: genome-wide agreement and extreme-tail differences among score, LRT, and Wald tests.
- `Score_Concordance`: Pearson and Spearman concordance of incidence and severity score-test association profiles.
- `Score_Enrichment`: overlap and enrichment among top-ranked incidence and severity SNPs.
- `ACAT_vs_Simes`: agreement between the two dual-trait combination methods.
- `Top_Ranked_Candidates`: exploratory genotype-equivalence-aware, LD-clumped score-test regions with trait-specific p/q values, effects, standard errors, and confidence intervals.
- `Candidate_Genes`: genes within 50 kb of each exploratory lead, or the nearest gene if none falls within that interval.
- `Power`: detectable standardized effects by sample size, MAF, power, and alpha.
- `Equivalence_Summary`: counts of SNP rows, chromosome-specific unique LOCO tests, and global exact or allele-complement genotype classes.
- `ACAT_Clump_Leads` and `SEV_Clump_Leads`: exploratory score-test clumping outputs.

### Additional sensitivity files

- `Model_Sensitivity.csv`: K-only versus K+PC3 and n=100 versus n=102 score-test comparisons.
- `Leave_One_Location_Out.csv`: genome-wide score-test comparisons after omitting each field location.
- `Leave_One_Location_Out_Candidates.csv`: candidate-level leave-one-location-out results.
- `Unique_Test_Classes_full.csv.gz`: full chromosome-specific unique-test table.
- `Global_Alias_Classes_full.csv.gz`: full global genotype-equivalence table.
- `Duplicate_Marker_Mapping_full.csv.gz`: SNP-to-equivalence-class mapping.

## Machine-readable manifest

`05_results/analysis_manifest.json` records the inferential hierarchy, analysis dimensions, principal result counts, region-selection rule, and interpretation note.
