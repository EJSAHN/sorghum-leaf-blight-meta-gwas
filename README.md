# Sorghum leaf blight multi-environment genomic analysis

Version **2.0.0** is the reproducible analysis workflow supporting the revised
sorghum leaf blight study based on West African field evaluations in Niger and
Senegal.

This release supersedes the earlier exploratory repository workflow. It
reconstructs accession phenotypes from the original plot-level observations
before transformation, preserves field/replicate structure, uses a prespecified
primary population of accessions evaluated in at least four locations, and uses
official GEMMA 0.98.5 mixed-model tests for inference.

## Inferential hierarchy

The code enforces the following hierarchy.

1. **Primary:** official GEMMA 0.98.5 leave-one-chromosome-out (LOCO) **score test**.
2. **Confirmatory sensitivity:** official GEMMA 0.98.5 LOCO **likelihood-ratio test**.
3. **Diagnostic only:** official GEMMA 0.98.5 LOCO **Wald test**.

The Wald results are retained to document small-sample variance-component
sensitivity; they are not used to define reported genome-wide discoveries.
Score-test p-values are integrated across incidence and severity with both ACAT
and Simes procedures. Benjamini-Hochberg FDR is calculated at 10%.

Top-ranked, LD-clumped regions are explicitly labelled **exploratory** unless
they pass the prespecified FDR criterion. The workflow does not perform or claim
Bayesian fine-mapping or credible sets.

## Main changes from the earlier exploratory workflow

- Builds a canonical 1,339-row plot-level leaf-blight dataset from five field
  locations.
- Distinguishes numeric zero from missing values.
- Adjusts raw incidence and severity for field stratum before constructing
  accession-level phenotypes.
- Applies rank-based inverse-normal transformation only after accession-level
  estimation.
- Uses 100 accessions represented in at least four locations as the primary
  population; all 102 genotype-matched accessions are retained as a sensitivity
  analysis.
- Computes standard severity, incidence-adjusted severity, and
  disease-positive-only severity.
- Builds LD-pruned, variance-standardized genomic structure markers and
  chromosome-specific LOCO kinship matrices.
- Runs official GEMMA score, LRT, and Wald tests with guaranteed-unique marker
  IDs and Linux-compatible SNP lists.
- Audits exact and allele-complement genotype aliases before candidate-region
  counting.
- Uses genotype-based LD clumping rather than fixed physical bins.
- Compares ACAT with Simes and reports effect estimates, standard errors, 95%
  confidence intervals, formal detectable-effect calculations, n=102
  sensitivity, K+PC3 sensitivity, and leave-one-location-out stability.

## Required inputs

Place the following files under a project directory. Raw inputs are deliberately
excluded from GitHub.

```text
PROJECT_ROOT/
  01_inputs/
    700k.vcf
    Phenotype.xlsx
    Sbicolor_454_v3.1.1.gene.gff3.gz
    raw_field/
      Maradi_Field_Niger_2022_MAY_2023.xlsx
      Bengou_Field_Niger_2022_MAY_18_2023.xlsx
      Field_data_all_locations_SEN_2022_MAY_18_2023.xlsx
```

`Phenotype.xlsx` is the corrected leaf-blight-only six-column workbook used in
the revision audit. The pipeline explicitly extracts only leaf blight incidence
and severity; values for other diseases are not used.

For the audited VCF, the expected SHA-256 is:

```text
8AE866DDBBF729A08EB53E577A94640C7968DEB5C5CE028780EA3C5DF647B723
```

The study-specific expected dimensions are encoded in
`config/final_config.json` and checked before analysis:

```text
1,339 plot-level observations
104 phenotype IDs
102 genotype-matched accessions
100 primary accessions
567,758 post-QC SNPs
```

## Windows one-command workflow

The intended local environment is Windows with a `D:` drive and WSL installed.
The local release package contains `SETUP_AND_RUN_LOCAL.ps1`, which:

- creates a clean project on `D:`;
- locates and copies or hard-links the audited inputs;
- creates a Python 3.11 environment on `D:`;
- downloads and verifies official GEMMA 0.98.5 inside WSL; and
- runs preparation, GEMMA, and final consolidation stages.

Generic direct invocation after the inputs and environment are prepared:

```powershell
python src\run_final_pipeline.py `
  --project-root "D:\projects\leaf_blight_PEIR1_final_v2" `
  --stage all `
  --gemma-wsl "/mnt/d/projects/leaf_blight_PEIR1_final_v2/01_tools/gemma-0.98.5"
```

Stages can be resumed separately:

```powershell
python src\run_final_pipeline.py --project-root "D:\path\to\project" --stage prepare
python src\run_final_pipeline.py --project-root "D:\path\to\project" --stage gemma --gemma-wsl "/mnt/d/path/to/gemma"
python src\run_final_pipeline.py --project-root "D:\path\to\project" --stage finalize
```

## Principal outputs

```text
05_results/
  gemma_primary_n100/
  gemma_primary_n100_pc3/
  gemma_full_n102/
  gemma_lolo_score/
  LeafBlight_PEIR1_Final_Analysis_Results.xlsx
  final_analysis_manifest.json
  READ_ME_FIRST_FINAL_ANALYSIS.txt

06_tables/
  Inference_Profile.csv
  Test_Sensitivity.csv
  Score_Concordance.csv
  Score_Enrichment.csv
  ACAT_vs_Simes.csv
  Top_Ranked_Candidates.csv
  Candidate_Genes.csv
  Power.csv
  Final_Model_Sensitivity.csv
  Final_Leave_One_Location_Out.csv
  Unique_Test_Classes_full.csv.gz
  Global_Alias_Classes_full.csv.gz
```

Large all-SNP outputs and binary caches remain local and are not committed to
GitHub.

## Reproducibility safeguards

- Input SHA-256 validation.
- Deterministic seed `4261991`.
- Explicit sample-order manifests.
- PLINK BED read-back checks.
- Synthetic unique GEMMA marker IDs mapped back to original SNP identifiers.
- LF-only chromosome SNP lists for WSL/Linux.
- Global and chromosome-specific kinship symmetry/eigenvalue audits.
- Exact/complement genotype-equivalence audit.
- Separate primary, confirmatory, and diagnostic test labels.
- Machine-readable final analysis manifest.

## Installation without the PowerShell helper

```bash
python -m venv .venv
.venv/Scripts/python -m pip install --upgrade pip
.venv/Scripts/python -m pip install -r requirements.txt
```

On Linux/macOS, use the platform-equivalent virtual-environment activation and
provide a native GEMMA 0.98.5 executable instead of the WSL path.

## Tests

```bash
python -m pytest -q
```

The test suite covers phenotype transformation, ACAT/Simes combination,
PLINK BED coding, genotype-equivalence classes, detectable-effect calculations,
and compact result consolidation. External GEMMA is not invoked in CI.

## Data availability

The code is public. The raw phenotype and genotype files are not redistributed
here. Genotype data should be obtained from the cited public data source, and
processed result tables should accompany the manuscript as Supplementary Data.

## License

MIT License. See `LICENSE`.

## Guarded GitHub publication helper

The release package includes `powershell/PUBLISH_GITHUB_v2.0.0.ps1`. It runs
compilation and unit tests, verifies the audited pre-update `main` SHA, creates a
backup tag for the exploratory repository state, refuses study-data files,
updates `main`, and pushes the annotated `v2.0.0` tag using the local Git
credential manager. No credentials are stored in the repository.

For the audited Windows system, `powershell/INSTALL_AND_PUBLISH_v2.0.0.ps1`
provides the safest release path: it installs the exact source on `D:`, runs the
complete study analysis and strict reference-reproduction audit, runs the unit
tests, and only then updates GitHub. A failed local analysis leaves GitHub
unchanged.
