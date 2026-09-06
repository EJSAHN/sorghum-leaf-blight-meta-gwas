# Sorghum leaf blight multi-environment genomic analysis

Version **2.0.2** provides a reproducible workflow for genomic analysis of sorghum leaf blight incidence and severity measured across field locations in Niger and Senegal.

## Analysis overview

The workflow:

- reconstructs a canonical 1,339-row plot-level dataset from five field locations;
- distinguishes numeric zero from missing observations;
- models field-stratum effects before constructing accession-level phenotypes;
- applies rank-based inverse-normal transformation after accession-level estimation;
- uses 100 accessions represented in at least four locations as the primary population;
- retains all 102 genotype-matched accessions as a sensitivity population;
- evaluates standard severity, incidence-adjusted severity, and disease-positive-only severity;
- builds LD-pruned genomic structure markers and chromosome-specific LOCO kinship matrices;
- combines incidence and severity p-values with ACAT and Simes procedures;
- groups exact and allele-complement genotype patterns before LD-based region summaries; and
- reports effect estimates, standard errors, 95% confidence intervals, detectable-effect calculations, and leave-one-location-out stability.

## Inferential hierarchy

1. **Primary:** official GEMMA 0.98.5 leave-one-chromosome-out (LOCO) score test.
2. **Confirmatory sensitivity:** official GEMMA 0.98.5 LOCO likelihood-ratio test.
3. **Diagnostic:** official GEMMA 0.98.5 LOCO Wald test.

Score-test p-values are combined across incidence and severity with both ACAT and Simes. Benjamini-Hochberg false-discovery-rate control is calculated at 10%. Top-ranked LD-clumped regions are labelled exploratory unless they satisfy the prespecified FDR criterion. The workflow does not perform Bayesian fine-mapping or construct credible sets.

## Required inputs

Arrange the input files as follows:

```text
<INPUT_ROOT>/
  700k.vcf
  Phenotype.xlsx
  Sbicolor_454_v3.1.1.gene.gff3.gz
  raw_field/
    Maradi_Field_Niger_2022_MAY_2023.xlsx
    Bengou_Field_Niger_2022_MAY_18_2023.xlsx
    Field_data_all_locations_SEN_2022_MAY_18_2023.xlsx
```

Only leaf blight incidence and severity are extracted from the field workbooks. Other disease measurements are not used.

The reference VCF used for the reported analysis has SHA-256:

```text
8AE866DDBBF729A08EB53E577A94640C7968DEB5C5CE028780EA3C5DF647B723
```

The expected study dimensions are stored in `config/final_config.json` and checked during execution:

```text
1,339 plot-level observations
104 phenotype identifiers
102 genotype-matched accessions
100 primary accessions
567,758 post-QC SNPs
```

## Windows workflow

Windows execution requires PowerShell, Python 3.11 or later, and WSL with a Linux distribution. From the repository root:

```powershell
powershell.exe -NoProfile -ExecutionPolicy Bypass `
  -File .\SETUP_AND_RUN_LOCAL.ps1 `
  -ProjectRoot "<PROJECT_ROOT>" `
  -InputRoot "<INPUT_ROOT>"
```

The helper creates the project structure, copies or hard-links the inputs, installs Python dependencies in the project directory, downloads and verifies official GEMMA 0.98.5 within WSL, and runs the selected stage. See `docs/WINDOWS_WORKFLOW.md`.

## Direct invocation

After inputs and dependencies are prepared:

```powershell
python src
un_final_pipeline.py `
  --project-root "<PROJECT_ROOT>" `
  --stage all `
  --gemma-wsl "<WSL_PATH_TO_GEMMA>"
```

Stages can be resumed independently:

```powershell
python src
un_final_pipeline.py --project-root "<PROJECT_ROOT>" --stage prepare
python src
un_final_pipeline.py --project-root "<PROJECT_ROOT>" --stage gemma --gemma-wsl "<WSL_PATH_TO_GEMMA>"
python src
un_final_pipeline.py --project-root "<PROJECT_ROOT>" --stage finalize
```

## Principal outputs

```text
05_results/
  gemma_primary_n100/
  gemma_primary_n100_pc3/
  gemma_full_n102/
  gemma_lolo_score/
  LeafBlight_MultiEnvironment_Analysis_Results.xlsx
  analysis_manifest.json
  ANALYSIS_SUMMARY.txt

06_tables/
  Inference_Profile.csv
  Test_Sensitivity.csv
  Score_Concordance.csv
  Score_Enrichment.csv
  ACAT_vs_Simes.csv
  Top_Ranked_Candidates.csv
  Candidate_Genes.csv
  Power.csv
  Model_Sensitivity.csv
  Leave_One_Location_Out.csv
  Leave_One_Location_Out_Candidates.csv
  Unique_Test_Classes_full.csv.gz
  Global_Alias_Classes_full.csv.gz
```

Large all-SNP outputs and binary caches remain in the analysis project and are not committed to the source repository.

## Reproducibility checks

- Input SHA-256 verification.
- Deterministic seed `4261991`.
- Explicit sample-order manifests.
- PLINK BED read-back verification.
- Unique GEMMA marker identifiers mapped back to the source SNP identifiers.
- LF-only chromosome SNP lists for Linux execution.
- Global and chromosome-specific kinship symmetry and eigenvalue checks.
- Exact and allele-complement genotype-equivalence analysis.
- Separate primary, confirmatory, and diagnostic test labels.
- Machine-readable analysis manifest.

## Installation without the PowerShell helper

```bash
python -m venv .venv
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
```

On Linux or macOS, provide a native GEMMA 0.98.5 executable instead of a WSL path.

## Tests

```bash
python -m pytest -q
```

The test suite covers phenotype transformation, ACAT and Simes combination, PLINK BED coding, genotype-equivalence classes, detectable-effect calculations, result consolidation, and public-distribution text checks. External GEMMA is not invoked in continuous integration.

## Data availability

The source code is public. Raw phenotype and genotype files are not redistributed in this repository. Genotype data should be obtained from the cited public data source, and processed result tables are distributed separately with the associated research article.

## Citation

Use the metadata in `CITATION.cff` or the DOI assigned to the corresponding Zenodo release.

## License

MIT License. See `LICENSE`.


## Auxiliary model diagnostics

OLS raw/INT-by-PC diagnostics, PC reconstruction, matched-mode-3 candidate
effect comparisons and residual-SD detectable effects are available through:

```text
python src/run_sensitivity.py --project-root PROJECT --output-dir NEW_OUTPUT --stage all
```

This command reads an existing analysis project and writes to a new, separate
folder. It does not repeat the full mixed-model scans or change study figures.
See [Model diagnostics](docs/DIAGNOSTICS.md) for inputs, effect units, estimator
labels and optional independently supplied reference tables.

A complete reproduction uses the main workflow followed by this auxiliary
command. Study data are supplied separately, not embedded in the source
archive. Source-file hashes can be checked with `python tools/verify_source.py`.
