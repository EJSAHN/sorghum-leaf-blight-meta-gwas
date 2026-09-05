# Windows workflow

## Requirements

- Windows PowerShell 5.1 or PowerShell 7
- Python 3.11 or later
- Windows Subsystem for Linux with a Linux distribution
- Sufficient free disk space for genotype caches and GEMMA outputs

## Input layout

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

## Run all stages

From the repository root:

```powershell
powershell.exe -NoProfile -ExecutionPolicy Bypass `
  -File .\SETUP_AND_RUN_LOCAL.ps1 `
  -ProjectRoot "<PROJECT_ROOT>" `
  -InputRoot "<INPUT_ROOT>"
```

The helper copies the source into the project, prepares an isolated Python environment, verifies the input checksum, downloads and verifies GEMMA 0.98.5 through WSL, and runs the analysis.

## Resume a stage

```powershell
powershell.exe -NoProfile -ExecutionPolicy Bypass `
  -File .\SETUP_AND_RUN_LOCAL.ps1 `
  -ProjectRoot "<PROJECT_ROOT>" `
  -Stage prepare

powershell.exe -NoProfile -ExecutionPolicy Bypass `
  -File .\SETUP_AND_RUN_LOCAL.ps1 `
  -ProjectRoot "<PROJECT_ROOT>" `
  -Stage gemma

powershell.exe -NoProfile -ExecutionPolicy Bypass `
  -File .\SETUP_AND_RUN_LOCAL.ps1 `
  -ProjectRoot "<PROJECT_ROOT>" `
  -Stage finalize
```

For `prepare`, provide `-InputRoot` unless the required files already exist under `<PROJECT_ROOT>/01_inputs`. Existing caches are reused unless `-Force` is supplied.
