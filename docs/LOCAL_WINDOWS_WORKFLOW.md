# Local Windows workflow

The local helper defaults to:

```text
Source analysis folder: C:\projects\leaf blight
Field-workbook search:  D:\2.Research\3. Prom project\Prom-African
New final project:      D:\projects\leaf_blight_PEIR1_final_v2
```

The source folder is treated as read-only. The final project is created on `D:`
and contains all environments, caches, logs, results, manuscript work, and
release materials.

WSL must be installed and able to run `sh`. The helper downloads the official
GEMMA 0.98.5 Linux static binary, verifies the published MD5 checksum, and
stores it under `01_tools` in the new project.

The analysis is resumable. Re-running without `-Force` reuses genotype,
structure, PLINK, GEMMA, and final-table caches when present.

## Installing the final code without running the analysis

From an extracted release package:

```powershell
powershell.exe -NoProfile -ExecutionPolicy Bypass `
  -File .\powershell\INSTALL_LOCAL_CODE_v2.0.0.ps1
```

The default destination is:

```text
D:\projects\sorghum-leaf-blight-meta-gwas_v2.0.0_local
```

## Publishing the audited source to GitHub

`powershell/INSTALL_AND_PUBLISH_v2.0.0.ps1` first runs the complete audited local analysis and unit tests, then calls `powershell/PUBLISH_GITHUB_v2.0.0.ps1` for guarded publication. GitHub is not modified if the full analysis or reference audit fails.

`powershell/PUBLISH_GITHUB_v2.0.0.ps1` performs the publication stage:

1. runs Python compilation and unit tests;
2. clones the current public repository;
3. verifies the pre-update `main` commit;
4. pushes an annotated backup tag for the exploratory workflow;
5. replaces the worktree with the audited v2.0.0 source;
6. refuses raw/derived study data and files larger than 25 MB;
7. commits and pushes `main`; and
8. pushes the annotated `v2.0.0` tag.

GitHub authentication is handled by the user's local Git credential manager.
The script never embeds or stores a token.
