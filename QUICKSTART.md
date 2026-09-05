# Quick start

## Windows

```powershell
powershell.exe -NoProfile -ExecutionPolicy Bypass `
  -File .\SETUP_AND_RUN_LOCAL.ps1 `
  -ProjectRoot "<PROJECT_ROOT>" `
  -InputRoot "<INPUT_ROOT>"
```

`<INPUT_ROOT>` must contain the files listed in `README.md`. The helper creates an isolated Python environment inside `<PROJECT_ROOT>` and prepares GEMMA 0.98.5 through WSL.

## Python

```bash
python -m venv .venv
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
python src/run_final_pipeline.py --project-root <PROJECT_ROOT> --stage prepare
python src/run_final_pipeline.py --project-root <PROJECT_ROOT> --stage gemma --gemma-wsl <PATH_TO_GEMMA>
python src/run_final_pipeline.py --project-root <PROJECT_ROOT> --stage finalize
```

The generated summary, workbook, and machine-readable manifest are written under `05_results`.
