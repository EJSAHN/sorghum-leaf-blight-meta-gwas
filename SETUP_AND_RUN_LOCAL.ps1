[CmdletBinding()]
param(
    [ValidateSet("all", "prepare", "gemma", "finalize")]
    [string]$Stage = "all",
    [Parameter(Mandatory = $true)]
    [string]$ProjectRoot,
    [string]$InputRoot = "",
    [string]$GemmaWslPath = "",
    [switch]$Force
)

$ErrorActionPreference = "Stop"
Set-StrictMode -Version Latest

function Banner([string]$Message) {
    Write-Host ""
    Write-Host ("=" * 78) -ForegroundColor Cyan
    Write-Host $Message -ForegroundColor Cyan
    Write-Host ("=" * 78) -ForegroundColor Cyan
}

function Copy-Or-Link {
    param([string]$Source, [string]$Destination)
    if (-not (Test-Path -LiteralPath $Source -PathType Leaf)) {
        throw "Source file not found: $Source"
    }
    $parent = Split-Path -Parent $Destination
    New-Item -ItemType Directory -Path $parent -Force | Out-Null
    if (Test-Path -LiteralPath $Destination -PathType Leaf) {
        $src = Get-Item -LiteralPath $Source
        $dst = Get-Item -LiteralPath $Destination
        if ($src.Length -eq $dst.Length) {
            Write-Host "Reusing existing file: $Destination"
            return
        }
        Remove-Item -LiteralPath $Destination -Force
    }
    if ([System.IO.Path]::GetPathRoot($Source) -eq [System.IO.Path]::GetPathRoot($Destination)) {
        try {
            New-Item -ItemType HardLink -Path $Destination -Target $Source -ErrorAction Stop | Out-Null
            Write-Host "Hard-linked: $Source -> $Destination"
            return
        }
        catch {
            Write-Host "Hard link unavailable; copying: $Source"
        }
    }
    Copy-Item -LiteralPath $Source -Destination $Destination -Force
}

function Get-PythonBase {
    if (Get-Command py.exe -ErrorAction SilentlyContinue) {
        & py.exe -3.11 -c "import sys; print(sys.executable)" *> $null
        if ($LASTEXITCODE -eq 0) { return @{Exe="py.exe"; Args=@("-3.11")} }
        & py.exe -3 -c "import sys; print(sys.executable)" *> $null
        if ($LASTEXITCODE -eq 0) { return @{Exe="py.exe"; Args=@("-3")} }
    }
    if (Get-Command python.exe -ErrorAction SilentlyContinue) {
        & python.exe -c "import sys; print(sys.executable)" *> $null
        if ($LASTEXITCODE -eq 0) { return @{Exe="python.exe"; Args=@()} }
    }
    throw "Python 3.11 or later was not found."
}

function To-WSLPath([string]$WindowsPath) {
    $full = [System.IO.Path]::GetFullPath($WindowsPath)
    if ($full.Length -lt 3 -or $full[1] -ne ':') {
        throw "A Windows drive path is required: $WindowsPath"
    }
    $drive = $full.Substring(0, 1).ToLower()
    $rest = $full.Substring(3).Replace("\", "/")
    return "/mnt/$drive/$rest"
}

Banner "Sorghum leaf blight analysis v2.0.1 - stage: $Stage"

$ProjectRoot = [System.IO.Path]::GetFullPath($ProjectRoot)
$dirs = @(
    $ProjectRoot,
    (Join-Path $ProjectRoot "01_inputs"),
    (Join-Path $ProjectRoot "01_inputs\raw_field"),
    (Join-Path $ProjectRoot "01_tools"),
    (Join-Path $ProjectRoot "02_pipeline"),
    (Join-Path $ProjectRoot "03_config"),
    (Join-Path $ProjectRoot "04_intermediate"),
    (Join-Path $ProjectRoot "05_results"),
    (Join-Path $ProjectRoot "06_tables"),
    (Join-Path $ProjectRoot "07_figures"),
    (Join-Path $ProjectRoot "08_logs"),
    (Join-Path $ProjectRoot "_temp"),
    (Join-Path $ProjectRoot "_pip_cache"),
    (Join-Path $ProjectRoot "_mplconfig")
)
foreach ($dir in $dirs) { New-Item -ItemType Directory -Path $dir -Force | Out-Null }

$pipelineDest = Join-Path $ProjectRoot "02_pipeline"
Banner "Installing the source into the analysis project"
Get-ChildItem -LiteralPath $pipelineDest -Force -ErrorAction SilentlyContinue |
    Remove-Item -Recurse -Force -ErrorAction SilentlyContinue
Get-ChildItem -LiteralPath $PSScriptRoot -Force |
    Where-Object { $_.Name -ne ".git" } |
    ForEach-Object {
        Copy-Item -LiteralPath $_.FullName -Destination $pipelineDest -Recurse -Force
    }
Copy-Item -LiteralPath (Join-Path $pipelineDest "config\final_config.json") `
    -Destination (Join-Path $ProjectRoot "03_config\final_config.json") -Force

if ($Stage -in @("all", "prepare")) {
    $inputs = Join-Path $ProjectRoot "01_inputs"
    if ($InputRoot) {
        $InputRoot = [System.IO.Path]::GetFullPath($InputRoot)
        Banner "Copying analysis inputs"
        Copy-Or-Link -Source (Join-Path $InputRoot "700k.vcf") `
            -Destination (Join-Path $inputs "700k.vcf")
        Copy-Or-Link -Source (Join-Path $InputRoot "Phenotype.xlsx") `
            -Destination (Join-Path $inputs "Phenotype.xlsx")
        Copy-Or-Link -Source (Join-Path $InputRoot "Sbicolor_454_v3.1.1.gene.gff3.gz") `
            -Destination (Join-Path $inputs "Sbicolor_454_v3.1.1.gene.gff3.gz")
        foreach ($name in @(
            "Maradi_Field_Niger_2022_MAY_2023.xlsx",
            "Bengou_Field_Niger_2022_MAY_18_2023.xlsx",
            "Field_data_all_locations_SEN_2022_MAY_18_2023.xlsx"
        )) {
            Copy-Or-Link -Source (Join-Path $InputRoot ("raw_field\" + $name)) `
                -Destination (Join-Path $inputs ("raw_field\" + $name))
        }
    }

    foreach ($required in @(
        "700k.vcf",
        "Phenotype.xlsx",
        "Sbicolor_454_v3.1.1.gene.gff3.gz",
        "raw_field\Maradi_Field_Niger_2022_MAY_2023.xlsx",
        "raw_field\Bengou_Field_Niger_2022_MAY_18_2023.xlsx",
        "raw_field\Field_data_all_locations_SEN_2022_MAY_18_2023.xlsx"
    )) {
        $path = Join-Path $inputs $required
        if (-not (Test-Path -LiteralPath $path -PathType Leaf)) {
            throw "Required input not found: $path"
        }
    }

    Banner "Validating input checksums"
    $hashRows = foreach ($file in Get-ChildItem -LiteralPath $inputs -Recurse -File) {
        $hash = Get-FileHash -LiteralPath $file.FullName -Algorithm SHA256
        [pscustomobject]@{
            Name = $file.Name
            RelativePath = $file.FullName.Substring($inputs.Length).TrimStart([char]92)
            SizeBytes = $file.Length
            SHA256 = $hash.Hash
        }
    }
    $hashRows | Export-Csv -LiteralPath (Join-Path $inputs "input_hashes.csv") -NoTypeInformation -Encoding UTF8
    $configData = Get-Content -LiteralPath (Join-Path $ProjectRoot "03_config\final_config.json") -Raw | ConvertFrom-Json
    $observed = ($hashRows | Where-Object Name -eq "700k.vcf" | Select-Object -First 1).SHA256
    $expected = [string]$configData.expected_vcf_sha256
    if ($expected -and $observed -ne $expected.ToUpper()) {
        throw "VCF SHA-256 mismatch: $observed"
    }
    Write-Host "Input checksum validation passed." -ForegroundColor Green
}

Banner "Creating or reusing the Python environment"
$env:TEMP = Join-Path $ProjectRoot "_temp"
$env:TMP = $env:TEMP
$env:PIP_CACHE_DIR = Join-Path $ProjectRoot "_pip_cache"
$env:MPLCONFIGDIR = Join-Path $ProjectRoot "_mplconfig"
$env:PYTHONUTF8 = "1"
$venv = Join-Path $ProjectRoot ".venv"
$python = Join-Path $venv "Scripts\python.exe"
if (-not (Test-Path -LiteralPath $python -PathType Leaf)) {
    $base = Get-PythonBase
    & $base.Exe @($base.Args) -m venv $venv
    if ($LASTEXITCODE -ne 0) { throw "Failed to create the Python environment." }
}
& $python -m pip install --upgrade pip
if ($LASTEXITCODE -ne 0) { throw "pip upgrade failed." }
& $python -m pip install -r (Join-Path $pipelineDest "requirements.txt")
if ($LASTEXITCODE -ne 0) { throw "Dependency installation failed." }

$gemmaWin = Join-Path $ProjectRoot "01_tools\gemma-0.98.5"
$gemmaGz = Join-Path $ProjectRoot "01_tools\gemma-0.98.5-linux-static-AMD64.gz"
if ($GemmaWslPath) {
    $gemmaWsl = $GemmaWslPath
} else {
    $gemmaWsl = To-WSLPath $gemmaWin
}

if ($Stage -in @("all", "gemma")) {
    Banner "Preparing official GEMMA 0.98.5"
    if (-not (Get-Command wsl.exe -ErrorAction SilentlyContinue)) {
        throw "WSL is required for this helper."
    }
    wsl.exe -e sh -lc "printf WSL_OK" *> $null
    if ($LASTEXITCODE -ne 0) { throw "No ready Linux distribution was found in WSL." }

    if (-not $GemmaWslPath -and -not (Test-Path -LiteralPath $gemmaWin -PathType Leaf)) {
        $url = "https://github.com/genetics-statistics/GEMMA/releases/download/v0.98.5/gemma-0.98.5-linux-static-AMD64.gz"
        Invoke-WebRequest -Uri $url -OutFile $gemmaGz -UseBasicParsing
        $gzWsl = To-WSLPath $gemmaGz
        wsl.exe -e sh -lc "gzip -dc '$gzWsl' > '$gemmaWsl' && chmod u+x '$gemmaWsl'"
        if ($LASTEXITCODE -ne 0) { throw "Failed to prepare GEMMA in WSL." }
    }

    $md5 = (wsl.exe -e sh -lc "md5sum '$gemmaWsl' | awk '{print `$1}'").Trim().ToLower()
    $expectedMd5 = "f5e90535ff6a36867dcb5f6b0fb24135"
    if ($md5 -ne $expectedMd5) { throw "GEMMA MD5 mismatch: $md5" }
    Write-Host "GEMMA 0.98.5 checksum validation passed." -ForegroundColor Green
}

Banner "Running stage: $Stage"
$runner = Join-Path $pipelineDest "src\run_final_pipeline.py"
$argsList = @(
    $runner,
    "--project-root", $ProjectRoot,
    "--stage", $Stage,
    "--config", (Join-Path $ProjectRoot "03_config\final_config.json")
)
if ($Stage -in @("all", "gemma")) { $argsList += @("--gemma-wsl", $gemmaWsl) }
if ($Force) { $argsList += "--force" }
$stamp = Get-Date -Format "yyyyMMdd_HHmmss"
$psLog = Join-Path $ProjectRoot "08_logs\powershell_run_$stamp.log"
$oldPreference = $ErrorActionPreference
$ErrorActionPreference = "Continue"
& $python @argsList 2>&1 | Tee-Object -FilePath $psLog
$exitCode = $LASTEXITCODE
$ErrorActionPreference = $oldPreference
if ($exitCode -ne 0) {
    throw "The analysis exited with code $exitCode. Complete output: $psLog"
}

Banner "ANALYSIS STAGE COMPLETE"
Write-Host "Project: $ProjectRoot" -ForegroundColor Yellow
$summaryPath = Join-Path $ProjectRoot "05_results\ANALYSIS_SUMMARY.txt"
$workbookPath = Join-Path $ProjectRoot "05_results\LeafBlight_MultiEnvironment_Analysis_Results.xlsx"
if (Test-Path -LiteralPath $summaryPath) {
    Write-Host "Summary: $summaryPath" -ForegroundColor Yellow
    Write-Host "Workbook: $workbookPath" -ForegroundColor Yellow
}
Write-Host "Log: $psLog"
