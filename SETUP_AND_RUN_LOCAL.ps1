[CmdletBinding()]
param(
    [ValidateSet("all", "prepare", "gemma", "finalize")]
    [string]$Stage = "all",
    [string]$ProjectRoot = "D:\projects\leaf_blight_PEIR1_final_v2",
    [string]$LegacyRoot = "C:\projects\leaf blight",
    [string]$PromRoot = "D:\2.Research\3. Prom project\Prom-African",
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

function Find-FirstFile {
    param(
        [string]$ExactName,
        [string[]]$SearchRoots,
        [string]$Pattern = ""
    )
    foreach ($root in $SearchRoots) {
        if (-not (Test-Path -LiteralPath $root)) { continue }
        $direct = Join-Path $root $ExactName
        if (Test-Path -LiteralPath $direct -PathType Leaf) {
            return (Get-Item -LiteralPath $direct).FullName
        }
    }
    foreach ($root in $SearchRoots) {
        if (-not (Test-Path -LiteralPath $root)) { continue }
        $hits = @(Get-ChildItem -LiteralPath $root -Recurse -File -Filter $ExactName -ErrorAction SilentlyContinue)
        if ($hits.Count -gt 0) {
            return ($hits | Sort-Object LastWriteTime -Descending | Select-Object -First 1).FullName
        }
    }
    if ($Pattern) {
        foreach ($root in $SearchRoots) {
            if (-not (Test-Path -LiteralPath $root)) { continue }
            $hits = @(Get-ChildItem -LiteralPath $root -Recurse -File -Filter $Pattern -ErrorAction SilentlyContinue)
            if ($hits.Count -gt 0) {
                return ($hits | Sort-Object LastWriteTime -Descending | Select-Object -First 1).FullName
            }
        }
    }
    return $null
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
    throw "Python 3.11 or the Windows py launcher was not found."
}

function To-WSLPath([string]$WindowsPath) {
    $drive = $WindowsPath.Substring(0,1).ToLower()
    $rest = $WindowsPath.Substring(3).Replace("\", "/")
    return "/mnt/$drive/$rest"
}

Banner "Sorghum leaf blight final revision pipeline v2.0.0 - stage: $Stage"

if (-not (Test-Path -LiteralPath "D:\")) { throw "D: drive was not found." }
$drive = New-Object System.IO.DriveInfo("D")
$freeGB = [math]::Round($drive.AvailableFreeSpace / 1GB, 2)
Write-Host "D: free space: $freeGB GB"
if ($freeGB -lt 8) { throw "At least 8 GB free space on D: is required." }

$dirs = @(
    $ProjectRoot,
    (Join-Path $ProjectRoot "00_submitted_snapshot"),
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
    (Join-Path $ProjectRoot "09_manuscript"),
    (Join-Path $ProjectRoot "10_response_letter"),
    (Join-Path $ProjectRoot "11_release"),
    (Join-Path $ProjectRoot "_temp"),
    (Join-Path $ProjectRoot "_pip_cache"),
    (Join-Path $ProjectRoot "_mplconfig")
)
foreach ($dir in $dirs) { New-Item -ItemType Directory -Path $dir -Force | Out-Null }

Banner "Installing versioned pipeline into the D: project"
$pipelineDest = Join-Path $ProjectRoot "02_pipeline"
Get-ChildItem -LiteralPath $pipelineDest -Force -ErrorAction SilentlyContinue |
    Remove-Item -Recurse -Force -ErrorAction SilentlyContinue
Copy-Item -Path (Join-Path $PSScriptRoot "*") -Destination $pipelineDest -Recurse -Force
Copy-Item -LiteralPath (Join-Path $pipelineDest "config\final_config.json") `
    -Destination (Join-Path $ProjectRoot "03_config\final_config.json") -Force

if ($Stage -in @("all", "prepare")) {
    Banner "Locating and copying audited inputs"
    $searchRoots = @(
        $LegacyRoot,
        $PromRoot,
        "$env:USERPROFILE\Downloads",
        "D:\projects\leaf_blight_PEIR1_4261991_v02_model_lock\01_inputs"
    )
    $vcf = Find-FirstFile -ExactName "700k.vcf" -SearchRoots $searchRoots
    $phenotype = Find-FirstFile -ExactName "Phenotype.xlsx" -SearchRoots $searchRoots
    $gff = Find-FirstFile -ExactName "Sbicolor_454_v3.1.1.gene.gff3.gz" -SearchRoots $searchRoots
    if (-not $vcf) { throw "700k.vcf was not found." }
    if (-not $phenotype) { throw "Phenotype.xlsx was not found." }
    if (-not $gff) { throw "Sbicolor_454_v3.1.1.gene.gff3.gz was not found." }

    $inputs = Join-Path $ProjectRoot "01_inputs"
    Copy-Or-Link -Source $vcf -Destination (Join-Path $inputs "700k.vcf")
    Copy-Or-Link -Source $phenotype -Destination (Join-Path $inputs "Phenotype.xlsx")
    Copy-Or-Link -Source $gff -Destination (Join-Path $inputs "Sbicolor_454_v3.1.1.gene.gff3.gz")

    $rawSpecs = @(
        @{Name="Maradi_Field_Niger_2022_MAY_2023.xlsx"; Pattern="*Maradi*Field*Niger*2022*.xlsx"},
        @{Name="Bengou_Field_Niger_2022_MAY_18_2023.xlsx"; Pattern="*Bengou*Field*Niger*2022*.xlsx"},
        @{Name="Field_data_all_locations_SEN_2022_MAY_18_2023.xlsx"; Pattern="*Field*data*all*locations*SEN*2022*.xlsx"}
    )
    foreach ($spec in $rawSpecs) {
        $source = Find-FirstFile -ExactName $spec.Name -SearchRoots $searchRoots -Pattern $spec.Pattern
        if (-not $source) { throw "Original field workbook not found: $($spec.Name)" }
        Copy-Or-Link -Source $source -Destination (Join-Path $inputs ("raw_field\" + $spec.Name))
    }

    Banner "Hashing inputs"
    $hashRows = foreach ($file in Get-ChildItem -LiteralPath $inputs -Recurse -File) {
        $hash = Get-FileHash -LiteralPath $file.FullName -Algorithm SHA256
        [pscustomobject]@{
            Name=$file.Name
            FullName=$file.FullName
            SizeBytes=$file.Length
            LastWriteTime=$file.LastWriteTime
            SHA256=$hash.Hash
        }
    }
    $hashRows | Export-Csv -LiteralPath (Join-Path $inputs "input_hashes.csv") -NoTypeInformation -Encoding UTF8
    $observed = ($hashRows | Where-Object Name -eq "700k.vcf" | Select-Object -First 1).SHA256
    $expected = "8AE866DDBBF729A08EB53E577A94640C7968DEB5C5CE028780EA3C5DF647B723"
    if ($observed -ne $expected) { throw "VCF SHA-256 mismatch: $observed" }
    Write-Host "VCF SHA-256 matches the audited baseline." -ForegroundColor Green
}

Banner "Creating or reusing the Python environment on D:"
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
    if ($LASTEXITCODE -ne 0) { throw "Failed to create Python environment." }
}
& $python -m pip install --upgrade pip
if ($LASTEXITCODE -ne 0) { throw "pip upgrade failed." }
& $python -m pip install -r (Join-Path $pipelineDest "requirements.txt")
if ($LASTEXITCODE -ne 0) { throw "Python dependency installation failed." }

$gemmaWin = Join-Path $ProjectRoot "01_tools\gemma-0.98.5"
$gemmaGz = Join-Path $ProjectRoot "01_tools\gemma-0.98.5-linux-static-AMD64.gz"
$gemmaWsl = To-WSLPath $gemmaWin
if ($Stage -in @("all", "gemma")) {
    Banner "Preparing official GEMMA 0.98.5 in WSL"
    if (-not (Get-Command wsl.exe -ErrorAction SilentlyContinue)) {
        throw "WSL is required for official GEMMA 0.98.5."
    }
    wsl.exe -e sh -lc "printf WSL_OK" *> $null
    if ($LASTEXITCODE -ne 0) { throw "WSL is installed but no Linux distribution is ready." }
    if (-not (Test-Path -LiteralPath $gemmaWin -PathType Leaf)) {
        $url = "https://github.com/genetics-statistics/GEMMA/releases/download/v0.98.5/gemma-0.98.5-linux-static-AMD64.gz"
        Write-Host "Downloading official GEMMA 0.98.5..."
        Invoke-WebRequest -Uri $url -OutFile $gemmaGz -UseBasicParsing
        $gzWsl = To-WSLPath $gemmaGz
        wsl.exe -e sh -lc "gzip -dc '$gzWsl' > '$gemmaWsl' && chmod u+x '$gemmaWsl'"
        if ($LASTEXITCODE -ne 0) { throw "Failed to decompress GEMMA in WSL." }
    }
    $md5 = (wsl.exe -e sh -lc "md5sum '$gemmaWsl' | awk '{print `$1}'").Trim().ToLower()
    $expectedMd5 = "f5e90535ff6a36867dcb5f6b0fb24135"
    if ($md5 -ne $expectedMd5) { throw "GEMMA MD5 mismatch: $md5" }
    Write-Host "Official GEMMA 0.98.5 verified." -ForegroundColor Green
}

Banner "Running pipeline stage: $Stage"
$runner = Join-Path $pipelineDest "src\run_final_pipeline.py"
$argsList = @($runner, "--project-root", $ProjectRoot, "--stage", $Stage, "--config", (Join-Path $ProjectRoot "03_config\final_config.json"))
if ($Stage -in @("all", "gemma")) { $argsList += @("--gemma-wsl", $gemmaWsl) }
if ($Force) { $argsList += "--force" }
$stamp = Get-Date -Format "yyyyMMdd_HHmmss"
$psLog = Join-Path $ProjectRoot "08_logs\powershell_final_v2_$stamp.log"
$oldPreference = $ErrorActionPreference
$ErrorActionPreference = "Continue"
& $python @argsList 2>&1 | Tee-Object -FilePath $psLog
$exitCode = $LASTEXITCODE
$ErrorActionPreference = $oldPreference
if ($exitCode -ne 0) {
    Write-Host "Pipeline failed. Complete output: $psLog" -ForegroundColor Red
    throw "Pipeline exited with code $exitCode."
}

Banner "FINAL PIPELINE STAGE COMPLETE"
Write-Host "Project: $ProjectRoot" -ForegroundColor Yellow
if (Test-Path -LiteralPath (Join-Path $ProjectRoot "05_results\READ_ME_FIRST_FINAL_ANALYSIS.txt")) {
    Write-Host "Summary: $(Join-Path $ProjectRoot '05_results\READ_ME_FIRST_FINAL_ANALYSIS.txt')" -ForegroundColor Yellow
    Write-Host "Workbook: $(Join-Path $ProjectRoot '05_results\LeafBlight_PEIR1_Final_Analysis_Results.xlsx')" -ForegroundColor Yellow
}
Write-Host "Log: $psLog"
