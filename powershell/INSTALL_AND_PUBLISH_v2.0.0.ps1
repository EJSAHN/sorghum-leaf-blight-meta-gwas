[CmdletBinding()]
param(
    [string]$LocalDestination = "D:\projects\sorghum-leaf-blight-meta-gwas_v2.0.0_local",
    [string]$AnalysisProject = "D:\projects\leaf_blight_PEIR1_final_v2",
    [switch]$ForceLocal,
    [switch]$SkipFullAnalysis,
    [switch]$SkipTests,
    [switch]$SkipRelease,
    [switch]$ForceRemoteMismatch
)

$ErrorActionPreference = "Stop"
Set-StrictMode -Version Latest

$Install = Join-Path $PSScriptRoot "INSTALL_LOCAL_CODE_v2.0.0.ps1"
$Publish = Join-Path $PSScriptRoot "PUBLISH_GITHUB_v2.0.0.ps1"
if (-not (Test-Path -LiteralPath $Install)) { throw "Installer not found: $Install" }
if (-not (Test-Path -LiteralPath $Publish)) { throw "Publisher not found: $Publish" }

$installArgs = @("-NoProfile", "-ExecutionPolicy", "Bypass", "-File", $Install, "-Destination", $LocalDestination)
if ($ForceLocal) { $installArgs += "-Force" }
& powershell.exe @installArgs
if ($LASTEXITCODE -ne 0) { throw "Local code installation failed." }

$AlreadyTested = $false
if (-not $SkipFullAnalysis) {
    $Runner = Join-Path $LocalDestination "SETUP_AND_RUN_LOCAL.ps1"
    if (-not (Test-Path -LiteralPath $Runner)) { throw "Final local runner not found: $Runner" }
    Write-Host ""
    Write-Host "Running the complete audited local analysis before publication." -ForegroundColor Cyan
    & powershell.exe -NoProfile -ExecutionPolicy Bypass -File $Runner -Stage all -ProjectRoot $AnalysisProject
    if ($LASTEXITCODE -ne 0) { throw "Full local analysis failed. GitHub was not modified." }

    $AnalysisPython = Join-Path $AnalysisProject ".venv\Scripts\python.exe"
    if ((-not $SkipTests) -and (Test-Path -LiteralPath $AnalysisPython)) {
        Push-Location -LiteralPath $LocalDestination
        try {
            & $AnalysisPython -m py_compile (Get-ChildItem (Join-Path $LocalDestination "src") -Filter "*.py" | ForEach-Object FullName)
            if ($LASTEXITCODE -ne 0) { throw "Python compilation failed after the full analysis." }
            & $AnalysisPython -m pytest -q
            if ($LASTEXITCODE -ne 0) { throw "Unit tests failed after the full analysis." }
            $AlreadyTested = $true
        }
        finally {
            Pop-Location
        }
    }
}

$publishArgs = @("-NoProfile", "-ExecutionPolicy", "Bypass", "-File", $Publish)
if ($SkipTests -or $AlreadyTested) { $publishArgs += "-SkipTests" }
if ($SkipRelease) { $publishArgs += "-SkipRelease" }
if ($ForceRemoteMismatch) { $publishArgs += "-ForceRemoteMismatch" }
& powershell.exe @publishArgs
if ($LASTEXITCODE -ne 0) { throw "GitHub publication failed." }

Write-Host ""
Write-Host "Local analysis validation and GitHub publication finished." -ForegroundColor Green
Write-Host "Local source: $LocalDestination" -ForegroundColor Yellow
if (-not $SkipFullAnalysis) {
    Write-Host "Final analysis project: $AnalysisProject" -ForegroundColor Yellow
}
