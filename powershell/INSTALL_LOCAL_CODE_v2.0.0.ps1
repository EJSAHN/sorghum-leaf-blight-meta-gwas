[CmdletBinding()]
param(
    [string]$Destination = "D:\projects\sorghum-leaf-blight-meta-gwas_v2.0.0_local",
    [switch]$Force
)

$ErrorActionPreference = "Stop"
Set-StrictMode -Version Latest

$SourceRoot = (Resolve-Path (Join-Path $PSScriptRoot "..")).Path
$VersionFile = Join-Path $SourceRoot "VERSION"
if (-not (Test-Path -LiteralPath $VersionFile)) {
    throw "VERSION file not found under source root: $SourceRoot"
}
$Version = (Get-Content -LiteralPath $VersionFile -Raw).Trim()
if ($Version -ne "2.0.0") {
    throw "Expected source version 2.0.0, found $Version"
}

if (Test-Path -LiteralPath $Destination) {
    if (-not $Force) {
        throw "Destination already exists: $Destination. Re-run with -Force only if replacement is intended."
    }
    Remove-Item -LiteralPath $Destination -Recurse -Force
}

New-Item -ItemType Directory -Path $Destination -Force | Out-Null
Get-ChildItem -LiteralPath $SourceRoot -Force |
    Where-Object { $_.Name -notin @(".git", ".venv", "__pycache__", ".pytest_cache") } |
    ForEach-Object {
        Copy-Item -LiteralPath $_.FullName -Destination $Destination -Recurse -Force
    }

Write-Host ""
Write-Host "Local final code installed:" -ForegroundColor Green
Write-Host $Destination -ForegroundColor Yellow
Write-Host ""
Write-Host "To run the full audited analysis:" -ForegroundColor Cyan
Write-Host ("powershell.exe -NoProfile -ExecutionPolicy Bypass -File `"{0}`" -Stage all" -f (Join-Path $Destination "SETUP_AND_RUN_LOCAL.ps1"))
