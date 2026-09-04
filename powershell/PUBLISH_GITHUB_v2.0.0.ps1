[CmdletBinding()]
param(
    [string]$RepositoryUrl = "https://github.com/EJSAHN/sorghum-leaf-blight-meta-gwas.git",
    [string]$RepositoryFullName = "EJSAHN/sorghum-leaf-blight-meta-gwas",
    [string]$WorkRoot = "D:\projects\sorghum-leaf-blight-meta-gwas_publish_v2.0.0",
    [string]$ExpectedRemoteMain = "37e9b2f22abb082d1efab2c58f4759767131b999",
    [string]$BackupTag = "pre-revision-exploratory-2026-04-27",
    [string]$ReleaseTag = "v2.0.0",
    [switch]$ForceRemoteMismatch,
    [switch]$SkipTests,
    [switch]$SkipRelease
)

$ErrorActionPreference = "Stop"
Set-StrictMode -Version Latest

function Banner([string]$Message) {
    Write-Host ""
    Write-Host ("=" * 78) -ForegroundColor Cyan
    Write-Host $Message -ForegroundColor Cyan
    Write-Host ("=" * 78) -ForegroundColor Cyan
}

function Get-GitExecutable {
    $command = Get-Command git.exe -ErrorAction SilentlyContinue
    if ($command) { return $command.Source }
    foreach ($candidate in @(
        "C:\Program Files\Git\cmd\git.exe",
        "C:\Program Files\Git\bin\git.exe",
        "$env:LOCALAPPDATA\Programs\Git\cmd\git.exe"
    )) {
        if (Test-Path -LiteralPath $candidate -PathType Leaf) { return $candidate }
    }
    throw "Git for Windows was not found. Install Git for Windows, then rerun this script."
}

function Run-Git {
    param([string[]]$Arguments, [string]$WorkingDirectory)
    Push-Location -LiteralPath $WorkingDirectory
    try {
        & $script:GitExe @Arguments
        if ($LASTEXITCODE -ne 0) {
            throw "git command failed: git $($Arguments -join ' ')"
        }
    }
    finally {
        Pop-Location
    }
}

function Get-PythonLauncher {
    if (Get-Command py.exe -ErrorAction SilentlyContinue) {
        & py.exe -3.11 -c "import sys; print(sys.version)" *> $null
        if ($LASTEXITCODE -eq 0) { return @{ Exe = "py.exe"; Prefix = @("-3.11") } }
        & py.exe -3 -c "import sys; print(sys.version)" *> $null
        if ($LASTEXITCODE -eq 0) { return @{ Exe = "py.exe"; Prefix = @("-3") } }
    }
    if (Get-Command python.exe -ErrorAction SilentlyContinue) {
        return @{ Exe = "python.exe"; Prefix = @() }
    }
    throw "Python was not found. Install Python 3.11 or the Windows py launcher."
}

$SourceRoot = (Resolve-Path (Join-Path $PSScriptRoot "..")).Path
$Version = (Get-Content -LiteralPath (Join-Path $SourceRoot "VERSION") -Raw).Trim()
if ($Version -ne "2.0.0") {
    throw "Expected source version 2.0.0, found $Version"
}
$script:GitExe = Get-GitExecutable

$GhForAuth = Get-Command gh.exe -ErrorAction SilentlyContinue
if ($GhForAuth) {
    & $GhForAuth.Source auth status *> $null
    if ($LASTEXITCODE -eq 0) {
        & $GhForAuth.Source auth setup-git
        if ($LASTEXITCODE -ne 0) {
            throw "GitHub CLI is authenticated, but gh auth setup-git failed."
        }
    }
}

$ManifestPath = Join-Path $SourceRoot "SOURCE_MANIFEST_SHA256.txt"
if (Test-Path -LiteralPath $ManifestPath) {
    Banner "Verifying source package SHA-256 manifest"
    $ManifestLines = Get-Content -LiteralPath $ManifestPath | Where-Object { $_ -and -not $_.StartsWith("#") }
    foreach ($line in $ManifestLines) {
        if ($line -notmatch "^([A-Fa-f0-9]{64})\s{2}(.+)$") {
            throw "Malformed source manifest line: $line"
        }
        $ExpectedHash = $Matches[1].ToUpper()
        $RelativePath = $Matches[2]
        $Target = Join-Path $SourceRoot $RelativePath
        if (-not (Test-Path -LiteralPath $Target -PathType Leaf)) {
            throw "Source manifest file is missing: $RelativePath"
        }
        $ObservedHash = (Get-FileHash -LiteralPath $Target -Algorithm SHA256).Hash.ToUpper()
        if ($ObservedHash -ne $ExpectedHash) {
            throw "Source manifest mismatch: $RelativePath"
        }
    }
    Write-Host "Source manifest verified." -ForegroundColor Green
}

Banner "Pre-publish source validation"
if (-not $SkipTests) {
    $PythonLauncher = Get-PythonLauncher
    $TestVenv = Join-Path $WorkRoot "_test_venv"
    New-Item -ItemType Directory -Path $WorkRoot -Force | Out-Null
    if (-not (Test-Path -LiteralPath (Join-Path $TestVenv "Scripts\python.exe"))) {
        & $PythonLauncher.Exe @($PythonLauncher.Prefix) -m venv $TestVenv
        if ($LASTEXITCODE -ne 0) { throw "Failed to create test environment." }
    }
    $TestPython = Join-Path $TestVenv "Scripts\python.exe"
    & $TestPython -m pip install --upgrade pip
    if ($LASTEXITCODE -ne 0) { throw "pip upgrade failed." }
    & $TestPython -m pip install -r (Join-Path $SourceRoot "requirements.txt")
    if ($LASTEXITCODE -ne 0) { throw "Dependency installation failed." }
    & $TestPython -m py_compile (Get-ChildItem (Join-Path $SourceRoot "src") -Filter "*.py" | ForEach-Object FullName)
    if ($LASTEXITCODE -ne 0) { throw "Python compilation failed." }
    Push-Location -LiteralPath $SourceRoot
    try {
        & $TestPython -m pytest -q
        if ($LASTEXITCODE -ne 0) { throw "Unit tests failed." }
    }
    finally {
        Pop-Location
    }
}
else {
    Write-Host "Tests skipped by explicit switch." -ForegroundColor DarkYellow
}

Banner "Cloning current GitHub repository"
$CloneDir = Join-Path $WorkRoot "repository"
if (Test-Path -LiteralPath $CloneDir) {
    Remove-Item -LiteralPath $CloneDir -Recurse -Force
}
New-Item -ItemType Directory -Path $WorkRoot -Force | Out-Null
& $script:GitExe clone $RepositoryUrl $CloneDir
if ($LASTEXITCODE -ne 0) {
    throw "Repository clone failed. Confirm network access and GitHub credentials."
}

# Git for Windows may reject repositories on filesystems that do not record
# Unix-style ownership. Trust only this exact fresh checkout, never a wildcard.
$SafeCloneDir = $CloneDir.Replace("\", "/")
$ConfiguredSafeDirs = @(& $script:GitExe config --global --get-all safe.directory 2>$null)
if ($ConfiguredSafeDirs -notcontains $SafeCloneDir) {
    & $script:GitExe config --global --add safe.directory $SafeCloneDir
    if ($LASTEXITCODE -ne 0) {
        throw "Could not register the fresh checkout as a Git safe.directory."
    }
}

Run-Git -Arguments @("fetch", "origin", "--tags", "--prune") -WorkingDirectory $CloneDir
$RemoteMain = (& $script:GitExe -C $CloneDir rev-parse "origin/main").Trim()
if ($LASTEXITCODE -ne 0) { throw "Could not resolve origin/main." }
Write-Host "Remote main before update: $RemoteMain"
if (($RemoteMain -ne $ExpectedRemoteMain) -and (-not $ForceRemoteMismatch)) {
    throw "Remote main changed from the audited SHA $ExpectedRemoteMain. No files were changed. Inspect the remote and rerun with -ForceRemoteMismatch only after confirming the difference."
}

Banner "Preserving the pre-revision repository state"
$RemoteTagLines = @(& $script:GitExe -C $CloneDir ls-remote --tags origin ("refs/tags/" + $BackupTag))
if ($LASTEXITCODE -ne 0) {
    throw "Could not query the remote backup tag."
}
$RemoteTag = ($RemoteTagLines -join "`n").Trim()
if (-not $RemoteTag) {
    Run-Git -Arguments @("tag", "-a", $BackupTag, $RemoteMain, "-m", "Archive of the exploratory workflow before the PEI revision") -WorkingDirectory $CloneDir
    Run-Git -Arguments @("push", "origin", $BackupTag) -WorkingDirectory $CloneDir
    Write-Host "Pushed backup tag: $BackupTag" -ForegroundColor Green
}
else {
    Write-Host "Backup tag already exists: $BackupTag"
}

Banner "Replacing the working tree with the audited v2.0.0 source"
Get-ChildItem -LiteralPath $CloneDir -Force |
    Where-Object { $_.Name -ne ".git" } |
    Remove-Item -Recurse -Force
Get-ChildItem -LiteralPath $SourceRoot -Force |
    Where-Object { $_.Name -notin @(".git", ".venv", "__pycache__", ".pytest_cache") } |
    ForEach-Object {
        Copy-Item -LiteralPath $_.FullName -Destination $CloneDir -Recurse -Force
    }

# Refuse to publish study data or unexpectedly large files.
$Forbidden = @(Get-ChildItem -LiteralPath $CloneDir -Recurse -File -Force | Where-Object {
    $_.FullName -notmatch "\\.git\\" -and (
        $_.Name -match "(?i)\\.(vcf|vcf\\.gz|bcf|bed|bim|fam|xlsx|xls|assoc\\.txt|cXX\\.txt|npy|npz|bin)$" -or
        $_.Length -gt 25MB
    )
})
if ($Forbidden.Count -gt 0) {
    $Forbidden | Select-Object FullName, Length | Format-Table -AutoSize
    throw "Refusing to publish raw/derived study data or files larger than 25 MB."
}

if (-not (& $script:GitExe -C $CloneDir config user.name)) {
    Run-Git -Arguments @("config", "user.name", "Ezekiel Ahn") -WorkingDirectory $CloneDir
}
if (-not (& $script:GitExe -C $CloneDir config user.email)) {
    Run-Git -Arguments @("config", "user.email", "ezekiel.ahn@usda.gov") -WorkingDirectory $CloneDir
}

Run-Git -Arguments @("add", "-A") -WorkingDirectory $CloneDir
Write-Host ""
Write-Host "Files staged for v2.0.0:" -ForegroundColor Yellow
& $script:GitExe -C $CloneDir status --short
if ($LASTEXITCODE -ne 0) { throw "Could not inspect staged files." }

$Staged = (& $script:GitExe -C $CloneDir diff --cached --name-only)
if (-not $Staged) {
    throw "No changes were staged; repository may already contain v2.0.0."
}

Banner "Committing and pushing main"
Run-Git -Arguments @(
    "commit", "-m",
    "Release v2.0.0: corrected multi-environment phenotype reconstruction and GEMMA LOCO inference"
) -WorkingDirectory $CloneDir
Run-Git -Arguments @("push", "origin", "main") -WorkingDirectory $CloneDir
$NewCommit = (& $script:GitExe -C $CloneDir rev-parse HEAD).Trim()
Write-Host "Pushed main commit: $NewCommit" -ForegroundColor Green

$ExistingReleaseTagLines = @(& $script:GitExe -C $CloneDir ls-remote --tags origin ("refs/tags/" + $ReleaseTag))
if ($LASTEXITCODE -ne 0) {
    throw "Could not query the remote release tag."
}
$ExistingReleaseTag = ($ExistingReleaseTagLines -join "`n").Trim()
if ($ExistingReleaseTag) {
    throw "Release tag already exists remotely: $ReleaseTag. Main was pushed, but the tag was not replaced."
}
Run-Git -Arguments @("tag", "-a", $ReleaseTag, "-m", "Sorghum leaf blight revised analysis workflow v2.0.0") -WorkingDirectory $CloneDir
Run-Git -Arguments @("push", "origin", $ReleaseTag) -WorkingDirectory $CloneDir
Write-Host "Pushed release tag: $ReleaseTag" -ForegroundColor Green

if (-not $SkipRelease) {
    $Gh = Get-Command gh.exe -ErrorAction SilentlyContinue
    if ($Gh) {
        & $Gh.Source auth status *> $null
        if ($LASTEXITCODE -eq 0) {
            & $Gh.Source release create $ReleaseTag --repo $RepositoryFullName --title "v2.0.0 - revised multi-environment analysis" --notes-file (Join-Path $CloneDir "RELEASE_NOTES_v2.0.0.md")
            if ($LASTEXITCODE -eq 0) {
                Write-Host "GitHub Release created." -ForegroundColor Green
            }
            else {
                Write-Host "Tag was pushed, but GitHub Release creation failed. Create the release manually from the tag." -ForegroundColor DarkYellow
            }
        }
        else {
            Write-Host "GitHub CLI is installed but not authenticated. The tag is pushed; create the release manually if desired." -ForegroundColor DarkYellow
        }
    }
    else {
        Write-Host "GitHub CLI is not installed. Main and tag are updated; create a GitHub Release from tag v2.0.0 manually if desired." -ForegroundColor DarkYellow
    }
}

Banner "GitHub update complete"
Write-Host "Repository: https://github.com/$RepositoryFullName" -ForegroundColor Yellow
Write-Host "New commit: $NewCommit" -ForegroundColor Yellow
Write-Host "Backup tag: $BackupTag"
Write-Host "Release tag: $ReleaseTag"
Write-Host "Local publish checkout: $CloneDir"
