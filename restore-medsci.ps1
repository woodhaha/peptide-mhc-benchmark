# MedSci Skill Recovery Tool
# Usage: pwsh -File restore-medsci.ps1
# One-click restore from multiple backup sources

param(
    [switch]$Force,
    [string]$TargetDir = "D:\Researching\Peptide epitope\.claude\skills\medsci"
)

$ErrorActionPreference = "Stop"
$BackupZip = "D:\Researching\Peptide epitope\.claude\skills\medsci-backup-20260620.zip"

Write-Host "╔══════════════════════════════════════════╗" -ForegroundColor Cyan
Write-Host "║   MedSci Skill Recovery Tool v1.0.0     ║" -ForegroundColor Cyan
Write-Host "╚══════════════════════════════════════════╝" -ForegroundColor Cyan
Write-Host ""

# Step 0: Check if already exists
if (Test-Path "$TargetDir\SKILL.md") {
    Write-Host "[CHECK] MedSci skill found at: $TargetDir" -ForegroundColor Green
    $fileCount = (Get-ChildItem -Recurse $TargetDir -File | Measure-Object).Count
    Write-Host "[CHECK] File count: $fileCount files" -ForegroundColor Green

    if (-not $Force) {
        Write-Host "[SKIP] Skill already exists. Use -Force to overwrite." -ForegroundColor Yellow
        exit 0
    }
    Write-Host "[FORCE] Will overwrite existing skill..." -ForegroundColor Yellow
}

# Step 1: Try ZIP backup
Write-Host ""
Write-Host "[STEP 1] Checking ZIP backup..." -ForegroundColor Cyan
if (Test-Path $BackupZip) {
    Write-Host "[FOUND] Backup ZIP: $BackupZip ($((Get-Item $BackupZip).Length) bytes)" -ForegroundColor Green
    try {
        Remove-Item -Recurse -Force $TargetDir -ErrorAction SilentlyContinue
        New-Item -ItemType Directory -Force -Path $TargetDir | Out-Null
        Expand-Archive -Path $BackupZip -DestinationPath (Split-Path $TargetDir -Parent) -Force
        Write-Host "[OK] Restored from ZIP backup!" -ForegroundColor Green

        # Verify
        if (Test-Path "$TargetDir\SKILL.md") {
            $files = (Get-ChildItem -Recurse $TargetDir -File | Measure-Object).Count
            Write-Host "[DONE] MedSci restored successfully: $files files" -ForegroundColor Green
            Get-ChildItem -Recurse $TargetDir -File | ForEach-Object { Write-Host "  - $($_.Name)" -ForegroundColor Gray }
            exit 0
        }
    } catch {
        Write-Host "[FAIL] ZIP restore failed: $_" -ForegroundColor Red
    }
} else {
    Write-Host "[MISS] No ZIP backup at: $BackupZip" -ForegroundColor Yellow
}

# Step 2: Try Yith memory restore (via Claude)
Write-Host ""
Write-Host "[STEP 2] Yith Memory Restore..." -ForegroundColor Cyan
Write-Host "[INFO] Run this in Claude Code to restore from Yith:" -ForegroundColor Yellow
Write-Host "  'Search Yith for MedSci Skill Complete Definition and restore all 14 skill files'" -ForegroundColor White
Write-Host "[INFO] Yith Memory ID: mem_mqlvduux_0e8089788964" -ForegroundColor Gray

# Step 3: Try Plan file restore
Write-Host ""
Write-Host "[STEP 3] Plan File Restore..." -ForegroundColor Cyan
$planFile = "$env:USERPROFILE\.claude\plans\stateless-toasting-hopper.md"
if (Test-Path $planFile) {
    Write-Host "[FOUND] Plan file: $planFile" -ForegroundColor Green
    Write-Host "[INFO] Run this in Claude Code:" -ForegroundColor Yellow
    Write-Host "  'Read $planFile and recreate the medsci skill from the plan'" -ForegroundColor White
} else {
    Write-Host "[MISS] No plan file at: $planFile" -ForegroundColor Yellow
}

# Step 4: Generate minimal skeleton
Write-Host ""
Write-Host "[STEP 4] Creating minimal recovery skeleton..." -ForegroundColor Cyan
$skeletonDir = "$TargetDir-recovery-skeleton"
New-Item -ItemType Directory -Force -Path "$skeletonDir\references" | Out-Null

@"
---
name: medsci
description: |
  Medical Science Research — unified literature research, review writing, and academic pipeline orchestration.
  Use when user asks to research biomedical topics, write reviews, plan papers, get peer review, or distill knowledge.
argument-hint: <topic or instruction>
allowed-tools: Bash, Read, Glob, Grep, Write, Edit, WebSearch, WebFetch
metadata:
  version: "1.0.0"
  author: "哈哥"
  prior_version_lessons: "Merged deep-research, medical-imaging-review v3.0.0, and ARS plugin v3.11.1. Recovery skeleton — use Yith memory mem_mqlvduux_0e8089788964 for full content."
---

# MedSci — Recovery Skeleton

Full skill content stored in Yith Archive (Memory ID: `mem_mqlvduux_0e8089788964`).
Run: "Search Yith for MedSci Skill Complete Definition and restore all skill files"
Or restore from: $planFile
"@ | Out-File "$skeletonDir\SKILL.md" -Encoding UTF8

@"
# Recovery Note
Full references/ files backed up in Yith memory.
Restore with: "Search Yith for MedSci Skill Complete Definition"
"@ | Out-File "$skeletonDir\references\README.md" -Encoding UTF8

Write-Host "[DONE] Skeleton created at: $skeletonDir" -ForegroundColor Green

# Summary
Write-Host ""
Write-Host "╔══════════════════════════════════════════╗" -ForegroundColor Cyan
Write-Host "║          Recovery Summary                ║" -ForegroundColor Cyan
Write-Host "╠══════════════════════════════════════════╣" -ForegroundColor Cyan
Write-Host "║  ZIP backup:      $((Test-Path $BackupZip) ? 'FOUND' : 'MISSING')                  ║" -ForegroundColor $(if (Test-Path $BackupZip) { 'Green' } else { 'Red' })
Write-Host "║  Yith memory:     mem_mqlvduux           ║" -ForegroundColor Green
Write-Host "║  Project memory:  .omc/project-memory    ║" -ForegroundColor Green
Write-Host "║  Plan file:       $((Test-Path $planFile) ? 'FOUND' : 'MISSING')                  ║" -ForegroundColor $(if (Test-Path $planFile) { 'Green' } else { 'Red' })
Write-Host "║  Session history: .omc/sessions/         ║" -ForegroundColor Green
Write-Host "╚══════════════════════════════════════════╝" -ForegroundColor Cyan
Write-Host ""
Write-Host "To restore: say 'restore medsci skill' in Claude Code" -ForegroundColor White
