#!/bin/bash
# MedSci Skill Recovery Tool (Bash version for Git Bash / Linux / macOS)
# Usage: bash restore-medsci.sh [--force]

set -e

TARGET_DIR="D:/Researching/Peptide epitope/.claude/skills/medsci"
BACKUP_ZIP="D:/Researching/Peptide epitope/.claude/skills/medsci-backup-20260620.zip"
FORCE=false

if [ "$1" = "--force" ]; then FORCE=true; fi

echo "=================================="
echo "  MedSci Skill Recovery Tool v1.0"
echo "=================================="
echo ""

# Check if already exists
if [ -f "$TARGET_DIR/SKILL.md" ]; then
    echo "[CHECK] MedSci skill found at: $TARGET_DIR"
    COUNT=$(find "$TARGET_DIR" -type f | wc -l)
    echo "[CHECK] File count: $COUNT files"
    if [ "$FORCE" = false ]; then
        echo "[SKIP] Skill already exists. Use --force to overwrite."
        exit 0
    fi
    echo "[FORCE] Will overwrite..."
fi

# Try ZIP backup
echo ""
echo "[STEP 1] Checking ZIP backup..."
if [ -f "$BACKUP_ZIP" ]; then
    echo "[FOUND] Backup ZIP exists"
    rm -rf "$TARGET_DIR"
    mkdir -p "$TARGET_DIR"
    unzip -o "$BACKUP_ZIP" -d "$(dirname "$TARGET_DIR")"
    if [ -f "$TARGET_DIR/SKILL.md" ]; then
        COUNT=$(find "$TARGET_DIR" -type f | wc -l)
        echo "[DONE] Restored $COUNT files from ZIP!"
        find "$TARGET_DIR" -type f | while read f; do echo "  - $(basename "$f")"; done
        exit 0
    fi
else
    echo "[MISS] No ZIP backup"
fi

# Try Yith / Plan
echo ""
echo "[STEP 2] Manual recovery via Claude Code:"
echo "  'Search Yith for MedSci Skill Complete Definition and restore all 14 skill files'"
echo "  Memory ID: mem_mqlvduux_0e8089788964"
echo ""

PLAN_FILE="$HOME/.claude/plans/stateless-toasting-hopper.md"
if [ -f "$PLAN_FILE" ]; then
    echo "[STEP 3] Plan file found. In Claude Code:"
    echo "  'Read $PLAN_FILE and recreate the medsci skill'"
fi

# Create minimal skeleton
echo ""
echo "[STEP 4] Creating minimal recovery skeleton..."
SKELETON="${TARGET_DIR}-recovery-skeleton"
mkdir -p "$SKELETON/references"
cat > "$SKELETON/SKILL.md" << 'SKILLEOF'
---
name: medsci
description: Medical Science Research skill. Full content in Yith memory mem_mqlvduux_0e8089788964.
metadata:
  version: "1.0.0-recovery"
  author: "哈哥"
---

# MedSci — Recovery Skeleton
Restore full skill: "Search Yith for MedSci Skill Complete Definition"
SKILLEOF
echo "[DONE] Skeleton at: $SKELETON"
echo ""
echo "To restore: say 'restore medsci skill' in Claude Code"
