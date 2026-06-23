# MedSci Skill Recovery Guide

## Backup Locations (in priority order)

| # | Storage | Location | Recovery Method |
|---|---------|----------|----------------|
| 1 | **Local files** | `.claude/skills/medsci/` (13 files) | Already active — no recovery needed |
| 2 | **Tarball backup** | `.claude/skills/medsci-backup-YYYYMMDD.tar.gz` | `tar -xzf <file>` |
| 3 | **Yith Archive** | Memory ID: `mem_mqlvduux_0e8089788964` | `yith_search "MedSci Skill Complete Definition"` |
| 4 | **Project Memory** | `.omc/project-memory.json` § architecture + directive | `project_memory_read` |
| 5 | **Plan file** | `~/.claude/plans/stateless-toasting-hopper.md` | Read and re-execute |
| 6 | **Session history** | `.omc/sessions/` | `session_search "medsci"` |

## Quick Recovery Commands

### Check if medsci exists
```bash
ls .claude/skills/medsci/SKILL.md
```

### Restore from tarball (if local backup exists)
```bash
cd .claude/skills && tar -xzf medsci-backup-*.tar.gz
```

### Restore from Yith memory
Ask Claude: "Search Yith for 'MedSci Skill Complete Definition' and restore the skill files"

### Restore from plan
Ask Claude: "Read ~/.claude/plans/stateless-toasting-hopper.md and recreate the medsci skill"

## Minimal Recovery (if all backups lost)

The essential structure:
```
.claude/skills/medsci/
├── SKILL.md           ← Router (can be minimal)
└── references/
    ├── WORKFLOW.md    ← Core workflow
    ├── DOMAINS.md     ← Domain knowledge
    └── TOOL_STRATEGY.md ← Tool fallback
```

The other 9 reference files can be re-derived from these 3 + re-distillation.

## Version Summary

- **v1.0.0** (2026-06-20): Initial merge of deep-research + medical-imaging-review + ARS plugin
- **Distilled sources**: QUEST, CCF-Figure, neoantigen-predictor, immune-pathway-analysis, biomarker-landscape-scanner, adme-property-predictor, sci-hub-search-skill/paper-search-mcp
- **Total knowledge fragments**: ~96 across 8 domains
