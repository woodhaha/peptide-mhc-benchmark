# Knowledge Distillation Protocol

## Concept

The MedSci skill can **absorb knowledge from other skills** over time, expanding its domain coverage and methodological toolkit. This protocol defines the mechanism: a structured 4-step cycle that extracts, validates, and merges knowledge fragments from external skill definitions into MedSci's domain registry.

---

## The Distill Loop

```
    ┌──────────┐
    │ OBSERVE  │  Identify knowledge-bearing content in a source skill
    └────┬─────┘
         ▼
    ┌──────────┐
    │ EXTRACT  │  Structure it as a Fragment: {domain, type, content, source, confidence, date}
    └────┬─────┘
         ▼
    ┌──────────┐
    │  MERGE   │  Hash-detect duplicates → auto-merge low-conflict → flag for human
    └────┬─────┘
         ▼
    ┌──────────┐
    │  APPLY   │  Write to DOMAINS.md Section C (Registry) or Section A (if promoted)
    └──────────┘
```

---

## Distillable Fragment Types

| Type | Description | Example |
|------|-------------|---------|
| `terminology` | Domain-specific term with definition | "Neoepitope: tumor-specific peptide from somatic mutation" |
| `taxonomy-axis` | Method classification axis | Axis 1: Architectural Priors, Axis 2: Data Priors, Axis 3: Validation Regime |
| `data-source` | Database, API, or repository | "IEDB: Immune Epitope Database, https://www.iedb.org/" |
| `benchmark` | Standard evaluation metric or dataset | "IEDB weekly benchmark: AUROC, PPV at top-N" |
| `method-family` | A family of related methods | "Pan-specific MHC binding predictors (NetMHCpan family)" |
| `key-author` | Leading researcher + their contribution area | "M. Nielsen (DTU) — NetMHCpan series, immunoinformatics" |
| `key-paper` | Foundational paper in the domain | "Jurtz et al. (2017) NetMHCpan 4.0, J Immunol" |
| `writing-convention` | Domain-specific writing rule | "Allele names: HLA-A*02:01 (4-digit resolution minimum)" |
| `quality-rule` | Domain-specific verification check | "Verify binding affinity units: IC50 (nM) vs % rank vs eluted ligand score" |
| `ars-agent-pattern` | ARS agent behavior worth noting | "synthesis_agent v3.7.3: Three-Layer Citation Emission (ref slug + anchor + kind)" |

---

## Fragment Schema

```json
{
  "id": "<content_hash_prefix_8>",
  "source_skill": "<skill-name>",
  "source_version": "<version-or-date>",
  "domain": "<domain-slug>",
  "type": "<fragment-type>",
  "content": "<the actual knowledge, markdown>",
  "confidence": "verified | extracted | tentative",
  "date_distilled": "<YYYY-MM-DD>",
  "last_verified": "<YYYY-MM-DD or null>",
  "related_fragments": ["<fragment-id>"],
  "notes": "<optional: extraction context, why this is useful>"
}
```

---

## Confidence Levels

| Level | Meaning | Merge Behavior |
|-------|---------|---------------|
| **verified** | Manually confirmed against primary source | Auto-merge, no flag |
| **extracted** | From skill text, structurally sound, plausible | Auto-merge with `[extracted]` annotation; flag if conflict |
| **tentative** | Inferred or uncertain; needs verification | Store in DOMAINS.md Section C with `[TENTATIVE]` tag; do NOT apply to Section A |

---

## Manual Distillation Command

**User trigger**: "distill from [skill-name]" / "absorb [skill-name]" / "learn from [skill-name]"

**Procedure**:
1. **Locate** the source skill's SKILL.md. Search paths: `~/.claude/skills/<name>/SKILL.md`, project `.claude/skills/<name>/SKILL.md`, ARS plugin skills directories.
2. **Read** SKILL.md + any domain-specific reference files (DOMAINS.md, references/, MODE_REGISTRY.md).
3. **Extract** all identifiable fragments of each type. For each fragment, assign a confidence level:
   - `verified`: Found the exact claim in a peer-reviewed source or official documentation
   - `extracted`: Plausible, structurally correct, consistent with domain knowledge
   - `tentative`: Inferred from context, ambiguous, or unverifiable from skill text alone
4. **Check** for duplicates by searching DOMAINS.md Section C for similar content.
5. **Merge**: New fragments appended to Section C. Conflicts (same term, different definition) logged to DISTILLATION_LOG.md with `[CONFLICT]` tag.
6. **Report**: Summarize what was extracted (N fragments of M types, K conflicts, L tentative).

---

## ARS Distillation Special Handling

The ARS plugin (academic-research-skills v3.11.1) has 42 agents across 4 skills. Distilling from ARS requires navigating this structure:

### Distilling from "ars:<skill-name>"

| Target | What to Extract |
|--------|----------------|
| `ars:deep-research` | MODE_REGISTRY entries for 7 modes, 14 agent role definitions, search strategy patterns |
| `ars:academic-paper` | 10 mode definitions, 12 agent role definitions, writing discipline patterns |
| `ars:academic-paper-reviewer` | 6 review mode definitions, 7 reviewer agent profiles, sprint contract protocol |
| `ars:academic-pipeline` | 10-stage pipeline state machine, integrity gate definitions, Material Passport schema references |

### ARS-Specific Fragment Types

| Type | What to Extract |
|------|----------------|
| `ars-agent-pattern` | Agent role + responsibility + trigger timing + phase boundary |
| `ars-mode-definition` | Mode name + spectrum position + outputs + oversight level |
| `ars-gate-definition` | Gate name + stage + check type (decision-heavy / integrity / machine-only) |
| `ars-schema-reference` | Schema name + purpose + key fields |

---

## Distillation Log

All distillation events are logged to `DISTILLATION_LOG.md` (created on first use):

```markdown
| Date | Source Skill | Fragments | Types | Conflicts | Status |
|------|-------------|-----------|-------|-----------|--------|
| 2026-06-20 | neoantigen-predictor | 8 | terminology:3, data-source:2, benchmark:1, key-paper:2 | 0 | merged |
| 2026-06-20 | ars:deep-research | 21 | ars-agent-pattern:14, ars-mode-definition:7 | 1 | merged (1 conflict: synthesis scope) |
```

---

## Merge Conflict Resolution

When two fragments define the same term differently:

1. **Both verified**: Keep both, note the disagreement in content as `[DISPUTED]`. Example: two skills give different AUC baselines for NetMHCpan — store both with source attribution.
2. **One verified, one extracted**: Verified wins. Extracted version stored as alternate with `[SUPERSEDED]` tag.
3. **Both extracted**: Flag for human review in DOMAINS.md Section C with `[CONFLICT: needs human review]`. Do not auto-resolve.
4. **Tentative + anything**: Tentative is demoted; the other version (verified or extracted) is applied.

---

## Auto-Distillation Hooks (Future)

When a skill in the `research-skills` family or a registered ARS skill is invoked in the same session:
- MedSci can optionally scan its output for domain patterns (opt-in, to avoid noise)
- Triggered by `MEDSCI_AUTO_DISTILL=1` environment variable
- Auto-distilled fragments are always `tentative` confidence; manual review required before promotion

---

## Related Files

| File | Role |
|------|------|
| `DOMAINS.md` Section C | Distilled knowledge registry |
| `DISTILLATION_LOG.md` | Distillation event log (created on first use) |
| `SKILL.md` § Distillation Quick Reference | User-facing distillation commands |
