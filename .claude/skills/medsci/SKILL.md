---
name: medsci
description: |
  Medical Science Research — unified literature research, review writing, and academic pipeline orchestration.
  
  Use when the user asks to:
  - Research biomedical/immunology topics: "research peptide-MHC binding", "investigate neoantigen prediction", "find papers on epitope discovery"
  - Write academic manuscripts: "write a review", "write a paper", "综述", "投稿", "draft a manuscript"
  - Plan or outline: "plan my paper", "outline", "大纲", "规划论文结构"
  - Get peer review: "review my paper", "审稿", "peer review"
  - Revise manuscripts: "revise", "modify manuscript", "修改稿", "respond to reviewers"
  - Check or convert: "check citations", "format convert", "abstract", "disclosure"
  - Absorb knowledge: "distill from <skill>", "absorb <skill>", "learn from <skill>"
  
  Do NOT use for:
  - Non-academic writing (blog posts, social media)
  - Non-biomedical topics without explicit domain adaptation
  - Simple fact lookup (use WebSearch directly)
  - Implementation/coding tasks (use appropriate coding agents)
  
  Default domain: biomedical, immunology, epitope prediction, cancer immunotherapy.
  Generalizable to any biomedical sub-domain via DOMAINS.md template system.
  
  Routes to ARS (academic-research-skills) plugin for heavy academic pipeline work.
  Executes autonomous 7-phase research pipeline for lighter literature tasks.
argument-hint: <topic or instruction>
allowed-tools: Bash, Read, Glob, Grep, Write, Edit, WebSearch, WebFetch
metadata:
  version: "1.0.0"
  author: "user"
  prior_version_lessons: "Merged deep-research (5-phase web-first), medical-imaging-review v3.0.0 (6-phase write-with-verify, citation integrity, hallucination patterns), and ARS plugin v3.11.1 (full academic pipeline: 14 commands, 42 agents, 10-stage workflow). Key innovations: 4-mode routing eliminates unnecessary pipeline stages for quick research, 3-tier tool fallback eliminates CLI dependency, domain template system allows adaptation beyond medical imaging, knowledge distillation protocol enables iterative expansion."
---

# MedSci — Medical Science Research Skill

Unified biomedical research skill merging fast web-first literature research with publication-quality academic pipeline orchestration.

---

## Quick Start

**Step 1 — Mode Detection (automatic):**

```
User intent includes "research"/"investigate"/"find papers" 
  → MODE 1: Autonomous Research (fast, no ARS dependency)
  
User intent includes "write a paper"/"review"/"综述"/"投稿"/"审稿"/"revise"
  → MODE 2: ARS Delegation (routes to appropriate ARS slash command)

User intent includes "comprehensive review"/"systematic review"/"投稿级综述"
  → MODE 3: Deep Review (ARS full pipeline + MedSci quality layers)

User intent includes "distill from"/"absorb"/"learn from"
  → MODE 4: Quick Distill (knowledge absorption from other skills)
```

**Step 2 — Execute the mode's workflow.** See `references/WORKFLOW.md` for complete phase-by-phase instructions.

**Step 3 — Deliver the output.** Research brief (Mode 1), manuscript (Mode 3), or distillation report (Mode 4).

---

## Core Principles

### 1. Voice to Evidence (from MIR)
Match language strength to evidence. When ≥2 independent groups confirm a finding, state it strongly. When single-source or contested, state it cautiously. When evidence is absent, say so.

**Avoid these LLM tells:**
- "has shown promising results" → cite the actual result
- "may suggest" → state what the evidence supports
- "interestingly," → let the finding speak
- "it is worth noting that" → just note it
- "in recent years," → give the year or drop it

### 2. Verify on the Way In (merged)
Metadata is verified BEFORE a paper enters the corpus — not after the draft is written. Every paper that survives filtering gets: DOI/PMID resolution (Rule 1), author verification (Rule 2), source type check (Rule 5).

### 3. Read First, Write After (from MIR)
Do not template-fill method descriptions. Read the actual paper → note actual modules/benchmarks/numbers → write from those notes → verify 1-2 numbers against the paper. If you cannot access the paper, do not write about its internal architecture or specific numbers.

### 4. Start Broad, Then Narrow (from DR)
Don't deep-read until you've scanned widely. The funnel: broad search (50+ results) → filter (>10 kept) → deep read (3-8 papers) → synthesize. Cross-reference web sources against academic papers.

### 5. Track What You've Read (from DR)
Deduplicate. Track search terms used, papers already read, and papers rejected with reasons. Don't re-search the same query; don't re-read the same paper.

### 6. Route, Don't Replicate (MedSci innovation)
ARS has 898 files and 42 agents. MedSci routes to it — it does NOT duplicate it. For full academic pipelines, delegate to ARS. MedSci adds: domain adaptation, MIR quality layers, tool fallback, and knowledge distillation.

---

## Mode Reference

| Mode | Trigger Keywords | Phases | Output | ARS Dependency |
|------|-----------------|--------|--------|---------------|
| **1: Autonomous Research** | research, investigate, find papers, what is known, explore, 调研 | 0-6 (lightweight) | Research Brief (.md) | None |
| **2: ARS Delegation** | write a paper, review, 综述, 投稿, 审稿, revise, plan, outline, abstract | Route to ARS command | ARS output | Required |
| **3: Deep Review** | comprehensive review, systematic review, 系统综述, 投稿级综述 | ARS pipeline + MedSci quality | Manuscript + Material Passport | Required |
| **4: Quick Distill** | distill from, absorb, learn from | Distill loop | Distillation Report + DOMAINS.md update | None |

### Mode 2 Routing Table

| User Intent | ARS Command | Model |
|------------|-------------|-------|
| Full paper pipeline | `/ars-full` | Opus |
| Chapter planning | `/ars-plan` | Sonnet |
| Detailed outline | `/ars-outline` | Sonnet |
| Literature review section | `/ars-lit-review` | Sonnet |
| Bilingual abstract | `/ars-abstract` | Sonnet |
| Peer review | `/ars-reviewer` | Opus |
| Revision coaching | `/ars-revision-coach` | Opus |
| Manuscript revision | `/ars-revision` | Sonnet |
| Citation check | `/ars-citation-check` | Sonnet |
| Format conversion | `/ars-format-convert` | Sonnet |
| AI disclosure | `/ars-disclosure` | Sonnet |
| Mark as read | `/ars-mark-read` | Sonnet |

---

## Heading Depth

Max 2 heading levels in body (H2 + H3). No H4 (`####`). Use `**Topic.**` bold lead-in for deeper grouping. No number prefixes on headings (1., 1.1, 1.2.3 are forbidden).

---

## Verdict Sentences

Each major section closes with 1 verdict sentence expressing authorial position. 3-5 per manuscript. Templates:
- "[Family] is currently the most effective approach for [problem]."
- "[Family] has yet to demonstrate clear advantage over [alternative]."
- "[Family] is best understood as complementary to [alternative], not a replacement."
- "The next [N] years will determine whether [family] becomes the standard or remains experimental."

Neutral catalogue is the LLM default — actively resist it.

---

## Equations & Vendor Names

- **Equations**: Display equations (`$$`) only in Boxes. Textbook formulas → describe in prose. Box 1 (metrics) is standard.
- **Vendor/Tool names**: Confined to Table 3 (Tools/Products). Body text uses category descriptors with table cross-reference.

---

## Mode 1: Standard Section Structure

### Research Brief
```
# Research Brief: <Topic>
## Key Takeaways
## Landscape Summary
## Key Papers
## Gaps & Opportunities
## References
```

### Review Manuscript
```
# <Title>: <Subtitle>
## Key Points (4-5 bullets)
## Abstract
## Introduction
## Datasets and Evaluation Metrics (Table 1, Box 1)
## Methods (3 H3 subsections — one per taxonomic axis)
## Applications
## Translation to Practice (Table 3)
## Outstanding Challenges
## Future Directions
## References
```

---

## Distillation Quick Reference

**Command**: "medsci, distill from `<skill-name>`"

**What happens**: MedSci reads the target skill's SKILL.md → extracts structured knowledge (terminology, data-sources, benchmarks, method-families) → merges into DOMAINS.md Section C → reports what was extracted.

**Confidence levels**:
- `verified` — manually confirmed against primary source → auto-merged
- `extracted` — from skill text, structurally sound → auto-merged with annotation
- `tentative` — inferred, needs verification → stored but not applied

**Target skills for distillation**: neoantigen-predictor, immune-pathway-analysis, biomarker-landscape-scanner, adme-property-predictor, and any of the 400+ registered biomedical skills.

---

## Reference Files

| File | Read When | Content |
|------|-----------|---------|
| `references/WORKFLOW.md` | Every session — workflow execution | Complete 4-mode phase-by-phase instructions |
| `references/TOOL_STRATEGY.md` | Phase 1 — tool detection | 3-tier fallback chain (MCP → CLI → WebSearch) |
| `references/DOMAINS.md` | Phase 0 (scope) + Phase 3 (taxonomy) | Default biomedical domain + templates + registry |
| `references/CITATION_INTEGRITY.md` | Phase 2, 4, 5, 6 (verification) | 5-rule citation verification protocol |
| `references/HALLUCINATION_PATTERNS.md` | Phase 4 (self-check) + Phase 6 (audit) | 11 hallucination patterns + ARS crosswalk |
| `references/PARADIGM.md` | Phase 0 (style anchor) | Exemplar capture for writing style |
| `references/TEMPLATES.md` | Phase 4 (write) | Research brief, manuscript, project file templates |
| `references/QUALITY_CHECKLIST.md` | Phase 6 (audit) | Mode-aware quality gates with self-check commands |
| `references/DISTILLATION_PROTOCOL.md` | Mode 4 (distill) | Knowledge absorption mechanism |
| `references/SKILL_CATALOG.md` | Mode 2 (route) + distillation planning | ARS 14 commands + related skills quick reference |
| `references/VERSIONING.md` | On iteration | Version history + bump checklist |

---

## Related Skills

| Skill | Relationship |
|-------|-------------|
| `deep-research` | Source: 5-phase web-first research workflow (merged into Mode 1) |
| `medical-imaging-review` | Source: citation integrity, hallucination patterns, 3-axis taxonomy, quality checklists (merged into reference files) |
| `academic-research-skills` (ARS) | Partner: full academic pipeline (Mode 2/3 delegate to it) |
| `ai-review-revision` | Complementary: revise-side for drafts with quality issues |
| Any registered biomedical skill | Distillation target: absorb via Mode 4 |
