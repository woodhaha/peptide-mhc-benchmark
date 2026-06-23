# Skill Catalog: ARS Commands & Related Skills Quick Reference

## ARS Plugin (academic-research-skills v3.11.1)

### 14 Slash Commands

| Command | Model | Purpose | Mode 2 Route |
|---------|-------|---------|-------------|
| `/ars-full` | Opus | Full pipeline: research → write → review → revise → finalize | "write a paper", "投稿", "publish" |
| `/ars-plan` | Sonnet | Socratic chapter-by-chapter planning | "plan my paper", "规划论文结构" |
| `/ars-outline` | Sonnet | Detailed outline + evidence map | "outline", "大纲" |
| `/ars-lit-review` | Sonnet | Annotated bibliography in paper format | "literature review for paper" |
| `/ars-abstract` | Sonnet | Bilingual abstract + keywords | "write abstract", "摘要" |
| `/ars-reviewer` | Opus | Simulated peer-review panel (EIC + 3 reviewers + DA) | "review my paper", "审稿" |
| `/ars-revision-coach` | Opus | Parse reviewer comments → Revision Roadmap + Response Letter | "respond to reviewers", "回复审稿人" |
| `/ars-revision` | Sonnet | Revised draft + point-by-point R&R responses | "revise manuscript", "修改稿" |
| `/ars-citation-check` | Sonnet | Citation error report | "check citations", "检查引文" |
| `/ars-format-convert` | Sonnet | Convert between LaTeX / DOCX / PDF / Markdown | "convert format", "格式转换" |
| `/ars-disclosure` | Sonnet | Venue-specific AI-usage disclosure statement | "AI disclosure", "AI声明" |
| `/ars-mark-read` | Sonnet | Record human-read signal for citation keys | "mark as read" |
| `/ars-unmark-read` | Sonnet | Rescind prior human-read mark | "unmark" |
| `/ars-cache-invalidate` | Sonnet | Drop cached verification entries | (internal) |

### 4 Core Skills

| Skill | Version | Agents | Modes | Layer |
|-------|---------|--------|-------|-------|
| **deep-research** | v2.9.4 | 14 | 7 (full, quick, socratic, review, lit-review, fact-check, systematic-review) | raw (unverified sources) |
| **academic-paper** | v3.2.0 | 12 | 10 (full, plan, outline-only, revision, revision-coach, abstract-only, lit-review, format-convert, citation-check, disclosure) | redacted (sanitized materials) |
| **academic-paper-reviewer** | v1.10.0 | 7 | 6 (full, re-review, quick, methodology-focus, guided, calibration) | verified_only |
| **academic-pipeline** | v3.11.1 | 5 | 2 (orchestrator, resume) | verified_only |

### 10-Stage Full Pipeline

```
Stage 1  RESEARCH     (deep-research, 14 agents)
Stage 2  WRITE        (academic-paper, 12 agents)
Stage 2.5 INTEGRITY   (integrity_verification_agent — MANDATORY)
Stage 3  REVIEW       (academic-paper-reviewer, 7 agents)
Stage 3→4 COACHING   (EIC Socratic, max 8 rounds)
Stage 4  REVISE       (academic-paper revision)
Stage 3' RE-REVIEW    (narrow reviewer team)
Stage 3'→4' COACHING (max 5 rounds)
Stage 4' RE-REVISE   (final, content frozen)
Stage 4.5 FINAL_INTEGRITY (100% claim verification)
── opt-in ──
Stage 4→5 CLAIM-AUDIT (claim_ref_alignment_audit_agent, ARS_CLAIM_AUDIT=1)
Stage 5  FINALIZE     (format-convert → MD/DOCX/LaTeX/PDF)
Stage 6  PROCESS_SUMMARY (orchestrator auto)
```

---

## MedSci Project Skills (D:\Researching\Peptide epitope\.claude\skills\)

| Skill | Path | Purpose |
|-------|------|---------|
| **medsci** | `.claude/skills/medsci/` | This skill — unified medical science research |

---

## User Global Skills (~/.claude/skills/)

| Skill | Purpose | Relevance to MedSci |
|-------|---------|-------------------|
| **research-skills** | Medical imaging review + paper slides + research proposal | medical-imaging-review integrated into MedSci references |
| **deep-research** | 5-phase literature research | Workflow merged into MedSci Mode 1 |
| **prompts/** | 14 research workflow templates (lit, review, draft, etc.) | Can be distilled into MedSci via Distillation Protocol |
| **paper-to-course** | Paper → interactive course conversion | Post-publication outreach |
| **notebooklm** | NotebookLM integration | AI-assisted literature digestion |

---

## Registered System Skills (Biomedical — Partial List)

The system has 400+ registered skills. Below are those most relevant to biomedical/immunology research:

### Epitope & Immunology
| Skill | Purpose |
|-------|---------|
| `neoantigen-predictor` | Neoantigen prediction from mutations |
| `adme-property-predictor` | ADME property prediction |
| `immune-pathway-analysis` | Immune pathway enrichment |
| `biomarker-landscape-scanner` | Biomarker landscape analysis |
| `biomed-outline-generator` | Biomedical paper outline generation |
| `biomedical-search-strategy-builder` | Search strategy construction |

### Clinical & Translational
| Skill | Purpose |
|-------|---------|
| `clinical-data-cleaner` | Clinical data preprocessing |
| `clinical-trial-info-extractor` | Trial information extraction |
| `real-world-evidence-study-designer` | RWE study design |
| `treatment-response-predictor-planner` | Treatment response prediction |

### Bioinformatics & Genomics
| Skill | Purpose |
|-------|---------|
| `differential-expression-analysis` | DEG analysis |
| `wgcna-analysis` | Co-expression network |
| `gsea` | Gene set enrichment |
| `variant-annotation` | Genetic variant annotation |
| `alphafold-db` | AlphaFold structure access |

### Literature & Writing
| Skill | Purpose |
|-------|---------|
| `literature-review` | Literature review generation |
| `systematic-review-screener` | Systematic review screening |
| `meta-analysis-methods-generator` | Meta-analysis methods |
| `journal-matchmaker` | Journal recommendation |

### External Research Tools & Benchmarks
| Resource | Type | Purpose | Status in MedSci |
|----------|------|---------|-----------------|
| **QUEST** (OSU-NLP) | Deep Research Agent | Multi-turn search+scholar+visit+python+memory agent with 4-stage training; 2B-35B models | ✅ Distilled 2026-06-20 (14 fragments → `deep-research-agents` domain) |
| BrowseComp | Benchmark | Fact-seeking via verifier-script objective eval | Tracked in DOMAINS.md §C |
| GAIA | Benchmark | Multi-step reasoning with ground-truth answers | Tracked in DOMAINS.md §C |
| DeepResearch Bench | Benchmark | Citation grounding assessment | Tracked in DOMAINS.md §C |

---

## Distillation Target Priority

These skills are prime candidates for distillation into MedSci's domain knowledge:

| Priority | Skill | Reason |
|----------|-------|--------|
| 🔴 High | `neoantigen-predictor` | Directly relevant to epitope prediction domain |
| 🔴 High | `immune-pathway-analysis` | Terminology + method taxonomies |
| 🔴 High | `biomarker-landscape-scanner` | Data source + benchmark knowledge |
| 🟡 Medium | `adme-property-predictor` | Adjacent domain — peptide drug properties |
| 🟡 Medium | `clinical-trial-info-extractor` | Validation regime knowledge |
| 🟢 Lower | `literature-review` | Method overlap with MedSci Mode 1/3 |
| 🟢 Lower | `journal-matchmaker` | Submission strategy knowledge |
| 🔵 **Done** | `QUEST` (OSU-NLP-Group) | ✅ Distilled 2026-06-20 — 14 fragments |
| 🔵 **Done** | `CCF-Figure` (Deepshare-Official) | ✅ Distilled 2026-06-20 — 26 fragments |
| 🔵 **Done** | `neoantigen-predictor` (AIPOCH) | ✅ Distilled 2026-06-20 — 17 fragments |
| 🔵 **Done** | `immune-pathway-analysis` (AIPOCH) | ✅ Distilled 2026-06-20 — 8 fragments |
| 🔵 **Done** | `biomarker-landscape-scanner` (AIPOCH) | ✅ Distilled 2026-06-20 — 6 fragments; 5-tier maturity framework + 17 hard rules |
| 🔵 **Done** | `adme-property-predictor` (AIPOCH) | ✅ Distilled 2026-06-20 — 10 fragments; 15+ ADME endpoints + Lipinski/QED/MPO |
| 🔵 **Done** | `sci-hub-search-skill` + `paper-search-mcp` | ✅ Integrated 2026-06-20 — 23 paper sources + OA-First Fallback Chain |

**How to distill**: `"medsci, distill from neoantigen-predictor"` → triggers DISTILLATION_PROTOCOL.md.

---

## Related Files

| File | Relationship |
|------|-------------|
| `DISTILLATION_PROTOCOL.md` | How to absorb knowledge from skills in this catalog |
| `DOMAINS.md` § C | Distilled knowledge is stored here |
| `WORKFLOW.md` | Mode 2 routes to ARS commands listed here |
| `SKILL.md` | Entry point that references this catalog for routing decisions |
