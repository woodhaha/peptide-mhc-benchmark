# Templates: Project Files & Output Structures

## 1. CLAUDE.md Template (Project Init)

```markdown
# <Project Name> — MedSci CLAUDE.md

## Project Info
- **Topic**: <topic>
- **Domain**: <from DOMAINS.md Section A or B>
- **Mode**: Research | Review
- **Target Journal** (review mode): <journal> | <tier>
- **Paradigm**: PARADIGM.md linked
- **Citation Integrity**: CITATION_INTEGRITY.md linked
- **Hallucination Self-Check**: HALLUCINATION_PATTERNS.md linked

## Terminology Standardization
| Term | Definition | Source |
|------|-----------|--------|
| <term> | <definition> | <from DOMAINS.md or project-specific> |

## Canonical Paper/Dataset Descriptions
(Prevent Pattern 9: Internal Inconsistency)
| Paper/Dataset | Canonical Description |
|--------------|----------------------|
| <paper> | <exact authors, year, journal, key finding> |
| <dataset> | <exact name, version, N, access> |

## Literature Inventory
Organized by taxonomic axes (from DOMAINS.md):
### Axis 1: <Algorithmic Priors>
- <family>: <papers> (N)
### Axis 2: <Data Priors>
- <family>: <papers> (N)
### Axis 3: <Validation Regime>
- <family>: <papers> (N)

## Verdict Positions (3-5)
1. <position>
2. <position>
...

## Writing-Time Guardrails
- Max heading depth: 2 (H2 + H3)
- No numbered headings
- Equations → Box 1 (not body)
- Tool/database names → Table (not body)
- Verdict sentences: ≥1 per major section
- Citation density: 1.5-2.5 refs per paragraph
```

---

## 2. IMPLEMENTATION_PLAN.md Template

```markdown
# Implementation Plan: <Project>

## Mode: Research | Review
## Date Started: <YYYY-MM-DD>

## Phase Checklist

### Phase 0: Scope & Paradigm
- [ ] Identify exemplar(s)
- [ ] Write paradigm spec
- [ ] Lock terminology table in CLAUDE.md

### Phase 1: Broad Discovery
- [ ] WebSearch: <query 1> (N results)
- [ ] WebSearch: <query 2> (N results)
- [ ] S2/PubMed/arXiv: <query> (N results)
- [ ] Build initial corpus (N papers)

### Phase 2: Collect, Filter & Verify
- [ ] Filter by relevance (N → M papers)
- [ ] Verify metadata (all M papers)
- [ ] Build literature matrix
- [ ] Gap analysis

### Phase 3: Deep Read & Taxonomy
- [ ] Deep read top 3-5 papers
- [ ] Classify all papers to 3 axes
- [ ] Draft taxonomy outline

### Phase 4: Synthesize & Write
- [ ] Write section: <section>
- [ ] Per-claim verification
- [ ] Verdict sentences placed
- [ ] Self-check every 5-6 paragraphs

### Phase 5: Citation Graph & Gaps (Review Mode)
- [ ] Forward citations for key papers
- [ ] Backward references for key papers
- [ ] Gap fill: negative results, recent preprints, diversity studies

### Phase 6: Quality Audit
- [ ] Citation integrity (all 5 rules)
- [ ] Hallucination self-check (11 patterns)
- [ ] Structural discipline (heading depth, verdicts, boxes)
- [ ] Voice and register
- [ ] Domain-specific checks

## Change Log
| Date | Phase | Change | Reason |
|------|-------|--------|--------|
```

---

## 3. Research Brief Template (Mode 1 Output)

```markdown
# Research Brief: <Topic>
**Date**: <YYYY-MM-DD> | **Mode**: Autonomous Research | **Sources**: <N> papers, <M> web

## Key Takeaways (2-4 bullets)
- <finding>
- <finding>

## Landscape Summary
### Major Approaches
| Family | Key Methods | Maturity |
|--------|------------|----------|
| <family> | <tools/papers> | mature / active / emerging |

### Consensus Findings
- <finding> — supported by [N1, N2, N3]

### Disagreements & Controversies
- <issue> — [N4] claims X, but [N5] finds Y

## Key Papers (3-5)
| # | Paper | Key Finding | Relevance |
|---|-------|-------------|-----------|
| 1 | <citation> | <one-sentence> | <why important> |

## Gaps & Opportunities
- <gap>

## References
[N1] <full citation>
[N2] <full citation>
...
```

---

## 4. Review Manuscript Template (Mode 3 Output)

```markdown
# <Title>: <Evocative Subtitle>

## Key Points
- <bullet 1>
- <bullet 2>
- <bullet 3>
- <bullet 4>

## Abstract
<structured abstract>

## Introduction
### Clinical / Biological Background
### Technical Challenge
### Scope and Contributions

## Datasets and Evaluation Metrics
(Table 1: public datasets)
(Box 1: evaluation metrics with formulas)

## Methods  ← 3-axis grouping (NOT flat list)
### <Axis 1: Algorithmic Priors>
**<Family A>.** ... **<Family B>.** ...
(Verdict sentence at end)

### <Axis 2: Data Priors>
**<Family A>.** ... **<Family B>.** ...
(Verdict sentence at end)

### <Axis 3: Validation Regime>
**<Family A>.** ... **<Family B>.** ...
(Verdict sentence at end)

(Table 2: representative methods)

## <Application / Downstream 1>
## <Application / Downstream 2>

## Translation to Practice
(Table 3: tools / products / regulatory status)

## Outstanding Challenges

## Future Directions

## References
```

---

## 5. Distillation Output Template (Mode 4 Output)

```markdown
# Distillation Report: <Source Skill>
**Date**: <YYYY-MM-DD> | **Source**: <skill-name> v<version>

## Extracted Fragments
| Type | Count | Confidence |
|------|-------|-----------|
| terminology | N | verified: N, extracted: N, tentative: N |
| data-source | N | ... |
| ... | ... | ... |

## Merge Results
- Auto-merged: N fragments → DOMAINS.md Section C
- Flagged for conflict: N fragments (see below)
- Tentative (stored, not applied): N fragments

## Conflicts
| Term | Source A | Source B | Resolution |
|------|---------|----------|------------|
| <term> | <def from DOMAINS> | <def from distill> | <kept / merged / needs review> |

## Next Steps
- [ ] Review tentative fragments
- [ ] Resolve conflicts
- [ ] Promote verified fragments to Section A (if domain is mature)
```

---

## Related Files

| File | Relationship |
|------|-------------|
| `WORKFLOW.md` | Phase structure that templates support |
| `DOMAINS.md` § A.1 | Terminology table to pull from for CLAUDE.md |
| `PARADIGM.md` | Style spec that drives writing templates |
| `CITATION_INTEGRITY.md` | Rules referenced in IMPLEMENTATION_PLAN.md |
