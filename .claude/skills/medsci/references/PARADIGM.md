# Paradigm Capture: Style Anchoring for Biomedical Writing

## Why Paradigm Capture

LLMs trained on the entire web absorb every register of biomedical writing and have no built-in preference between:
- A Nature Reviews Immunology flagship (2 heading levels, no numbers, verdict sentences, Boxes for equations)
- A mid-tier survey paper (numbered chapters, neutral catalogue, hedging everywhere, dense flat subsections)
- A student term paper (long introductions, excessive background, weak conclusions)

Without an exemplar anchor, the draft drifts toward "generic survey paper." Paradigm capture locks in the style BEFORE writing begins.

---

## Phase 0 Process (Research Mode: Quick; Review Mode: Full)

### Research Mode (Mode 1) — Quick Paradigm
1. Identify **1 exemplar** from the target tier (e.g., a Nature Briefings piece, a Cell Leading Edge preview)
2. Extract: heading structure, paragraph rhythm, citation density
3. Write 5-10 bullet spec to PARADIGM.md
4. Time: ~30 min

### Review Mode (Mode 3) — Full Paradigm
1. Identify **2-3 exemplars** from the target journal tier
2. Read carefully, extract full style spec
3. Write structured spec to PARADIGM.md
4. Time: 3-4 hours (highest-ROI time investment in the project)

---

## Exemplar Selection Criteria

- Same domain or same problem family as your topic
- Same journal tier as your target
- Published within the last 3 years
- Authored by recognized senior figures

---

## Suggested Exemplars by Journal Tier

### Tier 1: Flagship Reviews (Nature Reviews, Annual Reviews, Cell)

| Journal | Exemplar Suggestion |
|---------|-------------------|
| Nature Reviews Immunology | Any recent review on T-cell biology, cancer immunology, or computational immunology |
| Nature Reviews Cancer | Neoantigen / immunotherapy review |
| Annual Review of Immunology | MHC biology, antigen presentation |
| Nature Immunology | Perspective/review on computational methods in immunology |
| Cell | Leading Edge reviews on immunotherapy |

### Tier 2: Specialty Top-Tier

| Journal | Exemplar Suggestion |
|---------|-------------------|
| Cancer Discovery | Translational immunology reviews |
| Immunity | Resource articles + reviews |
| Journal of Experimental Medicine | Mechanistic reviews |
| Science Immunology | Focused reviews |
| Nature Machine Intelligence | Methods reviews (e.g., MUNIS paper) |

### Tier 3: Bioinformatics / Methods

| Journal | Exemplar Suggestion |
|---------|-------------------|
| Briefings in Bioinformatics | Survey-style reviews (note: more structured, often numbered) |
| Bioinformatics | Application notes + original articles |
| PLOS Computational Biology | Methods + benchmark papers |
| BMC Bioinformatics | Database + tool papers |

### Tier 4: Immunology Specialty

| Journal | Exemplar Suggestion |
|---------|-------------------|
| Journal of Immunology | Classic reviews |
| Frontiers in Immunology | Open-access reviews (variable quality — choose carefully) |
| Cancer Immunology Research | Translational immunology |

---

## What to Extract (Style Spec Template)

### Heading Structure
- Max depth: 2 levels (H2 + H3) or 3 levels?
- Numbered headings? (most flagship reviews: NO)
- Typical H2 count: 5-8
- Typical H3 count per H2: 2-5

### Paragraph Rhythm
- Opening pattern: topic sentence → context → specific
- Body pattern: claim → evidence → transition
- Closing pattern: summary → verdict → forward look
- Typical length: 4-8 sentences (Nature Reviews); longer for Annual Reviews

### Citation Density & Style
- Refs per paragraph: 1.5-2.5 (flagship); 3-5 (survey)
- Multi-citation cap: 3-4 max in one bracket
- Citation placement: claim-level (every claim attributed) or paragraph-level (attributed at end)?
- Reference count: cite what supports the argument; no artificial target

### Equation Handling
- Display equations in body? → usually NO (Boxes)
- Textbook formulas not displayed → described in prose
- Which formulas are novel enough to display? → only those with methodological insight

### Vendor / Database / Tool Names
- In body text? → usually NO (confined to tables)
- In tables? → YES with validation evidence
- How to reference tools in body → category descriptors ("pan-specific MHC-I predictors"), not brand names

### Authorial Voice
- Hedging frequency: low at flagship tier (state strongly when evidence permits)
- Verdict frequency: 1-2 per major section, clustered at section ends
- "We" vs passive: varies by journal (Nature Reviews: "we" uncommon; Frontiers: "we" common)
- Verdict templates to follow: "[Method family] is currently the most effective approach for [problem]." / "[Method] has yet to demonstrate clear advantage over [alternative]."

### Tables & Boxes
- Table count: 2-4 (rarely more)
- Box count: 1-3
- Key Points box: present? 4-5 bullets?
- Figure count: 3-5

---

## Anti-Patterns (Always Flag)

These patterns signal LLM default style, not flagship-review style. Flag them when reading exemplars so you know what to AVOID:

- Numbered section headings (1., 1.1, 1.2.3)
- H4 (`####`) subsubsections in body
- "has shown promising results" / "may suggest" / "interestingly" / "it is worth noting that"
- Neutral catalogue: 3+ paragraphs without a verdict or authorial position
- Equations inline in body paragraphs
- Vendor/tool names scattered in body text (vs confined to a table)
- >4 citations in a single bracket
- Every paragraph ending with a citation (some should end with verdict)

---

## Quick Spec (Research Mode Template)

```markdown
## Quick Paradigm: <exemplar>

**Source**: <journal> <year>, <first author> et al.
**Headings**: <max depth, numbered? yes/no>
**Citation density**: <refs per paragraph>
**Voice**: <hedging level: low/medium/high>
**Verdicts**: <present? where?>
**Key pattern to replicate**: <one thing to get right>
**Key pattern to avoid**: <one anti-pattern to watch for>
```

---

## Full Spec (Review Mode Template)

```markdown
## Full Paradigm: <2-3 exemplars>

### Exemplar 1: <citation>
### Exemplar 2: <citation>
### Exemplar 3: <citation> (optional)

### Heading Structure
### Paragraph Rhythm
### Citation Density & Style
### Equation Handling
### Tool/Database/Vendor Handling
### Authorial Voice
### Table & Box Patterns
### Figure Patterns
### Unique Conventions

### Composite Spec (merged from all exemplars)
### Anti-Pattern Checklist (from HALLUCINATION_PATTERNS.md § relevant patterns)
```

---

## Critical Practice

**Re-read PARADIGM.md before each Phase 4 (writing) session.** Phase 4 is multi-day; drift back to LLM-default style is the main risk. The spec you write here only helps if you actually consult it during writing.

---

## ARS Style Calibration Integration

ARS defines a `shared/style_calibration_protocol.md` with 6 dimensions:
1. **Register**: formal / semi-formal / accessible
2. **Hedging**: heavy / moderate / minimal
3. **Citation style**: dense / moderate / sparse
4. **Paragraph length**: long / medium / short
5. **Heading structure**: deep / shallow
6. **Narrative stance**: neutral / opinionated

When MedSci routes to ARS (Mode 2/3), the PARADIGM.md spec should be translated into these 6 dimensions for ARS's Style Calibration consumption.

---

## Related Files

| File | Relationship |
|------|-------------|
| `DOMAINS.md` § A.4 | Journal tier suggestions for exemplar selection |
| `HALLUCINATION_PATTERNS.md` | Anti-patterns to avoid while writing |
| `TEMPLATES.md` | Templates that embody the paradigm |
