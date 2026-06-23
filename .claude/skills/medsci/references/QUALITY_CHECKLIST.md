# Quality Checklist: Mode-Aware Gates

## Severity Levels

| Level | Label | Action |
|-------|-------|--------|
| **CRITICAL** | Placeholder DOI, wrong author, citation direction flipped | STOP and fix immediately |
| **HIGH** | Body-bibliography drift, universal scope claim, verdict absent in axis section | Fix before completing current section |
| **MEDIUM** | Numbered heading, equation inline, vendor name in body | Fix at section end |
| **LOW** | Multi-citation bracket >4 refs, >1 database version stat | Note; address during audit phase |

---

## Research Mode Checklist (Mode 1)

### Citation Integrity
- [ ] No placeholder DOIs in references
- [ ] Every reference's first/last author verified (spot-check 30%)
- [ ] No vendor white papers cited as peer-reviewed
- [ ] No fabricated performance numbers (spot-check 3 key papers)

### Structural Discipline
- [ ] Research Brief has all required sections (Key Takeaways, Landscape, Key Papers, Gaps)
- [ ] 3-5 key papers deep-read
- [ ] Gaps section is substantive, not generic

### Voice
- [ ] No "has shown promising results" / "interestingly" / "it is worth noting that"
- [ ] Claims match evidence strength

### Self-Check Commands
```bash
grep -i "promising\|interestingly\|worth noting" research_brief.md
grep -c "^###" research_brief.md  # Should have all required H3 sections
```

---

## Review Mode Checklist (Mode 3)

### Citation Integrity (Hard Gate — Non-Negotiable)
- [ ] No placeholder DOIs (`grep -i "xxx\|\[TBD\]"`)
- [ ] Every reference's first and last author verified against first-source
- [ ] Body-bibliography [N] reconciliation — spot-check 10 random per section
- [ ] Every directional claim verified against source
- [ ] No vendor/tool documentation cited as if peer-reviewed
- [ ] No duplicate references
- [ ] All database version statistics current (WebFetch verified)

### Structural Discipline
- [ ] Heading depth ≤ 2 levels (H2 + H3 only in body)
- [ ] No numbered headings
- [ ] No H4 (`####`) — use `**Topic.**` bold lead-in
- [ ] Methods section uses 3-axis grouping
- [ ] Verdict sentences present in 3-5 places
- [ ] Key Points box: 4-5 bullets, 1-3 sentences each
- [ ] Standard sections all present
- [ ] Tables 1-3 present and correctly populated
- [ ] Box 1 (metrics) present

### Voice and Register
- [ ] No LLM-tell phrases anywhere
- [ ] Hedging used only when evidence supports caution
- [ ] Strong findings stated strongly
- [ ] No neutral catalogue >3 paragraphs without verdict
- [ ] Domain-specific terms used correctly (cross-check DOMAINS.md § A.1)

### Equations and Boxes
- [ ] Display equations only in Boxes (not body)
- [ ] All equation definitions verified against original papers
- [ ] Metric formulas correct (Pattern 8)

### Tables
- [ ] Table 1: Public datasets with year, cases, annotation, access
- [ ] Table 2: Method comparison (12-20 papers)
- [ ] Table 3: Tools/products with validation evidence
- [ ] Total tables ≤ 4
- [ ] Each table has title, body, and footnote

### Domain-Specific (Biomedical/Immunology)
- [ ] All allele names use 4-digit resolution minimum (HLA-A*02:01, not HLA-A2)
- [ ] Binding affinity units consistent (IC50 nM or % rank, not mixed)
- [ ] Immunogenicity vs binding claims not conflated
- [ ] Universal scope claims checked (Pattern 11)
- [ ] Peptide lengths explicitly stated (8-11mer for MHC-I, 13-25mer for MHC-II)

### Self-Check Commands
```bash
# Structural
grep -E "^####" manuscript.md  # Should be EMPTY
grep -E "^#+ [0-9]+\." manuscript.md  # Should be EMPTY (no numbered headings)

# Citations
grep -i "xxx\|\[TBD\]\|placeholder" manuscript.md  # Should be EMPTY
grep -E "Zhang [A-Z], (Wang|Li|Liu) [A-Z]" manuscript.md  # Check for generic patterns

# Voice
grep -i "promising results\|may suggest\|interestingly\|worth noting\|in recent years" manuscript.md

# Domain
grep -E "HLA-[A-C][^0-9]" manuscript.md  # Check for non-4-digit allele names
grep -E "\b(all|every|any|universal) (peptide|MHC|epitope)" manuscript.md  # Pattern 11

# Verdict sentences
grep -E "is currently the|has yet to|is best understood as|will determine whether" manuscript.md
# Should find 3-5 matches
```

---

## ARS Integration: 7-Mode Failure Gate (Mode 2/3 only)

When routed through ARS pipeline, the `integrity_verification_agent` runs these checks at Stages 2.5 and 4.5. Verify:

- [ ] M1 (implementation bug): N/A for writing tasks
- [ ] M2 (hallucinated citation): Covered by Citation Integrity rules 1-5
- [ ] M3 (hallucinated result): Covered by HALLUCINATION_PATTERNS.md Patterns 2, 3
- [ ] M4 (shortcut reliance): Covered by Pattern 8 (metric errors)
- [ ] M5 (bug as insight): Covered by Patterns 3, 7, 10
- [ ] M6 (methodology fabrication): N/A for review writing
- [ ] M7 (frame-lock): Covered by Patterns 9, 11

---

## Distillation Quality Gate (Mode 4)

- [ ] All extracted fragments have confidence level assigned
- [ ] Verified fragments: source citation included
- [ ] Extracted fragments: structurally sound, no contradictions with existing domain knowledge
- [ ] Tentative fragments: explicitly tagged, not auto-merged
- [ ] Conflicts documented in DISTILLATION_LOG.md
- [ ] DOMAINS.md Section C updated

---

## What We Intentionally Removed

Following medical-imaging-review v3's lessons:
- ❌ "Hedging language used" — hedging-by-default is the LLM default, not a quality signal
- ❌ "80-120 references" — count targets drive fabrication
- ❌ "All major methods covered" — drives exhaustive enumeration over selective synthesis
- ❌ ">50% from last 3 years" — date-based filter has no relationship to quality
- ❌ "Performance metrics consistent" — checking format is not checking correctness

---

## Related Files

| File | Relationship |
|------|-------------|
| `CITATION_INTEGRITY.md` | Rules checked in Citation Integrity section |
| `HALLUCINATION_PATTERNS.md` | Patterns checked in Voice and Domain sections |
| `DOMAINS.md` § A.1 | Terminology table for domain checks |
| `PARADIGM.md` | Style spec for structural discipline checks |
