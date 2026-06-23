# Workflow: 4-Mode Unified Workflow

## Mode Selection (First Action)

Before executing any phase, detect the user's intent and select the mode:

```
User says "research X", "investigate", "find papers", "what is known" 
  → MODE 1: Autonomous Research (Phases 0-6, lightweight)
  
User says "write a paper", "write a review", "综述", "投稿", "revise", "审稿"
  → MODE 2: ARS Delegation (Route to appropriate ARS command)
  
User says "comprehensive review", "systematic review", "投稿级综述"
  → MODE 3: Deep Review (ARS pipeline + MIR quality layers)
  
User says "distill from X", "absorb X", "learn from X"
  → MODE 4: Quick Distill (DISTILLATION_PROTOCOL.md)
```

---

## Mode 1: Autonomous Research — 7-Phase Pipeline

For lightweight literature scans and background research. No ARS dependency.

### Phase 0: Scope & Paradigm
`[AUTOMATIC]` | Time: 5-30 min

**Goal**: Scope the problem and lock in terminology before searching.

**Actions**:
1. Parse the research question. Identify: domain, key terms, expected scope (broad/narrow), output format.
2. If topic maps to a domain in DOMAINS.md Section A or C → load terminology table.
3. If review-quality output expected → identify 1 exemplar (PARADIGM.md Quick Spec).
4. Note known unknowns (aspects you expect to NOT find).

**Deliverable**: 3-5 bullet scope + terminology loaded.

**Quality Gate**: User understands scope. If ambiguous, clarify with 1 question (not 3).

---

### Phase 1: Broad Discovery
`[AUTOMATIC]` | Time: 5-15 min

**Goal**: Map the landscape across multiple sources.

**Actions**:
1. Determine available tool tier (TOOL_STRATEGY.md).
2. Launch 2-3 parallel searches with complementary angles:
   - Angle 1: "X review 2024 2025" → find recent reviews
   - Angle 2: "X deep learning benchmark" → find methods papers
   - Angle 3: "X clinical application" → find translational work
3. If biomedical: add PubMed/E-utilities search.
4. Scan titles, snippets, citation counts. Identify:
   - Top 3-5 most cited papers
   - Most recent high-impact papers (last 12 months)
   - Existing reviews (mark as secondary sources — cite but don't treat as primary)
5. Build initial corpus: N papers (N = 15-30 for research, 50-100 for review).

**Quality Gates**:
- [ ] At least 2 search angles covered
- [ ] Both recent (2024-2025) and foundational papers identified
- [ ] Existing reviews marked as secondary

---

### Phase 2: Collect, Filter & Verify
`[AUTOMATIC]` | Time: 10-30 min (research), 1-3 hours (review)

**Goal**: Narrow the corpus and verify metadata on entry.

**Actions**:
1. **Filter**: Score each paper from Phase 1 on relevance (1-5). Keep ≥3. Drop duplicates.
2. **Verify metadata** (CITATION_INTEGRITY.md Rules 1, 2, 5):
   - Does the DOI/PMID resolve? (Rule 1)
   - Do first and last author names match the source? (Rule 2)
   - Is this a peer-reviewed source? (Rule 5)
3. **Categorize**: Tag each paper by method family, data type, and validation regime (from DOMAINS.md taxonomy).
4. **Gap analysis** (review mode only):
   - Any method families missing?
   - Negative results or null findings?
   - Recent preprints not yet captured?
   - Demographic/allele-diversity studies?
5. **Build literature matrix**: Table organized by taxonomic axes.

**Quality Gates**:
- [ ] All kept papers have verified metadata
- [ ] Filter criteria are explicit (not "I just picked these")
- [ ] Gap analysis done (review mode)

---

### Phase 3: Deep Read & Taxonomy
`[AUTOMATIC]` | Time: 15-45 min (research), 1-3 days (review)

**Goal**: Read deeply and organize findings thematically.

**Actions**:
1. **Select for deep read**: Top 3-5 papers (research mode) or top 8-15 papers (review mode). Selection criteria: citation count, recency, methodological novelty, direct relevance.
2. **Deep read each paper**:
   - Read abstract + results (and methods if methodological focus)
   - Extract: exact architecture name, exact performance numbers, exact dataset, exact validation scope
   - Note in working notes, NOT in manuscript yet (MIR discipline: read-first-write-after)
3. **Classify to taxonomic axes** (DOMAINS.md):
   - Each paper → primary axis + secondary axis (if applicable)
   - Identify method families within each axis
   - Note which families are well-represented vs sparse
4. **Identify verdict candidates**: Which families show clear advantage? Where is evidence contested?

**Quality Gates**:
- [ ] At least 3 papers deep-read (research) / 8 papers (review)
- [ ] All extracted numbers come from actual paper text, not memory
- [ ] Taxonomy is not flat — 3-axis structure applied
- [ ] Papers that can't be accessed are noted; no fabricated details about them

---

### Phase 4: Synthesize & Write
`[AUTOMATIC]` | Time: 5-15 min (research), 2-5 days (review)

**Goal**: Produce the output artifact with per-claim verification.

**Research Mode**:
1. Fill TEMPLATES.md § Research Brief Template.
2. Every claim attributed to a paper → verify the paper supports it.
3. Cross-reference web sources against academic sources.
4. Cite all claims (include paper titles and URLs).

**Review Mode** (follows MIR Phase 4 discipline):
1. Write section by section. For each section:
   - Write introduction paragraph
   - For each method family: re-read the cited paper → write 2-4 sentences with actual modules and numbers → verify citation (body-bibliography, number accuracy, direction)
   - Close with verdict sentence
2. **Per-paragraph self-check** (every 5-6 paragraphs): scan for HALLUCINATION_PATTERNS.md patterns.
3. Equations → Boxes, not body.
4. Tool/database names → Tables, not body.
5. Verdict sentences: 3-5 across the whole manuscript.
6. Update bibliography as you go (no deferred bibliography population).

**Quality Gates**:
- [ ] Every cited claim traceable to a source
- [ ] No template-filling (read-first-write-after discipline applied)
- [ ] Verdict sentences present (review mode)
- [ ] Self-check run every 5-6 paragraphs

---

### Phase 5: Citation Graph & Gaps
`[AUTOMATIC]` | Time: 5-15 min (research), 0.5-1 day (review)

**Research Mode** (lightweight):
1. For top 2-3 key papers: check who cites them (forward) and what they build on (backward).
2. Note any missed high-impact papers → add to references.

**Review Mode** (full):
1. **Forward citations**: For each key paper, find 5-10 citing papers. Any important follow-up work missed?
2. **Backward references**: For each key paper, check the reference list. Any foundational paper missed?
3. **Gap fill**: Targeted searches for: negative results, cross-site reproducibility, demographic-bias studies, recent 3-month preprints.
4. **Dedup**: Merge new findings into existing literature matrix.
5. **Re-verify cornerstone references** (CITATION_INTEGRITY.md Rules 1, 3) — citation graph exploration may surface errors.

**Quality Gates**:
- [ ] Citation graph explored for ≥3 key papers
- [ ] Gap fill searches executed
- [ ] Cornerstone references re-verified

---

### Phase 6: Quality Audit
`[AUTOMATIC]` | Time: 10-30 min (research), 4-8 hours (review)

**Goal**: Run the complete quality checklist before delivering.

**Actions**:
1. Run QUALITY_CHECKLIST.md for the active mode.
2. Execute all self-check commands.
3. Fix CRITICAL issues immediately, HIGH issues before delivery, MEDIUM issues if time permits.

**Research Mode**: Run research checklist (~10 min).
**Review Mode**: Run full checklist + ARS 7-mode failure crosswalk (~4-8 hours, can delegate to multi-agent review). If ≥5 hard factual errors or ≥10 citation drift instances → apply ARS revision workflow.

**Quality Gates**:
- [ ] No CRITICAL issues remaining
- [ ] All HIGH issues resolved or documented
- [ ] QUALITY_CHECKLIST.md self-check commands pass

---

### Deliver (End of Mode 1)

- **Research Mode**: Save to `outputs/<slug>_research_brief.md`
- **Review Mode**: Save to `manuscript_draft.md` + bibliography

---

## Mode 2: ARS Delegation

**Do not execute Mode 1 phases.** Instead, route directly to the appropriate ARS slash command.

### Routing Decision Tree

```
User intent contains:
├── "plan" / "规划" / "chapter plan" → /ars-plan
├── "outline" / "大纲" → /ars-outline
├── "literature review" / "lit review" (paper format) → /ars-lit-review
├── "abstract" / "摘要" → /ars-abstract
├── "write a paper" / "投稿" / "publish" / "full pipeline" → /ars-full
├── "review my paper" / "peer review" / "审稿" → /ars-reviewer
├── "respond to reviewers" / "回复审稿人" → /ars-revision-coach
├── "revise" / "修改稿" → /ars-revision
├── "check citations" / "检查引文" → /ars-citation-check
├── "convert format" / "格式转换" → /ars-format-convert
├── "AI disclosure" / "AI声明" → /ars-disclosure
└── "mark as read" → /ars-mark-read
```

### MedSci Domain Injection

Before routing to ARS, inject domain context:

1. Load terminology from DOMAINS.md § A.1
2. Load canonical paper descriptions from DOMAINS.md § A.6
3. Load paradigm from PARADIGM.md (if review mode)
4. Pass as context: "This is a biomedical/immunology paper. Use these standard terms: [list]. Use these canonical citations: [list]. Target journal style: [from PARADIGM.md]."

---

## Mode 3: Deep Review (ARS + MIR Quality)

**Goal**: Produce a publication-quality review with the full ARS pipeline PLUS MIR quality layers.

**Flow**:
1. **Execute ARS full pipeline** (`/ars-full`) — Stages 1-6
2. **At Stage 2.5 (integrity gate)**: Apply MedSci's CITATION_INTEGRITY.md 5-rule check IN ADDITION to ARS's 7-mode checklist
3. **At Stage 3 (review)**: MedSci provides domain-specific reviewer guidance: check for allele specificity overclaims (Pattern 11), database version drift (Pattern 10), binding vs immunogenicity conflation
4. **At Stage 4 (revise)**: Run HALLUCINATION_PATTERNS.md self-check on the revised draft
5. **At Stage 4.5 (final integrity)**: Run full QUALITY_CHECKLIST.md review mode
6. **Pre-submission**: Verify all 5 citation integrity rules, all 11 hallucination patterns, all structural discipline rules

**Quality Gates**: Both ARS gates (2.5, 4.5) AND MedSci gates (Phase 6) must pass.

---

## Mode 4: Quick Distill

**Goal**: Absorb knowledge from another skill.

**Flow**: Execute DISTILLATION_PROTOCOL.md distill loop:
1. Locate source skill SKILL.md
2. Extract fragments (terminology, data-sources, benchmarks, etc.)
3. Check duplicates against DOMAINS.md Section C
4. Merge: auto-merge low-conflict, flag conflicts
5. Report: summarize extracted fragments

---

## Related Files

| File | Referenced In Phase |
|------|-------------------|
| `TOOL_STRATEGY.md` | Phase 1 (tool detection) |
| `CITATION_INTEGRITY.md` | Phase 2 (verify), Phase 4 (verify), Phase 5 (re-verify), Phase 6 (audit) |
| `HALLUCINATION_PATTERNS.md` | Phase 4 (self-check), Phase 6 (audit) |
| `DOMAINS.md` | Phase 0 (scope), Phase 3 (taxonomy), Mode 2 (domain injection) |
| `PARADIGM.md` | Phase 0 (style anchor) |
| `TEMPLATES.md` | Phase 4 (output structure) |
| `QUALITY_CHECKLIST.md` | Phase 6 (audit) |
| `SKILL_CATALOG.md` | Mode 2 (ARS routing) |
| `DISTILLATION_PROTOCOL.md` | Mode 4 (distill) |
