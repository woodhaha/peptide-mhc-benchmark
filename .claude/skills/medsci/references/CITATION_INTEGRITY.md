# Citation Integrity: 5-Rule Verification Protocol

Adapted from `medical-imaging-review` v3.0.0 for the biomedical/immunology domain. Every citation must satisfy all 5 rules before it enters the manuscript bibliography.

---

## Rule 1: No Placeholder DOIs

**Check**: Any DOI containing `xxx`, `[TBD]`, `?`, stub patterns, or "to-be-filled" must be resolved before the manuscript advances.

**Domain note**: Biomedical papers often have PMID alternatives when DOIs are slow to register (preprints, early online). Always prefer DOI when available; fallback to PMID.

**Detection**:
```bash
grep -i "xxx\|\[TBD\]\|placeholder\|to.be.filled\|10.0000/\|xxxx" manuscript.md
```

**Resolution**:
1. WebFetch `https://pubmed.ncbi.nlm.nih.gov/?term=<title>` → extract PMID → extract DOI
2. WebFetch `https://api.crossref.org/works?query.bibliographic=<title>&rows=3` → extract DOI
3. Zotero MCP: search by title → extract DOI/PMID
4. Semantic Scholar API: `/paper/search?query=<title>` → extract externalIds

**Why this matters**: In biomedical publishing, placeholders signal workflow incompleteness. Editors at immunology journals (Immunity, J Immunol, Cancer Res) reject on sight.

---

## Rule 2: Author List Verification

**Check**: Every reference must have first and last author verified against the first source (PubMed, Crossref, or journal page).

**Domain notes**:
- Immunology papers often have 10-30 authors. For >6 authors: list first 6 + "et al."
- First and last author are the most important positions in biomedical publishing (first = did the work, last = PI/supervisor)
- Wrong PI (last author) is instantly detected by reviewers in the field

**Detection**: Look for 4-author patterns with generic surnames (see HALLUCINATION_PATTERNS.md Pattern 6).

**Resolution**:
1. WebFetch `https://pubmed.ncbi.nlm.nih.gov/?term=<title>` → extract full author list
2. WebFetch `https://api.crossref.org/works/<DOI>` → extract author array
3. Copy verbatim from source. Do not re-order, truncate creatively, or fill gaps.

---

## Rule 3: Body-Bibliography Reconciliation

**Check**: The `[N]` in the body must (a) exist in the bibliography and (b) be the paper the body sentence is actually attributing the claim to.

**Citation drift (the most common failure)**: When reorganizing sections, LLMs renumber or remap citations. The body says "Wohlwend et al. [43]" but bibliography [43] is Reynisson et al.

**Domain-specific drift patterns**:
- The same tool (NetMHCpan) cited from different papers in different sections without distinguishing them → body says [15] but the claim is from [23]
- pMHC binding predictors cited interchangeably (NetMHCpan vs MHCflurry vs MixMHCpred) → body attributes a feature to the wrong tool

**Resolution**:
- **Pattern A** (single misattribution): Find the correct paper, replace [N].
- **Pattern B** (numerical drift): Spot-check 10 random [N] per section. Correct all.

**Expected drift**: 5-15 instances per 100+ reference manuscript if careful, 30-40 if not.

---

## Rule 4: Conclusion-Direction Verification

**Check**: For every cited finding (AUC, IC50, PPV, HR, OR, p-value, effect size), the body sentence's directional claim must match the source's stated direction.

**Domain-specific failure modes**:
- **IC50 direction**: Lower IC50 = stronger binding. "IC50 of 500 nM" is weak binding (threshold typically <500 nM for "binder"). Verify the stated IC50 actually supports "strong binder" or "weak binder" claim.
- **% Rank direction**: Lower % rank = stronger binding. Often confused with raw IC50.
- **AUC direction**: Higher AUC = better. But AUC 0.92 vs 0.95 may not be statistically significant — check if the paper claims superiority.
- **Immunogenicity vs binding**: "Binds with high affinity" ≠ "is immunogenic." Verify whether the cited paper tested immunogenicity or only binding.

**Fix**: Quote directional claims verbatim from the source abstract. Do not paraphrase.

---

## Rule 5: First-Source over Vendor/Database Materials

**Check**: Preference order for clinical/scientific claims:

1. **Peer-reviewed journal article** (preferred source for scientific claims)
2. **Preprint with community validation** (acceptable for recent findings, must be labeled "preprint")
3. **Official database publication** (e.g., NAR database issue — preferred for database statistics)
4. **Regulatory document** (FDA 510(k), EMA guidance — only for regulatory facts)
5. **Software documentation** (only for implementation details, not for scientific claims)
6. **Company/vendor white paper** (DO NOT cite for clinical/scientific claims)

**Domain-specific patterns**:
- ✗ "NetMHCpan 4.1 User Manual" cited for binding prediction accuracy → cite Reynisson et al. (2020) NAR
- ✗ "IEDB Statistics Page" cited for database size → cite the NAR database issue paper
- ✗ "MHCflurry GitHub README" cited for model architecture → cite O'Donnell et al. (2018, 2020)
- ✗ "FDA 510(k) clearance letter" cited for clinical validation → cite the actual clinical trial publication

**When a primary source can't be found**: Label the claim with `[source: <type>, unverified]`. Do not invent a journal citation.

---

## Verification Workflow by Phase

| Phase | Rules Checked | How |
|-------|--------------|-----|
| Phase 2 (Collect + Filter) | Rules 1, 2, 5 | Verify at entry — DOI resolves, first/last author match, correct source type |
| Phase 4 (Synthesize + Write) | Rules 3, 4 | Per-paragraph — body-bibliography [N] match, directional claim matches source |
| Phase 5 (Citation Graph + Gaps) | Rules 1, 3 | Re-verify cornerstone references after graph exploration |
| Phase 6 (Quality Audit) | All 5 | Full audit — check every reference against every rule |

---

## When a Rule Failure is Found

1. **Stop** forward writing immediately
2. **Look up** correct metadata (WebFetch Crossref/PubMed/S2)
3. **Fix** citation in place (body + bibliography)
4. **Check** for related failures (same error often appears 2-3 more times)
5. **Log** the fix in IMPLEMENTATION_PLAN.md change log

**Time cost**: This protocol takes 10-20% of total writing time but eliminates >80% of credibility-killer issues.

---

## Domain-Specific Verification Commands

```bash
# Check for placeholder DOIs in biomedical formats
grep -E "(xxx|\[TBD\]|10\.0000/|10\.XXXX)" manuscript.md

# Check for generic Chinese 4-name author pattern
grep -E "Zhang [A-Z], (Wang|Li|Liu|Chen) [A-Z], " manuscript.md

# Check for vendor/corporate authors cited as peer-reviewed
grep -E "(IEDB|NetMHCpan|MHCflurry|IEDB|Immune Epitope).*Nature|Science|Cell" manuscript.md

# Check for IC50 direction errors (high IC50 = weak binding)
# Any IC50 > 5000 nM claimed as "strong" or "high affinity" = suspect
grep -E "IC50.*[5-9][0-9]{3,}|IC50.*strong|IC50.*high affinity" manuscript.md

# Check for universal scope claims (Pattern 11)
grep -E "\b(all|every|any|universal) (peptide|MHC|HLA|epitope|allele)" manuscript.md
```

---

## Related Files

| File | Relationship |
|------|-------------|
| `HALLUCINATION_PATTERNS.md` | The 11 patterns that citation verification addresses |
| `TOOL_STRATEGY.md` § Verification Mapping | How to verify at each tool tier |
| `DOMAINS.md` § A.6 Canonical Papers | Locked descriptions for cross-section consistency (Pattern 9 prevention) |
| `QUALITY_CHECKLIST.md` § Citation Integrity | Verification checkboxes for quality gates |
