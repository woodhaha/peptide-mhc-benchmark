# Hallucination Patterns: 11+7 Cross-Referenced Checklist

## Design Note

This file adapts the 9 hallucination patterns from `medical-imaging-review` v3.0.0 for the biomedical/immunology domain, adds 2 domain-specific patterns, and cross-references each against the ARS 7-Mode AI Research Failure Checklist (Lu 2026, integrated in ARS v3.11.1 via `integrity_verification_agent`).

---

## Patterns 1-9: Adapted from Medical-Imaging-Review

All examples have been replaced with biomedical/immunology equivalents. Cross-reference format: `→ ARS M<N>` maps to the corresponding ARS failure mode.

### Pattern 1: Real Paper, Fabricated Author List

**Description**: The cited paper exists but the author list is a generic pattern or mixes real with fabricated names.

**Domain examples**:
- ✗ "Zhang H, Liu W, Chen X, Li J. Nature Immunology. 2024." (generic 4-name Chinese pattern)
- ✗ "Smith R, Williams T, Brown A, Jones M. Immunity. 2025." (generic 4-name English pattern)
- ✓ "Wohlwend J, Nathan A, Shalon N, Crain CR, Tano-Menka R, Goldberg B, Richards E, Gaiha GD, Barzilay R. Nat Mach Intell. 2025." (actual MUNIS author list from source)

**Detection signals**:
- Exactly 4 authors — the most common fabrication pattern
- All 4 surnames are among the top 20 most common in a language
- Author initials all single-letter (A, B, C, D pattern)
- PI-level author count inconsistent with journal (e.g., 4 authors in Nature Immunology)

**Fix**: WebFetch PubMed/Crossref → copy entire author list verbatim from source.

→ **ARS M2**: Hallucinated citation

---

### Pattern 2: Inflated Performance Numbers

**Description**: The paper is real but specific AUC/IC50/PPV numbers are fabricated and inflated.

**Domain examples**:
- ✗ "NetMHCpan achieved AUC of 0.978 on IEDB benchmark" (fabricated 3-decimal number)
- ✗ "Model achieved IC50 prediction R² = 0.94" (suspiciously high)
- ✓ Quote the actual numbers from the paper's abstract or results section

**Detection signals**:
- Numbers ending in clean 3-decimal patterns (`.923`, `.978`, `.941`)
- Performance consistently above published field benchmarks
- Round numbers when papers typically report messy values (IC50 varies by orders of magnitude)
- AUC >0.98 claimed without specifying test set or allele restriction

**Fix**: Open the actual paper and quote the specific number from the results section. If paper is inaccessible, drop the specific numeric claim.

→ **ARS M2**: Hallucinated experimental result (direct mapping)

---

### Pattern 3: Conclusion-Direction Flip

**Description**: Real paper, real result — but the direction (higher/lower, binder/non-binder, immunogenic/non-immunogenic) is reversed.

**Domain examples**:
- ✗ "The G12V mutation decreased MHC binding affinity" — paper actually shows increased binding
- ✗ "ESM-2 embeddings underperformed BLOSUM62 features" — paper shows the opposite
- ✗ "The model was allele-specific" — paper demonstrates pan-allele generalization

**Detection**: Any sentence with directional language (increased/decreased, higher/lower, improved/degraded, specific/pan) + a citation bracket → verify direction against source.

**Fix**: Quote directional claims verbatim from source abstract. Do not paraphrase quantitative directional claims.

→ **ARS M2 + M5**: Hallucinated result + possibly reframed as insight (if the direction error benefits the narrative)

---

### Pattern 4: Vendor/Agency Material Cited as Peer-Reviewed

**Description**: Company white papers, FDA 510(k) letters, or agency reports cited with fabricated journal attribution.

**Domain examples**:
- ✗ "NetMHCpan 5.0 Technical Report. Nature Biotechnology. 2024." (software docs, not peer-reviewed)
- ✗ "FDA. Guidance on Neoantigen Prediction. N Engl J Med. 2025." (FDA doesn't publish in NEJM)
- ✗ "IEDB. Database Statistics. Bioinformatics. 2024." (website not a journal article)

**Detection signals**: The "author" is a company, database, or regulatory agency. The "journal" is top-tier but the claim is trivial. No PMID or DOI resolves.

**Fix**: Cite official publications only. For software: cite the method paper, not docs. For databases: cite the NAR database issue paper. For FDA: cite the actual guidance document ID, not a journal.

→ **ARS M2**: Hallucinated citation

---

### Pattern 5: Placeholder DOIs

**Description**: `xxx`, `[TBD]`, `x):xxx-xxx` stubs left in bibliography.

**Detection**: `grep -i "xxx\|\[TBD\]\|placeholder\|to.be.filled" manuscript.md`

**Fix**: WebFetch PubMed by PMID, or Crossref API by DOI, or Zotero MCP.

→ **ARS M2**: Hallucinated citation

---

### Pattern 6: Generic 4-Author Hallucination (Subset of Pattern 1)

Specifically the "4 common surnames + common initials" pattern, detected with a language-aware check.

**Detection heuristics by language**:
- Chinese: 4 of {Zhang, Wang, Li, Liu, Chen, Yang, Huang, Zhao, Wu, Zhou} + single initials
- English: 4 of {Smith, Johnson, Williams, Brown, Jones, Miller, Davis, Wilson} + single initials
- Indian: 4 of {Patel, Kumar, Singh, Sharma, Gupta, Reddy, Nair, Das} + single initials
- Japanese: 4 of {Sato, Suzuki, Takahashi, Tanaka, Watanabe, Yamamoto, Nakamura, Kobayashi}

**Fix**: Verify ALL 4 authors. At least 50% of the time, at least one is wrong.

→ **ARS M2**: Hallucinated citation

---

### Pattern 7: Citation Number Drift

**Description**: Body says "[43]" but bibliography [43] is a completely different paper.

**Domain example**: The body says "MUNIS achieved AUC 0.98 [43]" and reference [43] is "Reynisson et al. NetMHCpan 4.1. NAR. 2020." — the citation drifted during section reorganization.

**Detection**: Spot-check 10 random body-bibliography pairs per section. Expect 5-15 drift instances in a 100+ reference manuscript if careful, 30-40 if not.

**Fix**: Grep for the correct reference number and edit body in place.

→ **ARS M2 + M5**: Hallucinated citation + potentially reframed (if the drifted-to citation supports a different claim)

---

### Pattern 8: Metric Formula Errors

**Description**: Definitions of standard evaluation metrics are wrong.

**Domain examples**:
- ✗ "PPV = TP / (TP + FN)" — that's sensitivity/recall. PPV = TP / (TP + FP).
- ✗ "IC50 is the concentration at which 100% binding is inhibited" — it's 50%.
- ✗ "AUC is calculated as the area under the precision-recall curve" — that's AUC-PR, not AUC-ROC.
- ✗ "% Rank is the raw IC50 value" — % Rank compares IC50 to a background distribution of random peptides.

**Fix**: Verify every non-trivial metric formula against the original paper or a standard textbook.

→ **ARS M4**: Shortcut reliance (using memorized but wrong formula)

---

### Pattern 9: Internal Inconsistency Across Sections

**Description**: Same paper cited in different sections with different author lists, or same dataset appearing with different patient counts.

**Domain example**: The IEDB is cited as "687,924 entries" in Methods but "650,000 entries" in Discussion. MUNIS cited with 7 authors in Results but 9 authors in Introduction.

**Fix**: For each high-value paper/dataset, lock in a CANONICAL description in DOMAINS.md Section A.6. Every citation must match the canonical form.

→ **ARS M7**: Frame-lock (inconsistent frames adopted at different pipeline stages)

---

## Pattern 10 (NEW): Database Version Drift

**Description**: Citing a database statistic from a specific version when that number was from an older version, or fabricated entirely.

**Domain examples**:
- ✗ "IEDB 2024.02 contains 687,924 peptide-MHC binding entries" — the number might be from IEDB 2023 or fabricated
- ✗ "The latest NetMHCpan 4.1 was trained on 850,000 data points" — fabricated training set size
- ✗ "VDJdb contains 50,000+ TCR-pMHC pairs" — actual number is ~30,000

**Detection**: Any sentence with a database name + a version-specific statistic → verify by WebFetch on the database homepage or the latest NAR database issue paper.

**Fix**: Always WebFetch the database homepage for current statistics. Cite the NAR database issue paper for version-stable descriptions. Prefer qualitative descriptions ("contains hundreds of thousands of entries") when exact numbers are unnecessary.

→ **ARS M5**: Bug/error reframed as insight (if the version drift creates a false trend)

---

## Pattern 11 (NEW): Domain-Spanning Hyperbole

**Description**: Claiming a finding applies broadly when it was tested on a narrow subset. Particularly common in immunology due to HLA diversity.

**Domain examples**:
- ✗ "The model accurately predicts all peptide-MHC interactions" — tested only on HLA-A*02:01
- ✗ "This pan-specific predictor works for any HLA allele" — tested on 20 common alleles, none from underrepresented populations
- ✗ "MUNIS identifies all immunogenic CD8+ epitopes" — validated on 3 pathogens (influenza, HIV, EBV)
- ✗ "The method solves the TCR-pMHC prediction problem" — AUC drops from 0.9 to 0.6 on unseen epitopes

**Detection**: Any claim with universal quantifiers (all, any, every, universal, solves) + an epitope/MHC/TCR prediction method → verify the scope of validation.

**Fix**: Always note specific alleles tested, specific pathogens validated, specific tumor types where applicable. Replace "all" with the actual scope.

→ **ARS M7**: Frame-lock (early pipeline stage overgeneralizes the problem scope)

---

## Self-Check Workflow

### Per-Paragraph (every 5-6 paragraphs during writing)
- [ ] Scan for directional language → verify against source (Pattern 3)
- [ ] Scan for numeric claims + citation → verify number in paper (Pattern 2)
- [ ] Check new citations for 4-author generic patterns (Patterns 1, 6)
- [ ] Check for vendor/agency/site-name as "author" (Pattern 4)
- [ ] Check for universal quantifiers + narrow scope (Pattern 11)

### Per-Section
- [ ] Spot-check 10 body-bibliography matches (Pattern 7)
- [ ] Check for database version statistics (Pattern 10)
- [ ] Search top-5 cited papers for internal consistency (Pattern 9)

### Per-Manuscript
- [ ] `grep -i "xxx\|\[TBD\]\|placeholder"` (Pattern 5)
- [ ] Verify all metric formulas against original papers (Pattern 8)
- [ ] Check all database version citations (Pattern 10)
- [ ] Check all universal-scope claims (Pattern 11)

---

## ARS 7-Mode Failure Crosswalk

| ARS Mode | ARS Description | MedSci Pattern(s) | Severity |
|----------|----------------|-------------------|----------|
| M1 | Implementation bug passing AI self-review | (not a writing pattern) | — |
| M2 | Hallucinated citation | Patterns 1, 2, 4, 5, 6, 7 | **CRITICAL** |
| M3 | Hallucinated experimental result | Pattern 2 | **CRITICAL** |
| M4 | Shortcut reliance | Pattern 8 | HIGH |
| M5 | Bug reframed as novel insight | Patterns 3, 7, 10 | HIGH |
| M6 | Methodology fabrication | (not a writing pattern) | — |
| M7 | Frame-lock at early stage | Patterns 9, 11 | MEDIUM |

**Integration note**: When ARS pipeline is active (Mode 2/3), the `integrity_verification_agent` runs the 7-mode checklist automatically at Stages 2.5 and 4.5. MedSci's 11 patterns are supplemental — they catch writing-level hallucinations that the ARS machine-level checklist doesn't cover.
