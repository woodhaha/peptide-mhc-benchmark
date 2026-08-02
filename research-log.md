# Research Log — Multi-Allele Peptide-MHC Extension

*Decision timeline. Updated 2026-08-02.*

## 2026-08-02 — Session arc

### Phase 1: Resume + environment fix
- **Decision:** Resumed multi-allele extension (only unfinished of 4 prior tracks).
- **Fix:** Restored r-tensorflow venv base Python 3.10.20 (deleted by uv cleanup)
  via `uv python install` + pyvenv.cfg repoint. TF 2.17/keras 3.12 verified.
- mhcflurry labeling runs in base anaconda (venv torch broken); training in venv.

### Phase 2: Benchmark curation (the pivotal correction)
- **Discovery:** Prior run-001 "failed" alleles were curation-degenerate:
  - A*03:01 benchmark contaminated with 4 A*02:01 epitopes.
  - B*44:03 collapsed to 5 near-duplicate EE*IVVLF 9-mers (F1 capped).
- **Decision:** Curated evidence-ranked IEDB epitopes via API (25/allele), fixed
  seed-42 splits, 50/50 mixed biased+random training.
- **Result:** All 5 alleles passed (sens 100%, spec 97-99%, F1 0.94-0.98).

### Phase 3: Three-judge review + falsification
- **Decision:** Submitted claims to 3 adversarial reviewers (methodology /
  statistics / immunoinformatics).
- **Findings:** canonical F1 ~motif-tautology (anchor-checker = 1.0); pan collapse
  = data-size artifact; non-canonical partition partly motif-dict artifact.
- **Falsification experiments confirmed all 3.** Pan on full data recovered
  A*03:01 0.48->1.0. The naive "per-allele > pan" story was dead.

### Phase 4: Rebalancing → defensible positive
- **Decision:** Upsample B*07:02 (4x) in pan training.
- **Result:** Non-canonical 0.84-0.96 across all 5 alleles (5-seed SD <=0.05).
  Single rebalanced pan generalizes — one model replaces five.

### Phase 5: Literature reframe
- **Discovery:** Pan-specific is field consensus (NetMHCpan paradigm). Our
  pan-vs-per-allele result = replication, not novel.
- **Pivot:** Asked what the multi-allele work says about the manuscript's
  non-canonical claim → Type A/B distinction.

### Phase 6: Large-sample Type A/B revision
- **Decision:** Test Type A/B at scale (n=14,488, dual-predictor).
- **Result:** canonical ~0.98; Type B ~0.61-0.67 (moderate, not catastrophic);
  Type A 0.48 (NOT easy — n=25 was cherry-pick). "Systematically missed" reframed
  to "moderate penalty + extreme tail resistant."

### Phase 7: Structural validation (MUNIS)
- **Decision:** Installed MUNIS (jwohlwend/munis) on AutoDL RTX 6000D; scored
  full Type B population (n=3,112 + 3,112 canonical).
- **Result:** MUNIS Type B 0.647 = same plateau as sequence. Structure does NOT
  rescue non-canonical binders. Three-tier pipeline claim NOT supported.
- **Broad-motif retrain:** does NOT recover Type A (B*07:02 Type B hurt).
  Type A low sens is inherent, not a training gap.

### Phase 8: Consolidation (current)
- **Decision:** Wrote FINDINGS_SYNTHESIS.md (paper backbone) +
  manuscript_correction_checklist.md (6 items).
- **Status:** Research converged; findings complete and validated. Next: apply
  corrections to main manuscript OR write companion benchmarking note.

### Phase 9: Improvement round (5 directions)
- Type A evidence: mass-spec confirmed, 3-6% cross-allele => REAL coverage gap, not noise.
- Temporal holdout: post-2022 (unseen by MHCflurry) sensitivity preserved/improved =>
  BREAKS label-circularity concern. Strong positive.
- P2 virtual mutagenesis: P2CA->Glu63 = 3.83A; K/R/H all reach => salt-bridge is a
  general principle (charged P2), supports manuscript's "other charged P2" extension.
- NetMHCpan: BLOCKED by DTU license (no server install, IEDB API unreachable).
  MHCflurry + MUNIS suffice as independent predictors.
- Pan16: labeled 11 new alleles + 5 existing, rebalanced pan => pooled non-canonical
  0.819 (15-allele panel), 12/15 at 0.76-1.0. Weak tail = thin-training B*27:05/B*51:01.
  Scales "one model replaces the panel" to 15 alleles.

### Phase 10: CONCLUSION (current)
- 10 validated findings (F1-F10) in FINDINGS_SYNTHESIS.md. Research complete.
- Next is authoring: apply manuscript_correction_checklist.md to the main manuscript
  OR write a companion benchmarking note. Not more experiments.

### Phase 11: COMPANION NOTE (authoring, 2026-08-02)
- ml-paper-writing installed (from zechenzhangAGI/AI-research-SKILLs, v1.2.0).
- Full first draft written: 04_Manuscript/companion_multi_allele/companion_note_draft.md (~2,000 words, IMRaD, target Bioinformatics).
- All 15 numeric claims verified against source data (type_ab_dataset.json, summarize_type_ab.py re-run, munis csv recompute, multiseed_results.json, p2_virtual_scan.json).
- Figure 1 (4-panel findings panel) generated via plot_findings_panel.py, copied to companion dir.
- Citation discipline: reused verified refs from manuscript.md (54 refs); no hallucinated citations.
- Framing: positive (rebalanced pan = per-allele parity, 15-allele panel) + negative (Type B plateau ~0.65 cross-predictor) + structural contrast (P2-K/R/H hard tail).

### Phase 11b: COMPANION NOTE REVISIONS + FIG VECTORIZATION (2026-08-02)
- Draft review found 5 defects, all fixed (commit ccb9133):
  - [10] Robinson/IMGT-HLA verified via CrossRef (real, NAR 2015;43:D423-D431)
  - Metadata corrected (~2,000 words, 1 table, 2 figures)
  - [companion submission] self-refs -> formal [15]
  - "pre-registered" -> honest "fixed thresholds" (post-hoc 5th-pct canonical)
  - New Figure 2 (rebalanced pan 5-seed mean+/-SD); per-allele bars excluded
    (n=25 benchmark != large non-canonical set = apples-to-oranges)
- Figures 1-2 now have vector PDFs (commit 4d35a78).
- All matplotlib figures vectorized (commit ab8f141): md pic/ 15 + to_human/ 2.
  Excluded: docking/p53 PyMOL renders (2, raster by nature), analysis_100ns.png
  (no source script).
- Fixed falsified claim in plot_multi_allele_results.py (commit 17e8914):
  suptitle "per-allele > pan" + A*03:01 0.48 -> corrected to 5-seed means
  (0.984/0.920/0.976/0.920/1.000) and "parity" framing. Historical falsification
  records (falsify_pan_equalfull.py, to_human/progress_*) intentionally kept.

### Phase 11c: MERGE DECISION — companion note folded into main manuscript (2026-08-02)
- User decision: "不要创建姐妹篇，合并在一起" — do NOT split into a separate paper; merge the multi-allele extension into the main manuscript.
- Main manuscript already referenced the multi-allele findings in Discussion/Conclusions (4.3 MUNIS plateau, 4.6 Sixth, 4.7 pan16) but Results lacked the supporting evidence — a structural gap.
- Merge performed on 04_Manuscript/manuscript.md:
  - Methods: +2.8 Multi-Allele Data Curation, 2.9 Pan Model & Comparators, 2.10 Temporal Holdout
  - Results: +3.11 pan parity (5-seed), 3.12 canonical ~0.98, 3.13 Type B plateau (Table 7), 3.14 Type A inherent, 3.15 temporal, 3.16 pan16
  - +Table 7 (predictor plateau), +Figure 9 (rebalanced pan 5-seed) legend + image
  - Abstract conclusions + multi-allele sentence; Data/Code + pan models & type_ab dataset
  - Metadata: ~8,500 -> ~11,000 words, 6 -> 7 tables, 8 -> 9 figures
- Standalone companion note archived to 06_Logs/archived/companion_note_draft_superseded.md (content now in main manuscript).
- Verified: Tables 1-7 and Figures 1-9 + S1-S6 numbered continuously.

### Phase 11d: SYNC submitted version to merged manuscript (2026-08-02)
- User decision: one manuscript, not two. manuscript.md = single source of truth.
- manuscript_bioinformatics.md is the Bioinformatics-formatted rendering (Harvard refs),
  regenerated from manuscript.md via reformat_bioinformatics.py.
- Ran the converter: manuscript_bioinformatics.md now includes the merged multi-allele
  extension (Methods 2.8-2.10, Results 3.11-3.16, Table 7, Figure 9), Harvard citations
  intact, 10,507 words.
- NOTE: this regenerated version is NOT what the journal currently has (submitted June 30
  without multi-allele). It represents the revision-ready version if the manuscript is
  accepted with revisions. The originally-submitted content is preserved in git history.

### Phase 11e: FIGURE STYLE UNIFICATION (2026-08-02)
- Unified python-generated figure colors from Tailwind palette to Nature palette
  matching R figures (commit fa74606): blue #377EB8, red #E41A1C, orange #FF7F00,
  green #4DAF4A, gray #999999. Applied to plot_pan_vs_perallele (Figure 9),
  plot_findings_panel, plot_multi_allele_results, plot_type_ab. Verified Figure 9
  bar RGB~(48,112,176) ~ #377EB8.
- Unified fonts to Arial (commit 261f1d6): R figures use family='sans' = Arial on
  Windows; python figures had DejaVu Sans. Added plt.rcParams["font.family"]="Arial".
  Verified PDFs embed Arial, no DejaVu.
- All manuscript figures now consistent: Nature palette + Arial font.

### Phase 11f: CONSISTENCY AUDIT (2026-08-02)
- Verified merged manuscript.md internal consistency after all edits:
  - All 14 cited figures exist (Figure1-9 + S1-S5), paths resolve from manuscript location.
  - Tables 1-7 numbered continuously.
  - All 54 references present; max in-text citation 53 (ref 54 = Editorial Board bioRxiv,
    cited-only-in-reference, acceptable but flagged).
- Journal recommendation context: manuscript is with BiB (under review, NOT rejected).
  Backup options if needed: Bioinformatics (formatted version ready), NAR (needs repackaging).
- Manuscript is revision-ready and internally consistent.
