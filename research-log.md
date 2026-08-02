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
