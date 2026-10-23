# Three-Judge Review + Correction — 2026-10-23

Post-reframe adversarial review of `manuscript_bioinformatics_v2.md` (3 independent
agents: methodology/statistics, immunoinformatics, editor). All three judged the
draft reject/desk-reject. Corrections applied across this session.

## Verdicts

- **Judge 1 (methodology)**: reject — 5 FATAL + 17 MAJOR + 8 minor. Core evidence chains unclosed.
- **Judge 2 (immunoinformatics)**: desk-reject — went to the live IEDB API and found ground-truth contamination.
- **Judge 3 (editor)**: desk-reject ~60–65% — novelty pre-empted, anchor-rule arithmetic impossible, temporal holdout self-contradictory.

## Converged FATAL findings (≥2 judges)

1. anchor-rule "matched the trained model (sensitivity 1.0)" — numerically impossible
   (real sensitivity 0.786 under the paper's own motif).
2. 49-epitope IEDB ground truth contaminated (Judge 2, IEDB-verified).
3. 6,524-affinity validation is circular (same corpus), yet abstract claims "generalisation".
4. temporal holdout: model identity unspecified + training snapshot post-dates holdout + "breaks" vs "preserved" contradiction.
5. ~0.65 plateau built on 3 incompatible thresholds (MUNIS/BigMHC 0.950 canonical is a 5th-percentile arithmetic identity).
6. 62.24% over-attributed to label circularity.
7. Type A/B no operational definition; three incompatible sets (0.84/0.08/0.477).
8. 1:1:1 balanced accuracy without natural-prior metrics — the paper commits its own benchmark-composition sin.

## Ground-truth verification (IEDB API, 2026-10-23)

`02_Data/cleaned/verify_iedb_epitopes.py` queried every one of the 49 "A*02:01" positives:

- 42 legit A*02:01 Positive binders.
- 7 contaminated: `LLFDRFENL`, `LLHHAFDSL` (0 IEDB records); `YVKQNTLKL` (B*07:02);
  `ILRGSVAHK` (A*03:01, not the stated B*08:01); `QYDPVAALF` (A*24:02, the classic CMV pp65
  A24 epitope — mis-read as a non-canonical A*02:01 false negative); `CMTWNQMNL` (A*24:02);
  `VIYQYMDDL` (A2-negative).

Corrected: sensitivity **41/42 = 97.6%** (was 95.8% = 46/48), ROC AUC **0.970** (was 0.953),
precision 89.1%, F1 0.932. Single remaining false negative = `TLGIVCPIC` (p9=C).

## Re-run findings (transparently added to manuscript)

- **Threshold sensitivity** (`scripts/threshold_sensitivity.py`): reproduces manuscript Type B
  sensitivities exactly (DFFN 0.607 / MHCflurry 0.668 / MUNIS 0.647 / BigMHC 0.660), then shows
  the ~0.65 "plateau" is partly a threshold artifact — aligned to matched canonical sensitivity
  0.95, predictors diverge (0.454 / 0.497 / 0.647 / 0.653).
- **Multi-seed CI** (`02_Data/cleaned/multiseed_results.json`): CNN vs Deep FFN not separable at
  n=519 (95% CI [0.889,0.949] vs [0.872,0.940]); combined seed+binomial SE ~1.5 pp > the 1.35 pp gap.

Both findings support the paper's own thesis — its own positive claims also dissolve under
controlled re-analysis.

## Corrections applied (md + LaTeX + README + CLAUDE.md + audit_config)

- ground truth: 95.8%→97.6%, 0.953→0.970, precision/F1, 49→42 everywhere.
- anchor-rule: "matched (1.0)" → "underperforms (0.786)", model recovers 8/9 non-canonical.
- 62.24%: "is label circularity" → "label-definition divergence + circularity (unknown proportion)".
- Jessen → Nielsen (manuscript + multiseed JSON).
- temporal: "breaks the label-circularity concern" → "shows no detectable generalisation penalty".
- Type A: "inherently hard" → "under-represented in training" (Table 7 contradicts inherent).
- Box 1: added item (5) natural-prior reporting.
- plateau: "plateaued at ~0.65" → "reached ~0.65 at unmatched operating points" (abstract/Results/Discussion/Conclusion).
- architecture: added combined SE/CI caveat (CNN not separable at n=519).

## Commits (2026-10-23, master)

`470f965` verify IEDB · `a0b06ad` ground-truth correction · `b35f13f` Fig3 ROC ·
`8088134` anchor-rule script + Box1(5) · `eadd342` Type A reframe · `499fb7d` threshold
sensitivity · `a0f6201` transparent reporting.

## Remaining (multi-day / blocked / decision)

- temporal holdout leakage (re-run pan model with time-cutoff split).
- Type A/B operational definitions (per-allele motif).
- natural-prior actual metric (needs test-set per-class predictions).
- 6,524 de-circularize (re-scope validation set).
- NetMHCpan head-to-head (DTU license blocked).
- **Structural: novelty/no-tool → journal pivot (PLOS Comput Biol / Briefings).**

## Bottom line

The paper's honest value is now a benchmark-hygiene case study that caught its own
ground-truth contamination and threshold artifacts — not a novel non-canonical-frontier
finding. Bioinformatics (wants new tools) is a poor fit; PLOS Comput Biol or Briefings
fits the self-reflective negative-result narrative.
