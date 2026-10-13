# Manuscript Correction Checklist (2026-08-02)

Based on the multi-allele + Type A/B + cross-predictor + MUNIS (sequence) validation
work. These findings should be incorporated into the submitted manuscript (Briefings
in Bioinformatics) either as revisions or as honest caveats in the discussion.

## 1. CORE: the "systematically miss non-canonical binders" claim — weaken

**Current claim (Abstract/Results 3.9, 4.1):** non-canonical anchor binders are
"prediction-resistant configurations that sequence-only models systematically miss."

**Corrected finding (n=14,488, 3 predictors):**
- Canonical binders: caught at ~0.98 by Deep FFN, MHCflurry, AND MUNIS.
- Motif-failing (Type B) binders: caught at ~0.61-0.67 by ALL THREE (no predictor
  class rescues them).
- The manuscript's 3 FNs (ILRGSVAHK, QYDPVAALF, TLGIVCPIC) + KRAS P2-K are the
  HARD TAIL (MHCflurry: 3.8-13.8 uM; MUNIS: 3/4 low; Deep FFN: 0.000-0.254).

**Action:** Reframe from "systematically missed" → "the extreme energy-compensated
non-canonical tail is prediction-resistant (cross-predictor validated: Deep FFN,
MHCflurry, MUNIS all fail)."

## 2. CORE: the three-tier pipeline (sequence→AlphaFold→structure) — NOT supported

**Current claim (Discussion 4.5):** a three-tier screening pipeline (Tier 2 AlphaFold3,
Tier 3 energy minimization) would "capture these blind-spot configurations."

**Corrected finding:** MUNIS (sequence-based pMHC predictor, ESM-2+LSTM; NOT the structural Tier-2 equivalent)
catches Type B at 0.647 — statistically identical to sequence models (0.607-0.668).
A third sequence predictor (MUNIS) does NOT lift the non-canonical sensitivity ceiling at population scale; structure was not benchmarked.

**Action:** Re-frame the three-tier pipeline as a HYPOTHESIS needing experimental
validation, not an evidence-backed claim. Optionally note the single counterexample
(QYDPVAALF, MUNIS 0.972 vs sequence 0.000/6,199nM) as a promising but n=1 case.

## 3. Methods: benchmark composition — add the canonical-filter caveat

**Current claim:** "49 experimentally verified IEDB epitopes."

**Corrected finding:** the benchmark's canonical-motif positives overstate
generalization — a zero-parameter anchor-checker achieves sens 1.0/F1 0.85-1.0 on
the same canonical set. The real signal is the non-canonical/hard-tail.

**Action:** Add an anchor-checker baseline row; report non-canonical sensitivity as
the informative metric. → DONE in manuscript.md §3.5 (2026-08-02).

## 4. Methods: single-allele framing — the multi-allele result supports it

**Current:** HLA-A*02:01 single-allele focus framed as deliberate design choice.

**Corrected finding:** rebalanced pan model generalizes to 5 alleles (non-canonical
0.92-1.00, 5-seed SD ≤0.05). This SUPPORTS the single-allele depth argument — the
recipe is portable.

**Action:** Optionally add a one-line note that the recipe was validated to extend
across 5 common HLA-I alleles (as a strength, not a limitation).

## 5. Limitations (4.6): add the predictor-plateau fact

**Action:** Add a limitation: "Across two independent sequence predictors (our Deep
FFN, MHCflurry) and a second independent sequence predictor (MUNIS), non-canonical binders plateau
at ~0.65 sensitivity. The energy-compensated tail (e.g., P2-K salt bridge) remains
prediction-resistant to all current methods." → DONE in manuscript.md §4.6 (2026-08-02).

## 6. Positive to keep

- K-E63 salt-bridge structural hypothesis (two sequence predictors + structure all
  fail YKLVVVGAV, yet structure-informed analysis predicts binding — this contrast
  makes the structural story MORE compelling, not less).
- Label quality > architecture (the core thesis) — unaffected by this work.

## Summary table for the discussion

| Predictor class | Canonical sens | Type B sens |
|-----------------|:--------------:|:-----------:|
| Deep FFN (sequence) | 0.979 | 0.607 |
| MHCflurry (sequence) | 0.986 | 0.668 |
| MUNIS (sequence) | 0.950 | 0.647 |

All source data: `06_Logs/type_ab_*.json`, `05_Submission/munis_results/`.
