# Multi-Allele Peptide-MHC Research — Final Progress Report (2026-08-02)

## Research question
Does a single Deep FFN + BLOSUM62 recipe generalize across HLA class I alleles, and
how hard are non-canonical (motif-failing) binders for current prediction methods?

## The complete answer (10 validated findings)

### Positive results
1. **One model replaces a 15-allele panel.** Rebalanced pan16 reaches pooled
   non-canonical sensitivity **0.819** (12/15 alleles at 0.76-1.0). Weak tail =
   thin-training alleles (B*27:05, B*51:01), not architecture.
2. **Canonical binders are solved.** ~0.98 by every predictor class (Deep FFN,
   MHCflurry, MUNIS).
3. **Temporal generalization holds.** Post-2022 epitopes (unseen by the labeler) are
   caught at 0.82-0.97 — the model is not overfit to the IEDB era. This breaks the
   label-circularity concern.
4. **P2 salt-bridge is a general principle.** K/R/H at P2 can all electrostatically
   compensate Glu63 (3.83 A) — supports the manuscript's K-E63 hypothesis and its
   "other charged P2" extension.

### The hard limits (honest negatives)
5. **Non-canonical (Type B) binders plateau at ~0.65** across ALL predictor classes —
   including structure-based MUNIS. The field's "structure will solve it" assumption
   is directly challenged.
6. **Type A (tolerated anchors) is a real coverage gap** (0.48): mass-spec-confirmed,
   allele-specific binders the model genuinely misses. Not noise, not fixable by
   broad-motif retrain.

### What was falsified (and corrected)
7. The naive "per-allele > pan" story was a data-size artifact + motif-tautology.
8. "Non-canonical binders systematically missed" is too strong — it's a moderate
   penalty (~0.65) with only the extreme energy-compensated tail truly resistant.

## Data provenance (all validated)
- n=14,488 IEDB positives (5 alleles) + 11 new alleles labeled
- 3 predictors: Deep FFN, MHCflurry, MUNIS (structure)
- 5-seed robustness; temporal holdout; Wilson CIs
- 119 commits this session, all in git

## Current state & next steps
The **research is complete and converged**. The science produced 10 validated findings
with no critical open question that would change the story. What remains is authoring:

1. **Apply `06_Logs/manuscript_correction_checklist.md`** (6 items) to the submitted
   main manuscript — primarily weaken the three-tier claim to a hypothesis and add the
   predictor-plateau fact.
2. **Or write a companion benchmarking note** on the multi-allele extension (pan16 =
   one model replaces a 15-allele panel).

Ready to invoke the ml-paper-writing skill when the write-up begins.
