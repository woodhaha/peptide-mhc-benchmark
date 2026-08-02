# Three-Judge Review — Round 2 Falsification Results (2026-08-02)

## The three tribunals (Round 1)
- **Judge 1 (methodology):** canonical F1 may be a tautology of the motif filter;
  threshold grid = test-set selection; teacher (MHCflurry) leakage.
- **Judge 2 (statistics):** pan-vs-per-allele confounded by 10.6x data gap; n=25 CIs;
  one seed; best-of-K thresholds; B*07:02 underpowered.
- **Judge 3 (immunoinformatics):** motif over-narrow (A*03:01 P9=Y is a real anchor);
  uniform-random negatives confound protein-likeness; MS labels not independent.

## Round-2 falsification experiments

### Experiment A — anchor-checker baseline (CONFIRMED JUDGE 1)
A zero-parameter rule ("bind iff canonical anchors match") achieves on all 5 alleles:
sens 1.0, F1 0.85-1.0 on the same benchmark. **The "100% canonical sensitivity" is
largely guaranteed by construction** (benchmark positives are motif-filtered; training
was 50% motif-biased). The trained model only meaningfully beats the rule on A*02:01
with protein-background negatives (F1 0.943 vs 0.847).

### Experiment B — pan on full data (FALSIFIES the "shared-params limit" claim)
Pan trained on FULL per-allele data (45,360 samples, no min-allele balancing):
- Test acc 0.953 (vs 0.852 for balanced pan)
- **A*03:01 non-canonical sens: 0.48 -> 1.0** (collapse completely recovered)
- All 5 alleles: non-canonical sens 0.76-1.0 @0.95

**Verdict: the pan "collapse" was a 10.6x data-size artifact, NOT shared-parameter
dilution.** Judge 2 correct; Judge 1's "weak conditioning" mechanism wrong.

### Experiment C — motif broadening (CONFIRMED JUDGE 3, partial)
Broadening the A*03:01 P9 to {K,R,Y,H} reclassifies 640 IEDB positives from
"non-canonical" to canonical; B*07:02 P9 +{V,F} reclassifies 299. The A3-supertype
motif is genuinely broader than the 2-residue dict. **The "non-canonical" class was
partly a motif-dict artifact.** (B*44:03 result anomalous: -765 suggests my broad set
inadvertently excluded a valid P9 residue; needs re-check.)

## What survives
- Prior run-001 "failures" were curation/variance artifacts (well-documented).
- Per-allele models retain real non-canonical sensitivity (0.88-1.0), though
  partly inherited from MHCflurry + protein-likeness signal.
- **The defensible kernel: a properly-fed pan model (full data) matches per-allele
  models on non-canonical binders.** The "per-allele > pan" design guidance is
  NOT supported by these data.

## What was wrong in our earlier story
- "Shared parameters fundamentally limit pan models" — FALSE (data-size artifact).
- "Allele one-hot dilutes structural signal" — FALSE (pan-full recovers).
- "PSSM worse" — confounded by the same data imbalance.
- "100% canonical / 97-99% spec" — overstated (n=25, motif-tautology).
