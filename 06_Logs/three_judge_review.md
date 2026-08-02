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

## Round 3 — Corrected verdicts (all 3 judges converge)

### Dead (falsified by Round-2 experiments)
1. "100% / 97-99% / F1 0.94-0.98" canonical headline — n=25 Wilson CI floor 86.7%;
   motif tautology (anchor-checker = 1.0 F1); synthetic decoy negatives.
2. "Per-allele > pan for non-canonical" — was a 10.6x data-size artifact.
   Pan on full data recovers A*03:01 non-canonical 0.48->1.0, test acc 0.953.
3. "Shared parameters fundamentally limit pan" + "one-hot dilutes / PSSM worse"
   mechanism — falsified by Exp B; the comparison is now UNTESTED, not resolved.
4. "Canonical passes 5/5 alleles" as evidence of capability — matches a trivial rule.

### Survives (narrowly)
1. Prior run-001 "failures" were curation/variance artifacts (well-documented).
2. Real non-canonical sensitivity (0.88-1.0 per-allele, 0.76-1.0 pan-full) — but
   partly inherited from MHCflurry (teacher leakage) + protein-likeness signal.
3. A*02:01 protein-negative beat over the anchor rule (F1 0.943 vs 0.847) — the
   single non-tautological number, one allele, threshold-grid caveat.
4. Defensible negative-result claim: "with full data, a pan model reaches parity
   with per-allele models on non-canonical binders; per-allele superiority is not
   supported."

### Requirements before publishable (from all 3 judges)
- Protein-background negatives (not uniform-random).
- Motifs from independent validated sources (netMHCpan anchor preferences /
  supertype motifs), not hand-built dicts; recompute the canonical/non-canonical
  partition. B*44:03 anomalous -765 = hand-built broad set excluded a real anchor.
- Multi-seed (>=5) mean+/-SD; pre-register one threshold.
- External labels (e.g. affinity <500nM) never seen by MHCflurry to break
  IEDB/MHCflurry circularity.
- Report n=25 CIs; drop salt-bridge/structural-mechanism and per-allele-superiority
  from any abstract.
- Non-canonical must be split: Type A (dict-coverage bug) vs Type B (true
  energy-compensation binder). The manuscript's A*02:01 P2-K is Type B; the A*03:01
  P9-Y cases were Type A mislabeled by the narrow motif.

### Status
The multi-allele extension, as currently written, is NOT publishable. The falsified
claims must be removed or reframed. The single surviving kernel (pan-full = per-allele
parity on non-canonical) is a negative-result benchmarking finding, contingent on the
fixes above.
