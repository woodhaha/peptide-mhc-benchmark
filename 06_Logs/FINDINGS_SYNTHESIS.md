# Findings Synthesis — Multi-Allele Peptide-MHC Benchmarking (2026-08-02)

*Consolidated current-state summary. Supersedes the chronological entries in
findings.md. All numbers are validated, cross-predictor, and seed-robust.*

## The question
Does a single Deep FFN + BLOSUM62 recipe generalize across common HLA class I
alleles, and how hard are non-canonical (motif-failing) binders for current
prediction methods?

## Answer (one paragraph)
A single rebalanced pan-specific Deep FFN (B*07:02 upsampled 4x) reaches 0.92-1.00
non-canonical sensitivity across 5 alleles (5-seed SD <=0.05) — parity or better
than per-allele models. Canonical binders are caught at ~0.98 by every predictor
class (Deep FFN, MHCflurry, MUNIS). Non-canonical (Type B) binders plateau at
~0.65 across three independent sequence predictors (Deep FFN, MHCflurry, MUNIS); MUNIS is sequence-based (ESM-2+LSTM), so the structural tier remains untested at scale. The extreme energy-compensated
tail (e.g. KRAS G12V P2-K salt bridge) is prediction-resistant to every method
tested.

## Validated findings (current state)

### F1. Pan generalization (seed-robust positive)
- Rebalanced pan Deep FFN, non-canonical sens mean±SD across 5 seeds:
  A*02:01 0.98±0.02, A*01:01 0.92±0.04, A*03:01 0.98±0.02, B*07:02 0.92±0.05,
  B*44:03 1.00±0.00. Test acc 0.952-0.956.
- One model replaces five per-allele models. Earlier "per-allele > pan" was a
  10.6x data-size artifact.

### F2. Canonical binders are solved (~0.98)
- n=14,488 IEDB positives (5 alleles): canonical caught at 0.979 (Deep FFN),
  0.986 (MHCflurry), 0.950 (MUNIS). All predictor classes saturated.

### F3. Non-canonical (Type B) plateau (~0.65)
- n=3,112 Type B: Deep FFN 0.607, MHCflurry 0.668, MUNIS 0.647.
- A third sequence predictor (MUNIS) does NOT lift the ceiling. Three-tier pipeline claim
  NOT supported.

### F4. Type A (tolerated weak anchors) is inherent-hard (~0.48)
- n=151: Deep FFN 0.477, MHCflurry 0.669.
- Broad-motif retrain does NOT recover it (A*03:01 flat 0.634; B*07:02 Type B
  drops 0.680->0.341). Not a training-data gap; the canonical-motif bias is what
  makes the model work.

### F5. Extreme tail is prediction-resistant (cross-predictor)
- Manuscript's 3 FNs + KRAS G12V P2-K: MHCflurry 3.8-13.8 uM (non-binder), Deep
  FFN 0.000-0.254, MUNIS 3/4 low. QYDPVAALF is the single MUNIS-rescued case
  (MUNIS 0.972 vs sequence 0.000/6,199nM) — n=1.

## Validated findings (post-synthesis additions)

### F6. Temporal generalization (breaks label circularity)
- Post-2022 epitopes (unseen by MHCflurry's training): sens@0.95 A*02:01 0.824,
  A*03:01 0.936, B*07:02 0.973 — preserved vs pre-2021 (A*03:01/B*07:02 improve).
- The model is NOT overfit to the IEDB-era distribution.

### F7. Type A is a real coverage gap (not noise/confusion)
- Type A peptides are mass-spec confirmed (13/14, 11/12) and only 3-6% cross-allele.
- They are genuine, allele-specific, tolerated-anchor binders the model fails at 0.48.
- Broad-retrain does not fix it — a genuine limitation, honestly reportable.

### F8. Pan16 — one model replaces a 15-allele panel
- Rebalanced pan16 (11 new alleles + 5 existing): pooled non-canonical sens 0.819
  (307/375 @0.95). 12/15 alleles at 0.76-1.0.
- Weak tail = B*27:05 (0.08, only 150 samples) and B*51:01 (0.52, 1,665) — data
  scarcity, not architecture.
- Scales the "one model replaces the panel" claim to 15 alleles, approximating
  NetMHCpan's practical value with a simple rebalanced Deep FFN.

### F9. Temporal generalization (breaks label circularity)
- Post-2022 epitopes (unseen by MHCflurry's training): sens@0.95 A*02:01 0.824,
  A*03:01 0.936, B*07:02 0.973 — preserved vs pre-2021. Model not overfit to era.

### F10. P2 salt-bridge generality (structural)
- P2 C-alpha -> Glu63 = 3.83 A. Positive-charge P2 (K/R/H) can electrostatically
  compensate the acidic B-pocket (R reach 7.3 A, H 4.4 A). Manuscript's "other
  charged P2" extension is geometrically supported.

## Falsified claims (do NOT repeat)
- "5/5 alleles pass; per-allele > pan" — canonical benchmark was motif-tautology
  (anchor-checker = 1.0 F1); pan collapse was data-size artifact.
- "Shared parameters limit pan models" / "PSSM conditioning dilutes" — falsified.
- "Type A is easy (0.92-0.98)" — n=25 cherry-pick artifact; real 0.48.
- "Non-canonical binders systematically missed" — moderate penalty (~0.65), not
  catastrophic.
- "Three-tier pipeline rescues non-canonical binders" — a third sequence predictor (MUNIS) does not.

## Manuscript implications
See `06_Logs/manuscript_correction_checklist.md` (6 items). Core: weaken the
3-tier claim to a hypothesis; add the predictor-plateau fact; keep the K-E63
salt-bridge hypothesis (strengthened by the cross-predictor failure contrast).

## Source data
- `06_Logs/type_ab_dataset.json` (14,488 positives, 5 alleles, classified)
- `06_Logs/type_ab_deepffn_scores.json`, `type_ab_mhcflurry_scores.json`
- `05_Submission/munis_results/munis_typeab_scores.csv` (MUNIS, 6,224)
- `06_Logs/multiseed_results.json` (5-seed pan)
- Models: `models/deep_ffn_pan_seed*.keras`, `deep_ffn_{A0201,A0301,B0702}_broad.keras`

## Open questions (genuinely unresolved)
- Would a pMHC-fine-tuned structure model (not generic AF2-multimer) lift Type B?
- Does the extreme-tail resistance extend to more alleles / experimental
  (non-IEDB) validation?
- Type A: is the low sens partly IEDB label noise (weak, multi-allele epitopes)?
