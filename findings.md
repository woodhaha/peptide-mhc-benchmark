# Findings — Multi-Allele Extension

## Current Understanding (2026-08-02)

The single-allele manuscript (A*02:01, submitted to Bioinformatics) established:
**label quality > architecture** for peptide-MHC binding prediction. The multi-allele
extension asks whether the Deep FFN + BLOSUM62 recipe generalizes across HLA alleles.

**Key result so far:** The Deep FFN generalizes cleanly to other alleles when (a) the
evaluation benchmark is properly curated and (b) training uses fixed deterministic
splits with mixed biased+random data. 3/3 previously-"failed" alleles now pass.

## Patterns and Insights

1. **"Model failures" were benchmark artifacts.** A*03:01's benchmark was
   contaminated with A*02:01 epitopes; B*44:03's collapsed to 5 near-duplicate
   peptides (F1 math cap ~0.63). Neither reflected model capability. *Lesson: verify
   benchmark curation before trusting a negative result.*
2. **Training variance is a reproducibility trap.** Unseeded random data generation
   produced 22-89% sensitivity swings across runs for A*01:01. Seed-42 fixed splits
   eliminated it. *Lesson: always save training data to disk with a fixed seed.*
3. **Mixed training data matters.** Pure motif-biased data causes catastrophic
   overfitting (from prior run-001). 50/50 biased+random gives the model both the
   motif signal and the ability to discriminate.

## Lessons and Constraints

- BLOSUM62 canonical-motif epitopes are the reliable eval signal; non-canonical
  binders are a separate (harder) problem already documented for A*02:01.
- mhcflurry's Class1AffinityPredictor loads in base anaconda (torch OK) but not in
  the r-tensorflow venv (torch c10.dll broken). Split workflow: label in base, train
  in venv.
- Old `.h5`/`.keras` models saved under keras 2.x have a `quantization_config` kwarg
  that keras 3.12 rejects — retrain rather than patch-load.

## Open Questions

- Do the 2 remaining alleles (A*02:01, B*07:02) also pass with the same recipe?
- Does a pan-specific model (allele encoded in input) match per-allele models? (Step B)
- How does per-allele threshold behavior differ? (A*01:01/B*44:03 optimal ~0.97,
  A*03:01 ~0.95 — all above the original 0.95 binding-score rule.)
