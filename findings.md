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

## Update (2026-08-02, mission complete)

**5/5 alleles PASS** with the Deep FFN + BLOSUM62 recipe + curated IEDB benchmarks:
sens 100%, spec 97-99%, F1 0.94-0.98 for all. The prior run-001 "failures" were
benchmark-curation and training-variance artifacts.

**Non-canonical robustness confirmed:** models retain 0.88-1.0 sensitivity on IEDB
positives that violate canonical anchor motifs — they do not simply memorize anchors.
Precision 0.81-0.96 (lowest for A*02:01, its established hard case). This directly
supports the manuscript's claim that non-canonical binders (K-E63 salt bridge) are
the genuine hard case, not a failure of cross-allele generalization.

## Open Questions

- **Pan-specific model (Step B):** A single allele-aware model (BLOSUM62 + 5-dim
  allele one-hot) passes all 5 canonical benchmarks BUT loses non-canonical
  robustness — worst for A*03:01 (non-canonical sens 1.0 → 0.48). The allele
  one-hot dilutes allele-specific structural signal. Per-allele models remain
  better for non-canonical (neoantigen) discovery. **Resolved: per-allele wins.**
- Non-canonical precision is lowest for A*02:01 (0.81) — does a structural/PSSM
  hybrid recover it?
- Thresholds: A*03:01 optimal 0.95, others 0.97-0.99. Is there an allele-agnostic
  calibration?
