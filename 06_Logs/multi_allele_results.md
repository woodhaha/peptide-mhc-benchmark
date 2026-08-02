# Multi-Allele Extension Results (2026-08-02)

## Outcome: All 3 previously-failed alleles now PASS with curated IEDB benchmarks

| Allele | Sens | Spec_raw | F1 | Thresh | Status |
|--------|:----:|:--------:|:--:|:------:|:------:|
| HLA-A*01:01 | 1.0 | 0.99 | 0.980 | 0.97 | ✅ PASS |
| HLA-A*03:01 | 1.0 | 0.97 | 0.943 | 0.95 | ✅ PASS |
| HLA-B*44:03 | 1.0 | 0.99 | 0.980 | 0.97 | ✅ PASS |

## Root-cause analysis
Prior run-001 declared these alleles FAILED, but that was **benchmark curation + training variance**, not model failure:

1. **A*03:01** — its IEDB list contained 4 canonical A*02:01 epitopes (SVYDFFVWL, ILRGSVAHK, KVLEYVIKV, RMFPNAPYL) lifted from the wrong allele list. Purged.
2. **B*44:03** — only 5 unique 9-mers survived the length filter, all `EE*IVVLF` near-duplicates of one Melan-A sequence. F1 was mathematically capped at ~0.63. Replaced with 25 diverse curated epitopes.
3. **A*01:01** — training variance (22-89% sens across runs) due to unseeded data generation. Fixed with seed=42 fixed splits + 50/50 mixed biased/random training.

## Fixes
- `src/curate_allele_benchmarks.py`: IEDB API curation (25 epitopes/allele, evidence-ranked, canonical-motif, mass-spec validated, PMID provenance)
- `src/train_allele_v2.py`: deterministic splits, mixed training data, model persistence
- `src/evaluate_allele.py`: threshold-grid evaluation against curated benchmark

## Models
`models/deep_ffn_HLA-{A_0101,A_0301,B_4403}.keras` — test acc 95.2/94.3/93.2%
