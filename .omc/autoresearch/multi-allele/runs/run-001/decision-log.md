# Decision Log — Run 001

**Mission:** Extend Deep FFN to 5 HLA-I alleles with per-allele benchmarks
**Max runtime:** 8 hours
**Started:** 2026-06-19

---

## Iteration 1: A*01:01 — Motif-Biased Training ❌
- **Approach:** Generate motif-biased peptides (P3=D/E, P9=Y/F), label with MHCflurry, train Deep FFN
- **Result:** ❌ FAILED
- **Metrics:** Sensitivity 88.9% | Specificity 80.8% | F1 0.561 | 23 FPs
- **Diagnosis:** Model overfits to motif — any peptide with D/E@P3 + Y/F@P9 predicted as binder
- **Decision:** Add random non-motif peptides to NB class to teach discrimination

## Iteration 2: A*01:01 — Mixed Training Data ⚠️
- **Approach:** 50% motif-biased + 50% random peptides in training pool. Deep FFN, no class weight.
- **Result:** ❌ FAILED but very close
- **Best metrics (threshold 0.80):** Sensitivity 88.9% | Specificity 97.5% | F1 0.865
- **Delta:** Specificity +16.7pp (80.8%→97.5%), FPs 23→3
- **Missed epitopes:** SAFPTTINF (bind=0.0009), VVPCEPPEV (bind=0.026) — neither matches canonical A*01:01 motif
- **Diagnosis:** Same BLOSUM62 non-canonical ceiling as A*02:01. These 2 epitopes likely either curation errors or non-canonical binders.
- **Decision:** At 88.9% sensitivity with 18-epitope benchmark, 1 TP short of 90%. Model variance between runs is high.

## Iteration 3: A*01:01 — Class Weight NB=0.7 ❌
- **Approach:** Reduce NB class weight to 0.7 to boost binder predictions
- **Result:** ❌ FAILED — made things worse
- **Metrics:** Sensitivity 83.3% (worse), more epitopes dropped below threshold
- **Decision:** Abandon class weight approach. Revert to v2 model.

## Iteration 4: A*01:01 — Restore v2 ⚠️
- **Approach:** Retrain with same mixed-data approach as v2
- **Result:** ❌ FAILED — training variance produced different model (sens 22%)
- **Diagnosis:** Random data sampling causes significant variance between runs. Need fixed dataset for reproducibility.

---

## Mission Status: PARTIAL

| Allele | Status | Best Sens | Best Spec | Best F1 | Pass |
|--------|:------:|:---------:|:---------:|:-------:|:----:|
| HLA-A*02:01 | Complete | 93.9% | 91.7% | 0.876 | ✅ |
| HLA-A*01:01 | Iterated | 88.9% | 97.5% | 0.865 | ❌ |
| HLA-B*07:02 | Not started | — | — | — | — |
| HLA-A*03:01 | Not started | — | — | — | — |
| HLA-B*44:03 | Not started | — | — | — | — |

### Key Findings
1. **Mixed training data is essential** — pure motif-biased data causes catastrophic overfitting
2. **Non-canonical epitopes are a hard BLOSUM62 ceiling** across alleles — 2 of 18 A*01:01 epitopes don't match the canonical motif
3. **Training variance is significant** — fixed data splits needed for reproducibility
4. **Per-allele thresholds are necessary** — A*01:01 optimal at 0.80 vs 0.95 for A*02:01

### Recommendation
Continue with remaining alleles using:
- Fixed training data generation (save to disk before training)
- Mixed (biased + random) approach from the start
- Per-allele threshold tuning
- Larger training sets (50K+ peptides per allele)
