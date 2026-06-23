# Mission: Multi-Allele Extension — Step A (Independent Models)

## Goal
Extend Deep FFN peptide-MHC binding prediction from single allele (HLA-A*02:01) to 5 common HLA class I alleles with independent per-allele models and benchmarks.

## Target Alleles
| Allele | Frequency | Existing Model |
|--------|:---------:|:--------------:|
| HLA-A*02:01 | ~40% (Caucasian) | ✅ FFN_Deep.h5 |
| HLA-A*01:01 | ~15% | New |
| HLA-B*07:02 | ~12% | New |
| HLA-A*03:01 | ~12% | New |
| HLA-B*44:03 | ~5% | New |

## Target
| Metric | Per-Allele | Aggregate (mean) |
|--------|:----------:|:----------------:|
| Sensitivity | ≥ 90% | ≥ 90% |
| Specificity_raw | ≥ 80% | ≥ 80% |
| F1 | ≥ 0.85 | ≥ 0.85 |
| Classifier | binding_score ≥ 0.95 | same |

## Data Strategy
- Training data: MHCflurry-labeled random 9-mers per allele (~5,000-10,000 peptides per allele, balanced SB/WB/NB)
- Positive benchmark: IEDB-validated epitopes per allele (curated from IEDB database)
- Negative benchmark: 4-source diverse negatives + allele-specific negatives
- Reuse A*02:01 data from existing project

## Architecture
- Deep FFN: 360-BN-180-BN-90-45-3 (same as existing best model)
- BLOSUM62 encoding (180-dim)
- Training: Adam lr=0.001, batch_size=64, early_stopping patience=15, 150 epochs max

## Evaluator
**Path:** `.omc/autoresearch/multi-allele/evaluator.py`
**Command:** `python evaluator.py <model_dir> <allele>`
**Output:** `{"pass": bool, "allele": str, "sensitivity": ..., "specificity_raw": ..., "f1": ..., ...}`
**Pass condition:** sensitivity ≥ 0.90 AND specificity_raw ≥ 0.80 AND f1 ≥ 0.85

## Approach Priority
1. Generate training data per allele via MHCflurry
2. Train Deep FFN per allele (R keras pipeline)
3. Build per-allele IEDB positive benchmark
4. Evaluate with binding_score ≥ 0.95 classifier
5. Iterate per-allele models that don't pass

## Max Runtime
8 hours (5 alleles × ~1.5h training + evaluation per allele)

## Next Step (Step B)
Pan-specific model merging all alleles with allele encoding in input
