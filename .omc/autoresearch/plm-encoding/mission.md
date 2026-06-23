# Mission: PLM/Transformer Encoding for Peptide-MHC Binding Prediction

## Goal
Replace BLOSUM62 encoding with protein language model (PLM) embeddings (ESM-2) and evaluate whether learned evolutionary representations outperform hand-crafted substitution matrices for peptide-MHC binding prediction.

## Baseline (BLOSUM62 + Deep FFN)
| Metric | Value |
|--------|:-----:|
| Accuracy (A*02:01) | 91.9% |
| IEDB Sensitivity | 93.9% |
| IEDB Specificity (real negs) | 93.5% |
| IEDB F1 | 0.860 |

## Target
| Metric | Threshold |
|--------|:---------:|
| Accuracy | ≥ 91.9% (match or exceed BLOSUM62) |
| IEDB Sensitivity | ≥ 90% |
| IEDB Specificity_raw | ≥ 80% |
| IEDB F1 | ≥ 0.85 |
| Pass | All conditions met OR accuracy > 91.9% |

## Approach
1. **ESM-2 embedding extraction:** Use `facebook/esm2_t6_8M_UR50D` (8M params, 320-dim embeddings) to encode all peptides
2. **Train Deep FFN** with ESM-2 embeddings (320-dim input instead of 180-dim BLOSUM62)
3. **Compare** against BLOSUM62 baseline on:
   - Test accuracy
   - IEDB benchmark (49 positives + diverse negatives)
   - Per-class F1 scores
4. If t6 (8M) shows promise, scale to `esm2_t12_35M_UR50D` (480-dim) or `esm2_t30_150M_UR50D` (640-dim)

## Architecture
- Input: 320-dim ESM-2 mean pooling embedding
- Deep FFN: Dense(360)-BN-Drop(0.4)-Dense(180)-BN-Drop(0.4)-Dense(90)-Drop(0.3)-Dense(45)-Dense(3)
- Training: Adam lr=0.001, batch_size=64, early_stopping patience=15
- Data: Same 5,088 MHCflurry-labeled peptides used for BLOSUM62 baseline

## Evaluator
**Path:** `.omc/autoresearch/plm-encoding/evaluator.py`
**Command:** `python evaluator.py <model_path.h5> <encoding_type>`
**Encoding types:** `blosum62` | `esm2_t6`
**Output:** `{"pass": bool, "accuracy": X, "sensitivity": X, "specificity_raw": X, "f1": X, ...}`
**Pass condition:** sensitivity ≥ 0.90 AND specificity_raw ≥ 0.80 AND f1 ≥ 0.85

## Artifact Structure
```
.omc/autoresearch/plm-encoding/
  mission.md
  evaluator.py
  extract_esm_embeddings.py
  runs/run-001/
    evaluations/
    decision-log.md
```
