# Mission: Fix IEDB Specificity

## Goal
Improve IEDB benchmark specificity from **75% → ≥85%** while maintaining sensitivity ≥90%.

## Current Baseline
| Metric | Value |
|--------|:-----:|
| Sensitivity | 93.9% (46/49) |
| Specificity | 75.0% (15/20) |
| F1 | 0.920 |
| False Positives | 5 (Poly-L, Poly-M, Poly-V, Poly-I, Poly-F) |

## Target
| Metric | Threshold |
|--------|:---------:|
| Specificity | **≥ 85%** |
| Sensitivity | **≥ 90%** |
| Pass | Both conditions met |

## Approaches (in priority order)

### 1. Post-hoc Threshold Calibration
- Raise SB probability threshold from default (argmax)
- Class-prior weighted threshold (different thresholds per class)
- No retraining required — fastest path

### 2. Class-Weighted Training
- Penalize NB→SB misclassification in loss function
- Weight SB > WB > NB in inverse proportion to class frequency
- Modify `peptide_mhc_binding_study.R` training section

### 3. Data Rebalancing
- Oversample minority classes or undersample majority
- SMOTE for peptide feature space
- Combined with class weights if needed

## Evaluator
**Script:** `.omc/autoresearch/iedb-specificity/evaluator.R`
**Command:** `Rscript evaluator.R <model_path.h5>`
**Output:** `{"pass": true/false, "specificity": X, "sensitivity": Y, ...}`

## Artifact Structure
```
.omc/autoresearch/iedb-specificity/
  mission.md        ← this file
  evaluator.R       ← self-contained evaluator
  runs/
    <run-id>/
      evaluations/
        iteration-0001.json
        iteration-0002.json
      decision-log.md
```
