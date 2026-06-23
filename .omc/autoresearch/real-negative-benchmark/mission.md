# Mission: Real Negative Benchmark

## Goal
Build a diverse negative peptide benchmark (4 sources) and optimize the Deep FFN model to pass against real negatives — not just homopolymers.

## Baseline
| Metric | Current (homopolymer-only) | Realistic Expectation |
|--------|:--------------------------:|:---------------------:|
| Sensitivity | 93.9% (46/49) | Will hold or drop slightly |
| Specificity | 100% (20/20 homopolymers) | Likely 70-85% against real negatives |
| F1 | 0.968 | Likely 0.85-0.92 |
| False Positives | 0 | Expected: 30-60 against real negatives |

## Target
| Metric | Threshold |
|--------|:---------:|
| Sensitivity | **≥ 90%** |
| Specificity (raw, no filter) | **≥ 80%** |
| F1 | **≥ 0.85** |
| Pass | All three conditions met |

## Negative Set Design (4 sources, ~200 total)

| Source | Count | Description |
|--------|:-----:|-------------|
| Random 9-mers | 50 | Randomly generated, MHCflurry-verified NB (%rank > 2%) |
| IEDB non-binders | ~30 | Experimentally validated negatives from IEDB |
| MHCflurry NB | 50 | Predicted NB from MHCflurry (%rank > 2%), matched to positive label source |
| Natural protein windows | 50 | 9-mer sliding windows from non-immunogenic proteins (albumin, collagen), filtered to remove predicted binders |
| Homopolymers (control) | 20 | Retained from original benchmark, reported separately |

## Positive Set (unchanged)
49 IEDB-validated HLA-A*02:01 epitopes from viral, tumor, and melanoma antigens.

## Evaluator
**Path:** `.omc/autoresearch/real-negative-benchmark/evaluator.py`
**Command:** `python evaluator.py <model_path.h5>`
**Output:** `{"pass": true/false, "sensitivity": X, "specificity_raw": X, "specificity_filtered": X, "f1": X, "tp": N, "fp": N, "tn": N, "fn": N, "hp_overridden": N, "n_pos": N, "n_neg": N}`
**Pass condition:** sensitivity ≥ 0.90 AND specificity_raw ≥ 0.80 AND f1 ≥ 0.85

## Approaches (priority order)

### 1. Baseline measurement
- Run Deep FFN (FFN_Deep.h5) against new benchmark
- Record raw and filtered metrics
- Identify failure modes (which negative sources cause FPs?)

### 2. Post-hoc threshold calibration
- Grid search SB probability threshold
- Class-specific thresholds
- No retraining required

### 3. Class-weighted retraining
- Increase NB class weight to penalize false positives
- Vary weight from 2.0 to 10.0

### 4. Architecture modifications
- Add dropout, adjust layer sizes
- Test if ResNet or CNN handle diverse negatives better

### 5. Feature engineering
- Add physicochemical properties beyond BLOSUM62
- Position-specific features for anchor residues

## Artifact Structure
```
.omc/autoresearch/real-negative-benchmark/
  mission.md        ← this file
  evaluator.py      ← Python evaluator (created in iteration 1)
  runs/
    run-001/
      evaluations/
        iteration-0001.json
        iteration-0002.json
      decision-log.md
```

## Max Runtime
4 hours per run (covers baseline + up to 6 improvement iterations)
