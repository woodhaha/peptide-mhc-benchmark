# Deep Interview Spec: Fix IEDB Specificity

## Metadata
- Interview ID: f391c9c0-autoresearch-001
- Rounds: 4 (Topology + 3)
- Final Ambiguity Score: 16.9%
- Type: brownfield
- Generated: 2026-06-15
- Threshold: 20%
- Threshold Source: default
- Status: PASSED
- Evaluator Script: `.omc/autoresearch/iedb-specificity/evaluator.R`

## Clarity Breakdown
| Dimension | Score | Weight | Weighted |
|-----------|:-----:|:------:|:--------:|
| Goal Clarity | 0.88 | 35% | 0.308 |
| Constraint Clarity | 0.70 | 25% | 0.175 |
| Success Criteria | 0.85 | 25% | 0.213 |
| Context Clarity | 0.90 | 15% | 0.135 |
| **Total Clarity** | | | **0.831** |
| **Ambiguity** | | | **16.9%** |

## Topology
| Component | Status | Description | Coverage |
|-----------|--------|-------------|----------|
| Training modification | active | Class-weighted loss and/or data rebalancing in Deep FFN | Modify `peptide_mhc_binding_study.R` to penalize NB→SB |
| Threshold calibration | active | Post-hoc decision boundary tuning | Adjust SB/WB/NB probability thresholds |
| Re-evaluation | active | Self-contained R evaluator | Script at `.omc/autoresearch/iedb-specificity/evaluator.R` |

## Goal
Improve IEDB benchmark specificity from 75% → ≥85% while maintaining sensitivity ≥90%, using class-weighted training, data rebalancing, and/or threshold calibration on the Deep FFN model. The evaluator provides machine-readable pass/fail output for iterative refinement.

## Constraints
- Deep FFN architecture (360-180-90-45-3) preserved — no architecture changes
- BLOSUM62 encoding preserved
- IEDB benchmark: 49 validated HLA-A\*02:01 epitopes + 20 negative controls (9-mers only)
- Evaluator in R (keras, tensorflow, tidyverse, pROC)
- Windows 11, CPU only, R 4.6.0
- Model must be savable/loadable via HDF5

## Non-Goals
- Adding new benchmarking datasets
- Supporting additional HLA alleles
- Changing model architecture
- Modifying peptide encoding strategy
- Re-running the full training pipeline for comparison models

## Acceptance Criteria
- [ ] Evaluator script exists and runs successfully: `Rscript evaluator.R <model.h5>`
- [ ] Evaluator outputs valid JSON with keys: pass, specificity, sensitivity, f1, accuracy, tp, fp, tn, fn
- [ ] Pass condition: specificity ≥ 0.85 AND sensitivity ≥ 0.90
- [ ] Baseline (current Deep FFN) evaluator run confirms specificity 0.75, sensitivity 0.939
- [ ] At least one modified model achieves pass = true
- [ ] Each iteration logged in `.omc/autoresearch/iedb-specificity/runs/<run-id>/`

## Assumptions Exposed & Resolved
| Assumption | Challenge | Resolution |
|------------|-----------|------------|
| Evaluator design: CSV vs self-contained | Which approach? | Self-contained R script that loads model + runs benchmark |
| Pass condition definition | What thresholds? | Specificity ≥ 85% AND sensitivity ≥ 90% |
| Evaluator language | R or Python? | R (consistent with existing pipeline) |

## Technical Context
- **Main pipeline**: `Analysis/R_scripts/peptide_mhc_binding_study.R` (1,700+ lines)
- **Benchmark script**: `Analysis/R_scripts/benchmark_iedb_epitopes.R` (basis for evaluator)
- **Best model**: `Analysis/FFN_Deep.h5` (91.9% accuracy, 75% specificity on IEDB)
- **Benchmark data**: `Data/cleaned/iedb_benchmark_results.csv` (69 peptides)
- **Training data**: `Data/raw/real_peptides.csv` (5,088 peptides)
- **Environment**: R 4.6.0, keras 2.16.1, tensorflow 2.17.0, tidyverse 2.0.0

## Ontology (Key Entities)
| Entity | Type | Fields | Relationships |
|--------|------|--------|---------------|
| Deep FFN Model | core domain | 360-180-90-45-3, BLOSUM62 input, MHCflurry trained | produces predictions for benchmark peptides |
| IEDB Benchmark | core domain | 49 epitopes, 20 controls, 69 peptides total | evaluated by Deep FFN Model |
| Evaluator Script | supporting | R script, self-contained, accepts model path, JSON output | loads Model, reads benchmark, outputs pass/fail |
| Pass/Fail Threshold | supporting | Specificity ≥ 85%, Sensitivity ≥ 90% | defines success for Evaluator Script |
| Autoresearch Loop | external system | iterates modifications, calls evaluator, logs results | calls Evaluator Script, modifies Deep FFN Model |

## Interview Transcript
<details>
<summary>Full Q&A (4 rounds)</summary>

### Round 0 (Topology)
**Q:** Confirm 3 top-level components: Training modification, Threshold calibration, Re-evaluation
**A:** Looks right, proceed

### Round 1
**Q:** What specific, machine-readable pass/fail condition? (Component: Re-evaluation, Dimension: Success Criteria)
**A:** Specificity ≥ 85% AND sensitivity ≥ 90%

### Round 2
**Q:** Evaluator design — B (read CSV) vs A (wrap benchmark)? Language?
**A:** B — read the CSV, in R

### Round 3
**Q:** Self-contained evaluator (loads model + runs benchmark) or assumes CSV exists?
**A:** A — self-contained evaluator that takes the model path

</details>
