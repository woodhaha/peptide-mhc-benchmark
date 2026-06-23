# Deep Interview Spec: Real Negative Benchmark

## Metadata
- Interview ID: di-20260619-001
- Rounds: 4
- Final Ambiguity Score: 14%
- Type: brownfield
- Generated: 2026-06-19
- Threshold: 0.2
- Threshold Source: default
- Status: PASSED
- Mode: autoresearch

## Clarity Breakdown
| Dimension | Score | Weight | Weighted |
|-----------|-------|--------|----------|
| Goal Clarity | 0.95 | 0.35 | 0.333 |
| Constraint Clarity | 0.70 | 0.25 | 0.175 |
| Success Criteria | 0.90 | 0.25 | 0.225 |
| Context Clarity | 0.85 | 0.15 | 0.128 |
| **Total Clarity** | | | **0.860** |
| **Ambiguity** | | | **0.140** |

## Topology
| Component | Status | Description | Coverage |
|-----------|--------|-------------|----------|
| real-negative-benchmark | **active** | Replace homopolymer-only negatives with 4-source diverse negatives | ✅ This spec |
| multi-allele | active | Extend to multiple MHC class I alleles | → Next mission |
| plm-encoding | active | Replace BLOSUM62 with PLM embeddings | → Subsequent |
| p9-false-negatives | deferred | Non-canonical p9 FN fix | User deferred |
| submission-readiness | deferred | Formatting/finalization | Not auto-evaluable |

## Goal
Build a diverse negative peptide benchmark (4 sources: random 9-mers, IEDB non-binders, MHCflurry NB predictions, natural protein sliding windows) totaling ~200 negatives against 49 IEDB-validated positives. Run the Deep FFN model against this benchmark and optimize to achieve ≥90% sensitivity, ≥80% specificity (raw, without homopolymer filter), and ≥0.85 F1.

## Constraints
- Homopolymer negatives retained as a control subset (20 peptides), reported separately
- Homopolymer filter (≤2 unique AAs → NB) retained as post-hoc step; output BOTH raw and filtered metrics
- Pass condition based on RAW metrics (no filter)
- Training data must NOT include any benchmark negative peptides (data leakage prevention)
- Evaluator must be a standalone Python script runnable via `python evaluator.py <model_path.h5>`

## Non-Goals
- Does NOT retrain with new architectures (CNN/LSTM/ResNet) — use Deep FFN baseline
- Does NOT change the positive set (49 IEDB epitopes remain)
- Does NOT fix the 3 non-canonical p9 false negatives (deferred)
- Does NOT extend to multi-allele (separate mission)

## Acceptance Criteria
- [ ] Python evaluator script that: (a) loads a Keras .h5 model, (b) builds combined benchmark (49 POS + ~200 NEG from 4 sources), (c) predicts with BLOSUM62 encoding, (d) outputs JSON `{pass, sensitivity, specificity_raw, specificity_filtered, f1, tp, fp, tn, fn, hp_overridden, n_pos, n_neg}`
- [ ] Negative set contains peptides from all 4 sources (random, IEDB, MHCflurry, natural protein)
- [ ] No benchmark peptide appears in training data
- [ ] Baseline Deep FFN (FFN_Deep.h5) achieves sensitivity ≥ 90% on the new benchmark
- [ ] Baseline Deep FFN (FFN_Deep.h5) achieves specificity_raw ≥ 80% on the new benchmark
- [ ] Baseline Deep FFN (FFN_Deep.h5) achieves F1 ≥ 0.85 on the new benchmark
- [ ] Evaluator JSON output is machine-parseable (valid JSON on stdout, single line)

## Assumptions Exposed & Resolved
| Assumption | Challenge | Resolution |
|------------|-----------|------------|
| 100% specificity means model works | What about non-homopolymer negatives? | Need diverse negatives — 4 sources |
| R evaluator is fine | Python more flexible for MHCflurry integration? | Python evaluator chosen |
| Exploratory vs optimization | Benchmark artifact only, or performance target? | Optimization — model must hit thresholds |
| What thresholds? | 80% spec reasonable against real negatives? | ≥90% sens, ≥80% spec, ≥0.85 F1 |
| Homopolymer filter should stay? | Filter hides model weakness | Keep but report BOTH raw and filtered |

## Technical Context
- **Existing artifacts:** `evaluator_v2.R` (homopolymer filter logic, BLOSUM62 encoder, JSON output format), `benchmark_iedb_epitopes_v2.R` (benchmark construction)
- **Model:** `Analysis/FFN_Deep.h5` — Deep FFN, 152,823 params, 91.9% accuracy
- **Training data:** `Data/raw/real_peptides.csv` — 5,088 MHCflurry-labeled 9-mers (balanced SB/WB/NB)
- **BLOSUM62 matrix:** Hardcoded 20×20 in R scripts; must be replicated in Python
- **MHCflurry:** Available in Python environment (v2.2.0) for NB label verification
- **IEDB API:** Available for fetching known non-binders
- **Reference peptides:** 49 positive epitopes defined in `benchmark_iedb_epitopes.R`

## Ontology (Key Entities)
| Entity | Type | Fields | Relationships |
|--------|------|--------|---------------|
| NegativePeptide | core domain | sequence, source_type, true_label | belongs to Benchmark |
| Benchmark | core domain | positive_set, negative_set, metrics | evaluates Model |
| Model | core domain | architecture, encoding, weights_path | predicts BindingClass |
| Evaluator | supporting | script_path, threshold_config | runs Benchmark against Model |
| HomopolymerFilter | supporting | unique_aa_threshold, override_class | post-processes Model predictions |

## Evaluator
- **Language:** Python 3
- **Command:** `python .omc/autoresearch/real-negative-benchmark/evaluator.py <model_path.h5>`
- **Output:** Single-line JSON to stdout: `{"pass": bool, "sensitivity": float, "specificity_raw": float, "specificity_filtered": float, "f1": float, "tp": int, "fp": int, "tn": int, "fn": int, "hp_overridden": int, "n_pos": int, "n_neg": int}`
- **Pass condition:** sensitivity ≥ 0.90 AND specificity_raw ≥ 0.80 AND f1 ≥ 0.85
- **Exit code:** 0 on pass, 1 on fail

## Interview Transcript
<details>
<summary>Full Q&A (4 rounds)</summary>

### Round 1
**Q:** 负样本策略 — 随机9-mer / IEDB已知非结合肽 / MHCflurry负样本 / 天然蛋白滑窗负样本？
**A:** 全部四条路径
**Ambiguity:** 51% (Goal: 0.7, Constraints: 0.3, Criteria: 0.2, Context: 0.8)

### Round 2
**Q:** 探索型（构建+评估即可）还是优化型（有性能门槛）？
**A:** Model（优化型）
**Ambiguity:** 42% (Goal: 0.85, Constraints: 0.3, Criteria: 0.35, Context: 0.8)

### Round 3
**Q:** 特异性门槛推荐 ≥80%，敏感性 ≥90%，F1 ≥0.85？
**A:** OK
**Ambiguity:** 22% (Goal: 0.9, Constraints: 0.45, Criteria: 0.9, Context: 0.85)

### Round 4
**Q:** Evaluator — 扩展evaluator_v2.R还是Python？
**A:** Python evaluator
**Ambiguity:** 14% ✅ (Goal: 0.95, Constraints: 0.7, Criteria: 0.9, Context: 0.85)

</details>
