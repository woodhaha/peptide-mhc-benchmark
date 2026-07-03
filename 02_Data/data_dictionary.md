# Data Dictionary

## Raw Data (`02_Data/raw/`)

| File | Description | Source |
|------|-------------|--------|
| `feature.csv` | Original peptide features (BLOSUM encoding, physicochemical properties) | Custom feature engineering |
| `real_peptides.csv` | Real peptide dataset for prediction | Project dataset |

## Cleaned/Processed Data (`02_Data/cleaned/`)

| File | Description |
|------|-------------|
| `feature_resnet.csv` | Features processed for ResNet model input |
| `iedb_benchmark_results.csv` | IEDB benchmark evaluation results across models |
| `model_comparison.csv` | Performance comparison across all model architectures |
| `mutation_scan_results.csv` | Alanine scanning mutagenesis predictions |
| `new_peptide_predictions.csv` | Binding predictions for novel peptide candidates |
| `protein_epitope_scan.csv` | Full protein epitope scanning results |
| `protein_epitope_scan_extended.csv` | Extended protein epitope scan (additional proteins) |
| `cv_summary.csv` | Cross-validation summary statistics |
| `top20_epitope_candidates.csv` | Top 20 highest-confidence epitope candidates |

## Archives
- `data.tar.gz` — Backup archive of original data files
