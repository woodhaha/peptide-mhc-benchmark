# Results Summary — v2.1 (2026-06-23)

## Model Performance (6 architectures, MHCflurry labels)

| Model | Accuracy | Macro F1 | NB F1 | WB F1 | SB F1 | Params |
|-------|:--------:|:--------:|:-----:|:-----:|:-----:|:------:|
| **LSTM** | **91.5%** | **0.918** | 0.968 | 0.880 | 0.904 | 23,939 |
| Deep FFN | 90.8% | 0.910 | 0.962 | 0.869 | 0.899 | 152,823 |
| CNN | 90.2% | 0.904 | 0.956 | 0.862 | 0.895 | 1,053,863 |
| ResNet | 89.8% | 0.901 | 0.968 | 0.856 | 0.880 | 53,059 |
| FFN (Jessen) | 89.2% | 0.895 | 0.959 | 0.849 | 0.877 | 49,143 |
| Random Forest | 79.4% | 0.801 | 0.925 | 0.702 | 0.775 | -- |

## PSSM Comparison (Label-Model Interaction)

| Model | MHCflurry | PSSM | Gap |
|-------|:---------:|:----:|:---:|
| Random Forest | 79.4% | 85.3% | +5.9pp |
| FFN (Jessen) | 89.2% | 84.3% | -4.9pp |
| CNN | 90.2% | 84.3% | -5.9pp |
| Deep FFN | 90.8% | 81.4% | -9.4pp |
| ResNet | 89.8% | 75.5% | -14.3pp |
| LSTM | 91.5% | 69.6% | -21.9pp |

**Key finding**: PSSM ≠ universally better. Simple model (RF) improves; complex models degrade.
Correct interpretation: label-model interaction, not one-way "label quality."

## Cross-Validation

- 5-fold CV: **91.9% ± 0.7%**
- Range: 90.9% – 92.7%

## IEDB Benchmark

- Sensitivity: **93.9%** (46/49)
- Specificity: 75.0% (15/20)
- 3 FNs: ILRGSVAHK, QYDPVAALF, TLGIVCPIC (all non-canonical p9)
- 5 FPs: Poly-L/M/V/I/F (canonical anchors → bind MHC, not immunogenic)

## Mutation Scan

| Mutation | Peptide | Effect | ΔScore |
|----------|---------|--------|:------:|
| KRAS G12V | YKLVVVGAV | CREATED | +0.499 |
| p53 R248W | MNWRPILTI | CREATED | +0.429 |
| p53 R249S | NRSPILTII | CREATED | +0.303 |

## Docking

- Template: 1DUZ (HLA-A\*02:01 + LLFGYPVYV, 1.80Å)
- P2:CA → E63:CD = 4.7Å
- K sidechain reach = 7.6Å
- K-E63 salt bridge: STRUCTURALLY PLAUSIBLE

## Environment

- Python 3.13.9, TensorFlow 2.21, Keras 3
- OpenMM 8.5.2, Amber14 forcefield
- R 4.6.0, ggsci (Nature NPG palette)

## Output Files

| File | Contents |
|------|---------|
| `02_Data/cleaned/model_comparison.csv` | 6-model performance |
| `02_Data/cleaned/pssm_comparison.csv` | PSSM vs MHCflurry gap |
| `02_Data/cleaned/cv_summary.csv` | 5-fold CV |
| `02_Data/cleaned/iedb_benchmark_results_v2.csv` | IEDB benchmark |
| `02_Data/cleaned/protein_epitope_scan_extended.csv` | 3,357 9-mers |
| `02_Data/cleaned/mutation_scan_results.csv` | 144 mutation windows |
| `02_Data/cleaned/top20_epitope_candidates.csv` | Top 20 epitopes |
| `docking/energy_results.json` | Docking energy analysis |
| `03_Analysis/figures/` | 32 Nature figures (PNG+PDF) |
| `05_Submission/figures/` | Submission PDFs |
