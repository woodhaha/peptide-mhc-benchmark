# Peptide–MHC Binding Prediction with Deep Learning
## Project Summary — June 2026

---

### Project Overview

Deep learning models trained to classify 9-mer peptides as Strong Binder (SB), Weak Binder (WB), or Non-Binder (NB) to HLA-A\*02:01, reproducing and extending Jessen (2018) *"Deep Learning for Cancer Immunotherapy"*.

---

### Data Sources

| Source | Type | Peptides | Description |
|--------|------|----------|-------------|
| **MHCflurry 2.2.0** | ML-predicted binding | 5,088 | Real binding predictions (%rank < 0.5 → SB, 0.5–2.0 → WB, > 2.0 → NB) |
| **PSSM** | Biophysics-based | 5,007 | Position-specific scoring from crystallographic HLA-A\*02:01 motif |
| **Random Synthetic** | Baseline | 23,760 | Random labels with weak anchor bias |

---

### Model Performance

#### MHCflurry Data (Primary)

| Rank | Model | Accuracy | Macro F1 | NB F1 | WB F1 | SB F1 | Params |
|------|-------|----------|----------|-------|-------|-------|--------|
| 1 | **Deep FFN** | **91.9%** | **0.921** | 0.969 | 0.880 | 0.913 | 153K |
| 2 | FFN (Jessen) | 90.9% | 0.911 | 0.969 | 0.875 | 0.891 | 49K |
| 3 | CNN | 90.0% | 0.901 | 0.963 | 0.860 | 0.882 | 1.05M |
| 4 | ResNet-Style | 84.3% | 0.847 | 0.907 | 0.792 | 0.843 | 325K |
| 5 | LSTM | 83.3% | 0.836 | 0.935 | 0.752 | 0.822 | 24K |
| 6 | Random Forest | 81.1% | 0.814 | 0.927 | 0.723 | 0.793 | — |

#### 5-Fold Cross-Validation (MHCflurry)

| Fold | Accuracy | Macro F1 |
|------|----------|----------|
| 1 | 90.86% | 0.909 |
| 2 | 89.38% | 0.895 |
| 3 | 89.88% | 0.900 |
| 4 | 89.29% | 0.895 |
| 5 | 88.79% | 0.889 |
| **Mean ± SD** | **89.64% ± 0.79%** | **0.898 ± 0.008** |

---

### Data Source Comparison

| Data Source | Best Accuracy | Best Macro F1 | WB F1 | CV Mean |
|-------------|---------------|---------------|-------|---------|
| PSSM (biophysics) | 94.8% | 0.948 | 0.925 | 91.1% |
| MHCflurry (ML model) | 91.9% | 0.921 | 0.880 | 89.6% |
| Random Synthetic | 65.8% | 0.558 | 0.127 | 65.4% |

PSSM is easier because it uses linear additive positional scoring. MHCflurry captures non-linear inter-position interactions, making the WB class boundary genuinely fuzzy — matching real-world epitope prediction difficulty.

---

### New Peptide Predictions (Trained Model)

Model: **Deep FFN** (MHCflurry, 91.9%). Allele: **HLA-A\*02:01**.

| Peptide | p2-p9 | Pred | Confidence | P(SB) | P(WB) | P(NB) |
|---------|-------|------|------------|-------|-------|-------|
| `LMAFYLYEV` | M-V | **SB** | 100% | 100.0 | 0.00 | 0.00 |
| `FQTDRISYA` | Q-A | **SB** | 94.6% | 94.6 | 5.37 | 0.00 |
| `SLHLTNCFV` | L-V | **SB** | 97.6% | 97.6 | 2.42 | 0.00 |
| `YMFVILWVA` | M-A | **SB** | 100% | 100.0 | 0.04 | 0.00 |
| `LLLLLLLLL` | L-L | **SB** | 98.9% | 98.9 | 1.11 | 0.00 |
| `LLTDAQRIV` | L-V | **SB** | 72.1% | 72.1 | 27.8 | 0.03 |
| `AAAAAAAAA` | A-A | **NB** | 95.6% | 0.00 | 4.42 | 95.6 |
| `KGWGHSNGS` | G-S | **NB** | 100% | 0.00 | 0.01 | 100.0 |
| `DDDDDDDDD` | D-D | **NB** | 100% | 0.00 | 0.00 | 100.0 |
| `RRRRRRRRR` | R-R | **NB** | 100% | 0.00 | 0.00 | 100.0 |

---

### Model Architectures

| Model | Architecture | Description |
|-------|-------------|-------------|
| **FFN (Jessen)** | Dense(180)→Drop(0.4)→Dense(90)→Drop(0.3)→Dense(3) | Original Jessen 2018 |
| **CNN** | Conv2D(32,3×3)→Drop(0.25)→Flatten→FFN body | Jessen 2018 CNN |
| **LSTM** | LSTM(64)→Drop(0.3)→Dense(32)→Dense(3) | Sequential peptide model |
| **Deep FFN** | Dense(360)→BN→Drop(0.5)→Dense(180)→BN→Drop(0.4)→Dense(90)→Drop(0.3)→Dense(45)→Dense(3) | Best performer |
| **ResNet-Style** | Conv(32)→BN→ReLU→ResBlock×3(32→64→128, skip+projection)→GAP→Dense(128)→Dense(3) | Residual learning |
| **Random Forest** | 100 trees | Baseline classifier |

---

### Encoding

**BLOSUM62**: each amino acid → 20-dim substitution probability vector. Each 9-mer → 9×20×1 tensor. Input features: 180.

---

### Parallel Processing

- **doParallel**: 4-core cluster for cross-validation fold training
- **TensorFlow threading**: intra-op=4, inter-op=2
- **CV**: sequential fallback (robust)

---

### Output Files

```
models/
├── FFN_Deep.h5              ← Best model (91.9%)
├── FFN_Jessen2018.h5
├── CNN_Jessen2018.h5
├── LSTM_Extended.h5
├── ResNet_Style.h5
├── Random_Forest.rds
├── cv_FFN_CV_Fold_{1..5}.h5
└── model_manifest.csv

data/
├── feature.csv              ← 5,088 × 282 (FFN features)
├── feature_resnet.csv       ← 5,088 × 315 (ResNet features)
├── real_peptides.csv        ← MHCflurry-labeled data
├── new_peptide_predictions.csv
├── model_comparison.csv
└── cv_summary.csv

plots/
├── comparison_accuracy.png
├── comparison_per_class_f1.png
├── comparison_cv_folds.png
└── comparison_gap.png
```

---

### Key Findings

1. **Deep FFN is the best architecture** for peptide-MHC binding prediction (91.9% accuracy with MHCflurry labels)
2. **Simple architectures outperform complex ones** — extra depth helps, but ResNet/LSTM overcomplicate a positional task
3. **The WB (Weak Binder) class is the hardest** — F1 drops from 0.925 (PSSM) to 0.880 (MHCflurry), confirming the fuzzy boundary is the core challenge
4. **PSSM labels are ~3% easier** — linear biophysical model vs non-linear ML predictions
5. **Anchor residues (p2, p9) dominate** — poly-Leucine (L-L anchors) is confidently SB despite being repetitive
6. **Model generalizes correctly** — correctly identifies negative controls (poly-Asp, poly-Arg) and subtle cases (LLTDAQRIV at 72%/28% SB/WB)

---

### References

- Jessen, L.E. (2018). "Deep Learning for Cancer Immunotherapy." RStudio AI Blog.
- Reynisson, B. et al. (2020). "NetMHCpan-4.1 and NetMHCIIpan-4.0." *Nucleic Acids Research*, 48(W1), W449-W454.
- O'Donnell, T.J. et al. (2020). "MHCflurry 2.0: Improved Pan-Allele Prediction of MHC Class I-Presented Peptides." *Cell Systems*, 11(1), 42-48.
- Oluwagbemi et al. (2025). "AI-driven epitope prediction: a systematic review." *npj Vaccines*.

---

*Generated June 15, 2026 | Peptide Epitope Research Project*
