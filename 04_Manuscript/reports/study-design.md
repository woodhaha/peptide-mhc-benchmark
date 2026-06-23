# Study Design: Deep Learning for Peptide–MHCI Binding Prediction
## Mimicking & Extending Jessen (2018) for Epitope Research

---

### 1. Background & Rationale

Jessen (2018) demonstrated that a feed-forward neural network (FFN) could achieve ~95% accuracy classifying 9-mer peptides as Strong Binder (SB), Weak Binder (WB), or Non-Binder (NB) to HLA-A\*02:01, using netMHCpan-4.0 labels as ground truth. This study reproduces that pipeline and extends it with:

- 5-fold cross-validation (referenced in Jessen's conclusion as standard practice)
- Ensemble prediction (wisdom-of-the-crowd)
- Multiple encoding schemes (BLOSUM62 + one-hot + physicochemical properties)
- Performance comparison across FFN, CNN, LSTM, and Random Forest
- Generalization testing to additional HLA alleles

### 2. Research Questions

| RQ | Question | Metric |
|----|----------|--------|
| RQ1 | Can we reproduce Jessen's FFN ~95% accuracy on HLA-A\*02:01? | Accuracy, F1, AUC |
| RQ2 | Does 5-fold CV ensemble outperform single train/test split? | Accuracy, F1 |
| RQ3 | Which encoding (BLOSUM62 vs one-hot vs AAindex) performs best? | Accuracy across encodings |
| RQ4 | Does LSTM architecture improve over FFN for sequential peptide data? | Accuracy, F1 |
| RQ5 | How well does the model generalize to other common HLA alleles (A\*01:01, B\*07:02)? | Accuracy, F1 |

### 3. Methods Summary (adapted from Jessen 2018)

#### 3.1 Data Generation
- Generate N random 9-mer peptides (N = 100,000–1,000,000)
- Submit to netMHCpan-4.1 for binding prediction against target HLA alleles
- Labels: SB (%Rank < 0.5), WB (%Rank 0.5–2.0), NB (%Rank > 2.0)
- Balance by downsampling to equal class sizes
- Train/test split: 90/10 (RQ1) + 5-fold CV (RQ2)

#### 3.2 Peptide Encoding
- **BLOSUM62**: 9×20 matrix per peptide (original method)
- **One-hot**: 9×20 binary matrix
- **AAindex**: 9×N physicochemical property matrix (optional extension)

#### 3.3 Model Architectures

| Model | Architecture | Parameters |
|-------|-------------|-----------|
| **FFN** | Dense(180)→Dropout(0.4)→Dense(90)→Dropout(0.3)→Dense(3, softmax) | ~49K |
| **CNN** | Conv2D(32, 3×3)→Dropout(0.25)→Flatten→FFN | ~82K |
| **ResNet-Style** | Conv(32)→BN→ReLU → Residual Block×3 (32→64→128 filters, skip connections + 1×1 projections) → GlobalAvgPool → Dense(128)→Dropout → Softmax | ~150K |
| **LSTM** | LSTM(64)→Dropout(0.3)→Dense(32)→Dense(3, softmax) | ~45K |
| **RF** | 100 trees, default mtry | ~100 trees |

#### 3.4 Training
- Loss: categorical cross-entropy
- Optimizer: RMSprop (original) + Adam (extension)
- Epochs: 150, batch size: 50
- Validation split: 0.2 (internal holdout)
- Early stopping with patience = 10 (extension)

#### 3.5 Evaluation
- Primary: Accuracy on held-out 10% test set
- Secondary: F1-score (per class + macro), AUROC (one-vs-rest)
- Confusion matrix visualization
- Per-class precision/recall

### 4. R Implementation Plan

The companion script `peptide_mhc_binding_study.R` implements:

1. **Setup & Dependencies** — keras3/tensorflow, tidyverse, PepTools, randomForest
2. **Data Module** — peptide generation with unique `peptide_id`, netMHCpan prediction (or pre-computed data load), balancing
3. **Encoding Module** — BLOSUM62, one-hot encodings
4. **Model Module** — FFN, CNN, ResNet-Style, LSTM, Deep FFN, RF training functions
5. **Evaluation Module** — metrics, confusion matrices, ROC curves
6. **Feature Extraction** — raw BLOSUM features + penultimate-layer learned features → `feature.csv`
7. **Model Persistence** — all models saved to `models/*.h5` + `models/model_manifest.csv`
8. **Extension Module** — 5-fold CV, ensemble, multi-allele testing

### 5. Expected Outcomes

- Reproduction of ~95% FFN accuracy on HLA-A\*02:01 (RQ1)
- CV ensemble > single split performance (RQ2)
- BLOSUM62 > one-hot encoding (RQ3)
- LSTM ≈ FFN (sequential context helps marginally) (RQ4)
- Performance drop on alternate alleles proportional to training data availability (RQ5)

### 6. References

- Jessen, L.E. (2018). "Deep Learning for Cancer Immunotherapy." RStudio AI Blog.
- Reynisson, B. et al. (2020). "NetMHCpan-4.1 and NetMHCIIpan-4.0." *Nucleic Acids Research*.
- Oluwagbemi et al. (2025). "AI-driven epitope prediction: a systematic review." *npj Vaccines*.
