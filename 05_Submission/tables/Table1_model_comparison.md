# Table 1: Model Performance Comparison (MHCflurry-Labeled Data)

| Model | Architecture | Parameters | Accuracy | Macro F1 | NB F1 | WB F1 | SB F1 |
|-------|-------------|:----------:|:--------:|:--------:|:-----:|:-----:|:-----:|
| Deep FFN | 360→BN→180→BN→90→45→3 | 152,823 | **91.9%** | **0.921** | 0.969 | 0.880 | 0.913 |
| FFN (Jessen 2018) | 180→90→3 | 49,143 | 90.9% | 0.911 | 0.969 | 0.875 | 0.891 |
| CNN | Conv2D(32,3×3)→Flatten→FFN | 1,053,863 | 90.0% | 0.901 | 0.963 | 0.860 | 0.882 |
| ResNet-Style | Conv stem + 3 residual blocks (32/64/128) + GAP | 324,931 | 84.3% | 0.847 | 0.907 | 0.792 | 0.843 |
| LSTM | LSTM(64)→Dense(32)→3 | 23,939 | 83.3% | 0.836 | 0.935 | 0.752 | 0.822 |
| Random Forest | 100 trees, BLOSUM62 features | — | 81.1% | 0.814 | 0.927 | 0.723 | 0.793 |

**Note:** NB = Non-Binder, WB = Weak Binder, SB = Strong Binder. Training on 4,579 peptides; evaluation on 509 holdout peptides. All neural networks use BLOSUM62 encoding (180-dim input). Deep FFN is the best-performing model and was used for all downstream analyses.
