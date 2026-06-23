# Jessen (2018) — Deep Learning for Cancer Immunotherapy

## Metadata
- **Title:** Deep Learning for Cancer Immunotherapy
- **Author:** Leon Eyrich Jessen (Technical University of Denmark)
- **Published:** January 29, 2018
- **Venue:** RStudio AI Blog (`blogs.rstudio.com/ai/posts/2018-01-29-dl-for-cancer-immunotherapy`)
- **Citation:** Jessen, 2018
- **PDF:** `Literature/PDFs/Deep Learning for Cancer Immunotherapy.pdf`

## Status
This is the **direct methodological ancestor** of our project. Our manuscript's Jessen FFN baseline (citation [21]) is the 180→90→3 architecture described here.

## Key Findings

### Three-Model Comparison
| Model | Architecture | Test Accuracy | Notes |
|-------|-------------|:------------:|-------|
| **FFN** | 180→dropout(0.4)→90→dropout(0.3)→3 | **~95%** | 49,143 params; softmax output |
| **CNN** | Conv2D(32,3×3)→dropout(0.25)→FFN body | ~92% | 1M+ params; worse than FFN |
| **RF** | 100 trees on flattened BLOSUM62 | ~82% | Traditional ML baseline |

### Data Pipeline
```
1,000,000 random 9-mers (uniform AA sampling)
    → netMHCpan-4.0 binding prediction (HLA-A*02:01)
    → Classify: SB (rank<0.5%), WB (0.5%≤rank<2%), NB (rank≥2%)
    → Downsample to balanced: n(SB)=n(WB)=n(NB)=7,920 → 23,760 total
    → 90/10 train/test split
    → BLOSUM62 encoding: each peptide → 9×20 "image" matrix
    → Train FFN: 150 epochs, batch=50, 20% val split, RMSprop
```

### Encoded Peptide "Images"
Each 9-mer is encoded as a 9×20 matrix where:
- Row = position in peptide (1–9)
- Column = BLOSUM62 substitution score for that amino acid → each of 20 AAs
- Normalized from BLOSUM62 matrix (scores range roughly −4 to +11)
- Flattened to 180-dim vector for FFN input

### Sequence Logo Analysis
Strong binders show clear anchor motifs at p2 (L,M,I,V) and p9 (V,L) — matching established HLA-A*02:01 binding motif.

### Author's Key Observations
1. Deep learning captures non-linear pMHC interactions better than RF
2. CNN underperformance was surprising: peptide "images" lack spatial edges — p2 is always at column position 2, p9 at 9, so convolution across adjacent columns doesn't add information the way it does in MNIST
3. Data source matters: netMHCpan labels make this "a model of a model"
4. Ensemble predictions (5-fold CV × 5 models) recommended for production
5. Multi-source data (binding affinity + MS eluted ligands) improves generalization

## Relevance to Our Project

| Aspect | Jessen 2018 | Our Study (Zhou, Bai, Xu, Zhou, Han) |
|--------|-----------|--------------------------------------|
| Training peptides | 1,000,000 random → 23,760 balanced | 100,000 random → 5,088 balanced |
| Label source | netMHCpan-4.0 only | MHCflurry 2.2.0 + PSSM + Random |
| Architectures | 3 (FFN, CNN, RF) | 6 (+ Deep FFN, LSTM, ResNet) |
| Peptide encoding | BLOSUM62 only | BLOSUM62 + ESM-2 t6/t12/mean-pooled |
| External validation | None | IEDB 49 epitopes + 20 controls |
| Label quality analysis | Qualitative mention | Quantitative: PSSM 94.8% vs MHCflurry 91.9% |
| Protein scanning | None | 10 proteins + hotspot mutations |
| Structural biology | None | HDOCKlite docking + 10ns MD (OpenMM) |
| Code release | Blog post + GitHub gist | Full reproducible pipeline |

## Key Quotes

> "It is evident that the deep learning models capture the information in the system much better than the random forest model. However, the CNN model didn't perform as well as the straightforward FFN. This illustrates one of the pitfalls of deep learning - blind alleys."

> "There are a huge number of architectures available, and when combined with hyperparameter tuning the potential model space is breathtakingly large. To increase the likelihood of finding a good architecture and the right hyper-parameters it is important to know and understand the data you are modeling."

> "We are very careful to avoiding overfitting as this of course decreases the models extrapolation performance."

## Limitations (Author-Acknowledged)
1. Training labels are computational (netMHCpan), not experimental
2. Single HLA allele (A*02:01)
3. No external experimental validation
4. Random peptide generation (uniform AA frequencies) differs from natural proteome distribution
5. CNN architecture not optimized — single convolutional layer, no systematic hyperparameter search

## Our Study's Extensions Beyond Jessen
1. **Label quality quantification**: PSSM vs MHCflurry vs random — isolating the label effect
2. **ESM-2 embeddings**: Protein language model features vs hand-crafted BLOSUM62
3. **Additional architectures**: Deep FFN (batchnorm), LSTM, ResNet-style
4. **External validation**: IEDB benchmark (49 epitopes + 20 controls)
5. **Translational application**: Protein scanning + cancer mutation analysis (p53 R248W, KRAS G12V)
6. **Structural hypotheses**: HDOCKlite docking → K-E63 salt bridge hypothesis
7. **Systematic label quality analysis**: Quantifying the ~3pp performance gap

## R Code Reference

Key packages: `keras`, `tensorflow`, `tidyverse`, `ggseqlogo`, `PepTools`
Data URL: `https://git.io/vb3Xa`
Model: `keras_model_sequential()` + `layer_dense()` + `layer_dropout()` + `compile(loss='categorical_crossentropy', optimizer=optimizer_rmsprop())` + `fit(epochs=150, batch_size=50, validation_split=0.2)`

---

- [x] PDF fully read (14 pages)
- [x] Integrated into literature_matrix.md as Entry #0
- [x] Cross-referenced with manuscript ([21])