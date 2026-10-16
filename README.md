# Peptide-MHC Binding Prediction Benchmark

> **A systematic benchmark showing that benchmark composition can manufacture false conclusions — and that the non-canonical tail defines the real prediction frontier.**

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
![R 4.6.0+](https://img.shields.io/badge/R-4.6.0%2B-blue)
![Python 3.13+](https://img.shields.io/badge/Python-3.13%2B-blue)

Systematic benchmarking of six deep-learning architectures (with multi-seed estimates) for HLA-A\*02:01 peptide-MHC-I binding prediction, plus a multi-allele extension, external validation against 49 IEDB epitopes, and a data-integrity re-analysis that corrected three benchmark artifacts.

**Manuscript**: `06_Bioinformatics_Submission/manuscript_bioinformatics_v2.md` (16 references, ~4,400 words, 9 figures, 7 tables) · companion `06_Bioinformatics_Submission/companion_manuscript.md` (PLOS Comput Biol target)

---

## Core Findings

| Model | Accuracy (5-seed mean ± SD) |
|-------|:---------------------------:|
| CNN | **91.9% ± 0.9%** |
| Deep FFN | 90.6% ± 1.2% |
| FFN (Jessen) | 89.8% ± 1.1% |
| LSTM | 89.7% ± 0.9% |
| ResNet | 87.6% ± 1.8% |
| Random Forest | 79.7% ± 0.4% |

Architecture has a modest effect (top four within ~2 pp; CNN modestly best). Random labels fall to chance (~33%), confirming no noise overfitting. The genuinely hard problem is **non-canonical binding**: it plateaus at ~0.65 sensitivity across three independent sequence predictors (15-allele panel), defining the current frontier.

---

## Key Results

| Metric | Value |
|--------|:-----:|
| Best architecture | CNN (1.05M parameters) |
| IEDB sensitivity | **93.9%** (46/49) |
| IEDB ROC AUC | **0.947** |
| Anchor-rule baseline | sensitivity 1.0 on canonical peptides |
| Non-canonical plateau | ~0.65 across 3 sequence predictors |
| Cancer neoepitopes | p53 R248W (`MNWRPILTI`), KRAS G12V (`YKLVVVGAV`) |
| Salt bridge | K-E63, qualitative hypothesis (not quantified) |

---

## Repository Structure

```
├── 01_Literature/           # Literature review + notes
├── 02_Data/                 # Features, predictions, results
│   ├── raw/                 # Original feature matrices
│   └── cleaned/             # Processed benchmark results
├── 03_Analysis/             # ALL code
│   ├── R_scripts/           # R analysis (peptide_mhc_binding_study.R, 1,700+ lines)
│   ├── Python_scripts/      # Python benchmarks
│   ├── figures/             # Publication-quality figures (Nature palette)
│   └── models/              # Trained .keras models
├── 04_Manuscript/           # Manuscript + tables + reports
│   ├── manuscript.md        # Primary draft
│   ├── tables/              # Tables 1-5
│   └── reports/             # 10 analysis reports
├── 05_Submission/           # Cover letter + submission materials
│   └── cover_letter_bioinformatics.md
├── 06_Logs/                 # Decisions + change log
├── docking/                 # Structural analysis (separate workflow)
│   ├── salt_bridge_validate.py    # K-E63 validation
│   ├── md_salt_bridge.py          # MD simulation
│   ├── visualize_salt_bridge.py   # Publication figures
│   ├── 1DUZ.pdb                   # Template structure
│   ├── figures/                   # Docking figures
│   └── hdock_results/             # HDOCKlite outputs
├── CLAUDE.md                # Project constitution
└── RUN_PIPELINE.md          # Full pipeline docs
```

---

## Quick Start

### Prerequisites

```bash
# R 4.6.0+
install.packages(c("keras", "tensorflow", "pROC", "tidyverse", "caret", "ggsci", "gridExtra"))

# Python 3.13+
pip install openmm pdbfixer py3Dmol matplotlib numpy biopython tensorflow mhcflurry
```

### Reproduce All Results

```bash
# Step 1: Train models and generate benchmark results
cd 03_Analysis/R_scripts && Rscript run_main_pipeline.R

# Step 2: Regenerate all figures (Nature Publishing Group palette)
Rscript regenerate_all_figures_nature.R

# Step 3: Run structural analysis
cd ../../docking && python salt_bridge_validate.py
python visualize_salt_bridge.py

# Step 4: MHCflurry direct comparison
cd ../03_Analysis/Python_scripts && python mhcflurry_direct_compare.py
```

---

## Methodological Highlights

### Benchmark Artifacts (Methods Hygiene)

Three apparent findings dissolved under controlled re-analysis:

```
per-allele > pan         → 10.6× data-size artifact
label quality dominates  → degenerate random-null artifact
Deep FFN > CNN           → stale-checkpoint artifact
```

A true i.i.d. random-label control falls to chance (~33%), and multi-seed re-training shows a CNN modestly outperforms the Deep FFN.

### MHCflurry Label Circularity

When MHCflurry 2.2.0 is evaluated on the held-out test set (which it labelled), it achieves **100% accuracy** — quantifying the 8.1 percentage point gap as the information loss between a known decision boundary and a finite-sample approximation.

### Three-Tier Neoantigen Screening Pipeline

Proposed integration of structural information to capture non-canonical binding configurations:

| Tier | Method | Purpose |
|:----:|--------|---------|
| 1 | Sequence-based binding prediction | High-throughput initial filtering |
| 2 | AlphaFold3 pMHC structural scoring | Structure-conditioned assessment |
| 3 | Template-based energy minimisation | Structural plausibility validation |

### K-E63 Salt Bridge: A Predictor Blind Spot

The KRAS G12V neoepitope (YKLVVVGAV) carries a **charged lysine at P2** — a non-canonical anchor configuration. Template-based structural analysis and energy minimisation suggest a possible K-E63 salt bridge (qualitative, un-solvated electrostatic estimate only), a prediction-resistant configuration (cf. Mares et al., ICML 2025, AUROC 0.06–0.22).

---

## Figures

| Figure | Description |
|--------|-------------|
| **Figure 1** | Model performance comparison (7 models including MHCflurry baseline) |
| **Figure 2** | Architecture comparison + label-source controls — MHCflurry (90.6%), PSSM (82.6%), random (chance ~33%) |
| **Figure 3** | IEDB benchmark — ROC curve, confusion matrix, per-epitope scores |
| **Figure 4** | Protein epitope scanning — 10 proteins, 3,536 peptides |
| **Figure 5** | Cancer hotspot mutation analysis — p53 + KRAS |
| **Figure 6** | BLOSUM62 encoding heatmaps + learned feature representations |
| **Figure 7** | ESM-2 embedding comparison — 5 encoding strategies |
| **Figure 8** | K-E63 salt bridge structural hypothesis |

---

## Citation

```bibtex
@article{zhou2026benchmarking,
  title={The Non-Canonical Frontier in Peptide-MHC-I Binding Prediction: A Systematic Multi-Architecture and Multi-Allele Benchmark},
  author={Huang, Yunshi and Chen, Weisheng and Cai, Yiqi and Zhuang, Xiaopeng and Zhou, Zhuxian and Xue, Xiangyang and Zhou, Zhuha},
  journal={Submitted},
  year={2026}
}
```

---

## License

MIT License — see [LICENSE](LICENSE) for details.

## Contact

Zhuha Zhou — zhouzhuha@wmu.edu.cn — Department of Gastroenterology Surgery, The First Affiliated Hospital of Wenzhou Medical University
