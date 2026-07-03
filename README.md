# Peptide-MHC Binding Prediction Benchmark

> **Label quality, not architecture, defines the performance ceiling for peptide-MHC binding prediction.**

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
![R 4.6.0+](https://img.shields.io/badge/R-4.6.0%2B-blue)
![Python 3.13+](https://img.shields.io/badge/Python-3.13%2B-blue)

Systematic benchmarking of six deep learning architectures for HLA-A\*02:01 peptide-MHC binding prediction, with external validation against 49 experimentally verified IEDB epitopes, cancer hotspot mutation analysis (p53 R248W, KRAS G12V), and molecular docking-based structural hypotheses for neoepitope-MHC interactions.

**Manuscript**: `04_Manuscript/manuscript.md` (54 references, ~9,000 words, 8 figures, 6 tables)

---

## Core Finding

| Comparison | Accuracy | Macro F1 | WB F1 |
|-----------|:--------:|:--------:|:-----:|
| **PSSM labels** (deterministic) | **94.8%** | **0.948** | **0.925** |
| Deep FFN (MHCflurry labels) | 91.9% | 0.921 | 0.880 |
| Random Forest | 81.1% | 0.814 | 0.723 |
| Random labels | 65.8% | 0.558 | 0.000 |

Label source accounts for **~29 percentage points** of accuracy variation — exceeding the maximum inter-architecture difference of ~11 percentage points. The **weak binder (WB) class** is the primary locus of both label-induced and architecture-induced variation.

---

## Key Results

| Metric | Value |
|--------|:-----:|
| Best architecture | Deep FFN (152K parameters) |
| IEDB sensitivity | **93.9%** (46/49) |
| IEDB ROC AUC | **0.947** |
| ESM-2 t6 per-position accuracy | 93.3% |
| 10-protein scanning | **10/11 known epitopes** confirmed |
| Cancer neoepitopes | p53 R248W (`MNWRPILTI`), KRAS G12V (`YKLVVVGAV`) |
| Predicted salt bridge | K-E63, ~−30 kcal/mol (vs −3 kcal/mol canonical) |

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
├── 05_Submission/           # Cover letter + graphical abstract
│   ├── cover_letter_bioinformatics.md
│   └── graphical_abstract_concept.md
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

### Label Contamination Cascade

Three labelling strategies demonstrate how ML-derived training labels systematically compound performance ambiguity:

```
PSSM (deterministic)   → 94.8% accuracy  "Clean"
MHCflurry 2.2.0 (ML)  → 91.9% accuracy  "Complex boundary"
Random synthetic       → 65.8% accuracy  "Noise baseline"
```

**WB F1 drops from 0.925 → 0.880 → 0.000**, concentrating all performance variation in the ambiguous 50–500 nM affinity range.

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

The KRAS G12V neoepitope (YKLVVVGAV) carries a **charged lysine at P2** — a non-canonical anchor configuration. Template-based structural analysis and energy minimisation predict a K-E63 salt bridge with −30 kcal/mol electrostatic stabilisation, representing a **"prediction-resistant"** configuration that sequence-based models systematically miss (Mares et al., ICML 2025, AUROC 0.06–0.22).

---

## Figures

| Figure | Description |
|--------|-------------|
| **Figure 1** | Model performance comparison (7 models including MHCflurry baseline) |
| **Figure 2** | Label quality effect — PSSM vs MHCflurry vs Random |
| **Figure 3** | IEDB benchmark — ROC curve, confusion matrix, per-epitope scores |
| **Figure 4** | Protein epitope scanning — 10 proteins, 3,536 peptides |
| **Figure 5** | Cancer hotspot mutation analysis — p53 + KRAS |
| **Figure 6** | BLOSUM62 encoding heatmaps + learned feature representations |
| **Figure 7** | ESM-2 embedding comparison — 5 encoding strategies |
| **Figure 8** | K-E63 salt bridge structural hypothesis |

---

## Citation

```bibtex
@article{zhou2025benchmarking,
  title={Systematic Benchmarking of Deep Learning Architectures for Peptide-MHC Binding Prediction Reveals Label Quality, Not Architecture, Defines Performance},
  author={Zhou, Zhuha and Bai, Yongyu and Xu, Qigang and Zhou, Zhuxian and Han, Shaoliang},
  journal={Submitted to Briefings in Bioinformatics},
  year={2026}
}
```

---

## License

MIT License — see [LICENSE](LICENSE) for details.

## Contact

Zhuha Zhou — zhouzhuha@wmu.edu.cn — Department of Gastroenterology Surgery, The First Affiliated Hospital of Wenzhou Medical University
