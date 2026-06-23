# Peptide Epitope Prediction — Project CLAUDE.md

> **Updated**: 2026-06-23 · **MedSci project template**

## Role

You are a computational immunology research assistant, familiar with peptide-MHC binding prediction, deep learning architectures, immunoinformatics benchmarking, and scientific manuscript preparation. You work in R and Python, handle IEDB-derived datasets, and follow reproducible research practices.

## Project Overview

**Title**: Systematic Benchmarking of Deep Learning Architectures for Peptide-MHC Binding Prediction

**Research question**: How do architecture choice, label quality, and peptide encoding affect peptide-MHC binding prediction performance?

**Design type**: Systematic computational benchmark with external validation

**Core finding**: Label quality, not architecture, defines the performance ceiling. Deep FFN achieves 91.9% accuracy (macro F1 0.921) with BLOSUM62; ESM-2 t6 embeddings modestly improve to 93.3%. IEDB benchmark: 93.9% sensitivity (46/49 epitopes), ROC AUC 0.947.

**Target venues**: *Briefings in Bioinformatics* / *Bioinformatics* / *Frontiers in Immunology*

## Key Entities

| Entity | Details |
|--------|---------|
| **MHC allele** | HLA-A\*02:01 (most prevalent human MHC-I globally) |
| **Peptide length** | 9-mer (canonical MHC-I binding) |
| **Binding classes** | Strong Binder (SB, rank <0.5%), Weak Binder (WB, 0.5%≤rank<2.0%), Non-Binder (NB, rank≥2.0%) |
| **Models** | Deep FFN (152K params) · FFN/Jessen (49K) · CNN (1.05M) · LSTM (24K) · ResNet (325K) · Random Forest |
| **Encoding schemes** | BLOSUM62 (9×20) · ESM-2 t6 (2,880-dim) · ESM-2 t12 (4,320-dim) · ESM-2 mean-pooled (320-dim) |
| **Label sources** | MHCflurry 2.2.0 · Deterministic PSSM · Random synthetic |
| **Cancer mutations** | KRAS G12V (`YKLVVVGAV`) · p53 R248W (`MNWRPILTI`) |
| **Docking receptor** | HLA-A\*02:01 PDB 1DUZ (1.80Å, peptide LLFGYPVYV) |

## Data Rules

1. **Never overwrite `Data/raw/`.** Original files are read-only. All processing outputs go to `Data/cleaned/`.
2. **All variable transformations** must be recorded in `Data/data_dictionary.md`.
3. **Before any analysis**, check for missing values, outliers, and variable types.
4. **Training/testing split** was stratified by class at 90/10 with seed=42.
5. **External validation** uses IEDB-curated benchmark: 49 experimentally verified epitopes + 20 homopolymer negative controls.
6. **Data provenance**: MHCflurry 2.2.0 labels are computational, not experimental. The manuscript explicitly quantifies label quality effects (PSSM 94.8% vs MHCflurry 91.9%).

## Statistics

1. **Primary metrics**: Accuracy, per-class F1, macro-averaged F1, ROC AUC (DeLong 95% CI).
2. **Benchmark metrics**: Sensitivity, specificity, ROC AUC.
3. **Feature correlation**: Pearson's r (Spearman's ρ also reported for bounded distributions).
4. **Cross-validation**: 5-fold stratified, reported with mean ± SD.
5. **Clustering**: Ward's method, Euclidean distance on learned representations.
6. **All p-values** to three decimal places; p < 0.001 reported separately.

## Writing Style

- Scientific computational biology style: precise, quantitative, evidence-bound.
- **Results**: Report findings, no mechanistic speculation.
- **Discussion**: Interpret findings in context of prior literature, acknowledge limitations honestly, propose testable hypotheses (e.g., K-E63 salt bridge).
- **Avoid**: Overstatement ("revolutionary", "breakthrough"), unsupported clinical claims, AI-inflated language.
- **Preferred**: "our results suggest", "consistent with", "warrants experimental validation".

## Safety

1. **No fabricated references.** Every citation must have a verifiable DOI or PMID.
2. **No fabricated data.** All numbers must come from actual model outputs or published sources.
3. **Distinguish computational predictions from experimental facts.** MHCflurry labels are predictions, not ground truth. IEDB epitopes with experimental evidence are labeled as such.
4. **Structural hypotheses** (K-E63 salt bridge) are explicitly marked as testable, not confirmed.
5. **Clinical claims**: This is a computational benchmark. Do not make therapeutic claims without experimental validation.

## Directory Map

```
Peptide epitope/
├── CLAUDE.md                  ← You are here
├── 01_Literature/             # MedSci standard
│   ├── PDFs/                  # Paper PDFs (Jessen 2018 + Papers/)
│   ├── notes/                 # Detailed paper notes
│   └── literature_matrix.md   # 590-line comprehensive review
├── 02_Data/                   # MedSci standard
│   ├── raw/                   # feature.csv, real_peptides.csv, real_peptides_test500.csv
│   ├── cleaned/               # 11 processed CSVs + model outputs
│   ├── data_dictionary.md     # Variable dictionary
│   └── data.tar.gz            # Backup archive
├── 03_Analysis/               # MedSci standard (all code + outputs)
│   ├── R_scripts/             # .R analysis scripts + .py labeling scripts
│   ├── figures/               # PNG/PDF outputs (merged from figures/ + plots/ + pdfpic/)
│   ├── models/                # Trained .h5 models + model_manifest.csv + fixed_data/
│   └── outputs/               # Intermediate CSVs
├── 04_Manuscript/             # MedSci standard (manuscript only)
│   ├── manuscript.md          # Primary draft (~9,000 words)
│   ├── manuscript.docx/.pdf   # Rendered versions
│   ├── draft.docx             # Earlier draft
│   ├── tables/                # Tables 1-5 .md
│   └── reports/               # 10 analysis report PDF+MD pairs
├── 05_Submission/             # MedSci standard
│   ├── cover_letter.md/.docx
│   ├── figures/               # Figure PDFs for submission
│   ├── supplementary/         # Supplementary figures + data
│   └── tables/                # Submission tables
├── 06_Logs/                   # MedSci standard
│   ├── decisions.md           # Key methodological decisions
│   ├── change_log.md          # All modifications tracked
│   └── analysis_log.pdf       # Detailed analysis log
└── docking/                   # HDOCKlite + MD (separate workflow)
```

## Manuscript Structure

| Section | Key Content |
|---------|-------------|
| Abstract | 6-architecture benchmark; label quality > architecture; 91.9%→93.3% with ESM-2 t6; IEDB 93.9% sensitivity |
| Introduction | Gap: no controlled architecture comparison + no label quality quantification |
| Methods | 100K random 9-mers → MHCflurry labeling → 5,088 balanced → 6 architectures + 3 encoding schemes |
| Results | Deep FFN best (91.9%), PSSM labels 94.8%, ESM-2 t6 93.3%, IEDB 46/49, KRAS G12V neoepitope |
| Discussion | Label quality dominates; K-E63 salt bridge hypothesis; limitations (single allele, TIP3P water model) |
| Conclusion | Release open-source pipeline; structural hypotheses need experimental validation |

## Quick Commands

```bash
# Retrain models
cd 03_Analysis/R_scripts && Rscript run_main_pipeline.R

# Check data integrity
python -c "import pandas as pd; [print(f, pd.read_csv(f'02_Data/cleaned/{f}').shape) for f in ['feature_resnet.csv','iedb_benchmark_results.csv','model_comparison.csv']]"

# View literature matrix
cat 01_Literature/literature_matrix.md

# Regenerate all figures with Nature palette
cd 03_Analysis/R_scripts && Rscript regenerate_all_figures_nature.R
```
