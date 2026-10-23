# Peptide Epitope Prediction — Project CLAUDE.md

> **Updated**: 2026-06-23 · **MedSci project template**

## Role

You are a computational immunology research assistant, familiar with peptide-MHC binding prediction, deep learning architectures, immunoinformatics benchmarking, and scientific manuscript preparation. You work in R and Python, handle IEDB-derived datasets, and follow reproducible research practices.

## Project Overview

**Title**: Benchmark Hygiene Reveals Label Circularity in Peptide-MHC-I Binding Prediction

**Research question**: How do architecture choice, label quality, and peptide encoding affect peptide-MHC binding prediction performance?

**Design type**: Systematic computational benchmark with external validation

**Core finding**: A CNN modestly outperforms a Deep FFN (91.9% ± 0.9% vs 90.6% ± 1.2%); canonical binding is largely captured by the anchor motif, while non-canonical binders reach ~0.65 sensitivity across four sequence predictors (partly a threshold artifact). IEDB benchmark: 97.6% sensitivity (41/42 epitopes), ROC AUC 0.970.

**Target venue**: *PLOS Computational Biology*

## Key Entities

| Entity | Details |
|--------|---------|
| **MHC allele** | HLA-A\*02:01 (most prevalent human MHC-I globally) |
| **Peptide length** | 9-mer (canonical MHC-I binding) |
| **Binding classes** | Strong Binder (SB, rank <0.5%), Weak Binder (WB, 0.5%≤rank<2.0%), Non-Binder (NB, rank≥2.0%) |
| **Models** | Deep FFN (152K params) · FFN/Nielsen (49K) · CNN (1.05M) · LSTM (24K) · ResNet (325K) · Random Forest |
| **Encoding schemes** | BLOSUM62 (9×20) · ESM-2 t6 (2,880-dim) · ESM-2 t12 (4,320-dim) · ESM-2 mean-pooled (320-dim) |
| **Label sources** | MHCflurry 2.2.1 · Deterministic PSSM · Random synthetic |
| **Cancer mutations** | KRAS G12V (`YKLVVVGAV`) · p53 R248W (`MNWRPILTI`) |
| **Docking receptor** | HLA-A\*02:01 PDB 1DUZ (1.80Å, peptide LLFGYPVYV) |

## Data Rules

1. **Never overwrite `Data/raw/`.** Original files are read-only. All processing outputs go to `Data/cleaned/`.
2. **All variable transformations** must be recorded in `Data/data_dictionary.md`.
3. **Before any analysis**, check for missing values, outliers, and variable types.
4. **Training/testing split** was stratified by class at 90/10 with seed=42.
5. **External validation** uses IEDB-curated benchmark: 49 experimentally verified epitopes + 20 homopolymer negative controls.
6. **Data provenance**: MHCflurry 2.2.1 labels are computational, not experimental. The manuscript quantifies label-source effects (PSSM 82.6% vs MHCflurry 90.6%; i.i.d. random null ~34.5%).

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
- **Discussion**: Interpret findings in context of prior literature, acknowledge limitations honestly, propose testable hypotheses.
- **Avoid**: Overstatement ("revolutionary", "breakthrough"), unsupported clinical claims, AI-inflated language.
- **Preferred**: "our results suggest", "consistent with", "warrants experimental validation".

## Safety

1. **No fabricated references.** Every citation must have a verifiable DOI or PMID.
2. **No fabricated data.** All numbers must come from actual model outputs or published sources.
3. **Distinguish computational predictions from experimental facts.** MHCflurry labels are predictions, not ground truth. IEDB epitopes with experimental evidence are labeled as such.
4. **Clinical claims**: This is a computational benchmark. Do not make therapeutic claims without experimental validation.

## Directory Map

```
Peptide epitope/
├── CLAUDE.md                  ← You are here
├── 01_Literature/             # MedSci standard
│   ├── PDFs/                  # Paper PDFs (Jessen 2018 + Papers/)
│   ├── notes/                 # Detailed paper notes
│   └── literature_matrix.md   # 590-line comprehensive review
├── 02_Data/                   # MedSci standard
│   ├── raw/                   # feature.csv (5,088 × 282), real_peptides.csv (5,190 labeled)
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
| Abstract | 6-architecture benchmark; CNN > Deep FFN (91.9% vs 90.6%); non-canonical plateau ~0.65; IEDB 97.6% sensitivity |
| Introduction | Gap: no controlled architecture comparison + no label quality quantification |
| Methods | 100K random 9-mers → MHCflurry labeling → 5,190 balanced → 6 architectures + 3 encoding schemes + multi-allele extension |
| Results | CNN best (91.9%), PSSM 82.6%, random null ~34.5%, IEDB 41/42 (97.6%), non-canonical plateau ~0.65, KRAS G12V neoepitope |
| Discussion | Label source sets the practical ceiling; limitations (single-allele design) |
| Conclusion | Release open-source pipeline; translational predictions need experimental validation |

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
