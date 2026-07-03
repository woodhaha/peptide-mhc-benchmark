# Peptide Epitope Prediction — Complete Pipeline

> Last run: 2026-06-24 · Total runtime: ~3 min (excluding MD)

## Quick Start

```bash
cd "D:/Researching/Peptide epitope"

# Step 1: Retrain models (R)
cd 03_Analysis/R_scripts && Rscript run_main_pipeline.R && cd ../..

# Step 2: Run structural analysis
cd docking && python salt_bridge_validate.py && cd ..

# Step 3: Generate figures
cd docking && python visualize_salt_bridge.py && cd ..

# Step 4: Reformat manuscript
cd 04_Manuscript && python reformat_bioinformatics.py && cd ..
```

## Pipeline Architecture

```
Peptide epitope/
├── RUN_PIPELINE.md              ← You are here
│
├── 01_Literature/               # Literature review
│   ├── PDFs/                    # Source papers
│   ├── notes/                   # Detailed paper notes
│   └── literature_matrix.md     # 590-line comprehensive review
│
├── 02_Data/                     # Data layer (READ-ONLY raw/)
│   ├── raw/
│   │   ├── feature.csv          # 5,088 × 282 BLOSUM62 features
│   │   └── real_peptides.csv    # 5,190 labeled peptides (train/test split)
│   ├── cleaned/                 # 11 processed CSVs + model outputs
│   └── data_dictionary.md       # Variable definitions
│
├── 03_Analysis/                 # Analysis layer (ALL code lives here)
│   ├── R_scripts/
│   │   ├── run_main_pipeline.R          # Master R pipeline (1,700+ lines)
│   │   ├── run_optimized_pipeline.R     # Faster alternative
│   │   ├── peptide_mhc_binding_study.R  # Core ML training
│   │   ├── benchmark_iedb_epitopes.R    # IEDB validation
│   │   ├── scan_protein_epitopes.R      # Protein scanning
│   │   ├── mutation_scan.R              # Cancer mutation analysis
│   │   ├── plot_comparison.R            # Figure generation
│   │   ├── regenerate_all_figures_nature.R   # Nature-style figures
│   │   ├── label_peptides_mhcflurry.py  # MHCflurry label generation
│   │   ├── run_full_benchmark.py        # Full benchmark orchestrator
│   │   └── train_deep_ffn.py            # Python model training
│   ├── figures/                # All figures (PNG + PDF)
│   ├── models/                 # Trained .h5/.keras/.rds files
│   ├── outputs/                # Intermediate analysis outputs
│   └── model_manifest.csv      # Model registry
│
├── 04_Manuscript/              # Manuscript layer
│   ├── manuscript.md           # Original draft
│   ├── manuscript_bioinformatics.md  # Bioinformatics-formatted
│   ├── manuscript.docx/.pdf    # Rendered versions
│   ├── tables/                 # Tables 1-5 (.md)
│   ├── reports/                # 10 analysis reports (MD+PDF)
│   └── reformat_bioinformatics.py    # Format converter
│
├── 05_Submission/              # Submission layer
│   ├── cover_letter.md         # Generic cover letter
│   ├── cover_letter_bioinformatics.md  # Bioinformatics cover letter
│   ├── reviewer_response.md    # R&R response template
│   ├── journal_requirements.md # Target journal specs
│   ├── supplementary/          # Supplementary figures + data
│   ├── figures/                # Submission-ready figure PDFs
│   └── tables/                 # Submission tables
│
├── docking/                    # Structural analysis (separate workflow)
│   ├── salt_bridge_validate.py      # MAIN: K-E63 validation (8.1s)
│   ├── docking_pose_analysis.py     # HDOCK + B-pocket comparison
│   ├── docking_energy.py            # OpenMM binding energy
│   ├── md_salt_bridge.py            # MD simulation (experimental)
│   ├── visualize_salt_bridge.py     # Publication figure generator
│   ├── 1DUZ.pdb                     # Template structure
│   ├── receptor_1DUZ.pdb            # Prepared receptor
│   ├── peptide_*.pdb                # 10 peptide structures
│   ├── hdock_models/                # HDOCK output (PLACEHOLDER)
│   ├── references/                  # 7 reference HLA-A*02:01 structures
│   ├── figures/                     # Generated figures
│   │   ├── Figure_K63_salt_bridge.png   # 4-panel publication figure
│   │   ├── Figure_K63_salt_bridge.pdf   # Vector version
│   │   └── K63_salt_bridge_3D.html      # Interactive 3D view
│   ├── salt_bridge_validation.json   # Numerical results
│   ├── pose_analysis.json           # B-pocket comparison
│   ├── energy_results.json          # Anchor scoring
│   └── RESULTS.md                   # Results summary
│
├── 06_Logs/                    # Decision + change tracking
│   ├── decisions.md
│   └── change_log.md
│
├── CLAUDE.md                   # Project constitution
├── .gitignore                  # Git exclusions
└── .gitattributes              # Git LFS rules (.dcd files)
```

## Dependency Map

```
02_Data/raw/feature.csv
    │
    ▼
03_Analysis/R_scripts/run_main_pipeline.R
    │
    ├──→ 03_Analysis/models/*.h5, *.keras, *.rds
    ├──→ 03_Analysis/figures/Figure1-8.png
    ├──→ 03_Analysis/outputs/*.md, *.pdf
    │
    ▼
04_Manuscript/manuscript.md  ←──  manual writing
    │
    ▼
docking/salt_bridge_validate.py
    │
    ├──→ salt_bridge_validation.json
    ├──→ pose_analysis.json
    │
    ▼
docking/visualize_salt_bridge.py
    │
    ├──→ figures/Figure_K63_salt_bridge.png
    ├──→ figures/K63_salt_bridge_3D.html
    │
    ▼
04_Manuscript/reformat_bioinformatics.py
    │
    └──→ manuscript_bioinformatics.md

05_Submission/cover_letter_bioinformatics.md  ←── manual writing
```

## Key Results Summary

### ML Benchmark
| Metric | Value |
|--------|:-----:|
| Best architecture | Deep FFN (152K params) |
| Accuracy (MHCflurry labels) | 91.9% |
| Macro F1 | 0.921 |
| Accuracy (PSSM labels) | 94.8% |
| Label quality gap | 29 pp (PSSM vs random) |
| IEDB sensitivity | 93.9% (46/49) |
| IEDB ROC AUC | 0.947 |
| ESM-2 t6 accuracy | 93.3% |

### K-E63 Salt Bridge
| Metric | Value |
|--------|:-----:|
| P2 CA → E63 OE2 | 3.7 Å |
| K sidechain reach | 7.6 Å |
| K reach margin | +3.9 Å (overshoots) |
| Coulomb ΔE (salt bridge) | −29.6 kcal/mol |
| Canonical Leu burial ΔG | −3.0 kcal/mol |
| Salt bridge advantage | ~10× |
| B-pocket conservation | 11.4 ± 0.1 Å (8 structures) |
| Verdict | ✅ Structurally & energetically favorable |

## Environment

```bash
# Python
python 3.13+
pip install openmm pdbfixer py3Dmol matplotlib numpy biopython

# R
R 4.6.0+
install.packages(c("keras", "tensorflow", "pROC", "tidyverse", "caret"))

# System
git-lfs 3.7.1+
gh CLI (for GitHub push)
```

## Git History

```
940cd42 feat: Bioinformatics (Oxford) formatted manuscript
140c0c0 docs: add K-E63 salt bridge results summary
d5f32d2 fix: md_salt_bridge with pdbfixer + gradual restraint annealing
6eadd45 feat: add K-E63 salt bridge publication figures
4468954 fix: track DCD trajectory with Git LFS
f3a842f feat: docking code optimization — 4→3 consolidated scripts + K-E63 validated
ddedced feat: initialize Peptide Epitope Prediction project
```

## Reproducibility

All analysis is deterministic with seed=42. The complete pipeline from raw data to manuscript figures is executable via the scripts listed above. Trained models are serialized in HDF5/Keras/RDS format. Feature matrices are stored as CSV. The only non-deterministic component is the MD simulation (md_salt_bridge.py), which requires ≥100 ns equilibrated production run; the static analysis (salt_bridge_validate.py, 8.1s) provides all structural metrics reported in the manuscript.

## Next Steps for Submission

- [ ] Convert figures to TIFF (≥300 dpi)
- [ ] Package supplementary data
- [ ] Collect all author ORCIDs
- [ ] Final proofread of manuscript_bioinformatics.md
- [ ] Submit via ScholarOne: https://mc.manuscriptcentral.com/bioinformatics
