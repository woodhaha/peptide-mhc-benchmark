# Change Log

## 2026-06-23 — Comprehensive Results Analysis (v2.1)
- **LSTM is best model (91.5%)** — not Deep FFN. Uses only 24K params.
- **PSSM comparison REVERSED** from original manuscript:
  - Original: PSSM 94.8% > MHCflurry 91.9% (claimed label quality dominates)
  - Actual: PSSM worse for neural nets, RF only model that improves (+5.9pp)
  - Correct story: **label-model interaction**, not one-way label quality
- **CV: 91.9% ± 0.7%** (+2.3pp over original, tighter variance)
- **IEDB: 93.9% sensitivity, 75.0% specificity** — fully reproducible
  - 3 FNs: all non-canonical p9 anchors (K, F, C) — biologically justified
  - 5 FPs: all homopolymers with canonical anchors — binding ≠ immunogenicity
- **KRAS G12V neoepitope: YKLVVVGAV Δ=+0.499 (CREATED)** — reproducible
- **p53 R249S: NRSPILTII Δ=+0.303 (CREATED)** — NEW FINDING, was DESTROYED in original
- **Docking: K-E63 salt bridge structurally plausible** (P2:CA→E63 = 4.7Å, K reach=7.6Å)
- **Manuscript implications**:
  1. Replace Deep FFN with LSTM as best model
  2. Rewrite PSSM section: label-model interaction, not label quality
  3. Add p53 R249S as 3rd neoepitope
  4. Update CV numbers to 91.9%

---

## 2026-06-23 — Full Pipeline Re-run (v2.0)
- **6 models trained** on MHCflurry 2.2.0 labels (5,190 balanced peptides):
  - LSTM: 91.5% (best) | Deep FFN: 90.8% | CNN: 90.2% | ResNet: 89.8% | FFN: 89.2% | RF: 79.4%
- **PSSM comparison**: label quality gap quantified (MHCflurry vs PSSM → ~3pp)
- **5-Fold CV**: 91.9% ± 0.6%
- **IEDB benchmark**: Sensitivity 93.9%, Specificity 75.0% (46/49 TP)
- **Protein scan**: 10 proteins, 3,357 9-mers, 100 strong binders
- **Mutation scan**: 144 windows, KRAS G12V Δ=+0.499 (CREATED), p53 R248W Δ=+0.429 (CREATED)
- **32 Nature NPG figures**: 16 PNG (300dpi) + 16 PDF
- **Code optimized**: Python vectorized (15×), R encoding (8-20×), Keras 3 migration
- **Environment**: TF 2.21, Python 3.13.9, Keras 3
- **Skipped**: Docking, ESM-2 embeddings

---

## 2026-06-23 — MedSci Standard Directory Restructure
- Renamed 6 folders to MedSci numbered standard:
  - `Literature/` → `01_Literature/`
  - `Data/` → `02_Data/`
  - `Analysis/` → `03_Analysis/`
  - `Manuscript/` → `04_Manuscript/`
  - `Submission/` → `05_Submission/`
  - `Logs/` → `06_Logs/`
- Consolidated scattered directories into `03_Analysis/`:
  - `figures/` + `plots/` + `pdfpic/` → `03_Analysis/figures/`
  - `models/` → `03_Analysis/models/`
- Deduplicated 69 files between `04_Manuscript/` and `03_Analysis/`
  - Removed from `04_Manuscript/`: data/, figures/, models/, plots/, scripts/
  - Kept: tables/, reports/, manuscript.*, draft.docx
- Merged `Papers/` → `01_Literature/PDFs/`
- Updated `CLAUDE.md` directory map and quick commands
- `docking/` remains at top level (separate workflow)

## 2026-06-23 — Code Optimization
- **Python** (`label_peptides_mhcflurry_optimized.py`):
  - Vectorized peptide generation (numpy → 15× faster)
  - Pre-allocated float32 arrays for scores/labels
  - Single-pass balancing with numpy sort+slice
  - Fixed-width unicode dtype for labels
- **R** (`encoding_optimized.R` — drop-in replacements):
  - BLOSUM62 encoding: pre-computed index map → 8-12× faster
  - One-hot encoding: matrix indexing replaces for loops → 5× faster
  - PSSM scoring: matrix multiplication → 20× faster
  - Peptide generation: pre-allocated char matrix → 4× faster
  - AAindex encoding: pre-computed property matrix → 6× faster
  - Model metrics: pre-allocated arrays + vapply

## 2026-06-23 — Figure Renumbering & Pre-Submission Fixes
- Renamed 6 figure files to match final manuscript numbering
- Updated all 13 figure paths in `04_Manuscript/manuscript.md`
- Cleaned stale old-version files from `05_Submission/figures/`
- Added Jessen (2018) as Entry #0 in literature matrix
- Created `01_Literature/notes/jessen-2018-dl-cancer-immunotherapy.md`

---

- Renamed 6 figure files to match final manuscript numbering (Figures 1-8)
  - `Figure2_benchmark_roc` → `Figure3_iedb_benchmark` (Figure 3: IEDB benchmark)
  - `Figure3_blosum62_heatmap` → `Figure6_blosum62_heatmap` (Figure 6: Learned features)
  - `Figure4_mutation_delta` → `Figure5_mutation_delta` (Figure 5: Cancer mutations)
  - `Figure5_mutation_map` → `Figure8_mutation_map` (Figure 8: Docking)
  - `comparison_gap` → `Figure2_label_comparison` (Figure 2: Label comparison)
  - `benchmark_per_epitope` → `Figure4_protein_scan` (Figure 4: Protein scan)
- Updated all 13 figure path references in `Manuscript/manuscript.md`
- Cleaned up stale old-version files from `Submission/figures/`
- Added Jessen (2018) as Entry #0 in `Literature/literature_matrix.md`
- Created `Literature/notes/jessen-2018-dl-cancer-immunotherapy.md` (127-line detailed notes)
- Rewrote `CLAUDE.md` (31→129 lines, MedSci-grade project template)
- Orphaned file noted: `Figure6_pan_cancer_oncoplot.png` — not referenced in manuscript
- ⚠️ Note: `Submission/figures/` PDFs are PNG copies; proper PDF generation needed before final submission

---

- Reorganized project into 6-tier directory structure (01-06)
- Moved literature to `01_Literature/`
- Categorized data into raw/cleaned in `02_Data/`
- Consolidated scripts, outputs, and figures in `03_Analysis/`
- Organized manuscript files in `04_Manuscript/`
- Centralized submission materials in `05_Submission/`
- Archived logs and created decision tracking in `06_Logs/`
- Created `CLAUDE.md`, `data_dictionary.md`, and documentation placeholders

---

## Prior Work
See `06_Logs/decisions.md` for historical analysis decisions and iteration log.
