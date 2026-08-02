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

## 2026-08-02 — Conclusions update with multi-allele findings
- Updated 4.7 Conclusions: added pan16 (15-allele panel, 0.819 pooled non-canonical
  sensitivity), temporal holdout (post-2022 preserved, breaks label circularity),
  predictor-plateau (non-canonical ~0.65 across Deep FFN/MHCflurry/MUNIS), P2-K/R/H
  salt-bridge generality (3.83 Å), three-tier pipeline re-framed as hypothesis.
- Updated Key Points + Discussion three-tier (248): softened "systematically miss",
  added MUNIS-plateau caveat.
- Source: multi-allele extension (10 validated findings, FINDINGS_SYNTHESIS.md).

## 2026-08-02 — 研究进度记录（完整会话）
### 状态：研究已完成并收敛，全部提交，git clean（127 commits）

### 核心发现（10 个验证过的发现 F1-F10，见 FINDINGS_SYNTHESIS.md）
1. 单个重平衡 pan Deep FFN 泛化到 15 等位基因面板（非规范灵敏度 pooled 0.819，12/15 ≥0.76）
2. 规范结合剂被所有预测器类以 ~0.98 识别（已解决）
3. 非规范（Type B）结合剂在序列 + 结构预测器上停在 ~0.65 平台（MUNIS 不救）
4. Type A（容忍锚）是真实覆盖率缺口（0.48，质谱确认、非训练缺口）
5. 极端能量补偿尾巴（P2-K/R/H 盐桥）对所有预测器抵抗
6. 时间 holdout 验证：post-2022 表位灵敏度保持/提升（打破标签循环）
7. pan16 面板成立（15 等位基因，B*27:05/B*51:01 弱尾因数据稀缺）

### 关键决策/证伪
- "per-allele > pan" → 数据量假象（pan-full 恢复）
- "规范 100% 灵敏度" → motif 循环论证（anchor-checker = 1.0 F1）
- "三阶梯管线救非规范" → 假说（MUNIS 平台证据不支持）

### 主稿更新
- 4.7 Conclusions + Key Points + Discussion 三阶梯 → 融入多等位基因发现，三阶梯降级为假说

### 产物清单
- FINDINGS_SYNTHESIS.md（10 发现）、manuscript_correction_checklist.md（6 项修正）
- research-log.md（决策时间线，10 阶段）
- to_human/findings_panel.png（4 联图）+ multi_allele_results.png + type_ab_largesample.png
- 05_Submission/supplementary/findings_summary_tables.md（5 表）
- 04_Manuscript/manuscript.md（结论已更新）

### 下一步（仅写作，非研究）
- 将 manuscript_correction_checklist.md 应用到主稿其他章节
- 或写多等位基因扩展姊妹篇（可调 ml-paper-writing）
