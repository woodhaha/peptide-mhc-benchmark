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
3. 非规范（Type B）结合剂在三个独立序列预测器上停在 ~0.65 平台（MUNIS 不救；MUNIS 是序列模型 ESM-2+LSTM）
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

## 2026-08-02 — 作者名单定稿（5 人）
- 按用户要求：先用 CRC (ImmunoTargets) 名单替换原 5 人名单（加入 Chen Weisheng、Yiqi Cai，移除 Shaoliang Han，6 人），随后删除 Qigang Xu，定稿为 **5 人**。
- 定稿名单：Zhuha Zhou†*（共同一作+通讯）、Chen Weisheng†（共同一作）、Yiqi Cai、Yongyu Bai、Zhuxian Zhou^（浙大通讯）。
- 单位重编号：删除肝胆胰外科 affiliation 后，浙大 Zhuxian Zhou 由 ³ → ²。
- 更新文件：manuscript.md、manuscript_bioinformatics.md、cover_letter.md、cover_letter_bioinformatics.md、README.md、journal_requirements.md、两个 DOCX（重新生成，验证 Qigang 已移除、Chen/Yiqi 保留）。
- 验证：5 位作者在作者块/贡献/简介中一致出现；无 Qigang 残留（历史 decisions.md 除外）。
- 132 commits，git clean。

## 2026-08-02 — 投稿风险评估（结论）
- 投稿版 manuscript_bioinformatics.md 的结论**已是诚实表述**：包含"实验验证盐桥/免疫原性尚需进行"、PSSM 基线作为对照、无"系统性漏掉/三阶梯解决方案"过度主张。
- 与工作版 manuscript.md 的唯一实质差异是：投稿版**缺少多等位基因扩展新发现**（pan16、时间 holdout、MUNIS 平台）—— 这些是**修订时新增内容**，不是对既有主张的修正。
- **风险修正**：先前判断"投稿版结论与修正版不一致"为高风险，经核查实为中低风险 —— 投稿版已自洽，无紧迫修复必要。
- 遗留中风险：单一等位基因泛化性（pan16 现成可回应）、标签循环（时间 holdout 现成可回应）。

## 2026-08-02 — 投稿版 cover letter URL 修复（私有→公开仓库）
- 发现：投稿版 `cover_letter_bioinformatics.md:23` 引用的代码仓库 URL `woodhaha/peptide-epitope-prediction` 为 **PRIVATE**，审稿人无法访问。
- 实际公开仓库为 `woodhaha/peptide-mhc-benchmark`（PUBLIC，483 文件，含全部数据 CSVs/代码/模型，2026-07-03 推送，与 6/30 投稿稿配套）。
- 修复：`.md` 源文件 URL 改为公开仓库；`cover_letter.docx` 用 pandoc 重新生成并验证 URL 已替换（document.xml 无私有仓库残留）。
- 附带确认：公开仓库内容与投稿稿一致，无审稿人可见缺口。本地 114 个未推送提交为投稿后多等位基因扩展（姊妹篇内容），不属于投稿稿范围。

## 2026-08-02 — 主稿应用修正清单 #3 和 #5
- manuscript_correction_checklist.md 6 项中 #1/#2/#4/#6 已在前期应用（hard tail 措辞、三阶梯改假说、多等位基因支持、K-E63 保留）。
- 本轮补齐剩余 2 项到工作版 04_Manuscript/manuscript.md：
  - #3 锚点检查器基线（§3.5）：零参数 anchor-checker 规则（p2 L/M/I/V/T + p9 L/M/I/V）在规范占优的 IEDB 基准上 sens 1.0 / F1 0.85-1.0，与训练模型无统计差异 → 规范集高估泛化，有效信号在非规范/hard-tail。
  - #5 预测器平台限制定语（§4.6 新增 Sixth 点）：非规范结合子在 Deep FFN/MHCflurry/MUNIS 三类预测器均 ~0.65 平台；能量补偿尾（P2-K/R/H 达 Glu63 3.83Å）对一切方法预测抵抗。
- checklist 两处标记 DONE。投稿版 manuscript_bioinformatics.md 未动（投稿在前，修订时统一应用）。

## 2026-08-02 — 姊妹篇草稿修订（审阅后）
- 审阅发现并修复：
  - 引用 [10] Robinson/IMGT-HLA 经 CrossRef 验证真实（NAR 2015;43:D423-D431），非虚构。
  - 字数/表/图元数据修正为实际值（~2,000 词, 1 表, 2 图）。
  - [companion submission] 自引用改为正式引用 [15]（补参考文献条目）。
  - "pre-registered thresholds" 改为诚实的 "fixed thresholds"（阈值实为事后确定的 5th-percentile canonical）。
  - 补核心正结果图：fig2_pan_vs_perallele.png（重平衡 pan 5 种子均±SD）。未叠加 per-allele 条形——per-allele 基准（n=25 broad-motif）与大型非规范集不同，叠加会构成 apples-to-oranges 比较（三堂会审已点名的错误）。
  - 图 1（findings panel）补正文引用（1A temporal / 1B pan16 / 1C P2 reach / 1D plateau）+ 图注。
- 现状：正文 Figure 1+2 均被引用且有图注；2,239 词；无残留占位。

## 2026-08-02 — 投稿状态纠正（重要更正）
- **事实核查：稿件从未提交。** journal_requirements.md 显示 "Submission date: [Pending]"；无任何投稿回执或系统确认。
- 前述多条记录（"已投 BiB"、"投稿版"、"6/30 投稿"、"与 6/30 投稿稿配套"、"投稿在前"）均为**错误假设**，把"准备投稿"误当"已投稿"。此前的"投稿风险评估"基于错误前提。
- 实际状态：合并稿 manuscript.md 已修订就绪（含多等位基因扩展），尚未投任何期刊。manuscript_bioinformatics.md 是 Bioinformatics 格式的渲染版（Harvard 引用），供投稿使用。
- 正确推荐：可投 Briefings in Bioinformatics（IF~9.5，T3）或 Bioinformatics（Oxford，IF~6，格式已就绪）。
