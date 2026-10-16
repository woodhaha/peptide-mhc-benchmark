# Submission Package

Two sister papers on the non-canonical frontier of peptide-MHC-I binding prediction.

## Final deliverables (top level)

| File | Target journal | Description |
|------|---------------|-------------|
| `manuscript.pdf` | Bioinformatics | 主稿最终 PDF（7 页） |
| `companion_manuscript.pdf` | PLOS Comput Biol | 姐妹篇最终 PDF（11 页） |
| `supplementary.pdf` | (主稿附件) | 补充材料 PDF（9 页） |

## Source manuscripts

| File | Format |
|------|--------|
| `manuscript_bioinformatics_v2.md` | 主稿 Markdown（权威内容源） |
| `companion_manuscript.md` | 姐妹篇 Markdown（权威内容源） |
| `cover_letter.md` | 主稿投稿信（Bioinformatics） |
| `companion_cover_letter.md` | 姐妹篇投稿信（PLOS Comput Biol） |

## LaTeX source

| Directory | Template |
|-----------|----------|
| `latex_template/` | 主稿 LaTeX（官方 OUP `bioinfo.cls`）→ 生成 `manuscript.pdf` |
| `companion_latex/` | 姐妹篇 LaTeX（官方 PLOS `plos2015.bst`）→ 生成 `companion_manuscript.pdf` |

## Figures

| Location | Content |
|----------|---------|
| `../05_Submission/figures/` | 主稿 Fig 1–9 |
| `companion_figures/` | 姐妹篇 Fig 1–4 |
| `../05_Submission/supplementary/figures/` | 补充图 S1–S5 + S-MD1–4（9 张） |

## Data packages

| File | Content |
|------|---------|
| `manuscript_data_package.tar.gz` | 主稿数据包（12.5 MB） |
| `companion_data_package.tar.gz` | 姐妹篇数据包（10.6 MB） |

## Supplementary text

`../05_Submission/supplementary/` 含：
- `supplementary.tex` + `supplementary.pdf`（已编译）
- `MD_salt_bridge_supplementary.md`（MD 盐桥分析）
- `findings_summary_tables.md`（发现汇总表 S1–S5）
- `iedb_benchmark_results_v2.csv`（IEDB 基准数据）

## 投稿清单（提交顺序）

1. **Bioinformatics**：`manuscript.pdf` + `cover_letter.md` + 主稿 Fig 1–9
2. **PLOS Comput Biol**：`companion_manuscript.pdf` + `companion_cover_letter.md` + 姐妹篇 Fig 1–4（图需转 TIFF/EPS）
3. 附件：`supplementary.pdf`
