# Bioinformatics (Oxford) — Submission Requirements

> Target journal for the merged peptide-MHC manuscript. Original Paper type.
> Retrieved 2026-08-02 from OUP author guidelines.

## Journal Facts
- **IF**: ~5.4–5.5 (2024 data; 5-year ~7.1)
- **JCR**: Q1 (Biochemical Research Methods 11/85; Math & Comput Biology 7/65)
- **中科院分区**: 2025.3 升级版 生物学 3区 (小类均3区, 非Top); 2026.3 版或升 2区/小类1区 (待确认)
- **预警**: 非预警
- **Model**: 传统订阅刊 (非OA, 通常无版面费); CC BY OA option available
- **年发文**: ~649–736

## Article Types (word limits)
| Type | Limit |
|------|-------|
| **Original Papers** | ≤7 pages ≈ **5,000 words excl. figures** |
| Application Notes | ≤4 pages ≈ 2,600 words or 2,000+1 fig |
| Reviews | 3–8 pages ≈ 3,000–6,000 words |

⚠️ Hard limit: exceeding by ~20% → returned before review.

## Abstract (Original Papers)
**≤150 words**, structured: **Motivation / Results / Availability and Implementation / Contact / Supplementary Information**

## References
**OUP author-year style**: in-text Author (Year); alphabetical list; journal names abbreviated italic; volume bold.
(NOT Vancouver numbering.)

## Figures
| Item | Spec |
|------|------|
| Resolution | 350 DPI color/halftone; 1,200 DPI line |
| Formats | TIFF / EPS / PDF / PNG |
| Color | RGB |
| Column widths | 86 mm (single) / 178 mm (double) |
| Font | Arial or Helvetica, ≥7pt |
| Panel labels | lowercase (a) (b) (c) |

Our figures: Arial + vector PDF (close; verify dpi).

## ML Paper Special Requirements (critical)
- Dedicated subsection describing **training dataset composition** (cross-validation sets, independent test set, homology handling)
- Report **performance on independent test set** (in Methods for Original Papers)
- Leave-one-out → editorially rejected unless justified

Our manuscript has independent IEDB test set → compliant.

## Supplementary Material
- Separate file(s); all supplementary in **single file**
- Must be referenced in main text
- Cannot be altered after acceptance

## Submission Process
- Format-free first submission allowed (word/page limits still apply)
- Submit via OUP online system
- Required: ORCID, Data Availability statement, CRediT
- Single-anonymized review, ~3 reviewers, acceptance ~30%
- Code must be accessible & URL functional at review (GitHub ready)

## LaTeX Template
- OUP general template (Overleaf official)
- `\documentclass[numsec,webpdf,modern,large]`
- `\bibliographystyle{author-date}`
- GitHub improved: RagnarGrootKoerkamp/oxford-bioinformatics-template

## Key Gap
- Current manuscript.md: **11,314 words** → need ~5,000 (cut ~6,300)
- Word distribution: Abstract 377 | Key Points 264 | Intro 724 | Methods 1,474 |
  Results 3,060 | Discussion 2,824 | Data/Code 203 | Declarations 199 |
  Biographies 152 | References 1,228 | Figure Legends 630
