# Cover Letter — Submission to *PLOS Computational Biology*

**Date**: October 21, 2026

**To**: The Editors
*PLOS Computational Biology*

**Re**: Submission of research article — "Benchmark Hygiene and the Non-Canonical Frontier in Peptide-MHC-I Binding Prediction"

---

Dear Editors,

We are pleased to submit our manuscript for consideration in *PLOS Computational Biology*. This work addresses a methodological problem in computational immunology — that peptide-MHC-I binding-prediction benchmarks can manufacture false conclusions when confounders are uncontrolled — and extends it from a controlled single-allele study to a multi-allele benchmark of the non-canonical prediction frontier.

**Benchmark hygiene, applied to immunoinformatics.** Whereas the general data-leakage problem is well documented (Kapoor & Narayanan, *Patterns* 2023; Bernett et al., *Nat Methods* 2024), our distinct contribution is a **label-circularity** failure mode specific to peptide-MHC benchmarks — computational labels mistaken for experimental truth. Under a controlled re-analysis of our own six-architecture benchmark (5,190 MHCflurry-labelled HLA-A\*02:01 9-mers, multi-seed estimates with confidence intervals), three apparent findings dissolved: (i) a per-allele superiority over a pan model was a 10.6× training-data-size artifact; (ii) a claimed "label quality dominates architecture" result collapsed once the degenerate random-label null (which performs at chance, ~33%) was removed; and (iii) a "simple architecture beats complex ones" ranking was an artifact of stale checkpoints — re-training showed a convolutional network instead modestly outperforms a deep feed-forward network (91.9% ± 0.9% vs 90.6% ± 1.2%).

**Label circularity, quantified.** An independent predictor, NetMHCpan-4.1, agreed with the MHCflurry training labels on only 62.24% (323/519) of the held-out test peptides, versus 90.6% for a model trained directly on those labels. This ~28-point gap shows that reported accuracies measure fidelity to a particular labeler rather than true binding. A benchmark of 6,524 experimentally measured affinities (ROC AUC 0.942) independently confirms the quality of the training labels, while a temporal holdout of post-2022 epitopes — unseen by the labeler — establishes generalisation beyond the training distribution.

**A rebalanced pan model and the non-canonical frontier.** Extending the recipe to five common HLA-I alleles (n = 14,488 curated IEDB positives) and a 15-allele panel, a rebalanced pan-specific model matched per-allele specialists on four of five alleles and reached 0.819 pooled non-canonical sensitivity. Canonical binders were recognised at ~0.98 by every predictor class, while non-canonical (Type B) binders plateaued at ~0.65 sensitivity across three independent sequence predictors (our Deep FFN 0.607, MHCflurry 0.668, MUNIS 0.647) — a ceiling shared across sequence architectures that defines the current frontier.

All code, trained models, data, and re-analysis scripts are available at https://github.com/woodhaha/peptide-mhc-benchmark, and the predictor is released as the `pmhcbenchmark` Python package with a command-line tool and a standalone Windows GUI.

This manuscript has not been published previously and is not under consideration elsewhere. All authors have approved the manuscript and agree with its submission. The authors declare no competing interests. No specific grant funding was received.

**AI Tool Disclosure**: We disclose the use of AI-assisted tools (Claude Code, Anthropic, 2026) for code development and manuscript editing. All scientific content, interpretations, and conclusions were authored and approved by the human authors, who take full responsibility for the manuscript's content. No AI tools were listed as authors.

**Suggested Reviewers**:
1. **Professor Morten Nielsen** — Technical University of Denmark. Developer of NetMHCpan; expertise in peptide-MHC binding prediction and immunoinformatics benchmarking.
2. **Professor Bjoern Peters** — La Jolla Institute for Immunology. Director of the Immune Epitope Database (IEDB).
3. **Dr. Timothy O'Donnell** — Memorial Sloan Kettering Cancer Center. Developer of MHCflurry.

We appreciate your consideration and look forward to your response.

Respectfully,

**Zhuha Zhou, MD**<sup>*</sup>
Department of Gastroenterology Surgery
The First Affiliated Hospital of Wenzhou Medical University
Nanbaixiang Street, Ouhai District, 325000 Wenzhou, Zhejiang, China
E-mail: zhouzhuha@wmu.edu.cn
ORCID: 0009-0005-1818-4789

On behalf of all authors:
Yunshi Huang<sup>†</sup>, Weisheng Chen<sup>†</sup>, Yiqi Cai, Xiaopeng Zhuang, Zhuxian Zhou, Xiangyang Xue

<sup>†</sup> These authors contributed equally to this work.
<sup>*</sup> Corresponding Author

---

*Cover letter prepared for submission to PLOS Computational Biology. October 2026.*
