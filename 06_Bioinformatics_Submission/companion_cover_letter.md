# Cover Letter — Submission to *PLOS Computational Biology*

**Date**: August 3, 2026

**To**: The Editors
*PLOS Computational Biology*

**Re**: Submission of research article — "A Rebalanced Pan-Specific Model and the Non-Canonical Frontier: A Multi-Allele Benchmark of Peptide-MHC-I Binding Prediction"

---

Dear Editors,

We are pleased to submit our research article for consideration in *PLOS Computational Biology*. This work addresses two open questions in computational immunology: whether a single pan-specific sequence model can replace per-allele specialists, and how hard non-canonical peptide-MHC-I binders are for current prediction methods — including structure-based approaches.

**One rebalanced pan model replaces a per-allele panel.** We benchmarked a rebalanced pan-specific deep feed-forward network across five common HLA-I alleles (n = 14,488 curated IEDB positives) and a 15-allele panel, comparing against per-allele models, MHCflurry, and the structure-based predictor MUNIS. Across five seeds, the rebalanced pan model reached non-canonical sensitivity 0.92–1.00, matching per-allele specialists. We show that an apparent per-allele advantage reported in prior work was a 10.6× training-data-size artifact — a methodological warning that benchmark composition can manufacture apparent model differences.

**The non-canonical frontier.** Canonical binders were recognized at ~0.98 by every predictor class. In contrast, non-canonical (Type B) binders plateaued at ~0.65 sensitivity across DFFN (0.607), MHCflurry (0.668), and MUNIS (0.647) — a structure-based predictor achieves the same sensitivity ceiling as sequence-based methods, directly challenging the field's assumption that structural scoring will solve the non-canonical problem. Type A binders were inherently hard (0.477) and not recovered by broad-motif retraining. Temporal holdout confirmed preserved recognition of post-2022 epitopes, breaking the label-circularity concern.

We believe this work is well aligned with *PLOS Computational Biology*'s scope of rigorous computational methods research. The study provides: (i) a fully controlled multi-allele benchmark with seed-robust estimates and protein-background negatives; (ii) a direct test of whether structure-based prediction lifts the non-canonical sensitivity ceiling; and (iii) a reproducible pipeline with all code, models, and curated datasets publicly available at https://github.com/woodhaha/peptide-mhc-benchmark.

This manuscript has not been published previously and is not under consideration elsewhere. All authors have approved the manuscript and agree with its submission. The authors declare no competing interests. No specific grant funding was received.

**AI Tool Disclosure**: We disclose the use of AI-assisted tools (Claude Code, Anthropic, 2026) for code development and manuscript editing. All scientific content, interpretations, and conclusions were authored and approved by the human authors, who take full responsibility for the manuscript's content. No AI tools were listed as authors.

**Suggested Reviewers**:
1. **Professor Morten Nielsen** — Technical University of Denmark. Developer of NetMHCpan; expertise in peptide-MHC binding prediction and immunoinformatics benchmarking.
2. **Professor Bjoern Peters** — La Jolla Institute for Immunology. Director of the Immune Epitope Database (IEDB).
3. **Dr. Timothy O'Donnell** — Memorial Sloan Kettering Cancer Center. Developer of MHCflurry.

We appreciate your consideration and look forward to your response.

Respectfully,

**Zhuha Zhou, MD**<sup>†,*</sup>
Department of Gastroenterology Surgery
The First Affiliated Hospital of Wenzhou Medical University
Nanbaixiang Street, Ouhai District, 325000 Wenzhou, Zhejiang, China
E-mail: zhouzhuha@wmu.edu.cn
ORCID: 0009-0005-1818-4789

On behalf of all authors:
Weisheng Chen<sup>†</sup>, Yiqi Cai, Xiaopeng Zhuang, Zhuxian Zhou

<sup>†</sup> These authors contributed equally to this work.
<sup>*</sup> Corresponding Author

---

*Cover letter prepared for submission to PLOS Computational Biology. August 2026.*
