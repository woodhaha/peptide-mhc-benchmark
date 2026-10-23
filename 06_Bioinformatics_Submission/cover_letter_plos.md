# Cover Letter — Submission to *PLOS Computational Biology*

**Date**: October 23, 2026

**To**: The Editors
*PLOS Computational Biology*

**Re**: Submission of a Research Article — "Label-Fidelity Decomposition Reveals Label Circularity in Peptide-MHC-I Benchmarking"

---

Dear Editors,

We are pleased to submit our manuscript for consideration as a Research Article in *PLOS Computational Biology*.

Our study is a **benchmark-hygiene case study** in immunoinformatics that names a failure mode — *label circularity* — distinct from the data leakage now widely documented, and demonstrates its countermeasures by turning them on our own prior analysis. Three points summarise the contribution:

1. **We name label circularity as a failure mode distinct from data leakage.** Leakage is the same data reappearing in train and test; label circularity is subtler — peptide-MHC-I predictors are trained on computationally labelled data, so reported accuracy can reproduce the labeler's decision boundary rather than measure true binding. We separate a labelling scheme's two orthogonal properties — *label self-consistency* versus *cross-labeler transferability* (an independent state-of-the-art predictor, NetMHCpan-4.1, agrees on only 62.24% of the same test peptides) — and distil the countermeasures into a concrete, reusable benchmark-hygiene checklist (Box 1): equalised training data, a true i.i.d. random-label null, frozen snapshots with multi-seed confidence intervals, label provenance with a temporal holdout and a pinned state-of-the-art comparator, and reporting under the natural class prior.

2. **The framework caught its own ground truth.** During re-analysis we found that seven of our 49 "HLA-A\*02:01" benchmark epitopes were mis-annotated — two with no IEDB record, four restricted to other alleles, and one with no positive A\*02:01 record — including the classic CMV pp65 A\*24:02 epitope `QYDPVAALF`, which we had mis-read as a non-canonical A\*02:01 false negative. Correcting this changed the headline sensitivity from 95.8% to 97.6%, and removed a spurious "non-canonical" data point from our argument.

3. **A threshold-sensitivity analysis qualified our own headline.** The apparent ~0.65 "non-canonical plateau" across four predictors is partly a threshold artifact: at matched operating points the predictors diverge (Type B sensitivity 0.45–0.65), rather than converging on a shared ceiling.

The self-corrections are the point. Prior work (Kapoor & Narayanan, *Patterns* 2023; Bernett et al., *Nat Methods* 2024) established that training-data contamination drives the reproducibility crisis; we extend this taxonomy to the subtler, peptide-MHC-specific failure of label circularity. Rather than asserting a checklist in the abstract, we demonstrate that it works by letting it catch real errors in our own prior analysis — a benchmark-hygiene framework that detects its own artifacts is one others can adopt.

All code, trained models, data, and re-analysis scripts are available at https://github.com/woodhaha/peptide-mhc-benchmark, and the predictor is released as the pip-installable `pmhcbenchmark` package. A data availability statement, CRediT author contributions, and competing-interests declaration are included.

This manuscript has not been published previously and is not under consideration elsewhere. All authors have approved the manuscript and agree with its submission. The authors declare no competing interests. No specific grant funding was received.

**AI Tool Disclosure**: We disclose the use of AI-assisted tools (Claude Code, Anthropic, 2026) for code development and manuscript editing. All scientific content, interpretations, and conclusions were authored and approved by the human authors, who take full responsibility for the manuscript's content. No AI tools were listed as authors.

**Suggested Reviewers**:
1. **Professor Morten Nielsen** — Technical University of Denmark. Developer of NetMHCpan; peptide-MHC binding prediction and immunoinformatics benchmarking.
2. **Professor Bjoern Peters** — La Jolla Institute for Immunology. Director of the Immune Epitope Database (IEDB).
3. **Professor Sayash Kapoor** — Princeton University. Reproducibility and data-leakage methodology in machine learning.

We appreciate your consideration and look forward to your response.

Respectfully,

**Xiangyang Xue, Professor**<sup>*</sup>
Zhejiang Key Laboratory of Intelligent Cancer Biomarker Discovery and Translation
The First Affiliated Hospital of Wenzhou Medical University, Wenzhou, China
E-mail: wzxxy@wmu.edu.cn · ORCID: 0000-0002-8670-2016

**Zhuha Zhou, MD**<sup>*</sup>
Department of Gastroenterology Surgery
The First Affiliated Hospital of Wenzhou Medical University
Nanbaixiang Street, Ouhai District, 325000 Wenzhou, Zhejiang, China
E-mail: zhouzhuha@wmu.edu.cn · ORCID: 0009-0005-1818-4789

On behalf of all authors: Yunshi Huang<sup>†</sup>, Weisheng Chen<sup>†</sup>, Yiqi Cai, Xiaopeng Zhuang, Zhuxian Zhou

<sup>†</sup> These authors contributed equally to this work.
<sup>*</sup> Corresponding authors
