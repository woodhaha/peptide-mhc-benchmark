# Cover Letter — Submission to *Bioinformatics* (Oxford)

**Date**: October 21, 2026

**To**: The Editors
*Bioinformatics*
Oxford University Press

**Re**: Submission of original paper — "The Non-Canonical Frontier in Peptide-MHC-I Binding Prediction: A Benchmark-Hygiene Case Study"

---

Dear Editors,

We are pleased to submit our manuscript for consideration as an Original Paper in *Bioinformatics*. Our study makes a methodological contribution to immunoinformatics benchmarking: we show that benchmark composition can manufacture false conclusions, and we quantify the resulting label-circularity gap.

**Benchmark artifacts are common and consequential.** During a controlled re-analysis of our own six-architecture benchmark (5,190 MHCflurry-labelled HLA-A\*02:01 9-mers, multi-seed estimates with confidence intervals), three apparent findings dissolved: (i) an apparent per-allele superiority over a pan model was a 10.6× training-data-size artifact; (ii) a claimed "label quality dominates architecture" result collapsed once the degenerate random-label null control (which performs at chance, ~33%) was removed; and (iii) a reported "simple architecture beats complex ones" ranking was an artifact of stale model checkpoints — re-training with fixed seeds showed a convolutional network instead modestly outperforms a deep feed-forward network (91.9% ± 0.9% vs 90.6% ± 1.2%, the top four architectures within ~2 pp).

**Label circularity, quantified.** An independent state-of-the-art predictor, NetMHCpan-4.1, agreed with the MHCflurry training labels on only 62.24% (323/519) of the same held-out test peptides, versus 90.6% for a model trained directly on those labels. This ~28-point gap shows that reported accuracies measure fidelity to a particular labeler rather than true binding — a caveat with direct implications for how peptide-MHC benchmarks should be reported.

**Experimental label validation.** A benchmark of 6,524 experimentally measured binding affinities (IC50) confirms the quality of the training labels (the MHCflurry-trained model retains ROC AUC 0.942 on this set), while a temporal holdout of post-2022 epitopes — unseen by the labeler — establishes generalisation beyond the training distribution.

**The non-canonical frontier.** A zero-parameter anchor rule matches the trained model on canonical IEDB peptides, while non-canonical binders plateau at ~0.65 sensitivity across three independent sequence predictors (a 15-allele panel) — identifying where the genuinely hard problem lies.

We believe this work aligns with *Bioinformatics*' scope of methodological benchmarks and reproducible computational pipelines. All code, trained models, data, and re-analysis scripts are available at https://github.com/woodhaha/peptide-mhc-benchmark, and the predictor is released as the `pmhcbenchmark` Python package (pip-installable) with a command-line tool and a standalone Windows GUI.

This manuscript has not been published previously and is not under consideration elsewhere. All authors have approved the manuscript and agree with its submission. The authors declare no competing interests. No specific grant funding was received.

**AI Tool Disclosure**: In accordance with journal guidelines, we disclose the use of AI-assisted tools (Claude Code, Anthropic, 2026) for code development and manuscript editing. All scientific content, interpretations, and conclusions were authored and approved by the human authors, who take full responsibility for the manuscript's content. No AI tools were listed as authors.

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

*Cover letter prepared for submission to Bioinformatics. October 2026.*
