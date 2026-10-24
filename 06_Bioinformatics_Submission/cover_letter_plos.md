# Cover Letter — Submission to *PLOS Computational Biology*

**Date**: October 24, 2026

**To**: The Editors
*PLOS Computational Biology*

**Re**: Submission of a Benchmark — "Benchmark Hygiene Reveals the Non-Canonical Frontier of Peptide-MHC-I Binding Prediction"

---

Dear Editors,

We are pleased to submit our manuscript for consideration as a **Benchmark article** in *PLOS Computational Biology*, following the benchmark category you established (Peters et al., *PLoS Comput Biol* 2018;14:e1006494).

Our study is a **benchmark-hygiene case study** in immunoinformatics that names a failure mode — *cutpoint sensitivity* in rank-based labels — and demonstrates its countermeasures through a controlled single-allele benchmark. Four points summarise the contribution:

1. **We quantify cutpoint sensitivity as a benchmark-hygiene failure mode.** Peptide-MHC-I predictors are trained on computationally labelled data, so reported accuracy can reproduce a labeler's decision boundary rather than measure true binding. MHCflurry and NetMHCpan both reach ~0.95 AUC on 6,524 measured affinities and rank peptides almost identically (Spearman ρ = 0.96), yet agree on only 62% of hard labels — a gap that rises to 96% under the natural class prior and reflects stricter percentile cutpoints in NetMHCpan, not label definition or fidelity. We distil the countermeasures into a concrete, reusable benchmark-hygiene checklist (Box 1): a true i.i.d. random-label null, frozen snapshots with multi-seed confidence intervals, continuous-score rank correlation against a pinned state-of-the-art comparator, and reporting under a realistic class prior.

2. **The framework caught its own ground truth.** During re-analysis we found that seven of our 49 "HLA-A\*02:01" benchmark epitopes were mis-annotated — two with no IEDB record, four restricted to other alleles, and one with no positive A\*02:01 record — including the classic CMV pp65 A\*24:02 epitope `QYDPVAALF`, which we had mis-read as a non-canonical A\*02:01 false negative. Correcting this changed the headline sensitivity from 95.8% to 97.6%, and removed a spurious "non-canonical" data point from our argument.

3. **The 62% disagreement is largely a sampling and cutpoint artifact.** Under the natural class prior, cross-labeler agreement is 96% rather than 62%; the gap reflects class-balanced sampling and stricter percentile cutpoints in one labeler, not a fundamental disagreement or label error.

4. **The single-allele scope is a methodological choice, not a limitation.** We held the allele constant (HLA-A\*02:01) to isolate label-source variation from allele-driven confounds. The cutpoint offset we quantify is itself allele-specific: re-scoring the same peptides for A\*01:01 and A\*03:01 drops rank agreement to ρ = 0.73 and 0.82 (from 0.96) precisely because those peptides are A\*02:01-selected and mostly non-binders elsewhere — evidence that a per-allele balanced benchmark, not a pooled multi-allele one, is the correct unit of analysis. We supply the A\*02:01 benchmark, the Box 1 checklist, and the code so others can build the remaining per-allele benchmarks.

**Fit to the benchmark criteria.** Our manuscript meets the seven criteria for a good benchmark (Peters et al., 2018, Box 1): (1) peptide-MHC binding prediction is an active research area; (2) we benchmark the leading publicly available predictors spanning the field's two dominant objectives — MHCflurry and NetMHCpan (binding affinity) and MUNIS and BigMHC (presentation/eluted-ligand) — alongside a student model; (3) the input and expected-output datasets (the 6,524-affinity benchmark and the re-verified 42-epitope set) are freely available and reusable; (4) results are fully reproducible, with published code and a frozen manifest of file hashes; (5) we report a small set of orthogonal, practically relevant metrics (AUC, sensitivity, and a cutpoint-offset diagnostic); (6) all benchmarked methods are public; and (7) our own Deep FFN is used as a probe of label calibration rather than a novel tool, a caveat we state prominently.

The benchmark-hygiene findings are the point. Prior work (Kapoor & Narayanan, *Patterns* 2023; Bernett et al., *Nat Methods* 2024) established that training-data contamination drives the reproducibility crisis; we extend this taxonomy to the subtler, peptide-MHC-specific failure of cutpoint sensitivity in rank-based labels. We distil the countermeasures into a reusable benchmark-hygiene checklist (Box 1) that others can apply to their own benchmarks.

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
