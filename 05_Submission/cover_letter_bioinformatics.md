# Cover Letter — Submission to *Bioinformatics*

---

**Date**: June 24, 2026

**To**: The Editors  
*Bioinformatics*  
Oxford University Press

**Re**: Submission of original research article — "Systematic Benchmarking of Deep Learning Architectures for Peptide-MHC Binding Prediction Reveals Label Quality, Not Architecture, Defines Performance"

---

Dear Editors,

We wish to submit an original research article entitled **"Systematic Benchmarking of Deep Learning Architectures for Peptide-MHC Binding Prediction Reveals Label Quality, Not Architecture, Defines Performance"** for consideration for publication in *Bioinformatics*.

### Summary of the Work

Computational prediction of peptide-MHC binding is a cornerstone of epitope-based vaccine design and cancer immunotherapy target identification. Despite rapid proliferation of deep learning tools in this space, systematic architecture comparisons under controlled conditions and quantification of label quality effects remain limited.

This study makes five contributions:

1. **Label quality quantitation.** Across three labelling strategies (MHCflurry 2.2.0, deterministic PSSM, random synthetic), label source accounted for ~29 percentage points of accuracy variation — exceeding the maximum inter-architecture difference (~11 pp). The performance gap is concentrated in the weak binder class (F1: 0.925 [PSSM] vs. 0.880 [MHCflurry]). Given that over 55% of IEDB entries are computationally labelled, we recommend explicit reporting of label provenance in future benchmarking.

2. **Six-architecture controlled comparison.** A deep feed-forward network (152K parameters) achieved 91.9% accuracy, outperforming CNN (1.05M), ResNet-style (325K), and LSTM (24K) architectures — confirming that architectural complexity beyond depth with batch normalisation provides diminishing returns for this well-characterised anchor-motif task.

3. **ESM-2 encoding evaluation.** ESM-2 t6 per-position embeddings (2,880-dim) modestly improved accuracy to 93.3%, while mean-pooled embeddings collapsed to 65.9%, providing direct evidence that residue-position-specific features are essential for MHC-I binding.

4. **IEDB external validation.** The model achieved 93.9% sensitivity (46/49 experimentally verified epitopes, ROC AUC 0.947), and correctly identified 10 of 11 known epitopes across 10 scanned proteins.

5. **K-E63 salt bridge structural hypothesis.** For the KRAS G12V neoepitope, template-based structural analysis and energy minimisation identify a predicted lysine-Glu63 salt bridge (estimated Coulombic stabilisation −30 kcal/mol vs. −3 kcal/mol for canonical hydrophobic P2 burial), representing — to our knowledge — the first quantitative structural evidence for a charged P2 anchor compensation mechanism in HLA-A*02:01.

All code, trained models, feature matrices, and prediction results are provided as an open-source reproducible pipeline. The analysis is implemented in R and Python, with Jupyter notebooks documenting each analytical step.

### Fit with *Bioinformatics*

This manuscript aligns with *Bioinformatics*' scope in several respects:

- **Original methods**: Controlled architecture comparison + label provenance analysis under a unified computational framework
- **Reproducibility**: Complete open-source pipeline with documented code, trained models, and all intermediate data
- **Biological relevance**: Direct application to cancer neoepitope discovery (p53, KRAS) with structural validation
- **Broad interest**: The label quality finding has implications beyond immunoinformatics — it addresses a systemic concern in computational biology where machine-learning-derived labels are used as ground truth

### Declarations

- This manuscript has not been published previously and is not under consideration elsewhere.
- All authors have read and approved the manuscript and agree with its submission to *Bioinformatics*.
- The authors declare no competing interests.
- No specific funding was received for this work.
- All code and data are publicly available at: https://github.com/woodhaha/peptide-epitope-prediction (DOI to be assigned upon acceptance).

### Suggested Reviewers

We suggest the following potential reviewers with relevant expertise:

1. **Dr. Morten Nielsen** — Technical University of Denmark — developer of NetMHCpan; expertise in peptide-MHC binding prediction
2. **Dr. Tim O'Donnell** — Memorial Sloan Kettering Cancer Center — developer of MHCflurry; expertise in computational immunology
3. **Dr. John Sidney** — La Jolla Institute for Immunology — expertise in MHC binding assays and epitope validation
4. **Dr. Sebastian B. Mohr** — ETH Zurich — expertise in structural bioinformatics and deep learning for molecular interactions
5. **Dr. Yuan Luo** — Northwestern University — expertise in machine learning for biomedical applications

We would be happy to provide any additional information required during the review process.

Respectfully,

**Zhuha Zhou** (Corresponding author)  
Department of Gastroenterology Surgery  
The First Affiliated Hospital of Wenzhou Medical University  
Nanbaixiang Street, Ouhai District  
325000 Wenzhou, Zhejiang, China  
E-mail: zhouzhuha@wmu.edu.cn  
ORCID: 0009-0005-1818-4789

On behalf of all authors:  
Yongyu Bai, Qigang Xu, Zhuxian Zhou, Shaoliang Han

---
*Cover letter prepared for submission to Bioinformatics. June 2026.*
