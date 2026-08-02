# Cover Letter — Submission to *Bioinformatics* (Oxford)

**Date**: August 3, 2026

**To**: The Editors
*Bioinformatics*
Oxford University Press

**Re**: Submission of original paper — "Label Quality Defines the Performance Ceiling for Peptide-MHC Binding Prediction: A Six-Architecture Benchmark with Multi-Allele Generalization"

---

Dear Editors,

We are pleased to submit our manuscript for consideration as an Original Paper in *Bioinformatics*. Our study addresses a foundational question in immunoinformatics: what actually defines the performance ceiling for peptide-MHC-I binding prediction?

**Label quality, not architecture, defines the performance ceiling.** Across three labelling strategies (deterministic PSSM, MHCflurry 2.2.0, random synthetic), label source accounted for ~29 percentage points of accuracy variation — exceeding the maximum inter-architecture difference (~11 pp) by nearly threefold — and the variation concentrated exclusively in the weak-binder class (F1 0.723–0.925). We establish WB F1 as a quantitative sentinel for label contamination: ML-derived training labels, themselves trained on >55% computationally labelled IEDB data, compound the ambiguous SB–WB–NB boundary across predictor generations. We recommend reporting WB F1 and label provenance as standard benchmark practice.

We systematically benchmarked six architectures under fully controlled conditions (identical training data, encoding, evaluation): a deep feed-forward network achieved 91.9% accuracy, outperforming CNN (90.0%), ResNet (84.3%), LSTM (83.3%), and Random Forest (81.1%). BLOSUM62 generalized more reliably than ESM-2 per-position embeddings under cross-validation (92.2% vs. 82.2%), despite ESM-2's single-split advantage — a practical finding for groups applying protein language models to modestly sized peptide datasets.

External validation against 49 experimentally verified IEDB epitopes yielded 93.9% sensitivity (ROC AUC 0.947). A multi-allele extension showed that one rebalanced pan-specific model matches per-allele specialists across five common HLA-I alleles (0.819 pooled sensitivity on a 15-allele panel), while non-canonical binders plateau at ~0.65 sensitivity across sequence- and structure-based predictors alike — the energy-compensated tail (e.g., the KRAS G12V P2-Lys salt bridge) is genuinely prediction-resistant to all current methods.

We believe this work aligns well with *Bioinformatics*' scope of methodological benchmarks and reproducible computational pipelines. All code, trained models, feature matrices, and curated datasets are available at https://github.com/woodhaha/peptide-mhc-benchmark.

This manuscript has not been published previously and is not under consideration elsewhere. All authors have approved the manuscript and agree with its submission. The authors declare no competing interests. No specific grant funding was received.

**AI Tool Disclosure**: In accordance with journal guidelines, we disclose the use of AI-assisted tools (Claude Code, Anthropic, 2026) for code development and manuscript editing. All scientific content, interpretations, and conclusions were authored and approved by the human authors, who take full responsibility for the manuscript's content. No AI tools were listed as authors.

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
Weisheng Chen<sup>†</sup>, Yiqi Cai, Yongyu Bai, Zhuxian Zhou

<sup>†</sup> These authors contributed equally to this work.
<sup>*</sup> Corresponding Author

---

*Cover letter prepared for submission to Bioinformatics. August 2026.*
