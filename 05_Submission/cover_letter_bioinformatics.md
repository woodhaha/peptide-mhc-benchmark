# Cover Letter — Submission to *Briefings in Bioinformatics*

**Date**: July 3, 2026

**To**: The Editors
*Briefings in Bioinformatics*
Oxford University Press

**Re**: Submission of original research article -- "Systematic Benchmarking of Deep Learning Architectures for Peptide-MHC Binding Prediction Reveals Label Quality, Not Architecture, Defines Performance"

---

Dear Editors,

**Label quality, not architecture, defines the performance ceiling.** We submit a manuscript that demonstrates, through systematic benchmarking, that the immunoinformatics field should invest in better training labels rather than more complex neural architectures. Across three labelling strategies, label source accounted for ~29 percentage points of accuracy variation, exceeding the maximum inter-architecture difference of ~11 percentage points. A simple position-specific scoring matrix (PSSM) baseline achieved 94.8% accuracy -- outperforming all six deep learning architectures tested.

We systematically compared six architectures under fully controlled conditions (identical training data, encoding, and evaluation protocols): a deep feed-forward network (FFN, 152K parameters) achieved 91.9% accuracy, outperforming CNN (1.05M), ResNet-style (325K), and LSTM (24K) architectures, confirming that architectural complexity beyond depth with batch normalisation provides diminishing returns for this well-characterised anchor-motif task. ESM-2 t6 per-position embeddings modestly improved accuracy to 93.3%, while mean-pooled embeddings collapsed to 65.9%, providing direct evidence that residue-position-specific features are essential for MHC-I binding.

External validation against 49 experimentally verified IEDB epitopes yielded 93.9% sensitivity (ROC AUC 0.947). In a translational application, cancer hotspot mutation scanning across p53 and KRAS identified neoepitope candidates including p53 R248W and KRAS G12V. Molecular docking of the KRAS G12V neoepitope (YKLVVVGAV) predicts a K-E63 salt bridge with HLA-A*02:01 (estimated Coulombic stabilisation of approximately -30 kcal/mol vs. approximately -3 kcal/mol for canonical hydrophobic P2 burial), representing a testable structural hypothesis for a charged P2 anchor compensation mechanism.

These findings connect directly to Preibisch et al. (2026), who report that over 55% of IEDB entries are computationally rather than experimentally labelled. Our results quantify the downstream consequences of this label provenance problem: models trained on derivative labels inherit the biases of their teachers, and the weak binder class is the primary locus of this variation (F1: 0.925 [PSSM] vs. 0.880 [MHCflurry]). We recommend that future benchmarking studies explicitly report label provenance and include biophysical baselines as comparators.

All code, trained models, feature matrices, and prediction results are provided as a fully reproducible open-source pipeline (R and Python, with Jupyter notebooks documenting each analytical step), publicly available at https://github.com/woodhaha/peptide-epitope-prediction.

We believe this work aligns well with *Briefings in Bioinformatics*, which has a strong tradition of publishing methodological benchmark studies and reproducible computational pipelines of broad interest to the immunoinformatics and computational biology communities.

This manuscript has not been published previously and is not under consideration elsewhere. All authors have approved the manuscript and agree with its submission. The authors declare no competing interests. No specific grant funding was received for this work.

**Suggested Reviewers**

1. **Professor Morten Nielsen** -- Department of Health Technology, Technical University of Denmark (DTU), Kongens Lyngby, Denmark. Developer of NetMHCpan; leading expertise in peptide-MHC binding prediction and immunoinformatics benchmarking.

2. **Professor Bjoern Peters** -- La Jolla Institute for Immunology (LJI), La Jolla, California, USA. Director of the Immune Epitope Database (IEDB); expertise in computational immunology and epitope prediction benchmarking standards.

3. **Dr. Timothy O'Donnell** -- Memorial Sloan Kettering Cancer Center, New York, USA (formerly OpenVax). Developer of MHCflurry; deep expertise in deep learning architectures for peptide-MHC binding prediction.

We appreciate your consideration and look forward to your response.

Respectfully,

**Zhuha Zhou, MD**<sup>dagger,\*</sup>
Department of Gastroenterology Surgery
The First Affiliated Hospital of Wenzhou Medical University
Nanbaixiang Street, Ouhai District, 325000 Wenzhou, Zhejiang, China
E-mail: zhouzhuha@wmu.edu.cn
ORCID: 0009-0005-1818-4789

On behalf of all authors:
Yongyu Bai<sup>dagger</sup>, Qigang Xu, Zhuxian Zhou, Shaoliang Han

<sup>dagger</sup> These authors contributed equally to this work.
<sup>\*</sup> Corresponding Author

---

*Cover letter prepared for submission to Briefings in Bioinformatics. July 2026.*
