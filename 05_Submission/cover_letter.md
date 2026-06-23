# Cover Letter

**Date**: June 19, 2026

**To**: The Editors  
*Briefings in Bioinformatics*  
Oxford University Press

**Re**: Submission of manuscript — "Systematic Benchmarking of Deep Learning Architectures for Peptide-MHC Binding Prediction Reveals Label Quality, Not Architecture, Defines Performance"

---

Dear Editors,

We are pleased to submit our manuscript entitled "Systematic Benchmarking of Deep Learning Architectures for Peptide-MHC Binding Prediction Reveals Label Quality, Not Architecture, Defines Performance" for consideration for publication in *Briefings in Bioinformatics*.

### Significance and Fit

Computational prediction of peptide-MHC binding is a foundational task in immunoinformatics with direct applications in cancer immunotherapy and neoantigen discovery. Our manuscript addresses three gaps we believe are of broad interest to the *Briefings in Bioinformatics* readership.

First, and most importantly, we quantify the impact of training label quality on downstream model performance — a question that is widely acknowledged but rarely measured systematically. Across three labelling strategies (MHCflurry 2.2.0, a position-specific scoring matrix based on crystallographic anchor preferences, and random synthetic labels), label source accounted for ~29 percentage points of accuracy variation, exceeding the maximum inter-architecture difference of 10.8 percentage points. A simple PSSM baseline achieved 94.8% accuracy, outperforming all neural network architectures. These findings carry direct implications for how the field benchmarks and reports model performance, particularly given emerging evidence that over 55% of IEDB entries are computationally rather than experimentally labelled (Preibisch et al., 2026). We recommend that future benchmarking studies explicitly report label provenance.

Second, we systematically benchmark six deep learning architectures under fully controlled conditions (identical training data, encoding, and evaluation protocols) and evaluate ESM-2 protein language model embeddings against the standard BLOSUM62 encoding. ESM-2 t6 per-position embeddings achieved the highest neural network accuracy (93.3%), modestly exceeding BLOSUM62 (91.9%), while a mean-pooled ablation collapsed to 65.9%, confirming that positional information is essential for MHC-I binding prediction. The shallower t6 embeddings (8M parameters, 2,880 dimensions) outperformed the deeper t12 variant (35M parameters, 4,320 dimensions), revealing a data-to-dimensionality matching constraint with practical implications for groups working with modestly sized training sets.

Third, we demonstrate a complete translational pipeline: protein epitope scanning (10 proteins, 3,536 peptides), cancer hotspot mutation analysis (16 mutations across p53 and KRAS), and external validation against 49 experimentally verified IEDB epitopes (93.9% sensitivity, 75.0% specificity, ROC AUC 0.947). Molecular docking generates the testable structural hypothesis that a K-E63 salt bridge may compensate for the non-canonical lysine at P2 of the KRAS G12V neoepitope (YKLVVVGAV). We present this as a docking-derived hypothesis requiring experimental validation by X-ray crystallography or cryo-EM.

### Novel Contributions

- **Label quality, not architecture, defines the performance ceiling**: Systematic quantification across three labelling strategies, demonstrating that label provenance accounts for more performance variation (~29 pp) than architecture choice (~11 pp). A simple PSSM baseline (94.8%) outperforms all neural networks, underscoring that biophysical knowledge should be included as a comparator in computational benchmarking.
- **ESM-2 embeddings offer modest gains over BLOSUM62**: ESM-2 t6 per-position embeddings (93.3%) modestly exceed BLOSUM62 (91.9%), while mean-pooled embeddings collapse to 65.9%, confirming positional information is essential. The shallower t6 variant outperforms the deeper t12 variant due to data-to-dimensionality matching — a practical finding for groups applying protein language models to peptide-MHC tasks.
- **Systematic six-architecture comparison under controlled conditions**: Confirms that architectural complexity beyond additional depth with batch normalisation provides diminishing returns for positional binding tasks; the Deep FFN (152,823 parameters) outperforms CNN (1,053,863), LSTM (23,939), and ResNet-style (324,931) architectures.
- **K-E63 salt bridge as a testable structural hypothesis**: Molecular docking predicts electrostatic compensation between a charged P2 lysine and the conserved Glu63 in the HLA-A\*02:01 B-pocket as a potential non-canonical anchor compensation mechanism. This is presented as a docking-derived hypothesis requiring experimental validation.
- **Fully reproducible open-source pipeline**: All code (1,700+ lines, R/Python), trained models, feature matrices, and all prediction results are provided in a single data package with a documented decision log.

### Editorial Considerations

This manuscript has not been published previously and is not under consideration elsewhere. All authors have approved the manuscript and agree with its submission to *Briefings in Bioinformatics*.

We have no conflicts of interest to declare. This work did not receive specific grant funding.

### Suggested Reviewers

1. **Professor Morten Nielsen** — Department of Health Technology, Technical University of Denmark (DTU), Kongens Lyngby, Denmark  
   *Expertise*: Immunoinformatics, peptide-MHC binding prediction, NetMHCpan development

2. **Associate Professor Jiangning Song** — Biomedicine Discovery Institute, Monash University, Melbourne, Australia  
   *Expertise*: Systematic benchmarking of bioinformatics methods, HLA-I peptide binding prediction tools

3. **Professor Bjoern Peters** — La Jolla Institute for Immunology (LJI), La Jolla, California, USA  
   *Expertise*: Immune Epitope Database (IEDB), computational immunology, epitope prediction benchmarking

### AI Tool Disclosure

In accordance with *Briefings in Bioinformatics* author guidelines, we disclose the use of AI-assisted tools (Claude Code, Anthropic, June 2026) for code development, literature search, and manuscript editing. All scientific content, interpretations, and conclusions were authored and approved by the human authors, who take full responsibility for the manuscript's content. No AI tools were listed as authors.

---

We appreciate your consideration of our manuscript and look forward to your response.

Sincerely,

**Zhuha Zhou, MD**<sup>†,\*</sup>  
Department of Gastroenterology Surgery  
The First Affiliated Hospital of Wenzhou Medical University  
Nanbaixiang Street, Ouhai District, 325000 Wenzhou, Zhejiang, China  
E-mail: zhouzhuha@wmu.edu.cn  
ORCID: 0009-0005-1818-4789

**Yongyu Bai, MD**<sup>†</sup>  
Department of Gastroenterology Surgery  
The First Affiliated Hospital of Wenzhou Medical University

**Qigang Xu, MD**  
Department of Hepatobiliary and Pancreatic Surgery  
The First Affiliated Hospital of Wenzhou Medical University

**Zhuxian Zhou, PhD**  
Center for Bionanoengineering, College of Chemical and Biological Engineering  
Zhejiang University, Hangzhou 310027, China

**Shaoliang Han, MD**  
Department of Gastroenterology Surgery  
The First Affiliated Hospital of Wenzhou Medical University

<sup>†</sup> These authors contributed equally to this work.  
<sup>\*</sup> Corresponding Author
