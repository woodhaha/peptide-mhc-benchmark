# Deep Learning for Epitope Prediction: A Literature Map

> Comprehensive review of the computational epitope prediction landscape, prepared for manuscript writing.
> Last updated: 2026-06-23

---

## Table of Contents

1. [Scope and Definitions](#1-scope-and-definitions)
2. [Landscape Overview](#2-landscape-overview)
3. [Major Methodological Paradigms](#3-major-methodological-paradigms)
4. [Tools & Models: Comparative Analysis](#4-tools--models-comparative-analysis)
5. [Datasets and Benchmarking](#5-datasets-and-benchmarking)
6. [TCRpMHC Interaction Prediction](#6-tcrpmhc-interaction-prediction)
7. [Generative Epitope Design](#7-generative-epitope-design)
8. [Clinical Translation & Validation](#8-clinical-translation--validation)
9. [Knowledge Gaps & Open Problems](#9-knowledge-gaps--open-problems)
10. [Annotated Bibliography of Key Papers](#10-annotated-bibliography-of-key-papers)
11. [Tool Availability & Code Index](#11-tool-availability--code-index)

---

## 1. Scope and Definitions

### Epitope Types

| Type | Definition | Prediction Difficulty |
|------|-----------|----------------------|
| **Linear B-cell epitope** | Continuous stretch of amino acids recognized by antibodies | Moderate |
| **Conformational B-cell epitope** | Discontinuous residues brought together by 3D folding, recognized by antibodies | **Hard** (AUC typically 0.580.80) |
| **T-cell epitope (MHC-I)** | 811mer peptide bound to HLA class I, recognized by CD8+ T cells | Relatively mature (AUC 0.920.98) |
| **T-cell epitope (MHC-II)** | 1325mer peptide bound to HLA class II, recognized by CD4+ T cells | Moderate (AUC 0.820.90) |

### Key Terminology
- **Immunogenicity**  the ability to trigger an immune response (broader than binding)
- **Antigenicity**  the ability to be recognized by antibodies or TCRs (binding alone)
- **pMHC**  peptideMHC complex; the ligand recognized by T-cell receptors
- **Neoantigen / Neoepitope**  tumor-specific peptide arising from somatic mutation


---

## 2. Landscape Overview

### Historical Trajectory

| Era | Dominant Methods | Representative Tools |
|-----|-----------------|---------------------|
| **Pre-2015** | Sequence motifs, position-specific scoring matrices (PSSM) | SYFPEITHI, BIMAS, ProPred |
| **20152019** | Shallow ML: SVM, Random Forest, early ANNs | NetMHCpan 3.0/4.0, MHCflurry 1.0, BepiPred-2.0 |
| **20202022** | Deep neural networks, protein language models emerge | NetMHCpan 4.1, MHCflurry 2.0, BepiPred-3.0, DeepImmuno |
| **20232025 (current)** | Transformers + GNNs + structural data, hybrid architectures | MUNIS, BigMHC, GraphEPN, DiscoTope 3.0, EpiGraph |

### Key Shift (20232025)

The field has undergone a paradigm shift driven by three converging forces:

1. **Protein language models** (ESM-2, ProtBERT, ProtTrans) replaced hand-crafted sequence features as the dominant encoding strategy.
2. **AlphaFold2/3** made structural data widely available, enabling structure-conditioned models for conformational epitopes.
3. **Graph neural networks** emerged as the architecture of choice for capturing spatial clustering of epitope residues on protein surfaces.

### Current State Summary

- **T-cell epitope prediction (MHC-I binding):** Mature field; top tools achieve AUROC >0.95. The challenge has shifted from binding prediction to *immunogenicity* prediction.
- **B-cell epitope prediction (conformational):** Active frontier; best tools achieve AUCPR 0.230.25, reflecting the inherent difficulty of the problem.
- **TCRpMHC prediction:** Exploding new subfield; 30+ methods in 20242025 alone. Generalization to unseen epitopes remains poor.
- **Generative design:** Nascent but rapidly growing; diffusion models show particular promise for structure-guided epitope library design.

### Top Venues Publishing in This Area (20242025)

- *Briefings in Bioinformatics* (2 major reviews)
- *npj Vaccines* (systematic review with comprehensive benchmark table)
- *Frontiers in Immunology* (multiple methods papers + reviews)
- *PROTEOMICS* (CD8+ T-cell epitope review)
- *Nature Machine Intelligence* (UniPMT framework)
- *NeurIPS* (AsEP benchmark)
- *bioRxiv* (fast-moving preprints dominate new method releases)

---

## 3. Major Methodological Paradigms

### 3.1 Protein Language Models (PLMs) as Encoders

PLMs pre-trained on massive protein sequence databases (UniRef, BFD) now serve as the foundation for most state-of-the-art epitope predictors.

| PLM | Size | Training Data | Used By |
|-----|------|--------------|---------|
| **ESM-2 (650M15B params)** | Up to 15B | UniRef90 | MUNIS, BepiPred-3.0, EpiGraph, EpitopeTransfer, TransHLA |
| **ESM-1v** | 650M | UniRef90 (five variants) | SEMA 1D |
| **ESM-IF1** (inverse folding) | ~100M | CATH structures | SEMA 3D, DiscoTope 3.0 |
| **ProtBERT** | 420M | BFD + UniRef100 | Various research pipelines |
| **ProtTrans (T5-XL)** | 3B | BFD + UniRef100 | Early PLM-based methods |

**Key insight from EpitopeTransfer (2025):** Fine-tuning PLMs on phylogenetically related epitope data before downstream classification significantly outperforms using frozen embeddings or general fine-tuning (AUC gains of +0.090.13 over BepiPred-3.0, EpiDope, EpitopeVec, and epitope1D).

### 3.2 Graph Neural Networks for Conformational Epitopes

GNNs exploit the fact that conformational B-cell epitopes form spatially clustered patches on the protein surface. Key architectures:

| Model | Graph Construction | GNN Type | Key Finding |
|-------|-------------------|----------|-------------|
| **GraphBepi (2023)** | AlphaFold2 structure  residue proximity graph | Graph Attention + BiLSTM | +5.5% AUC over sequence-only |
| **EpiGraph (2024)** | ESM-2 + ESM-IF1 embeddings  surface graph | GAT with residual connections | AUCPR 0.24 vs 0.20 BepiPred-3.0 |
| **GraphEPN (2025)** | VQ-VAE + surface meshing | Graph Transformer | Highest accuracy/MCC on multiple benchmarks |
| **UniPMT (2025)** | Heterogeneous graph (peptide+MHC+TCR) | GraphSAGE + deep matrix factorization | +15% PR-AUC over pMTnet |

**Key finding:** Graph homophily (the tendency of epitope residues to cluster) provides a strong inductive bias that GNNs exploit effectively.

### 3.3 Transformer-Based Architectures

Self-attention mechanisms capture long-range dependencies in both protein sequences and MHC-peptide interactions.

| Model | Architecture | Application | Performance |
|-------|-------------|-------------|-------------|
| **MUNIS (2025)** | ESM-2  BiLSTM with attention | MHC-I binding + immunogenicity | AUROC 0.980 |
| **BERTMHC (2021)** | Fine-tuned BERT | MHC-II binding | AUC 0.882 vs 0.877 NetMHCIIpan |
| **TransHLA (2024)** | Transformer + Residue CNN | HLA-I & II epitope detection | 91.95% AUC (class I), 88.14% (class II) |
| **STMHCpan (2025)** | Self-attention MHC pan-allelic | MHC-I binding | AUROC 0.961 (top in 17-predictor benchmark) |
| **DiscoTope 3.0 (2024)** | Inverse-folding transformer + CNN | Conformational B-cell | Best F1/MCC among structure-based |

### 3.4 Hybrid and Multi-Task Architectures

The strongest 20242025 models increasingly combine multiple paradigms:

- **MUNIS**  Transformer embeddings + recurrent processing  the synergy of PLM representations with sequential modeling of peptide-MHC interaction
- **GraphBepi / EpiGraph**  PLM features  graph construction  GNN classification  leveraging both evolutionary sequence signals and 3D proximity
- **UniPMT**  Heterogeneous GNN with multi-task learning  simultaneously predicts PMT, PM, and PT binding within one framework

### 3.5 Traditional Machine Learning (Still Relevant)

Despite the deep learning revolution, traditional ML retains niches:

- **SVM/Random Forest:** TCR-H (2024) uses SVM with physicochemical features, achieving strong generalization (AUC 0.870.89 on strict epitope splitting)
- **Gradient boosting:** epitope3D (2022) uses XGBoost with graph features from 3D structures
- **Ensemble methods:** IEDB consensus approach, BigMHC's 7-network ensemble

Traditional methods remain valuable when interpretability is paramount (e.g., SHAP analysis in TCR-H) or when labeled data is extremely scarce.

---

## 4. Tools & Models: Comparative Analysis

### 4.1 T-Cell Epitope Prediction Tools

#### MHC-I Binding Prediction

| Tool | Year | Architecture | Training Data | AUROC | Code Available |
|------|------|-------------|---------------|-------|---------------|
| **MUNIS** | 2025 | ESM-2 + BiLSTM | >650K peptide-HLA pairs | ~0.980 | Not yet public |
| **STMHCpan** | 2025 | Self-attention transformer | Eluted ligand MS data | 0.961 | Preprint only |
| **BigMHC** | 2023 | Ensemble of 7 DNNs + LSTM | MS eluted ligands + BA | 0.973 | Yes |
| **NetMHCpan 4.1** | 2020 | ANN (pan-specific) | IEDB binding affinities | ~0.925 | Yes (web + download) |
| **MHCflurry 2.0** | 2020 | ANN ensemble | IEDB + MS | ~0.920.93 | Yes (GitHub) |
| **TransHLA** | 2024 | Transformer + ResCNN | IEDB + CEDAR | 0.9195 | Yes (GitHub) |
| **DeepHLApan** | 2019 | RNN multi-task | IEDB | >0.90 | Yes |
| **CapsNet-MHC_AN** | 2023 | Capsule network | IEDB | Competitive | Yes |
| **MHCnuggets** | 2020 | LSTM | IEDB | ~0.82 | Yes (GitHub) |

**Critical finding (2025 bioRxiv benchmark):** Models trained on mass spectrometry-eluted ligand data *consistently outperform* those trained on binding affinity data alone. Training data quality matters more than architectural choices.

#### MHC-II Binding Prediction

| Tool | Year | AUROC | Notes |
|------|------|-------|-------|
| **MixMHC2pred** | 2022 | Best in 2024 benchmark | Also predicts peptide motifs |
| **NetMHCIIpan 4.0** | 2020 | ~0.877 | Pan-specific, widely used |
| **BERTMHC** | 2021 | ~0.882 | Transformer outperforms ANN slightly |
| **TransHLA** | 2024 | 0.8814 | Covers both class I and II |
| **DeepMHCII** | 2022 | Competitive | CNN-based |

MHC-II prediction consistently lags behind MHC-I because the open binding groove introduces greater structural complexity and peptide length variation.

### 4.2 B-Cell Epitope Prediction Tools

#### Linear B-Cell Epitope Prediction

| Tool | Year | Architecture | Performance |
|------|------|-------------|-------------|
| **DeepBCE** | 2023 | CNN | Acc 87.8%, AUC 0.945, F1 0.871 |
| **DeepLBCEPred** | 2023 | BiLSTM + multi-scale CNN + attention | Acc ~66%, MCC ~0.35 |
| **BepiPred-2.0** | 2019 | Random Forest (sequence features) | AUC ~0.580.62 |
| **EpitopeTransfer** | 2025 | Fine-tuned ESM-2 + RF | AUC +0.09 over BepiPred-3.0 |

#### Conformational B-Cell Epitope Prediction (Active Frontier)

| Tool | Year | Approach | AUC | AUCPR | Code |
|------|------|----------|-----|--------|------|
| **GraphEPN** | 2025 | VQ-VAE + Graph Transformer |  | Highest on multiple benchmarks |  |
| **EpiGraph** | 2024 | ESM-2/IF1 + GAT |  | 0.24 | Yes (GitHub) |
| **DiscoTope 3.0** | 2024 | Inverse-folding transformer + CNN |  | 0.232 | Yes (web server) |
| **GraphBepi** | 2023 | GNN + BiLSTM | +5.5% over SOTA |  | Yes |
| **SEMA 2.0** | 2024 | Ensemble (sequence + 3D) | 0.760.80 |  | Yes |
| **BepiPred-3.0** | 2022 | ESM-2 + ANN | ~0.700.71 | 0.177 | Yes (web server) |
| **epitope3D** | 2022 | XGBoost + graph features |  | Competitive | Yes |
| **EpiDope** | 2021 | BiLSTM | ~0.67 |  | Yes |

**Important caveat:** AUCPR values in the 0.200.25 range appear low, but they reflect the genuine difficulty of conformational epitope prediction on highly imbalanced benchmark sets where epitope residues are <10% of surface residues. The field trend is positive (0.177 BepiPred-3.0  0.24 EpiGraph  GraphEPN improvements).

### 4.3 Antibody-Specific Epitope Prediction

A distinct sub-problem: given a known antibody sequence/structure, predict where it binds on the antigen.

| Tool | Year | Key Innovation |
|------|------|---------------|
| **WALLE (NeurIPS 2024)** | 2024 | 310 improvement over general binding-site methods; PLM + bipartite GNN |
| **ISPIPab** | 2024 | Ensemble of 2 feature-based + docking + hierarchical clustering |
| **AlphaFold3** | 2024 | Closed-source but predicts antibody-antigen complexes with improved accuracy |
| **ImmuneFold** | 2025 | ESMFold + LoRA fine-tuning for TCRs/antibodies; zero-shot binding prediction |

### 4.4 Integrated Vaccine Design Frameworks

| Framework | Year | Capabilities |
|-----------|------|-------------|
| **VenusVaccine** | 2025 | Dual-attention integrating sequence + structural embeddings |
| **DeepVacPred** | 20212023 | End-to-end from epitope selection to multi-epitope vaccine construct design |
| **IntegralVac** | 2023 | Integrated immunoinformatics + DL framework |

---

## 5. Datasets and Benchmarking

### 5.1 Core Databases

| Database | Content | Size | URL |
|----------|---------|------|-----|
| **IEDB** | T-cell & B-cell epitopes | ~1.61M entries | iedb.org |
| **CEDAR** | Cancer-associated epitopes | ~224K entries | cedar.iedb.org |
| **BCiPep** | Linear B-cell epitopes | ~3,031 entries |  |
| **CEDb** | Conformational B-cell epitopes | ~225 entries |  |
| **AntiJen** | T-cell & B-cell epitopes | ~24,000 entries |  |
| **VDJdb** | TCR sequences with known epitope specificities | Growing | vdjdb.cdr3.net |
| **IEDB Benchmark** | Curated MHC binding benchmarks (BD2009, BD2013) | 135K176K data points | tools.iedb.org/benchmark |

### 5.2 The Data Leakage Crisis

A **2026 bioRxiv preprint by Preibisch et al.** revealed a systemic problem:

> **55.8% of assessable IEDB entries are labeled by computational models, not verified experimentally.**

This creates a vicious cycle:
1. NetMHCpan predicts which peptides bind to which allele in multi-allelic MS data
2. Those predicted labels enter IEDB as "experimental" data
3. The next generation of predictors trains on this contaminated data
4. **AUROC remains deceptively high (>0.89) while true discovery rate collapses** (Sensitivity@Top2% drops from 0.58 apparent to 0.23 true)

**Recommendation for manuscript:** Audit IEDB-derived training/test data using Preibisch et al.'s four-category framework (Clean / Biased / Multi-allelic / Insufficient metadata). Preference for mono-allelic, experimentally resolved data.

### 5.3 Best Practices for Benchmark Splitting

| Concern | Recommendation | Reference |
|---------|---------------|-----------|
| Peptide-level leakage | Split by protein/antigen source | Kim et al. (2014) |
| Near-identical proteins | BLAST 80% identity threshold | fcampelo/epitopes R package |
| Allele-level leakage | Hold out entire alleles for testing | DeepMHCflare approach |
| Confirmation bias | Use mono-allelic, experimentally resolved data only | Preibisch et al. (2026) |
| Metric insufficiency | Report AUROC + Precision@K + Sensitivity@Top% |  |

### 5.4 Existing Standardized Benchmarks

| Benchmark | Year | Scope | Key Finding |
|-----------|------|-------|-------------|
| **IEDB BD2009 / BD2013** | 2009/2013 | MHC-I binding | Foundation for most tool development |
| **17-predictor benchmark** | 2025 (bioRxiv) | MHC-I (44 alleles, >290K peptides) | STMHCpan #1 (0.961); eluted ligand training critical |
| **MHC-II benchmark (11 tools)** | 2024 (*Frontiers*) | MHC-II (20 allotypes) | MixMHC2pred + NetMHCIIpan-4.1 top |
| **AsEP (NeurIPS 2024)** | 2024 | Antibody-specific epitope prediction | Clustered epitope groups; WALLE 310 better |
| **epitope3D benchmark set** | 2022 | Conformational B-cell | Used by EpiGraph, GraphBepi, DiscoTope 3.0 |
| **DiscoTope 3.0 test set** | 2024 | Conformational B-cell | New standard for conformational epitope evaluation |

---

## 6. TCRpMHC Interaction Prediction

### 6.1 Why This Matters

MHC binding is necessary but **not sufficient** for T-cell immunogenicity. A peptide-MHC complex must also be recognized by a TCR to trigger an immune response. This is the critical missing link in current epitope prediction pipelines.

### 6.2 Key Reviews

- **Qi et al. (2025)**  *"A Roadmap for TCRpMHC Binding Prediction"*  *Briefings in Bioinformatics*  Most comprehensive survey; compared 22 pMHC + 30 TCR specificity methods
- **Zeng et al. (2025)**  *"Leveraging AI for Neoantigen Prediction"*  *Cancer Research*  Benchmarked antigen presentation + TCR binding predictors
- **Cheng et al. (2025)**  *"Progress of Deep Learning Prediction of CD8+ T-Cell Epitopes"*  *PROTEOMICS*  PLM-based encoding schemes

### 6.3 Major Methods

| Method | Year | Approach | Key Innovation | Performance |
|--------|------|----------|---------------|-------------|
| **UniPMT** | 2025 | Heterogeneous GNN + contrastive learning | Unified P-M-T, P-M, P-T prediction | +15% PR-AUC over pMTnet |
| **TCR-H** | 2024 | SVM + physicochemical features + SHAP | Strong generalization to unseen epitopes | AUC 0.870.89 |
| **TCRen** | 2024 | Structure-based | TCR recognition of unseen epitopes |  |
| **pMTnet** | 2021 | Siamese neural network | Widely used baseline |  |
| **PanPep** | 2023 | Memory neural network | Pan-peptide TCR recognition |  |

### 6.4 Critical Challenges

1. **Data imbalance:** >50% of training peptides from flu/SARS-CoV-2 studies (*Ly, Bonn & Prinz rebuttal, 2025*)
2. **Poor generalization:** Most models degrade severely under strict out-of-distribution evaluation (novel epitope variants)
3. **CDR3-only limitation:** Methods using full TCR sequence (+ chain) outperform CDR3-only but are data-limited
4. **Cross-reactivity:** TCRs predicted as neoepitope-specific may cross-react with viral epitopes (Zaghla 2025 thesis)

---

## 7. Generative Epitope Design

### 7.1 Emerging Paradigm

Instead of screening natural epitopes, generative models design *novel* peptides with desired properties (MHC binding, immunogenicity, stability). This is the frontier of computational epitope research.

### 7.2 Methods Landscape

| Approach | Example | Year | Key Contribution |
|----------|---------|------|-----------------|
| **Diffusion models** | RFdiffusion for pMHC-I | 2025 (ICML) | Structure-conditioned; exposes predictor blind spots |
| **MCTS + fine-tuning** | AlphaVacc | 2025 (bioRxiv) | Iterative arena competition; 11/12 candidates validated in vitro |
| **GAN** | Dual-projection WGAN (IBM) | 2023/2024 | Entropy-regularized Softmax for discrete AA sampling |
| **Conditional diffusion** | CDiffusion-AMP | 2024 | Outperforms GAN/VAE on AMP generation |

### 7.3 Critical Finding: Predictor Blind Spots

The 2025 ICML study (Mares et al.) demonstrated that **state-of-the-art sequence-based binding predictors largely fail to recognize structurally valid, novel peptides designed by RFdiffusion**  AUROCs dropped as low as 0.060.22. This exposes a fundamental limitation: current tools are good at recognizing natural epitope patterns but fail on peptides outside their training distribution.

### 7.4 Shift from GANs to Diffusion Models

The field is moving rapidly toward diffusion models for peptide design, mirroring trends in protein design more broadly. Reviewer critiques from eLife (2024) explicitly note that *"GANs are known for their training difficulty and mode collapse"* while recommending diffusion alternatives.

---

## 8. Clinical Translation & Validation

### 8.1 The False-Positive Problem

The gap between computational prediction and experimental validation remains stark:

- **TESLA Consortium (2020):** ~6% of top-predicted peptides confirmed as immunogenic (still the benchmark finding)
- **Zaghla et al. (2025):** 21.6% validation rate  78% false-positive rate
- **Cross-reactivity risk:** TCRs presumed neoepitope-specific may target CMV epitopes
- **HLA-LOH escape:** Neoepitopes frequently unpresented when the cognate HLA allele is lost (77% of cases)

### 8.2 Next-Generation Prediction Tools

| Tool | Year | Innovation | PPV |
|------|------|-----------|-----|
| **imNEO** | 2025 | 30 immunogenicity factors, 7 ML algorithms | 85% (top 10) |
| **neoIM** | 2025 | RF trained on MHC-presented non-self peptides | >30% improvement over existing tools |
| **KAIST B-cell neoantigen model** | 2024 | Dual prediction of T-cell + B-cell reactivity | FDA IND planned 2027 |

### 8.3 Current Clinical Trials (20242025)

| Trial | Platform | Key Outcome |
|-------|----------|-------------|
| BioNTech (mRNA) | Personalized mRNA | Pre-existing NeoAg clones expand post-vaccine |
| Geneos Therapeutics | DNA plasmid + IL-12 + pembrolizumab | 30.6% ORR in HCC |
| Nouscom | Adenovirus + MVA (up to 60 NeoAg) | ~700 SFU in NSCLC/melanoma |
| Gritstone GRANITE | Adenovirus + saRNA | Failed ctDNA primary endpoint |
| Genocea/ATLAS | Functionally validated NeoAg | Feasible in low-TMB; limited efficacy |

### 8.4 Validation Methodology Hierarchy

| Tier | Method | Cost/Time | Reliability |
|------|--------|-----------|-------------|
| 1 | Competitive peptide exchange assay | Low | Moderate |
| 2 | T2 binding assay (HLA-A*02:01) | Low | Moderate |
| 3 | ELISPOT / FluoroSpot (IFN-) | Medium | High (functional T-cell response) |
| 4 | IP-MS (immunopeptidomics) | High | High (direct presentation evidence) |
| 5 | TCR clonotype tracking (scTCRseq) | High | Highest (clonal expansion evidence) |

---

## 9. Knowledge Gaps & Open Problems

### 9.1 Critical Gaps

1. **Immunogenicity vs. Binding**: Current tools predict binding well but immunogenicity poorly. The key missing ingredients are: TCR recognition, antigen processing efficiency, and host immune context.

2. **Conformational B-cell epitopes**: Despite progress with GNNs and structural data, AUCPR values remain at 0.230.25. This is partly an intrinsic difficulty problem (epitope residues are <10% of surface) and partly a data problem (<250 confirmed conformational epitopes in CEDb).

3. **MHC-II prediction lags**: The open binding groove creates greater peptide length variation and binding mode flexibility than MHC-I, making accurate prediction harder.

4. **Data contamination**: The Preibisch et al. (2026) finding that 55.8% of IEDB is prediction-labeled represents a crisis for the field that has not yet been addressed by most tool developers.

5. **Generalization to unseen alleles and pathogens**: Most models overfit to common alleles (HLA-A*02:01) and well-studied pathogens (influenza, SARS-CoV-2). Performance on rare alleles and neglected pathogens is largely unknown.

6. **Non-canonical epitopes**: Spliced peptides, post-translationally modified peptides, and non-canonical reading frames remain poorly served by current tools.

7. **CD4+ T cell epitopes undervalued**: Most pipelines focus on CD8+ (MHC-I), but evidence increasingly shows CD4+ T cells are critical for durable antitumor immunity.

8. **Interpretability**: Most deep learning models are black boxes. SHAP/LIME analysis is possible (TCR-H demonstrates this) but not standard practice.

### 9.2 Opportunities for Novel Contribution

| Gap | Potential Contribution |
|-----|----------------------|
| Standardized, contamination-free benchmark | Curate a "Clean IEDB" benchmark using Preibisch et al. framework |
| Phylogeny-aware transfer learning | Apply EpitopeTransfer approach to neglected pathogens |
| Multi-task immunogenicity prediction | Combine MHC binding + TCR recognition + processing in one model |
| B-cell+T-cell unified prediction | Extend KAIST's 2024 approach; currently only one published model |
| Structure-guided generative design | Build on RFdiffusion approach with experimental validation loop |
| Interpretable epitope predictors | SHAP-based feature attribution as standard output |

### 9.3 Recommended Standard Practices (for manuscript methodology)

1. **Data splitting:** Split by protein source, not by peptide. Use BLAST at >80% identity to prevent leakage.
2. **Benchmark metrics:** Report AUROC, AUPR, Precision@K, and Sensitivity@Top%  not AUROC alone.
3. **External validation:** Test on at least one completely independent dataset with no sequence overlap to training.
4. **Data provenance audit:** Report the fraction of training/test data that is experimentally confirmed vs. computationally labeled.
5. **Negative data handling:** Specify how negative examples were defined (experimental non-binders vs. inferred negatives).
6. **Allele coverage:** Report performance disaggregated by HLA allele, not just aggregate.
7. **Model availability:** Release code + trained weights; provide web server for experimentalists.

---

## 10. Annotated Bibliography of Key Papers

### Foundation Paper — Direct Methodological Ancestor

**0. Jessen (2018)** — *"Deep Learning for Cancer Immunotherapy"*
- RStudio AI Blog, Published Jan. 29, 2018 | URL: `blogs.rstudio.com/ai/posts/2018-01-29-dl-for-cancer-immunotherapy`
- **Author:** Leon Eyrich Jessen (Technical University of Denmark)
- **Why read:** This is the **direct methodological ancestor of our project**. Jessen demonstrated that a simple feed-forward neural network (180→90→3, 49,143 params) could classify peptides into SB/WB/NB with ~95% accuracy using netMHCpan-4.0 labels and BLOSUM62 encoding. Our manuscript systematically extends this work with six architectures, ESM-2 embeddings, IEDB external validation, and structural docking.
- **Key contributions:**
  - **FFN model**: 180→dropout(0.4)→90→dropout(0.3)→3 (softmax), achieving 94.6% accuracy on 10% held-out test data (n=23,760 balanced: 7,920 per class)
  - **CNN model**: Conv2D(32, 3×3)→dropout(0.25)→flatten→FFN body, achieving 92% accuracy — notably *worse* than the simpler FFN, establishing that peptide "images" lack the spatial edge/continuity structure that CNNs exploit
  - **Random Forest**: 100 trees on flattened BLOSUM62 features, 82% accuracy — establishing baseline for traditional ML
  - **Data pipeline**: 1,000,000 random 9-mers → netMHCpan-4.0 labeling → downsampled to balanced 23,760 (SB/WB/NB) → 90/10 train/test split → BLOSUM62 encoding (9×20 matrix per peptide)
  - **Biological primer**: Comprehensive tutorial on adoptive T-cell therapy, pMHC biology, cancer immunotherapy rationale — the bridge between immunology and deep learning
  - **Reproducibility**: All R code released via RStudio blog + GitHub; Keras/TensorFlow R interface demonstrated
- **Relation to our project:**
  | Aspect | Jessen 2018 | Our Study |
  |--------|-----------|-----------|
  | Training data | 1M peptides, balanced to 23,760 | 100K peptides, balanced to 5,088 |
  | Label source | netMHCpan-4.0 | MHCflurry 2.2.0 (+ PSSM + random controls) |
  | Architectures | FFN, CNN, RF (3) | Deep FFN, FFN, CNN, LSTM, ResNet, RF (6) |
  | Encoding | BLOSUM62 only | BLOSUM62 + ESM-2 t6/t12/mean-pooled |
  | External validation | None | IEDB 49 epitopes + 20 controls |
  | Docking | None | HDOCKlite (5 peptides × 1DUZ) + MD |
  | Mutation analysis | None | p53 + KRAS hotspot scanning |
  | Key finding | Deep learning captures non-linear pMHC interactions | Label quality > architecture; ESM-2 offers modest gains |
- **Limitations noted in the original post**: (1) Labels are model-derived (netMHCpan), making the FFN "a model of a model"; (2) No external experimental validation; (3) Single HLA allele (A\*02:01); (4) CNN underperformance vs FFN not fully explained — speculated to reflect lack of spatial edge structure in peptide encoding
- **Code**: R scripts at `https://git.io/vb3Xa` (training data) + RStudio blog repository
- **Status**: [x] PDF read and distilled · [x] Cross-referenced with our manuscript · [x] Notes at `Literature/notes/jessen-2018-dl-cancer-immunotherapy.md`

---

### Core Review Papers (Start Here)

**1. Oluwagbemi et al. (2025)**  *"AI-driven epitope prediction: a systematic review, comparative analysis, and practical guide for vaccine development"*
- *npj Vaccines*, DOI: 10.1038/s41541-025-01258-y
- **Why read:** Most comprehensive 2025 review; benchmark table covering B-cell + T-cell tools; practical recommendations.
- **Key table:** Table 1 (all tools compared with architectures, metrics, availability).

**2. Qi et al. (2025)**  *"A roadmap for T cell receptor-peptide-bound major histocompatibility complex binding prediction by machine learning: glimpse and foresight"*
- *Briefings in Bioinformatics*, 26(4)
- **Why read:** Deepest survey of TCRpMHC prediction; benchmarks 22 pMHC + 30 TCR methods; proposes a structured roadmap.

**3. Hotop et al. (2024)**  *"Integrating machine learning to advance epitope mapping"*
- *Frontiers in Immunology*, DOI: 10.3389/fimmu.2024.1463931
- **Why read:** Bridges experimental epitope mapping (peptide arrays, cryo-EM, crystallography) with ML/DL methods; discusses feature engineering and dataset constraints.

**4. Cheng et al. (2025)**  *"Progress of Deep Learning Prediction of CD8+ T-Cell Epitopes"*
- *PROTEOMICS*, DOI: 10.1002/pmic.70101
- **Why read:** Focuses specifically on CD8+ epitopes; PLM-based encoding schemes and their comparative advantages.

### Key Method Papers

**5. Mares et al. (2025)**  *"Generation of structure-guided pMHC-I libraries using Diffusion Models"*
- ICML 2025, arXiv:2507.08902
- **Why read:** Demonstrates that current MHC binding predictors fail on structurally valid novel peptides (AUROCs 0.060.22); structure-guided generative design.

**6. Zhao et al. (2025)**  *"UniPMT: A unified deep learning framework for peptide-MHC-TCR binding prediction"*
- *Nature Machine Intelligence*
- **Why read:** State-of-the-art for joint pMHC-TCR prediction; heterogeneous GNN + contrastive learning; +15% PR-AUC over pMTnet.

**7. Preibisch et al. (2026)**  *"Resolution of recursive data corruption to transform T-cell epitope discovery"*
- *bioRxiv*, DOI: 10.1101/2026.03.30.710191
- **Why read:** Reveals 55.8% IEDB data contamination; proposes audit framework; **essential reading for anyone training epitope models**.

**8. Rajeshwar et al. (2024)**  *"TCR-H: Explainable machine learning prediction of TCRpMHC binding"*
- *Frontiers in Immunology*
- **Why read:** Best example of interpretable epitope prediction using SHAP; strong generalization on strict splitting.

### Structural/Conformational Epitope Papers

**9. Del Vecchio et al. (2024)**  *"EpiGraph: B cell epitope prediction by capturing spatial clustering property of the epitopes using graph attention network"*
- *Scientific Reports*
- **Why read:** State-of-the-art conformational epitope prediction; GAT with PLM embeddings; AUCPR 0.24.

**10. Risum et al. (2024)**  *"DiscoTope 3.0: Improved B-cell epitope prediction using AlphaFold2 modeling and inverse folding latent representations"*
- **Why read:** Inverse-folding transformer approach; most widely adopted recent method; strong benchmark.

### Benchmark & Dataset Papers

**11. Kim et al. (2014)**  *"Dataset size and composition impact the reliability of performance benchmarks for peptide-MHC binding predictions"*
- *BMC Bioinformatics*
- **Why read:** Classic paper on why cross-validation overestimates performance; still the go-to reference for benchmark methodology.

**12. bioRxiv benchmark (2025)**  *"Comprehensive Evaluation and Interpretative Insights of Peptide-HLA Binding Prediction Tools"*
- *bioRxiv*, DOI: 10.1101/2025.04.10.648169
- **Why read:** Benchmarked 17 MHC-I predictors across 44 alleles (>290K peptides); best head-to-head comparison available.

**13. Zaghla et al. (2025)**  *"Systematic evaluation of (neo)epitope predictions and their correlation with clinically observed T-cell responses"*
- FU Berlin thesis
- **Why read:** Reveals 78% false-positive rate in neoepitope prediction; cross-reactivity with CMV epitopes found.

### Transfer Learning & PLM Papers

**14. EpitopeTransfer (2025)**  *"Phylogeny-Aware Transfer Learning for Linear B-Cell Epitopes"*
- *bioRxiv*, DOI: 10.1101/2025.04.17.649425
- **Why read:** Demonstrates phylogeny-aware fine-tuning of ESM models; transfer learning blueprint for data-scarce pathogens.

**15. Shashkova et al. (2022)**  *"SEMA: Antigen B-cell conformational epitope prediction using deep transfer learning"*
- *Frontiers in Immunology*
- **Why read:** Early adoption of ESM-1v and ESM-IF1 for epitope prediction; contact number regression approach; AUC 0.76.

### Clinical Translation Papers

**16. neoIM (2025)**  *"Accelerating Neoantigen Discovery: A High-Throughput Approach to Immunogenic Target Identification"*
- *Vaccines*, 13(8), 865
- **Why read:** >30% improvement over existing tools on ELISpot benchmarks; shows utility as checkpoint inhibitor biomarker.

**17. imNEO (2025)**  *"High-accuracy neoantigen prediction through AI-based integrated analysis"*
- *J Clin Oncol*, ASCO 2025
- **Why read:** Integrates 30 immunogenicity factors; 85% PPV for top 10 ranked neoantigens.

### KRAS-Targeted Therapeutics

**18. Saetang et al. (2025)**  *"Identification and characterization of oncogenic KRAS G12V inhibitory peptides by phage display, molecular docking and molecular dynamic simulation"*
- *Computers in Biology and Medicine*, 192:110272 | PMID: 40300294
- **Why read:** First phage-display-derived peptide inhibitors selective for KRAS G12V over WT. Linear 23-mer peptides (Pep I, Pep II) achieve 70-75% CRC cell viability reduction (SW620) at 400μM with minimal WT toxicity. Pep II binding: G12V −35.96 vs WT −18.06 kcal/mol. Demonstrates subtractive biopanning strategy for mutation-selective peptide therapeutics. **Relevance to our project:** Validates KRAS G12V as a peptide-accessible target; complementary to our neoepitope approach (YKLVVVGAV for immunotherapy vs Pep I/II for direct inhibition). See `Literature/notes/saetang-2025-kras-g12v-peptide.md`.

---

## 11. Tool Availability & Code Index

### Web Servers (Ready to Use)

| Tool | URL | Epitope Type |
|------|-----|-------------|
| IEDB Analysis Resource | tools.iedb.org | T-cell (MHC-I, MHC-II) |
| NetMHCpan 4.1 | services.healthtech.dtu.dk | MHC-I binding |
| NetMHCIIpan 4.0 | services.healthtech.dtu.dk | MHC-II binding |
| BepiPred 3.0 | services.healthtech.dtu.dk | B-cell (linear + conformational) |
| DiscoTope 3.0 | services.healthtech.dtu.dk | B-cell conformational |
| MHCflurry | github.com/openvax/mhcflurry | MHC-I binding |
| EpiGraph | epigraph.kaist.ac.kr | B-cell conformational |

### GitHub Repositories (Open Source)

| Tool | Repository | Language |
|------|-----------|----------|
| MHCflurry | github.com/openvax/mhcflurry | Python |
| BigMHC | github.com/yourh/bigmhc | Python |
| TransHLA | github.com/SkywalkerLuke/TransHLA | Python |
| EpiGraph | github.com/sj584/EpiGraph | Python |
| DeepImmuno | github.com/yourh/deepimmuno | Python |
| WALLE / AsEP | github.com/biochunan/AsEP-dataset | Python |
| ISPIPab | github.com/mcarroll8/ISPIPab | Python |
| MHCnuggets | github.com/KarchinLab/mhcnuggets | Python |
| CapsNet-MHC | github.com/yourh/capsnet-mhc | Python |
| RFdiffusion pMHC-I | github.com/sermare/struct-mhc-dev | Python |

### Key Python Libraries for Epitope Research

```python
# Core immunoinformatics
pip install mhcflurry mhcnuggets

# Protein language models
pip install fair-esm  # ESM-2 from Meta

# Structure handling
pip install biopython pdb-tools

# Benchmark utilities
# R package: fcampelo/epitopes (train/test split with BLAST deduplication)
```

---

*This literature map was compiled on 2026-06-14 via systematic web search of PubMed, Semantic Scholar, bioRxiv, and Google Scholar, with deep reading of key papers. It should be updated before manuscript submission to capture new 2026 publications.*

### Search Methodology

- **Sources:** PubMed, Semantic Scholar, Google Scholar, bioRxiv, Nature, Frontiers, ScienceDirect
- **Coverage period:** Emphasis on 20232026, with key earlier papers included where foundational
- **Search terms:** "deep learning epitope prediction," "B-cell/T-cell epitope," "MHC binding prediction," "TCRpMHC," "conformational epitope GNN," "protein language model epitope," "neoantigen prediction," "generative epitope design"
- **Inclusion criteria:** Peer-reviewed or high-citation preprint; reports quantitative benchmark performance; architectural novelty; open-source availability; clinical validation
