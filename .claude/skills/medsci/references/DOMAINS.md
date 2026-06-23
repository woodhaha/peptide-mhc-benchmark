# Domain Knowledge Registry

## Section A: Default Domain — Biomedical / Immunology / Epitope Prediction

### A.1 Terminology Standardization

| Unified Term | Definition | Avoid Using |
|-------------|------------|-------------|
| **Epitope** | The specific peptide fragment recognized by the immune system (antibody or TCR) | "antigenic determinant" (ambiguous), "target peptide" |
| **pMHC** | peptide-MHC complex; the ligand presented on cell surfaces | "MHC-peptide", "HLA-ligand" (less precise) |
| **MHC-I** | Major Histocompatibility Complex class I; presents 8-11mer peptides to CD8+ T cells | "HLA class I" (HLA is the human gene; MHC is the molecule) |
| **MHC-II** | MHC class II; presents 13-25mer peptides to CD4+ T cells | — |
| **HLA allele** | Specific MHC variant, e.g., HLA-A*02:01 | "MHC type", "HLA type" |
| **Binding affinity** | Strength of peptide-MHC binding, measured as IC50 (nM) or % rank | "binding score" (ambiguous) |
| **Immunogenicity** | Ability to trigger an immune response (T-cell activation); broader than binding | "immunoreactivity" |
| **Neoantigen / Neoepitope** | Tumor-specific peptide from somatic mutation | "cancer epitope" (imprecise), "mutant peptide" |
| **Pan-specific** | Model that generalizes to unseen HLA alleles | "pan-allele", "allele-independent" |
| **Eluted ligand (EL)** | Peptide physically extracted from MHC molecules, identified by MS | "MS ligand", "MHC ligandome" |
| **Binding affinity (BA)** | In vitro measured peptide-MHC binding strength | "IC50 data", "affinity data" |
| **BLOSUM62** | Blocks Substitution Matrix; 20×20 amino acid substitution scores | — |
| **ESM-2** | Evolutionary Scale Modeling; protein language model from Meta FAIR | "ESM", "ESM-1b" (different version) |
| **AUC / AUROC** | Area Under ROC Curve; standard binary classification metric | "AUC" alone (ambiguous: AUC-PR also exists) |
| **PPV** | Positive Predictive Value at given rank (e.g., PPV@top100) | "precision" (PPV is precision at a specific cutoff) |

### A.2 3-Axis Method Taxonomy

For classifying peptide-MHC binding / epitope prediction methods:

#### Axis 1: Algorithmic Priors (What kind of model?)
- **Sequence-based (shallow ML)**: SVM, Random Forest, early ANNs (pre-2020)
- **Deep feed-forward (FFN/CNN)**: MHCflurry 1.0/2.0, DeepImmuno, ConvNeXt-MHC
- **Recurrent (LSTM/RNN)**: DeepHLApan, earlier NetMHCpan variants
- **Transformer / Attention**: MUNIS (ESM-2 + BiLSTM), TransHLA, HLApollo, MHCRoBERTa
- **Graph Neural Networks (GNN)**: GraphMHC, GraphEPN, EpiGraph
- **Protein Language Model (PLM) transfer**: ESM-2-based, ProtBERT-based, ProtTrans-based
- **Ensemble / Hybrid**: NetMHCpan 4.1 (ANN ensemble), MHCflurry 2.0 (ensemble FNN), BigMHC (7 NN ensemble)

#### Axis 2: Data Priors (What data and how is it used?)
- **Binding affinity (BA) data**: Trained on IC50 measurements (IEDB BA data)
- **Eluted ligand (EL) data**: Trained on mass-spec identified peptides (IEDB EL data)
- **Multi-task (BA + EL)**: Joint training on both data types
- **Transfer learning chain**: BA → EL → stability → immunogenicity (Genesis/ImmugenX)
- **Self-supervised pre-training**: Protein language model pre-training on UniRef/UniProt
- **Negative set construction**: Random natural peptides, UniProt decoys, negative-set switching (HLApollo)

#### Axis 3: Validation Regime (How was it tested?)
- **Cross-validation (CV)**: k-fold within the same dataset
- **Cross-allele**: Test on HLA alleles not seen in training
- **Cross-species**: Test on non-human MHC (e.g., H-2 in mice)
- **External benchmark**: IEDB weekly benchmark, independent test set
- **Prospective experimental validation**: New epitopes tested in vitro (ELISpot, tetramer)
- **Clinical validation**: Correlated with patient outcomes, vaccine trial results

### A.3 Core Databases

| Database | URL | Description | Data Type |
|----------|-----|-------------|-----------|
| **IEDB** | https://www.iedb.org/ | Immune Epitope Database; gold standard for epitope data | BA, EL, T-cell, B-cell |
| **NetMHCpan** | https://services.healthtech.dtu.dk/services/NetMHCpan-4.1/ | Pan-specific MHC-I binding predictor; benchmark baseline | Predictions |
| **MHCflurry** | https://github.com/openvax/mhcflurry | Open-source MHC-I binding predictor with antigen processing | Predictions + training data |
| **SYFPEITHI** | http://www.syfpeithi.de/ | Classical MHC ligand database; motif-based scoring | Sequence motifs |
| **TCRdb** | http://bioinfo.life.hust.edu.cn/TCRdb/ | T-cell receptor sequence database | TCR sequences |
| **VDJdb** | https://vdjdb.cdr3.net/ | TCR-pMHC interaction database | TCR-pMHC pairs |
| **CEDAR** | https://cedar.iedb.org/ | Cancer Epitope Database and Analysis Resource | Cancer epitopes |
| **TumorAgDB** | (literature) | Tumor antigen database v2.0; 187K+ entries | Tumor antigens |
| **PDB** | https://www.rcsb.org/ | Protein Data Bank; pMHC structures | 3D structures |
| **AlphaFold DB** | https://alphafold.ebi.ac.uk/ | Predicted protein structures | 3D predictions |
| **UniProt** | https://www.uniprot.org/ | Universal protein database; negative set construction | Protein sequences |
| **dbSNP / COSMIC** | https://cancer.sanger.ac.uk/cosmic | Somatic mutation catalogs for neoantigen identification | Mutations |

### A.4 Key Journals & Venues

| Tier | Journals |
|------|----------|
| **Top-tier immunology** | Immunity, Nature Immunology, Science Immunology, Journal of Experimental Medicine |
| **Top-tier methods** | Nature Methods, Nature Biotechnology, Nature Machine Intelligence |
| **Bioinformatics** | Bioinformatics, Briefings in Bioinformatics, PLOS Computational Biology |
| **Immunology** | Journal of Immunology, Frontiers in Immunology, European Journal of Immunology |
| **Cancer** | Cancer Research, Cancer Immunology Research, Journal for ImmunoTherapy of Cancer |
| **Proteomics** | Molecular & Cellular Proteomics, PROTEOMICS, Journal of Proteome Research |
| **Preprints** | bioRxiv, Research Square |

### A.5 Standard Benchmarks

| Benchmark | Metric | Target | Notes |
|-----------|--------|--------|-------|
| IEDB weekly benchmark | AUROC, PPV | MHC-I binding | Gold standard; new data added weekly |
| IEDB MHC-II benchmark | AUROC | MHC-II binding | Harder than MHC-I |
| CEDAR neoepitope | PPV@top-N | Immunogenicity | Cancer-specific |
| VDJdb TCR-pMHC | AUROC, PPV | TCR recognition | Sparse; mostly well-studied epitopes |
| Wells et al. CD8+ | AP, AUROC | CD8+ T-cell response | Used by MUNIS, HLApollo |
| TESLA (NEJM 2020) | Clinical response | Neoantigen vaccine | 6 patients; limited but clinical gold |

### A.6 Canonical Papers (for Pattern 9 consistency)

When citing the following papers, use ONLY these verified descriptions:

| Paper | Canonical Description |
|-------|----------------------|
| **Jurtz et al. (2017) — NetMHCpan 4.0** | ANN-based pan-specific MHC-I binding predictor trained on BA + EL data; introduced NNAlign framework. J Immunol. |
| **Reynisson et al. (2020) — NetMHCpan 4.1** | Updated NetMHCpan with expanded training data and improved pan-specific performance. Nucleic Acids Res. |
| **O'Donnell et al. (2018) — MHCflurry 1.0** | Open-source FNN-based MHC-I binding predictor; allele-specific models. Cell Syst. |
| **O'Donnell et al. (2020) — MHCflurry 2.0** | Pan-allele ensemble FNN with antigen processing integration. Cell Syst. |
| **Wohlwend et al. (2025) — MUNIS** | ESM-2 + BiLSTM hybrid; CD8+ T-cell epitope immunogenicity; ROC-AUC ~0.98, AP ~26% better than next-best. Nat Mach Intell. |
| **Liu et al. (2024) — HLApollo** | Transformer with ESM-1b embeddings; negative-set switching; 78.7% AP for presentation. Nat Commun. |
| **Chu et al. (2024) — TransHLA** | Hybrid Transformer + Residual CNN with ESM-2; 91.95% AUC for HLA-I. GigaScience. |
| **Bravi et al. (2024) — Genesis/ImmugenX** | Modular encoder-only transformer; transfer learning chain BA→EL→stability→immunogenicity. PLOS Comp Biol. |

---

## Section B: Domain Template

Copy this template to adapt MedSci to a new biomedical sub-domain:

```markdown
### B.X Terminology Standardization
| Unified Term | Definition | Avoid Using |
|-------------|------------|-------------|
| <term> | <precise definition> | <imprecise alternatives> |

### B.X N-Axis Method Taxonomy
#### Axis 1: <Architectural/Algorithmic dimension>
- <family 1>
- <family 2>

#### Axis 2: <Data/Feature dimension>
- <family 1>
- <family 2>

#### Axis 3: <Validation/Evaluation dimension>
- <family 1>
- <family 2>

### B.X Core Databases
| Database | URL | Description | Data Type |

### B.X Key Journals
| Tier | Journals |

### B.X Standard Benchmarks
| Benchmark | Metric | Target | Notes |

### B.X Canonical Papers
| Paper | Canonical Description |
```

**Instructions**:
1. Fill in the domain name and replace X with the next available number
2. Define at least 10 terminology entries
3. Define 3 taxonomic axes appropriate for the domain
4. List at least 5 core databases
5. List at least 5 canonical papers with verified descriptions
6. Submit for human review before promoting to Section A

---

## Section C: Domain Registry (Grows via Distillation)

| Domain | Date Added | Source | Status | Fragments |
|--------|-----------|--------|--------|-----------|
| `epitope-prediction` | 2026-06-20 | Project bootstrap | **active** (Section A) | terminology:15, taxonomy:3, databases:12, journals:11, benchmarks:6, papers:8 |
| `deep-research-agents` | 2026-06-20 | QUEST (OSU-NLP-Group) | **extracted** (Section C) | terminology:9, taxonomy:1, benchmarks:8, data-source:3, method-family:5, key-paper:1, quality-rule:3 |
| `scientific-figure-generation` | 2026-06-20 | CCF-Figure (Deepshare-Official) | **extracted** (Section C) | terminology:18, method-family:4, quality-rule:5, writing-convention:6, paper-classification:7, diagram-types:11 |
| `neoantigen-prediction` | 2026-06-20 | neoantigen-predictor (AIPOCH) | **extracted** (Section C) | terminology:12, data-source:10, method-family:4, benchmark:2, quality-rule:4, key-paper:3, chemistry:2 |
| `immune-pathway-analysis` | 2026-06-20 | immune-pathway-analysis (AIPOCH) | **extracted** (Section C) | terminology:8, method-family:3, data-source:2, quality-rule:4, cli-params:30+ |
| `biomarker-landscape-mapping` | 2026-06-20 | biomarker-landscape-scanner (AIPOCH) | **extracted** (Section C) | terminology:12, method-family:1, quality-rule:17, maturity-framework:5, output-structure:10 |
| `adme-prediction` | 2026-06-20 | adme-property-predictor (AIPOCH) | **extracted** (Section C) | terminology:20, method-family:3, data-source:6, quality-rule:8, adme-endpoints:15+, druglikeness:3 |
| `paper-access` | 2026-06-20 | sci-hub-search-skill + paper-search-mcp | **integrated** (TOOL_STRATEGY.md) | terminology:5, data-source:23, tool-chain:4, fallback-strategy:1 |
| `immunopeptidomics` | — | — | planned | — |
| `tcr-pmhc` | — | — | planned | — |
| `neoantigen-design` | — | — | planned | — |
| `cancer-vaccine` | — | — | planned | — |

### Fragment Index

For distilled fragments that don't yet form a complete domain section, this table tracks them:

| ID | Domain | Type | Content (first 80 chars) | Source | Confidence | Date |
|----|--------|------|--------------------------|--------|-----------|------|
| `qst-001` | `deep-research-agents` | `terminology` | Deep Research Agent: multi-turn LLM system orchestrating search, scholar loo... | QUEST | extracted | 2026-06-20 |
| `qst-002` | `deep-research-agents` | `benchmark` | BrowseComp: fact-seeking benchmark; verifier-script-based objective evaluation | QUEST | extracted | 2026-06-20 |
| `qst-003` | `deep-research-agents` | `benchmark` | GAIA: fact-seeking benchmark; multi-step reasoning questions with ground-truth... | QUEST | extracted | 2026-06-20 |
| `qst-004` | `deep-research-agents` | `benchmark` | DeepResearch Bench + LiveResearchBench: citation grounding benchmarks; assess ... | QUEST | extracted | 2026-06-20 |
| `qst-005` | `deep-research-agents` | `method-family` | Verifiable Rubric-Tree Pipeline: synthetic task gen → rubric tree → auto-verif... | QUEST | extracted | 2026-06-20 |
| `qst-006` | `deep-research-agents` | `method-family` | Multi-turn Agent Rollout (search→scholar→visit→python→memory→synthesize): core... | QUEST | extracted | 2026-06-20 |
| `qst-007` | `deep-research-agents` | `method-family` | 4-Stage Training Pipeline: Vanilla → SFT → Mid-Training (MT) → RL (verifiable... | QUEST | extracted | 2026-06-20 |
| `qst-008` | `deep-research-agents` | `quality-rule` | Citation Grounding Verification: every factual claim in synthesized report mus... | QUEST | extracted | 2026-06-20 |
| `qst-009` | `deep-research-agents` | `quality-rule` | Verifier-Script Objective Reward: auto-generated Python verifier per objective... | QUEST | extracted | 2026-06-20 |
| `qst-010` | `deep-research-agents` | `quality-rule` | Rubric-Based Open-Ended Evaluation: LLM-as-judge comparing on comprehensivenes... | QUEST | extracted | 2026-06-20 |
| `qst-011` | `deep-research-agents` | `data-source` | QUEST Hugging Face collection: RL Data, SFT Objective Data, SFT Open-Ended Da... | QUEST | extracted | 2026-06-20 |
| `qst-012` | `deep-research-agents` | `data-source` | FAISS-indexed search and scholar caches for agent retrieval during rollouts | QUEST | extracted | 2026-06-20 |
| `qst-013` | `deep-research-agents` | `key-paper` | Xie J, Lin T, Wang Z et al. (2026). QUEST: Training Frontier Deep Research Age... | QUEST | verified | 2026-06-20 |
| `qst-014` | `deep-research-agents` | `data-source` | Serper API (web search) + Jina API (page reading/summarization) — external ser... | QUEST | extracted | 2026-06-20 |
| `ccf-001` | `scientific-figure-generation` | `terminology` | CCF Figure: intelligent assistant for auto-generating top-conference-level sci... | CCF-Figure | extracted | 2026-06-20 |
| `ccf-002` | `scientific-figure-generation` | `method-family` | Classification-Before-Generation: classify paper type first → select optimal d... | CCF-Figure | extracted | 2026-06-20 |
| `ccf-003` | `scientific-figure-generation` | `paper-classification` | Paper Type A: Method Paper → Method Overview / Model Architecture diagrams | CCF-Figure | extracted | 2026-06-20 |
| `ccf-004` | `scientific-figure-generation` | `paper-classification` | Paper Type B: Mechanism/Analysis Paper → Mechanism Zoom-in / Contrastive Diagr... | CCF-Figure | extracted | 2026-06-20 |
| `ccf-005` | `scientific-figure-generation` | `paper-classification` | Paper Type C: Benchmark/Evaluation → Evaluation Pipeline / Data Construction Pi... | CCF-Figure | extracted | 2026-06-20 |
| `ccf-006` | `scientific-figure-generation` | `paper-classification` | Paper Type D: Scaling Law/Empirical → Trend Summary Panel / Experimental Matri... | CCF-Figure | extracted | 2026-06-20 |
| `ccf-007` | `scientific-figure-generation` | `paper-classification` | Paper Type E: Robot/Embodied AI → Closed-Loop Feedback Diagram (Perception→Pla... | CCF-Figure | extracted | 2026-06-20 |
| `ccf-008` | `scientific-figure-generation` | `paper-classification` | Paper Type F: Interdisciplinary/AI for X → Cross-domain Application Framework | CCF-Figure | extracted | 2026-06-20 |
| `ccf-009` | `scientific-figure-generation` | `paper-classification` | Paper Type G: Survey/Position Paper → Taxonomy Tree / Timeline / Comparison Ma... | CCF-Figure | extracted | 2026-06-20 |
| `ccf-010` | `scientific-figure-generation` | `diagram-type` | Method Overview: high-level system workflow showing components and data flow; ... | CCF-Figure | extracted | 2026-06-20 |
| `ccf-011` | `scientific-figure-generation` | `diagram-type` | Model Architecture: detailed neural network structure with layers, dimensions,... | CCF-Figure | extracted | 2026-06-20 |
| `ccf-012` | `scientific-figure-generation` | `diagram-type` | Mechanism Zoom-in: local detail callout diagram magnifying one critical sub-co... | CCF-Figure | extracted | 2026-06-20 |
| `ccf-013` | `scientific-figure-generation` | `diagram-type` | Contrastive Diagram: side-by-side comparison of methods or before/after effect... | CCF-Figure | extracted | 2026-06-20 |
| `ccf-014` | `scientific-figure-generation` | `diagram-type` | Evaluation Pipeline: data flow from raw input → processing stages → metrics → ... | CCF-Figure | extracted | 2026-06-20 |
| `ccf-015` | `scientific-figure-generation` | `writing-convention` | Background Rule: pure white #FFFFFF or very light gray #F5F5F5 — no dark/gradi... | CCF-Figure | extracted | 2026-06-20 |
| `ccf-016` | `scientific-figure-generation` | `writing-convention` | Style Rule: 2D flat vector ONLY — explicitly prohibits 3D perspective, neon gr... | CCF-Figure | extracted | 2026-06-20 |
| `ccf-017` | `scientific-figure-generation` | `writing-convention` | Color Rule: functional coloring only; max 3 color groups; accent colors restri... | CCF-Figure | extracted | 2026-06-20 |
| `ccf-018` | `scientific-figure-generation` | `writing-convention` | Typography Rule: horizontal text only; max 2 font sizes; Chinese labels with E... | CCF-Figure | extracted | 2026-06-20 |
| `ccf-019` | `scientific-figure-generation` | `writing-convention` | Target Venues: NeurIPS, ICML, ICLR, CVPR, ACL, Nature Machine Intelligence, IE... | CCF-Figure | extracted | 2026-06-20 |
| `ccf-020` | `scientific-figure-generation` | `quality-rule` | 5-Failure-Mode Prevention: built-in self-check catching (1) wrong diagram type... | CCF-Figure | extracted | 2026-06-20 |
| `ccf-021` | `scientific-figure-generation` | `quality-rule` | 3-Round Iteration Protocol: max 3 revision rounds for figure refinement; each ... | CCF-Figure | extracted | 2026-06-20 |
| `ccf-022` | `scientific-figure-generation` | `quality-rule` | Constraint-Heavy Prompting: embed rigid visual rules directly in generation pro... | CCF-Figure | extracted | 2026-06-20 |
| `ccf-023` | `scientific-figure-generation` | `method-family` | Bilingual Prompt Template Library: Chinese+English copy-ready prompts for each... | CCF-Figure | extracted | 2026-06-20 |
| `ccf-024` | `scientific-figure-generation` | `data-source` | DeepShare article (June 2026): "先判断论文类型，再选择最合适的图示结构" — fo... | CCF-Figure | extracted | 2026-06-20 |
| `ccf-025` | `scientific-figure-generation` | `key-paper` | DeepShare-Official (2026). CCF-Figure: 为 AI / 计算机科学论文自动生成顶会级科研... | CCF-Figure | verified | 2026-06-20 |
| `ccf-026` | `scientific-figure-generation` | `diagram-type` | Closed-Loop Feedback Diagram: Perception→Planning→Execution→Feedback loop; spe... | CCF-Figure | extracted | 2026-06-20 |
| `neo-001` | `neoantigen-prediction` | `terminology` | Neoantigen: tumor-specific variant peptide from non-synonymous somatic mutatio... | neoantigen-predictor | extracted | 2026-06-20 |
| `neo-002` | `neoantigen-prediction` | `method-family` | 4-Stage Neoantigen Pipeline: (1) Mutant Peptide Generation (8-11mer from VCF) ... | neoantigen-predictor | extracted | 2026-06-20 |
| `neo-003` | `neoantigen-prediction` | `benchmark` | MHC Binding Thresholds: Strong binder = Rank<0.5% or IC50<50nM; Weak binder = ... | neoantigen-predictor | verified | 2026-06-20 |
| `neo-004` | `neoantigen-prediction` | `data-source` | NetMHCpan 4.1: pan-specific MHC-I binding predictor; gold-standard baseline; R... | neoantigen-predictor | verified | 2026-06-20 |
| `neo-005` | `neoantigen-prediction` | `data-source` | IEDB (Immune Epitope Database): validated HLA-peptide binding data; immunogeni... | neoantigen-predictor | verified | 2026-06-20 |
| `neo-006` | `neoantigen-prediction` | `key-paper` | Reynisson et al. (2020). NetMHCpan-4.1 and NetMHCIIpan-4.0: improved predictio... | neoantigen-predictor | verified | 2026-06-20 |
| `neo-007` | `neoantigen-prediction` | `key-paper` | Wells et al. (2022). Neoantigen prediction: perspectives on the present and fut... | neoantigen-predictor | verified | 2026-06-20 |
| `neo-008` | `neoantigen-prediction` | `method-family` | Immunogenicity Scoring (5-component weighted): Foreignness(0.30) + Anchor Mutat... | neoantigen-predictor | extracted | 2026-06-20 |
| `neo-009` | `neoantigen-prediction` | `method-family` | Priority Ranking (3-weight): MHC binding(0.40) + Immunogenicity(0.35) + Clinic... | neoantigen-predictor | extracted | 2026-06-20 |
| `neo-010` | `neoantigen-prediction` | `quality-rule` | Experimental Validation Required: MHC binding prediction accuracy ~85%; immunog... | neoantigen-predictor | verified | 2026-06-20 |
| `neo-011` | `neoantigen-prediction` | `quality-rule` | Clinical Research Only: prediction results should not be sole basis for clinica... | neoantigen-predictor | extracted | 2026-06-20 |
| `neo-012` | `neoantigen-prediction` | `data-source` | TCGA (The Cancer Genome Atlas): tumor mutation data; COSMIC: somatic mutation c... | neoantigen-predictor | extracted | 2026-06-20 |
| `neo-013` | `neoantigen-prediction` | `chemistry` | Kyte-Doolittle Hydrophobicity Scale: I(4.5), V(4.2), L(3.8), F(2.8), C(2.5), ... | neoantigen-predictor | verified | 2026-06-20 |
| `neo-014` | `neoantigen-prediction` | `chemistry` | Amino Acid Molecular Weights: G(75), A(89), S(105), P(115), V(117), T(119), I(... | neoantigen-predictor | verified | 2026-06-20 |
| `neo-015` | `neoantigen-prediction` | `data-source` | Allele Frequency Net Database (allelefrequencies.net): population HLA allele fr... | neoantigen-predictor | extracted | 2026-06-20 |
| `neo-016` | `neoantigen-prediction` | `data-source` | IMGT/HLA (ebi.ac.uk/ipd/imgt/hla): official HLA sequence database | neoantigen-predictor | extracted | 2026-06-20 |
| `neo-017` | `neoantigen-prediction` | `quality-rule` | VCF/FASTA Input Validation: verify CHROM/POS/REF/ALT columns; HLA nomenclature ... | neoantigen-predictor | extracted | 2026-06-20 |
| `imp-001` | `immune-pathway-analysis` | `terminology` | GSVA (Gene Set Variation Analysis): transforms gene-by-sample matrix into pathw... | immune-pathway-analysis | extracted | 2026-06-20 |
| `imp-002` | `immune-pathway-analysis` | `method-family` | GSVA/ssGSEA + limma Pipeline: (1) Read expression matrix + gene-set table → (2... | immune-pathway-analysis | extracted | 2026-06-20 |
| `imp-003` | `immune-pathway-analysis` | `quality-rule` | FDR Fallback Protocol: when no pathways pass adj.P.Val ≤ threshold, fall back t... | immune-pathway-analysis | extracted | 2026-06-20 |
| `imp-004` | `immune-pathway-analysis` | `quality-rule` | Append-Only Provenance: output_manifest.txt + run_record.txt append-only for re... | immune-pathway-analysis | extracted | 2026-06-20 |
| `imp-005` | `immune-pathway-analysis` | `data-source` | MSigDB Reactome: human Reactome pathways via msigdbr R package; filter immune-r... | immune-pathway-analysis | extracted | 2026-06-20 |
| `imp-006` | `immune-pathway-analysis` | `quality-rule` | Smoke-Test Validation: bundled minimal test data validates execution + fallback ... | immune-pathway-analysis | extracted | 2026-06-20 |
| `imp-007` | `immune-pathway-analysis` | `method-family` | Heatmap Visualization: diverging 3-color gradient; case_group first→control_gro... | immune-pathway-analysis | extracted | 2026-06-20 |
| `imp-008` | `immune-pathway-analysis` | `data-source` | R/Bioconductor: GSVA package, limma package, ComplexHeatmap — standard immunoin... | immune-pathway-analysis | extracted | 2026-06-20 |
| `bio-001` | `biomarker-landscape-mapping` | `terminology` | Biomarker Landscape Scan: field-level evidence map distinguishing use case, bio... | biomarker-landscape-scanner | extracted | 2026-06-20 |
| `bio-002` | `biomarker-landscape-mapping` | `maturity-framework` | 5-Tier Biomarker Maturity: Tier1(Exploratory signal, discovery-only) → Tier2(E... | biomarker-landscape-scanner | extracted | 2026-06-20 |
| `bio-003` | `biomarker-landscape-mapping` | `quality-rule` | 17 Hard Rules: never present exploratory association as maturity; always separa... | biomarker-landscape-scanner | extracted | 2026-06-20 |
| `bio-004` | `biomarker-landscape-mapping` | `quality-rule` | Validation ≠ Maturity: a biomarker may have external validation yet remain Tie... | biomarker-landscape-scanner | extracted | 2026-06-20 |
| `bio-005` | `biomarker-landscape-mapping` | `output-structure` | Mandatory 10-Section Output: A(Topic Framing)→B(Retrieval Audit)→C(Structured ... | biomarker-landscape-scanner | extracted | 2026-06-20 |
| `bio-006` | `biomarker-landscape-mapping` | `terminology` | Targeted vs Full Field Mode: narrow subdomain scan → Sections A,C(partial),D,H... | biomarker-landscape-scanner | extracted | 2026-06-20 |
| `adm-001` | `adme-prediction` | `terminology` | ADME: Absorption, Distribution, Metabolism, Excretion — 4 pillars of drug phar... | adme-property-predictor | extracted | 2026-06-20 |
| `adm-002` | `adme-prediction` | `adme-endpoints` | Absorption Panel: HIA(>80% good), Caco-2(>70 high permeability), Solubility(>0... | adme-property-predictor | extracted | 2026-06-20 |
| `adm-003` | `adme-prediction` | `adme-endpoints` | Distribution Panel: Vd(0.1-10 L/kg typical), PPB(>90% high binding), BBB(LogBB... | adme-property-predictor | extracted | 2026-06-20 |
| `adm-004` | `adme-prediction` | `adme-endpoints` | Metabolism Panel: CYP Inhibition(IC50, <1μM high DDI risk), CYP Substrate(clas... | adme-property-predictor | extracted | 2026-06-20 |
| `adm-005` | `adme-prediction` | `adme-endpoints` | Excretion Panel: CL(<5 low, 5-15 moderate, >15 high mL/min/kg), T1/2(2-8h typi... | adme-property-predictor | extracted | 2026-06-20 |
| `adm-006` | `adme-prediction` | `druglikeness` | Integrated Scores: QED(0-1, >0.6 good), Muegge Bioavailability(0-6, >4), MPO M... | adme-property-predictor | extracted | 2026-06-20 |
| `adm-007` | `adme-prediction` | `quality-rule` | CRITICAL: predictions are computational estimates only — do NOT replace experim... | adme-property-predictor | verified | 2026-06-20 |
| `adm-008` | `adme-prediction` | `quality-rule` | Model Accuracy: LogP R²=0.85-0.95; Solubility R²=0.65-0.80; HIA Acc=75-85%; BB... | adme-property-predictor | verified | 2026-06-20 |
| `adm-009` | `adme-prediction` | `method-family` | ADME Workflow Chain: Chemical Structure Converter → Lipinski Rule Filter → ADM... | adme-property-predictor | extracted | 2026-06-20 |
| `adm-010` | `adme-prediction` | `data-source` | Dependencies: Python 3.10+, RDKit (cheminformatics), dataclasses; models CPU-on... | adme-property-predictor | extracted | 2026-06-20 |

---

## Related Files

| File | Relationship |
|------|-------------|
| `DISTILLATION_PROTOCOL.md` | How fragments get added to Section C |
| `DISTILLATION_LOG.md` | Log of all distillation events |
| `HALLUCINATION_PATTERNS.md` | Domain-specific pattern detection (Patterns 10, 11) |
| `CITATION_INTEGRITY.md` | Verification rules reference canonical papers from A.6 |
| `PARADIGM.md` | Exemplar suggestions drawn from A.4 (Key Journals) |
