# Distillation Log

| Date | Source Skill | Fragments | Types | Conflicts | Status |
|------|-------------|-----------|-------|-----------|--------|
| 2026-06-20 | QUEST (OSU-NLP-Group) | 14 | terminology:9, benchmark:3, method-family:3, quality-rule:3, data-source:3, key-paper:1 | 0 | merged → DOMAINS.md §C `deep-research-agents` |
| 2026-06-20 | CCF-Figure (Deepshare-Official) | 26 | terminology:1, method-family:2, paper-classification:7, diagram-types:6, writing-convention:5, quality-rule:3, data-source:1, key-paper:1 | 0 | merged → DOMAINS.md §C `scientific-figure-generation` |
| 2026-06-20 | neoantigen-predictor (AIPOCH) | 17 | terminology:1, method-family:3, data-source:5, benchmark:1, key-paper:2, quality-rule:3, chemistry:2 | 0 | merged → DOMAINS.md §C `neoantigen-prediction` |
| 2026-06-20 | immune-pathway-analysis (AIPOCH) | 8 | terminology:1, method-family:2, data-source:2, quality-rule:3 | 0 | merged → DOMAINS.md §C `immune-pathway-analysis` |
| 2026-06-20 | biomarker-landscape-scanner (AIPOCH) | 6 | terminology:2, maturity-framework:1, quality-rule:2, output-structure:1 | 0 | merged → DOMAINS.md §C `biomarker-landscape-mapping` |
| 2026-06-20 | adme-property-predictor (AIPOCH) | 10 | terminology:1, adme-endpoints:4, druglikeness:1, quality-rule:2, method-family:1, data-source:1 | 0 | merged → DOMAINS.md §C `adme-prediction` |
| 2026-06-20 | sci-hub-search-skill + paper-search-mcp | — | integrated into TOOL_STRATEGY.md + SKILL_CATALOG.md | — | integrated → `paper-access` domain |

---

## 2026-06-20: QUEST Distillation Detail

### Source
- **Repository**: https://github.com/OSU-NLP-Group/QUEST
- **Paper**: Xie J, Lin T, Wang Z et al. (2026). "QUEST: Training Frontier Deep Research Agents with Fully Synthetic Tasks." arXiv:2605.24218.
- **License**: MIT
- **Accessed**: 2026-06-20

### Extraction Method
WebFetch on GitHub repo page + raw README.md. Direct source reading — no secondary summarization.

### Fragments Extracted (14 total)

| ID | Type | Content | Confidence |
|----|------|---------|-----------|
| qst-001 | terminology | Deep Research Agent definition | extracted |
| qst-002 | benchmark | BrowseComp | extracted |
| qst-003 | benchmark | GAIA | extracted |
| qst-004 | benchmark | DeepResearch Bench + LiveResearchBench | extracted |
| qst-005 | method-family | Verifiable Rubric-Tree Pipeline | extracted |
| qst-006 | method-family | Multi-turn Agent Rollout loop | extracted |
| qst-007 | method-family | 4-Stage Training Pipeline (SFT→MT→RL) | extracted |
| qst-008 | quality-rule | Citation Grounding Verification | extracted |
| qst-009 | quality-rule | Verifier-Script Objective Reward | extracted |
| qst-010 | quality-rule | Rubric-Based Open-Ended Evaluation | extracted |
| qst-011 | data-source | HuggingFace collection | extracted |
| qst-012 | data-source | FAISS-indexed caches | extracted |
| qst-013 | key-paper | QUEST paper (arXiv:2605.24218) | verified |
| qst-014 | data-source | Serper + Jina API services | extracted |

### Merge Results
- **Auto-merged**: 14 fragments → DOMAINS.md Section C
- **Conflicts**: 0 (new domain `deep-research-agents`, no existing fragments to conflict with)
- **Tentative**: 0

### Relevance to MedSci
QUEST's methodology is directly relevant to MedSci's Mode 1 (Autonomous Research). The verifiable rubric-tree pipeline concept can strengthen how MedSci self-verifies research claims. The citation grounding verification protocol aligns with CITATION_INTEGRITY.md Rule 4 (conclusion-direction verification). The multi-turn agent rollout loop informs potential future enhancement of the Mode 1 autonomous pipeline.

---

## 2026-06-20: CCF-Figure Distillation Detail

### Source
- **Repository**: https://github.com/Deepshare-Official/CCF-Figure
- **Organization**: DeepShare (深度之眼)
- **License**: MIT
- **Stars**: 33
- **Skill Name**: ccf-figure
- **Accessed**: 2026-06-20

### Extraction Method
WebFetch on GitHub repo page + raw README.md. Direct source reading.

### Fragments Extracted (26 total)

| ID | Type | Content | Confidence |
|----|------|---------|-----------|
| ccf-001 | terminology | CCF Figure definition | extracted |
| ccf-002 | method-family | Classification-Before-Generation methodology | extracted |
| ccf-003~009 | paper-classification | 7 paper types: A(Method)→G(Survey/Position) | extracted |
| ccf-010~014,026 | diagram-type | 6 diagram structures (Method Overview, Model Arch, Mechanism Zoom-in, Contrastive, Evaluation Pipeline, Closed-Loop Feedback) | extracted |
| ccf-015~019 | writing-convention | 5 visual spec rules: background, style, color, typography, target venues | extracted |
| ccf-020~022 | quality-rule | 5-Failure-Mode Prevention, 3-Round Iteration Protocol, Constraint-Heavy Prompting | extracted |
| ccf-023 | method-family | Bilingual Prompt Template Library | extracted |
| ccf-024 | data-source | DeepShare June 2026 article (foundational design philosophy) | extracted |
| ccf-025 | key-paper | CCF-Figure repository (MIT, 2026) | verified |

### Merge Results
- **Auto-merged**: 26 fragments → DOMAINS.md Section C
- **Conflicts**: 0 (new domain `scientific-figure-generation`, no existing fragments)
- **Tentative**: 0

### Relevance to MedSci
CCF-Figure's classification-before-generation methodology bridges manuscript writing (MedSci Mode 3) and figure preparation. The 7 paper-type classification system can help MedSci identify which diagram structures best communicate epitope prediction research. The strict visual specification (flat 2D, max 3 colors, horizontal typography, white background) serves as a quality gate against AI-generated figures that would be rejected by journal editors. The 5 failure modes complement HALLUCINATION_PATTERNS.md — they're the visual equivalent of text hallucination patterns.

---

## 2026-06-20: neoantigen-predictor Distillation Detail

### Source
- **Skill**: neoantigen-predictor (AIPOCH)
- **Source repo**: https://github.com/aipoch/medical-research-skills
- **License**: MIT
- **Accessed**: 2026-06-20

### Fragments Extracted (17 total)

| ID | Type | Content | Confidence |
|----|------|---------|-----------|
| neo-001 | terminology | Neoantigen definition | extracted |
| neo-002 | method-family | 4-Stage Neoantigen Pipeline | extracted |
| neo-003 | benchmark | MHC Binding Thresholds (Strong/Weak/Non-binder) | verified |
| neo-004 | data-source | NetMHCpan 4.1 reference | verified |
| neo-005 | data-source | IEDB reference | verified |
| neo-006 | key-paper | Reynisson et al. 2020 | verified |
| neo-007 | key-paper | Wells et al. 2022 | verified |
| neo-008 | method-family | 5-Component Immunogenicity Scoring | extracted |
| neo-009 | method-family | 3-Weight Priority Ranking | extracted |
| neo-010 | quality-rule | Experimental Validation Required | verified |
| neo-011 | quality-rule | Clinical Research Only | extracted |
| neo-012 | data-source | TCGA, COSMIC | extracted |
| neo-013 | chemistry | Kyte-Doolittle Hydrophobicity Scale (20 AA) | verified |
| neo-014 | chemistry | Amino Acid Molecular Weights (20 AA) | verified |
| neo-015 | data-source | Allele Frequency Net Database | extracted |
| neo-016 | data-source | IMGT/HLA database | extracted |
| neo-017 | quality-rule | VCF/FASTA Input Validation | extracted |

### Key Knowledge Gems
- **MHC Binding Quantification**: Rank% <0.5 = Strong; 0.5-2% = Weak; >2% = Non-binder — the standard for neoantigen calling
- **5-Component Immunogenicity Formula**: Foreignness(30%) + AnchorMutation(25%) + SelfSimilarity(20%) + Hydrophobicity(15%) + Clonality(10%)
- **3-Weight Priority Formula**: MHC(40%) + Immunogenicity(35%) + Clinical(25%)
- **Kyte-Doolittle Scale**: Complete 20 amino acid hydrophobicity reference table
- **Clinical Caveat**: MHC binding ~85% accurate; immunogenicity correlation only ~60-70%

---

## 2026-06-20: immune-pathway-analysis Distillation Detail

### Source
- **Skill**: immune-pathway-analysis (AIPOCH)
- **Source repo**: https://github.com/aipoch/medical-research-skills
- **License**: MIT
- **Accessed**: 2026-06-20

### Fragments Extracted (8 total)

| ID | Type | Content | Confidence |
|----|------|---------|-----------|
| imp-001 | terminology | GSVA definition | extracted |
| imp-002 | method-family | GSVA/ssGSEA + limma Pipeline | extracted |
| imp-003 | quality-rule | FDR Fallback Protocol | extracted |
| imp-004 | quality-rule | Append-Only Provenance | extracted |
| imp-005 | data-source | MSigDB Reactome via msigdbr | extracted |
| imp-006 | quality-rule | Smoke-Test Validation | extracted |
| imp-007 | method-family | Heatmap Visualization | extracted |
| imp-008 | data-source | R/Bioconductor stack (GSVA, limma, ComplexHeatmap) | extracted |

### Key Knowledge Gems
- **GSVA Algorithm**: Enrichment of empirical expression distribution of genes in pathway S relative to genes not in S
- **limma Contrast**: case_group - control_group; positive logFC = more active in case
- **FDR Fallback**: When no pathways pass FDR threshold, rank by |t| statistic instead
- **Provenance Design**: append-only manifest + run_record for reproducibility
