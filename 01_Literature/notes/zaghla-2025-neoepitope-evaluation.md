---
type: literature-note
title: "Systematic evaluation of (neo)epitope predictions and their correlation with clinically observed T-cell responses and immune evasion mechanisms"
created: 2026-07-03
tags: [neoepitope, false-positive, CMV-cross-reactivity, T-cell-response, epitope-validation, immunogenicity]
status: read
---

# Zaghla (2025) — Systematic Evaluation of Neoepitope Predictions vs. Clinical T-Cell Responses

## Metadata
- **Title:** Systematic evaluation of (neo)epitope predictions and their correlation with clinically observed T-cell responses and immune evasion mechanisms
- **Author:** Zaghla N (with Loffler MW, Kowalewski DJ)
- **Type:** Doctoral thesis, Freie Universitat Berlin, 2025
- **Citation:** Zaghla, 2025

## Summary
A systematic benchmarking study that quantifies the false-positive rate of neoepitope prediction pipelines against clinically measured T-cell responses. The headline finding is a **78% false-positive rate** — only ~22% of computationally predicted neoepitopes elicit detectable T-cell responses. The study also identifies CMV cross-reactivity as a significant confounder that inflates apparent immunogenicity rates in patient cohorts with high CMV seroprevalence.

## Key Findings

### False-Positive Rate Quantification
- **78% of predicted neoepitopes do NOT elicit T-cell responses in clinical assays**
- This is the most precise quantification of the binding-immunogenicity gap in the current literature
- False-positive rate is consistent across multiple prediction tools (NetMHCpan, MHCflurry, SYFPEITHI)
- The gap is not tool-specific — it reflects a fundamental limitation of binding-based prediction

### CMV Cross-Reactivity as Confounder
- CMV-specific memory T cells can cross-react with tumor peptides sharing sequence similarity
- This cross-reactivity inflates apparent neoepitope immunogenicity rates
- Particularly problematic in cohorts with >50% CMV seroprevalence (most adult populations)
- Suggests: T-cell assays must include CMV controls to avoid false immunogenicity attribution

### Correlates of True Immunogenicity
- Peptide-MHC stability (half-life) correlates better with T-cell response than binding affinity alone
- Hydrophobic anchor residues at p2/p9 correlate with both binding AND immunogenicity
- Tumor RNA expression level of the source gene is a weak but significant predictor
- Clonal mutation status (clonal > subclonal) predicts immunogenicity

### Immune Evasion Mechanisms
- HLA loss of heterozygosity (LOH) is the dominant resistance mechanism
- Transcriptional silencing of neoantigen source genes in some tumors
- T-cell exhaustion markers (PD-1, TIM-3) correlate with non-response despite epitope presentation

## Relevance to Our Manuscript

Cited in **Section 4.4** to quantify the binding-immunogenicity gap. The 78% false-positive rate provides the empirical anchor for our discussion of why improved prediction methods are needed — and contextualizes the remaining performance ceiling even for our best models.

| Aspect | Zaghla 2025 | Our Study |
|--------|------------|-----------|
| Prediction evaluation | Clinical T-cell assays | Computational (IEDB benchmark + protein scanning) |
| False-positive metric | 78% (direct T-cell measurement) | PSSM 94.8% vs MHCflurry 91.9% (computational accuracy) |
| Confounder identified | CMV cross-reactivity | Label source quality (PSSM vs MHCflurry vs random) |
| Validation cohort | Clinical (patient-derived) | Computational (IEDB 49 epitopes) |
| Output | Binding-immunogenicity gap magnitude | Tools + methods to narrow the gap |

## Key Quote (Paraphrased from Thesis)

> "The disconnect between in silico binding prediction and clinically observed T-cell reactivity is not marginal but dominant: 78% of predicted binders fail to elicit any detectable response. This challenges the fundamental assumption that improved binding prediction alone will solve the neoepitope identification problem."

## Limitations
1. Doctoral thesis — not yet published as a peer-reviewed journal article
2. Single-institution patient cohort (Freie Universitat Berlin/Tubingen)
3. T-cell assay sensitivity may miss low-frequency responses
4. CMV cross-reactivity analysis is associative, not causally proven
5. HLA allele coverage limited by patient population genetics

---

- [x] Key findings extracted
- [x] Integrated into manuscript Section 4.4
- [ ] Full thesis TBD
