---
type: literature-note
title: "Challenges in developing personalized neoantigen cancer vaccines"
created: 2026-07-03
tags: [neoantigen, cancer-vaccine, immunotherapy, epitope-prediction, immune-tolerance, review]
status: read
---

# Katsikis et al. (2023) — Challenges in Developing Personalized Neoantigen Cancer Vaccines

## Metadata
- **Title:** Challenges in developing personalized neoantigen cancer vaccines
- **Authors:** Katsikis PD, Ishii KJ, Schliehe C
- **Journal:** Nature Reviews Immunology, 2023, 23:745-760
- **DOI:** [10.1038/s41577-023-00937-y](https://doi.org/10.1038/s41577-023-00937-y)
- **Citation:** Katsikis et al., 2023
- **Type:** Review

## Summary
Comprehensive review of the major hurdles facing personalized neoantigen cancer vaccines, spanning the full pipeline from antigen prediction through delivery platforms to clinical translation. The authors systematically catalog failure points at each stage: neoantigen identification (both computational and mass-spectrometry-based), vaccine platform selection (mRNA, peptide, DC, viral vector), and biological barriers (immune tolerance, tumor heterogeneity, immunosuppressive microenvironment).

## Key Findings

### Antigen Prediction Bottlenecks
- Current MHC-binding prediction algorithms produce high false-positive rates, especially for MHC class II
- Only a fraction of predicted binders elicit T-cell responses in vivo (<1-5% in most studies)
- The "prediction-immunogenicity gap" remains the central unsolved problem
- Multi-omics integration (genomics + transcriptomics + immunopeptidomics) improves prediction but is resource-intensive

### Delivery Platform Trade-offs
| Platform | Advantages | Disadvantages |
|----------|-----------|---------------|
| mRNA/LNP | Rapid manufacturing, modular | Stability, cold-chain logistics |
| Synthetic long peptides | Proven safety, off-the-shelf components | HLA restriction, weak immunogenicity alone |
| Dendritic cell | Personalized cellular therapy | Cost, manufacturing complexity |
| Viral vector | Strong innate immune activation | Anti-vector immunity, limited payload |

### Biological Barriers
- Central and peripheral immune tolerance mechanisms limit neoantigen immunogenicity
- Tumor heterogeneity means single-target vaccines risk immune escape
- Immunosuppressive tumor microenvironment (Tregs, MDSCs, checkpoint ligands) dampens vaccine responses
- Clonal vs. subclonal neoantigen selection determines vaccine efficacy

### Clinical Translation Gaps
- Patient selection remains empirical (tumor mutational burden is a crude proxy)
- No validated biomarker to predict vaccine responders
- Manufacturing turnaround time (8-12 weeks) is too slow for many metastatic patients
- Combination strategies (vaccine + checkpoint inhibitor) show promise but add complexity

## Relevance to Our Manuscript

Cited in the **Introduction** to connect epitope prediction methodology to clinical vaccine development. Our manuscript's focus on quantifying the binding-immunogenicity gap directly addresses one of the three central challenges identified by Katsikis et al. The review provides authoritative context for why improved prediction methods matter translationally.

| Aspect | Katsikis 2023 | Our Study |
|--------|-------------|-----------|
| Scope | Clinical vaccine pipeline review | Computational prediction benchmarking |
| Gap identified | Prediction-immunogenicity gap (qualitative) | Prediction-immunogenicity gap (quantified: 78% false-positive rate per Zaghla 2025) |
| MHC class | Both I and II | HLA-A*02:01 (class I) |
| Validation type | Literature synthesis | Computational + IEDB external + structural |

## Key Quotes

> "The identification of immunogenic neoantigens remains the rate-limiting step in the development of personalized cancer vaccines, with current computational methods achieving only modest predictive accuracy."

> "Bridging the gap between MHC-binding prediction and T-cell recognition represents perhaps the single most important challenge for the field."

## Limitations
1. Review article — no primary experimental data
2. Field moving rapidly — some delivery platform data predates mRNA COVID-19 vaccine era
3. MHC class II prediction discussed but less actionable for our HLA class I-focused study

---

- [x] Abstract and key sections read
- [x] Integrated as citation in manuscript Introduction
- [ ] Full text TBD if needed for revision
