# Saetang et al. (2025) — KRAS G12V Inhibitory Peptides

## Metadata
- **Title:** Identification and characterization of oncogenic KRAS G12V inhibitory peptides by phage display, molecular docking and molecular dynamic simulation
- **Authors:** Jirakrit Saetang, Montarop Yamabhai, Kuntalee Rangnoi, Napat Prompat, Thaiyawat Haewphet, Surasak Sangkhathat, Varomyalin Tipmanee, Soottawat Benjakul
- **Journal:** *Computers in Biology and Medicine*, Vol. 192 (Pt A), 110272
- **Date:** June 2025 (Epub 2025-04-28)
- **DOI:** `10.1016/j.compbiomed.2025.110272`
- **PMID:** 40300294

## Key Findings

### Peptide Discovery
- **Method:** Phage display with subtractive bio-panning (linear + cyclic libraries, counter-selected against WT KRAS)
- **Result:** Two 23-mer linear peptides — **Pep I** and **Pep II**
- **Surprise:** Linear peptides outperformed cyclic peptides in selectivity

### Binding (MD Simulation)
| Peptide | KRAS G12V ΔG (kcal/mol) | KRAS WT ΔG (kcal/mol) | Selectivity |
|---------|:------------------------:|:----------------------:|:-----------:|
| Pep I | Binds functional regions | Binds non-functional regions | Moderate |
| **Pep II** | **−35.96** | −18.06 | **2×** |

### Cellular Activity
| Cell Line | KRAS Status | Viability Reduction |
|-----------|:-----------:|:-------------------:|
| SW620 (CRC) | G12V | **70–75%** @400μM, 48h |
| NCI-H2444 (lung) | G12V | ~70–75% |
| Caco-2 (CRC) | WT | 20–30% (DMSO-level) |

### Mechanism
- ERK1/2 phosphorylation ↓
- p21 upregulation
- G0/G1 cell cycle accumulation
- Reduced efficacy upon KRAS G12V knockdown → target specificity confirmed

## Relevance to Our Project

This paper validates KRAS G12V as a **peptide-accessible target** through a completely different mechanism from our neoepitope approach:

| | Saetang 2025 | Our YKLVVVGAV |
|---|---|---|
| **Role** | Direct KRAS inhibitor | T-cell epitope (pMHC) |
| **Length** | 23-mer | 9-mer |
| **Target** | KRAS G12V protein | HLA-A*02:01/peptide complex |
| **Mechanism** | Block KRAS signaling | TCR recognition → CTL killing |
| **CRC model** | SW620 | Same (G12V mutant) |

**Complementary:** Inhibitor peptide + neoepitope vaccine could simultaneously block KRAS signaling AND mark G12V cells for immune clearance.

**Key takeaway for our work:** The 23-mer peptides bind the KRAS protein surface. The K-E63 salt bridge mechanism we identified for YKLVVVGAV is about peptide-MHC anchoring — a different interface. Both demonstrate that G12V creates novel peptide-accessible surfaces, just at different molecular interfaces.

## Status
- [x] Read and distilled
- [x] Added to literature_matrix.md
- [x] Cross-referenced with docking/KRAS G12V analysis
