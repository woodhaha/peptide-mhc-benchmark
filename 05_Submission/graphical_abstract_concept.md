# Graphical Abstract Concept

**Manuscript**: "Systematic Benchmarking of Deep Learning Architectures for Peptide-MHC Binding Prediction Reveals Label Quality, Not Architecture, Defines Performance"

**Target Journal**: *Briefings in Bioinformatics*

**Design**: Three-panel horizontal layout, left-to-right narrative flow

**Software**: BioRender (preferred for biological accuracy of MHC structure) or Adobe Illustrator

**Color Palette**: Nature-style -- muted blue (#4A90D9), teal-green (#2E8B85), warm amber (#D4952A), with neutral grey (#6B7280) for supporting elements and white background.

---

## Panel 1 (Left): Architecture Benchmark -- "Architecture Matters Less Than You Think"

**Narrative**: Six deep learning architectures compete, but the simple one wins.

**Visual Elements**:
- Six stylised neural network icons arranged in a row, varying in complexity (LSTM: simple rectangle, FFN: 3-layer stacked, CNN: grid pattern, ResNet: skip-connection diagram, and two ESM-2 variants: small vs. large protein-language-model icons)
- Each architecture labelled with parameter count (e.g., LSTM 24K, FFN 152K, CNN 1.05M, ResNet 325K, ESM-2 t6 8M, ESM-2 t12 35M)
- A podium/ranking visual: Deep FFN on top (91.9%), with a surprise element -- a simple PSSM matrix icon floating above all six at 94.8%, visually breaking the "more complex = better" expectation
- Color coding: neural architectures in blue gradient; PSSM in amber (signaling "non-neural baseline wins")

**Key Annotation**: "PSSM 94.8% > Best Neural Network 91.9%"

---

## Panel 2 (Center): Label Quality Waterfall -- "The Real Performance Ceiling"

**Narrative**: Label quality, not architecture, defines how high you can go.

**Visual Elements**:
- A descending waterfall/bar chart with three bars:
  - **PSSM labels**: 94.8% (amber, tallest bar)
  - **MHCflurry labels**: 91.9% (blue, middle bar)
  - **Random labels**: 65.8% (grey, shortest bar)
- The gap between PSSM and MHCflurry (~2.9 pp) is shaded and labelled "Label quality gap"
- The gap between MHCflurry and Random (~26.1 pp) is shaded darker and labelled "Contamination ceiling"
- A magnified inset or callout box zooms into the **weak binder (WB) class**, showing F1 scores: PSSM 0.925 vs. MHCflurry 0.880 -- highlighting the WB class as the primary locus of label-driven variation
- A small annotation connects to "55.8% of IEDB entries computationally labelled (Preibisch et al., 2026)" as a footnote-style ribbon at the bottom

**Key Annotation**: "Label source accounts for ~29 pp variation -- exceeding maximum inter-architecture difference (~11 pp)"

---

## Panel 3 (Right): Translational Pipeline -- "From Benchmark to Bedside"

**Narrative**: The pipeline works end-to-end: scan proteins, find mutations, propose mechanism.

**Visual Elements**:
- **Step 1 (top)**: Epitope scanning schematic -- 10 protein icons (shapes of varying sizes) with peptide fragments (small coloured segments) being screened; annotation: "3,536 peptides scanned, 10/11 known epitopes identified"
- **Step 2 (middle)**: Cancer hotspot mutation analysis -- stylised p53 and KRAS protein ribbons with mutation sites marked (p53 R248W, KRAS G12V highlighted in amber); neoepitope candidate peptides displayed below as sequence tiles
- **Step 3 (bottom)**: Structural hypothesis visualisation -- a simplified molecular docking inset showing HLA-A*02:01 binding groove (blue ribbon/cartoon), the KRAS G12V peptide (YKLVVVGAV, amber stick), and a dashed line or electrostatic halo between K (P2 lysine) and Glu63 in the B-pocket, labelled "Predicted K-E63 salt bridge (~-30 kcal/mol)"
- A small "Experimental validation needed" badge or asterisk noting this is a docking-derived hypothesis
- Far right: a verification icon (checkmark with "93.9% IEDB sensitivity, AUC 0.947") and an open-source badge (GitHub icon)

**Key Annotation**: "93.9% sensitivity on 49 experimentally verified IEDB epitopes"

---

## Layout & Design Notes

- **Overall aspect ratio**: 16:4 (wide format, suitable for journal website banner)
- **Flow**: Left-to-right with subtle connecting arrows or gradient ribbon between panels
- **Typography**: Sans-serif (suggest Helvetica or Arial for journal compatibility); panel titles in bold 14pt, annotations in 9pt regular
- **Data precision**: Report all percentages to one decimal place; F1 scores to three decimal places
- **White space**: Generous margins between panels; avoid crowding
- **File format for submission**: Vector PDF (primary), with a high-resolution PNG (300 dpi, minimum 1200px wide) as fallback

## Suggested Implementation Workflow

1. **BioRender**: Create MHC structural elements (HLA-A*02:01 binding groove) and protein ribbon diagrams using BioRender's immunology library
2. **Adobe Illustrator** (or Inkscape as free alternative): Assemble the three-panel layout, add bar charts, network icons, annotations, and typography
3. **Export**: PDF (vector) for manuscript submission; PNG at 300 dpi for online preview

---

*Concept prepared July 2026. To be implemented by a scientific illustrator or the corresponding author using BioRender/Illustrator.*
