# Cancer Hotspot Mutation Scan -- Neoantigen Epitope Prediction

> **Model**: Deep FFN (MHCflurry-trained, 91.9% accuracy, 89.6% CV)  
> **Allele**: HLA-A*02:01  
> **Scope**: 16 hotspot mutations across p53 (7) and KRAS (9)  
> **Date**: June 15, 2026

---

## Summary

| Category | Count | Clinical Impact |
|----------|:-----:|-----------------|
| Neoepitopes CREATED | 2 | New immunotherapeutic targets |
| Epitopes ENHANCED | 2 | Stronger existing binders |
| Epitopes DESTROYED | 3 | Immune evasion mechanism |
| No epitope change | 9 | Silent for MHC-I binding |
| **Total** | **16** | |

---

## Neoepitopes Created (NB -> SB/WB)

These mutations create entirely new HLA-A*02:01 epitopes -- candidates for personalized cancer vaccines.

### 1. p53 R248W -- CREATED neoepitope

| Property | Value |
|----------|-------|
| **Mutation** | R > W at position 248 |
| **Clinical context** | Most frequent p53 mutation across all cancers |
| **WT peptide** | MNRRPILTI (NB, score 0.002) |
| **Mutant peptide** | **MNWRPILTI** (WB, score 0.409) |
| **Effect** | Non-binder -> Weak Binder (+0.407) |
| **p2 anchor** | N (weak) |
| **p9 anchor** | I (canonical) |
| **Mechanism** | R>W at p3 introduces bulky hydrophobic residue, improving MHC groove fit |

### 2. KRAS G12V -- CREATED neoepitope

| Property | Value |
|----------|-------|
| **Mutation** | G > V at position 12 |
| **Clinical context** | Pancreatic and lung adenocarcinoma driver |
| **WT peptide** | YKLVVVGAG (NB, score 0.000) |
| **Mutant peptide** | **YKLVVVGAV** (WB, score 0.475) |
| **Effect** | Non-binder -> Weak Binder (+0.475) |
| **p2 anchor** | K (non-canonical) |
| **p9 anchor** | G -> **V** (becomes canonical!) |
| **Mechanism** | G>V at p9 replaces glycine with valine -- perfect C-terminal anchor |

---

## Epitopes Enhanced (WB -> SB)

These mutations substantially boost existing weak binders to strong binders.

### 3. KRAS G12V -- ENHANCED epitope

| Property | Value |
|----------|-------|
| **WT peptide** | LVVVGAGGV (WB, score 0.509) |
| **Mutant peptide** | **LVVVGAVGV** (SB, score 0.816) |
| **Effect** | WB -> SB (+0.307) |
| **p9 anchor** | G -> **V** (canonical anchor created) |

### 4. KRAS G12C -- ENHANCED epitope

| Property | Value |
|----------|-------|
| **Mutation** | G > C at position 12 |
| **Clinical context** | Lung adenocarcinoma -- druggable (sotorasib, adagrasib) |
| **WT peptide** | LVVVGAGGV (WB, score 0.509) |
| **Mutant peptide** | **LVVVGACGV** (SB, score 0.816) |
| **Effect** | WB -> SB (+0.307) |
| **Clinical implication** | Enhanced SB may synergize with KRAS G12C inhibitors for combination immunotherapy |

---

## Epitopes Destroyed (SB/WB -> NB)

These mutations eliminate existing epitopes -- potential immune evasion mechanisms.

### 5. KRAS G13D -- DESTROYED epitope

| Property | Value |
|----------|-------|
| **Mutation** | G > D at position 13 |
| **Clinical context** | Colon cancer |
| **WT peptide** | GVGKSALTI (WB, score 0.471) |
| **Mutant peptide** | DVDKSALTI (NB, score 0.001) |
| **Effect** | WB -> NB (-0.470) |
| **Mechanism** | G>D at p2 introduces negatively charged aspartate -- destroys p2 anchor |

### 6. KRAS A146T -- DESTROYED epitope

| Property | Value |
|----------|-------|
| **Mutation** | A > T at position 146 |
| **Clinical context** | Colon cancer |
| **WT peptide** | GIPFIETSA (WB, score 0.436) |
| **Mutant peptide** | GIPFIETST (NB, score 0.020) |
| **Effect** | WB -> NB (-0.416) |
| **Mechanism** | A>T at p9 -- threonine is a weaker anchor than alanine |

### 7. p53 R249S -- DESTROYED epitope

| Property | Value |
|----------|-------|
| **Mutation** | R > S at position 249 |
| **Clinical context** | Aflatoxin-associated hepatocellular carcinoma |
| **WT peptide** | GMNRRPILT (WB, score 0.401) |
| **Mutant peptide** | GMNRSPILT (NB, score 0.245) |
| **Effect** | WB -> NB (-0.156) |
| **Mechanism** | R>S at p5 -- polar to small polar change affects TCR-facing residue |

---

## Silent Mutations (No Epitope Change)

| Mutation | Reason |
|----------|--------|
| p53 R175H | Mutation outside 9-mer windows with binding potential |
| p53 G245S | Mutation in non-anchor position of non-binding windows |
| p53 R273H | DNA contact mutant -- no A*02:01 epitope overlap |
| p53 R282W | LFS hotspot -- no A*02:01 epitope overlap |
| p53 Y220C | Structural cavity mutant -- no epitope overlap |
| KRAS G12D | G>D at p9 destroys anchor (like G13D but different window) |
| KRAS G12R | G>R at p9 introduces charged residue -- destroys anchor |
| KRAS Q61H/L/R | Q61 mutations in non-binding windows |

---

## Clinical Recommendations

### Tier 1 -- Validate Immediately

| Rank | Neoepitope | Mutation | Rationale |
|:----:|------------|----------|-----------|
| 1 | **MNWRPILTI** | p53 R248W | Most common p53 mutant; pan-cancer target |
| 2 | **LVVVGAVGV** | KRAS G12V | Pancreatic/lung cancer; enhanced to SB |
| 3 | **LVVVGACGV** | KRAS G12C | Druggable; combination therapy potential |

### Tier 2 -- Mechanistic Investigation

| Rank | Finding | Mutation | Question |
|:----:|---------|----------|----------|
| 4 | GVGKSALTI destroyed | KRAS G13D | Does epitope loss contribute to immune evasion in MSI-H colon cancer? |
| 5 | GIPFIETSA destroyed | KRAS A146T | Same question for A146T colon cancers |

### Validation Protocol

1. **In silico**: Confirm with netMHCpan-4.1 and MHCflurry consensus
2. **In vitro**: T2 binding assay (HLA-A*02:01 stabilization)
3. **Ex vivo**: IFN-gamma ELISpot with patient PBMCs
4. **In vivo**: HLA-A*02:01 transgenic mouse immunization

---

## Methods

- **Model**: Deep FFN (360-180-90-45-3), MHCflurry 2.2.0 trained, 91.9% accuracy
- **Scan method**: 9 sliding 9-mer windows spanning each mutation position
- **Comparison**: WT vs mutant binding score (SB*1.0 + WB*0.5)
- **Significance threshold**: delta score > 0.1 AND class change required

---

## Output Files

| File | Contents |
|------|----------|
| `data/mutation_scan_results.csv` | Full 144-row results with all windows |
| `mutation_scan_report.md` | This report |

---

*Generated June 15, 2026 | Peptide Epitope Research Project*
