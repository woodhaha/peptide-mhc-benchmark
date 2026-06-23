# HLA-A*02:01 Peptide-MHC Docking — Neoepitope Validation

> Generated: 2026-06-16 · Proteindocker Protocol · Phase A

## Project Background

From the peptide-MHC deep learning pipeline (Deep FFN 91.9% accuracy, MHCflurry 2.2.0):
- **p53 R248W** and **KRAS G12V** produce predicted neoepitopes
- Goal: Validate these predictions with structure-based peptide-MHC docking

## Key Neoepitope Candidates

| Priority | Peptide | Source | Effect | P2 | P2 Canonical? | Key Question |
|:---:|--------|--------|--------|:---:|:---:|------|
| **P1** | `YKLVVVGAV` | KRAS G12V | CREATED (NB→WB, Δ+0.475) | K | ❌ NON-CANONICAL | Can K form salt bridge with E63? |
| **P2** | `MNWRPILTI` | p53 R248W | CREATED (NB→WB, Δ+0.407) | M | ✅ Canonical | Does W-at-P7 improve groove burial? |
| **P3** | `LVVVGAVGV` | KRAS G12V | ENHANCED (WB→SB, Δ+0.307) | V | ✅ Canonical | Best anchor profile — expected strong binding |

## Key Structural Finding

**KRAS G12V neoepitope has K at P2** — a positively charged Lys in the normally hydrophobic HLA-A*02:01 B-pocket. However, **GLU63 is 2.9Å from P2** in the template structure (1DUZ). This creates a potential K-E63 salt bridge, which would represent a **non-canonical anchor compensation mechanism**.

## Files Prepared

```
docking/
├── 1DUZ.pdb                          # Full template (HLA-A2 + B2M + peptide)
├── receptor_1DUZ.pdb                 # Receptor only (chains A+B)
├── peptide_p53_R248W_neo.pdb         # MNWRPILTI
├── peptide_p53_R248W_wt.pdb          # MNRRPILTI
├── peptide_KRAS_G12V_neo.pdb         # YKLVVVGAV
├── peptide_KRAS_G12V_wt.pdb          # YKLVVVGAG
├── peptide_KRAS_G12V_enh.pdb         # LVVVGAVGV
├── peptide_KRAS_G12V_enh_wt.pdb      # LVVVGAGGV
├── peptide_KRAS_G12C_enh.pdb         # LVVVGACGV
├── peptide_KRAS_G12C_enh_wt.pdb      # LVVVGAGGV
├── analyze_neoepitopes.py            # Anchor analysis script
└── DOCKING_README.md                 # This file
```

## Template Structure: 1DUZ

| Property | Value |
|----------|-------|
| PDB ID | 1DUZ |
| Resolution | 1.80 Å (X-ray) |
| Allele | HLA-A*02:01 |
| Bound peptide | LLFGYPVYV (HTLV-1 Tax 11-19) |
| Chains | A=heavy, B=β2m, C=peptide |
| P2 anchor contacts | GLU63, LYS66, TYR99, MET45, TYR7, VAL67, TYR159, PHE9 |
| P9 anchor contacts | THR143, LYS146, TYR84, ASP77, THR80, LEU81, TYR116 |

## Docking Protocol (HDOCK Web Server)

### Step 1: Access
Go to http://hdock.phys.hust.edu.cn/

### Step 2: Submit each pair
- **Receptor**: Upload `receptor_1DUZ.pdb`
- **Ligand**: Upload each `peptide_*.pdb`
- Wait for results (~15-30 min per job)

### Step 3: Analyze
Compare docking scores and poses for each mutant/WT pair:
- Mutant docking score > WT? → neoepitope structurally validated
- Mutant pose similar to WT? → same binding mode, check anchor contacts
- Mutant pose different? → conformational rearrangement mechanism

### Alternative: Local HADDOCK3 (WSL required)
```bash
# In WSL Ubuntu:
pip install haddock3
haddock3 docking-peptide.cfg
```

## PDB ID: 1DUZ Reference Peptide

The template peptide LLFGYPVYV has anchors:
- P2=L (canonical Leu in B-pocket)
- P9=V (canonical Val in F-pocket)

Our candidate peptides should be compared against this binding mode.

## Next Steps After Docking

1. Compare WT vs mutant poses for each pair
2. Calculate anchor RMSD relative to template
3. Identify key hydrogen bond/salt bridge differences
4. Cross-validate with MHCflurry binding scores
5. Select top candidates for MD refinement (~100 ns)
