
# Phase A: Neoepitope Peptide-MHC Docking Validation
## Structural Analysis Report

**Date**: 2026-06-16
**Project**: HLA-A*02:01 peptide-MHC docking for KRAS G12V & p53 R248W neoepitopes
**Template**: 1DUZ (HLA-A*02:01 + LLFGYPVYV, 1.80A)
**Method**: Template-based anchor scoring + comparative B-pocket analysis

---

## 1. Neoepitope Candidates

| Peptide | Source | Effect | ΔScore | P2 | P9 | P2 Canon? | Anchor Score |
|---------|--------|--------|:---:|:---:|:---:|:---:|:---:|
| `YKLVVVGAV` | KRAS G12V | CREATED (NB→WB) | +0.475 | K | V | ❌ | 2.3 |
| `MNWRPILTI` | p53 R248W | CREATED (NB→WB) | +0.407 | N | I | ❌ | 1.0 |
| `LVVVGAVGV` | KRAS G12V | ENHANCED (WB→SB) | +0.307 | V | V | ✅ | 8.4 |
| `YKLVVVGAG` | KRAS WT | Non-binder | -- | K | G | ❌ | -2.3 |

## 2. Key Discovery: K-E63 Salt Bridge Hypothesis

**Observation**: KRAS G12V neoepitope has Lys at P2 anchor position, which is
non-canonical for HLA-A*02:01 (canonical P2 = L/M/I/V). However, GLU63 in the
B-pocket is positioned 4.7A from the P2 CA in the template structure.

**Comparative B-pocket analysis of 8 HLA-A*02:01 structures**:
- B-pocket volume: 11.3-11.5A (highly conserved, stdev < 0.1A)
- E63-P2 distance: 4.5-4.9A (consistent across all structures)
- P2 residue: 100% Leu in deposited structures (no K ever observed)
- E63 orientation: carboxyl group points toward P2 side chain

**Hypothesis**: Lys at P2 can form a stabilizing salt bridge with GLU63 in the
B-pocket. The K side chain length (4 methylene + NH3+) is sufficient to reach
the E63 carboxyl group (estimated K NZ to E63 OE1/OE2 distance: 2.5-3.5A).

This would represent a **novel non-canonical anchor compensation mechanism**:
```
Canonical:  P2(L) → hydrophobic burial in B-pocket (ΔG ≈ -3 kcal/mol)
Proposed:   P2(K) → K-NZ:H3+  ---  E63-OE1/OE2 (salt bridge) + partial burial
```

## 3. p53 R248W: P7 Conformational Switch

**Observation**: Both WT (MNRRPILTI) and neo (MNWRPILTI) share the same P2=N
anchor. The CREATED signal comes entirely from the P7 R→W mutation.

**Structural interpretation**:
- R at P7: positively charged, solvent-exposed, no groove contact
- W at P7: bulky hydrophobic, potential groove surface burial
- The mutation does NOT affect the canonical anchor positions
- The +0.407 ΔScore suggests improved groove complementarity via W burial

## 4. KRAS G12V Enhanced: Canonical Control

`LVVVGAVGV` has P2=V, P9=V — the optimal canonical anchor pair for HLA-A*02:01.
Anchor score (8.4) matches known immunodominant epitopes (GIL: 8.3, NLV: 8.0).
The G→A at P6 is a minimal change that slightly increases hydrophobicity.

**Prediction**: This peptide should show the strongest actual binding among all
candidates. Serves as a positive control for the docking validation.

## 5. Docking Validation Plan

### Priority 1: KRAS G12V neo vs WT
- Receptor: receptor_1DUZ.pdb (HLA-A*02:01 + B2M)
- Ligands: peptide_KRAS_G12V_neo.pdb vs peptide_KRAS_G12V_wt.pdb
- Submit to: http://hdock.phys.hust.edu.cn/
- Key metric: K-E63 salt bridge formation in top poses

### Priority 2: p53 R248W neo vs WT
- Ligands: peptide_p53_R248W_neo.pdb vs peptide_p53_R248W_wt.pdb
- Key metric: P7 conformational difference, groove burial

### Priority 3: KRAS G12V enh vs WT
- Ligands: peptide_KRAS_G12V_enh.pdb vs peptide_KRAS_G12V_enh_wt.pdb
- Expected: strong binding for both, minimal difference

### Post-docking MD validation
- 100 ns MD simulation for top complexes
- K-E63 distance stability over trajectory
- MM/PBSA binding free energy calculation

## 6. Significance

If validated, the K-E63 salt bridge mechanism would be:
1. The first documented case of a charged P2 anchor in HLA-A*02:01
2. A novel neoantigen presentation mechanism for KRAS G12V
3. Clinically relevant (KRAS G12V in 30% of pancreatic, 20% of lung cancers)
4. Potentially generalizable to other K-containing P2 peptides

## 7. Files Prepared

```
docking/
├── 1DUZ.pdb                        # Template structure
├── receptor_1DUZ.pdb               # Prepared receptor
├── peptide_*.pdb                    # 10 peptide structures
├── analyze_neoepitopes.py          # Initial anchor analysis
├── template_docking.py             # Template-based scoring
├── compare_pockets.py              # Comparative B-pocket analysis
├── DOCKING_README.md               # Protocol for HDOCK submission
└── references/                     # 7 reference HLA-A*02:01 structures
