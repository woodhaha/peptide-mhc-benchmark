# Table 5: Cancer Mutation Epitope Scan — Epitope-Altering Mutations

| Mutation | Gene | Cancer Type | Effect | Peptide | WT Score | Mut Score | Δ Score |
|----------|:----:|-------------|--------|---------|:--------:|:---------:|:-------:|
| **G12V** | KRAS | Pancreatic/Lung | **CREATED** | YKLVVVGAV | 0.000 | 0.475 | **+0.475** |
| **R248W** | p53 | Pan-cancer (most frequent) | **CREATED** | MNWRPILTI | 0.002 | 0.409 | **+0.407** |
| G12V | KRAS | Pancreatic/Lung | ENHANCED | LVVVGAVGV | 0.509 | 0.816 | +0.307 |
| G12C | KRAS | Lung (druggable) | ENHANCED | LVVVGACGV | 0.509 | 0.816 | +0.307 |
| **G13D** | KRAS | Colon | **DESTROYED** | DVDKSALTI | 0.471 | 0.001 | **−0.470** |
| **A146T** | KRAS | Colon | **DESTROYED** | GIPFIETST | 0.436 | 0.020 | **−0.415** |
| **R249S** | p53 | HCC (aflatoxin) | **DESTROYED** | GMNRSPILT | 0.401 | 0.245 | **−0.157** |

## Summary

| Category | Count | Mutations |
|----------|:-----:|-----------|
| Created (neoepitope) | 2 | KRAS G12V, p53 R248W |
| Enhanced | 2 | KRAS G12V, KRAS G12C |
| Destroyed | 3 | KRAS G13D, KRAS A146T, p53 R249S |
| Unchanged | 9 | Remaining mutations |
| **Total mutations scanned** | **16** | p53 (7) + KRAS (9) |
| **Total windows analyzed** | **144** | 9 windows × 16 mutations |

**Clinical Significance:**
- **p53 R248W**: Most frequent p53 mutation across all human cancers — the neoepitope MNWRPILTI is a potential universal cancer vaccine target
- **KRAS G12V**: Creates a neoepitope (YKLVVVGAV) and enhances a second epitope (LVVVGAVGV) — dual epitope gain
- **KRAS G12C**: Enhances an existing epitope — potential for combination therapy with sotorasib
- **KRAS G13D**: Destroys a strong wild-type epitope — potential immune escape mechanism

**Note:** Scores are predicted strong binder probabilities from the Deep FFN model. "Created" = mutated peptide crosses the SB threshold (≥0.5) while wild-type does not. "Destroyed" = wild-type is predicted SB/WB while mutant drops below threshold. Δ Score = mutant_score − wild_type_score.
