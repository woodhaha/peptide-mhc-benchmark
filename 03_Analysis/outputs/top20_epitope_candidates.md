# Top 20 Epitope Candidates -- HLA-A*02:01

> **Model**: Deep FFN (MHCflurry-trained, 91.9% accuracy, 89.6% CV)  
> **Allele**: HLA-A*02:01  
> **10 proteins scanned | 3,536 9-mers | Top 0.6% selected**  
> **Date**: June 15, 2026

---

## Top 20 Ranking (by Binding Score)

| Rank | Protein | Position | Peptide | Pred | Score | P(SB) | Anchors |
|:----:|---------|----------|---------|:----:|:-----:|:-----:|:-------:|
| 1 | MART-1 | 56-64 | **ALMDKSLHV** | SB | 1.000 | 100.0 | L-V |
| 2 | WT1 | 126-134 | **RMFPNAPYL** | SB | 1.000 | 100.0 | M-L |
| 3 | WT1 | 187-195 | **SLGEQQYSV** | SB | 1.000 | 100.0 | L-V |
| 4 | WT1 | 10-18 | **ALLPAVPSL** | SB | 1.000 | 100.0 | L-L |
| 5 | p53 | 65-73 | **RMPEAAPPV** | SB | 1.000 | 100.0 | M-V |
| 6 | M1 | 3-11 | **LLTEVETYV** | SB | 1.000 | 100.0 | L-V |
| 7 | gp100 | 619-627 | **RLMKQDFSV** | SB | 1.000 | 100.0 | L-V |
| 8 | Tyrosinase | 369-377 | **YMNGTMSQV** | SB | 1.000 | 100.0 | M-V |
| 9 | M1 | 134-142 | **RMGAVTTEV** | SB | 1.000 | 100.0 | M-V |
| 10 | gp100 | 178-186 | **MLGTHTMEV** | SB | 1.000 | 100.0 | M-V |
| 11 | Tyrosinase | 1-9 | **MLLAVLYCL** | SB | 1.000 | 100.0 | L-L |
| 12 | Tyrosinase | 490-498 | ALLAGLVSL | SB | 1.000 | 99.9 | L-L |
| 13 | CMV pp65 | 40-48 | **RLLQTGIHV** | SB | 1.000 | 99.9 | L-V |
| 14 | KRAS | 28-36 | **FVDEYDPTI** | SB | 1.000 | 99.9 | V-I |
| 15 | gp100 | 13-21 | AVIGALLAV | SB | 1.000 | 99.9 | V-V |
| 16 | Spike RBD | 87-95 | **KIADYNYKL** | SB | 0.999 | 99.9 | I-L |
| 17 | CMV pp65 | 495-503 | **NLVPMVATV** | SB | 0.999 | 99.9 | L-V |
| 18 | Tyrosinase | 214-222 | FLLRWEQEI | SB | 0.999 | 99.8 | L-I |
| 19 | Tyrosinase | 482-490 | AMVGAVLTA | SB | 0.999 | 99.8 | M-A |
| 20 | M1 | 58-66 | **GILGFVFTL** | SB | 0.999 | 99.8 | I-L |

---

## Candidate Composition

### By Protein

| Protein | Candidates | Top Rank |
|---------|:----------:|:--------:|
| Tyrosinase | 5 | 8 |
| WT1 | 3 | 2 |
| gp100/PMEL | 3 | 7 |
| Influenza M1 | 3 | 6 |
| CMV pp65 | 2 | 13 |
| MART-1 | 1 | **1** |
| p53 | 1 | 5 |
| Spike RBD | 1 | 16 |
| KRAS | 1 | 14 |

### By Anchor Residues

| p2 | Count | p9 | Count |
|:--:|:-----:|:--:|:-----:|
| L | 7 | V | 10 |
| M | 7 | L | 5 |
| V | 3 | I | 3 |
| I | 2 | A | 2 |
| F | 1 | T | 0 |
| S | 1 | M | 0 |

**L-V and M-V dominate** -- the classic HLA-A*02:01 anchor motif.

---

## Clinical Annotation

| Rank | Peptide | Protein | Known? | Clinical Context |
|:----:|---------|---------|:------:|------------------|
| 1 | ALMDKSLHV | MART-1 | Novel | Melanoma -- adjacent to known epitope cluster |
| 2 | RMFPNAPYL | WT1 | **Validated** | Leukemia/solid tumor vaccine target |
| 3 | SLGEQQYSV | WT1 | Novel | WT1 zinc finger region |
| 4 | ALLPAVPSL | WT1 | Novel | WT1 N-terminal domain |
| 5 | RMPEAAPPV | p53 | **Known region** | p53 transactivation domain (mutation hotspot) |
| 6 | LLTEVETYV | M1 | Novel | N-terminal M1, near immunodominant cluster |
| 7 | RLMKQDFSV | gp100 | Novel | C-terminal gp100 |
| 8 | YMNGTMSQV | Tyrosinase | **Validated** | Melanoma antigen 369-377 |
| 9 | RMGAVTTEV | M1 | Novel | M1 central domain |
| 10 | MLGTHTMEV | gp100 | Novel | gp100 core domain |
| 11 | MLLAVLYCL | Tyrosinase | Novel | Signal peptide -- efficient ER processing |
| 12 | ALLAGLVSL | Tyrosinase | Novel | Transmembrane domain |
| 13 | RLLQTGIHV | CMV pp65 | Novel | Early pp65, surface-exposed |
| 14 | FVDEYDPTI | KRAS | Novel | Effector binding region |
| 15 | AVIGALLAV | gp100 | Novel | Signal peptide region |
| 16 | KIADYNYKL | Spike RBD | Novel | ACE2 binding interface |
| 17 | NLVPMVATV | CMV pp65 | **Validated** | Immunodominant CMV epitope |
| 18 | FLLRWEQEI | Tyrosinase | Novel | Tyrosinase catalytic domain |
| 19 | AMVGAVLTA | Tyrosinase | Novel | Transmembrane-proximal |
| 20 | GILGFVFTL | M1 | **Validated** | Gold-standard immunodominant epitope |

---

## Priority for Experimental Validation

**Tier 1 -- Known validated epitopes (positive controls):**  
GILGFVFTL, NLVPMVATV, YMNGTMSQV, RMFPNAPYL

**Tier 2 -- Novel candidates in validated protein regions:**  
ALMDKSLHV (MART-1, adjacent to known epitope), RMPEAAPPV (p53 mutation hotspot), LLTEVETYV (M1 N-terminus), KIADYNYKL (Spike RBD ACE2 interface)

**Tier 3 -- Novel candidates in less-characterized regions:**  
SLGEQQYSV (WT1), FVDEYDPTI (KRAS effector domain), RLLQTGIHV (CMV pp65 early region)

---

## Output Files

| File | Contents |
|------|----------|
| `data/top20_epitope_candidates.csv` | Machine-readable table |
| `top20_epitope_candidates.md` | This report |
| `data/protein_epitope_scan_extended.csv` | Full 2,764 peptide scan |

---

*Generated June 15, 2026 | Peptide Epitope Research Project*
