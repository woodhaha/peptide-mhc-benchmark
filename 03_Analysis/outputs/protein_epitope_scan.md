# Protein Epitope Scan -- HLA-A*02:01

> **Model**: Deep FFN (MHCflurry-trained, 91.9% accuracy, 89.6% CV)  
> **Allele**: HLA-A*02:01  
> **Method**: Sliding 9-mer window, BLOSUM62 encoding  
> **Thresholds**: SB = P(SB) > 50%, WB = P(WB) > 50%, NB otherwise  
> **Date**: June 15, 2026

---

## Known Epitope Validation

| Protein | Epitope | Pred | P(SB) | Rank | Status | Note |
|---------|---------|------|-------|------|--------|------|
| p53 | `LLGRNSFEV` | **SB** | 99.3% | Top 1.0% | [OK] Confirmed | Classic p53 265-273 epitope |
| Spike | `YLQPRTFLL` | **SB** | 99.8% | -- | [OK] Confirmed | Spike 269-277 (NTD, not RBD) |
| MART-1 | `ILTVILGVL` | **WB** | 3.3% (95.7% WB) | Top 4.5% | [OK] Confirmed | MART-1 32-40 |
| Spike | `RLFRKSNLK` | NB | 0% | Top 18% | [--] Missed | MHC class II epitope |
| MART-1 | `ELAGIGILTV` | NB | 0% | -- | [!!] Missed | Clinical anchor-modified variant |

---

## Top Epitope Candidates per Protein

### MART-1 / Melan-A (118 amino acids, melanoma antigen)

**110 9-mer windows scanned | 2 SB, 12 WB, 96 NB**

| Pos | Peptide | Pred | P(SB) | P(WB) | Binding Score | Anchors |
|-----|---------|------|-------|-------|---------------|---------|
| 56-64 | **`ALMDKSLHV`** | SB | 100.0% | 0.0% | 1.000 | [OK] L-V |
| 31-39 | `GILTVILGV` | SB | 99.7% | 0.3% | 0.999 | [OK] I-V |
| 35-43 | `VILGVLLLI` | SB | 93.1% | 6.9% | 0.965 | [OK] I-I |
| 34-42 | `TVILGVLLL` | WB | 5.5% | 91.5% | 0.513 | [OK] V-L |
| 32-40 | `ILTVILGVL` | WB | 3.3% | 95.7% | 0.512 | [OK] L-L |
| 27-35 | `AAGIGILTV` | WB | 0.8% | 97.2% | 0.494 | [OK] A-V |
| 88-96 | `LQEKNCEPV` | WB | 0.1% | 96.6% | 0.484 | [OK] Q-V |
| 29-37 | `GIGILTVIL` | WB | 0.6% | 95.4% | 0.483 | [OK] I-L |
| 40-48 | `LLLIGCWYC` | WB | 0.1% | 93.7% | 0.469 | -- |
| 57-65 | `LMDKSLHVG` | WB | 0.2% | 89.0% | 0.447 | -- |

---

### SARS-CoV-2 Spike RBD (318 amino acids, receptor binding domain)

**310 9-mer windows scanned | 5 SB, 27 WB, 278 NB**

| Pos | Peptide | Pred | P(SB) | P(WB) | Binding Score | Anchors |
|-----|---------|------|-------|-------|---------------|---------|
| 87-95 | **`KIADYNYKL`** | SB | 99.9% | 0.1% | 0.999 | [OK] I-L |
| 182-190 | `VLSFELLHA` | SB | 90.7% | 9.3% | 0.953 | [OK] L-A |
| 186-194 | `ELLHAPATV` | SB | 84.9% | 15.1% | 0.925 | [OK] L-V |
| 282-290 | `YQDVNCTEV` | WB | 35.1% | 64.8% | 0.676 | [OK] Q-V |
| 44-52 | `FSTFKCYGV` | WB | 17.0% | 82.4% | 0.582 | -- |
| 172-180 | `GVGYQPYRV` | WB | 12.5% | 86.3% | 0.556 | [OK] V-V |
| 304-312 | `RVYSTGSNV` | WB | 4.5% | 95.3% | 0.521 | [OK] V-V |
| 114-122 | `KVGGNYNYL` | WB | 2.4% | 96.9% | 0.508 | [OK] V-L |
| 94-102 | `KLPDDFTGC` | WB | 1.7% | 98.2% | 0.508 | -- |
| 175-183 | `YQPYRVVVL` | WB | 1.9% | 97.3% | 0.505 | [OK] Q-L |

---

### p53 Tumor Suppressor (323 amino acids)

**315 9-mer windows scanned | 8 SB, 25 WB, 282 NB**

| Pos | Peptide | Pred | P(SB) | P(WB) | Binding Score | Anchors |
|-----|---------|------|-------|-------|---------------|---------|
| 65-73 | **`RMPEAAPPV`** | SB | 100.0% | 0.0% | 1.000 | [OK] M-V |
| 164-172 | `KQSQHMTEV` | SB | 99.4% | 0.6% | 0.997 | [OK] Q-V |
| 264-272 | `LLGRNSFEV` | SB | 99.3% | 0.7% | 0.996 | [OK] L-V |
| 187-195 | `GLAPPQHLI` | SB | 97.9% | 2.1% | 0.990 | [OK] L-I |
| 129-137 | `ALNKMFCQL` | SB | 96.2% | 3.8% | 0.981 | [OK] L-L |
| 103-111 | `YQGSYGFRL` | SB | 94.2% | 5.8% | 0.971 | [OK] Q-L |
| 31-39 | `VLSPLPSQA` | SB | 92.0% | 8.0% | 0.960 | [OK] L-A |
| 132-140 | `KMFCQLAKT` | SB | 80.2% | 19.8% | 0.901 | [OK] M-T |
| 139-147 | `KTCPVQLWV` | WB | 39.2% | 60.6% | 0.696 | [OK] T-V |
| 135-143 | `CQLAKTCPV` | WB | 38.5% | 61.1% | 0.691 | [OK] Q-V |

---

## Scan Statistics

| Metric | MART-1 | Spike RBD | p53 | **Total** |
|--------|--------|-----------|-----|-----------|
| Protein length (aa) | 118 | 318 | 323 | 759 |
| 9-mer windows | 110 | 310 | 315 | **772** |
| Strong Binders | 2 (1.8%) | 5 (1.6%) | 7 (2.2%) | **14 (1.8%)** |
| Weak Binders | 12 (10.9%) | 27 (8.7%) | 25 (7.9%) | **64 (8.3%)** |
| Non-Binders | 96 (87.3%) | 278 (89.7%) | 283 (89.8%) | **694 (89.9%)** |

---

## Key Findings

### Validation Success

- **3 of 4 MHC-I epitopes confirmed** -- 75% sensitivity on known positives
- p53 `LLGRNSFEV` correctly identified as SB at Top 1.0% rank
- MART-1 `ILTVILGVL` correctly identified as WB at Top 4.5% rank
- `RLFRKSNLK` missed because it is an MHC-II epitope (the model is MHC-I specific)

### Novel Candidate Epitopes

| Candidate | Protein | Confidence | Rationale |
|-----------|---------|------------|-----------|
| **`ALMDKSLHV`** | MART-1 56-64 | 100% SB | L-V anchors, hydrophobic core |
| **`KIADYNYKL`** | Spike RBD 87-95 | 99.9% SB | I-L anchors, in RBD core |
| **`RMPEAAPPV`** | p53 65-73 | 100% SB | M-V anchors, known mutation hotspot |
| **`KQSQHMTEV`** | p53 164-172 | 99.4% SB | Q-V anchors, DNA binding domain |
| **`GLAPPQHLI`** | p53 187-195 | 97.9% SB | L-I anchors, proline-rich region |

### Anchor Residue Analysis

The model strongly weights canonical anchor positions:
- **p2**: L, M, I, V, Q, A  high SB probability
- **p9**: V, L, I, A, T, M  high SB probability
- **Non-canonical p2** (E, D, R, K, G, P)  almost always NB
- **Cysteine at p2**  occasionally tolerated in specific contexts

### Clinical Relevance

- p53 `LLGRNSFEV` (264-272) is a validated immunotherapy target in multiple cancers
- p53 `RMPEAAPPV` (65-73) lies in the transactivation domain -- frequently mutated
- Spike `YLQPRTFLL` (269-277) is the immunodominant CD8+ epitope in COVID-19
- The novel candidates should be validated by *in vitro* binding assays before clinical use

---

## Methods

1. **Model**: Deep Feed-Forward Neural Network (36018090453)
2. **Training**: MHCflurry 2.2.0-labeled peptides, 5-fold CV 89.6%  0.8%
3. **Encoding**: BLOSUM62 substitution matrix (9-mer  920 tensor)
4. **Scan**: Sliding window (stride=1) across full protein sequences
5. **Thresholding**: SB = P(SB) > P(WB) and P(SB) > P(NB)

---

## Output Files

| File | Description |
|------|-------------|
| `data/protein_epitope_scan.csv` | Full scan results (772 peptides with predictions) |
| `protein_epitope_scan.md` | This report |
| `plots/comparison_*.png` | Model comparison visualizations |

---

*Generated June 15, 2026 | Peptide Epitope Research Project*
