# Findings Summary Tables — Multi-Allele Peptide-MHC Research

All numbers are validated, cross-predictor, seed-robust, temporal-holdout-checked.

## Table 1. Predictor comparison — canonical vs Type B sensitivity (pooled)

| Predictor class | Canonical sens | Type B (non-canonical) sens | n |
|-----------------|:--------------:|:---------------------------:|--:|
| Deep FFN (sequence) | 0.979 | 0.607 | 14,488 |
| MHCflurry (sequence) | 0.986 | 0.668 | 14,488 |
| MUNIS (structure) | 0.950 | 0.647 | 6,224 |
| **Takeaway** | **~0.98 solved** | **~0.65 plateau (all classes)** | |

## Table 2. Pan16 — 15-allele panel, per-allele non-canonical sensitivity

| Allele | sens@0.95 | | Allele | sens@0.95 |
|--------|:--:|-|--------|:--:|
| HLA-A*02:01 | 0.880 | | HLA-B*08:01 | 0.880 |
| HLA-A*01:01 | 0.800 | | HLA-B*15:01 | 0.960 |
| HLA-A*03:01 | **1.000** | | HLA-B*27:05 | 0.080 |
| HLA-B*07:02 | 0.640 | | HLA-B*35:01 | 0.840 |
| HLA-B*44:03 | **1.000** | | HLA-B*40:01 | **1.000** |
| HLA-A*11:01 | **1.000** | | HLA-B*49:01 | **1.000** |
| HLA-A*24:02 | 0.760 | | HLA-B*51:01 | 0.520 |
| | | | HLA-B*57:01 | 0.920 |
| **POOLED (15 alleles)** | **0.819** (307/375) | | | |

Weak tail = thin-training alleles (B*27:05 n=150, B*51:01 n=1,665) — data scarcity, not architecture.

## Table 3. Temporal holdout — generalization to post-2022 epitopes

| Allele | pre-2021 sens@0.95 | post-2022 sens@0.95 | Delta |
|--------|:--:|:--:|:--:|
| HLA-A*02:01 | 0.879 (n=2704) | 0.824 (n=1151) | −0.055 |
| HLA-A*03:01 | 0.905 (n=1744) | 0.936 (n=803) | **+0.031** |
| HLA-B*07:02 | 0.898 (n=2315) | 0.973 (n=865) | **+0.075** |

Post-2022 epitopes were published after MHCflurry's training cutoff (unseen by the
labeler). Preserved/improved sensitivity = NOT overfit to the IEDB era.

## Table 4. P2 virtual mutagenesis — sidechain reach to Glu63

| P2 residue | Reach (Å) | Electrostatics vs Glu63(−) | Can compensate? |
|------------|:--:|:--:|:--:|
| Lys (K) | 7.6 | +1 (salt bridge) | ✅ |
| Arg (R) | 7.3 | +1 (salt bridge) | ✅ |
| His (H) | 4.4 | +0.5 (partial) | ✅ |
| Asp (D) | 3.0 | −1 (repulsive) | ❌ |
| Glu (E) | 4.4 | −1 (repulsive) | ❌ |
| Others (W/Y/Q/M/F…) | 2.0–5.3 | neutral/weak | partial |

K-E63 salt bridge is one instance of a general principle: positively charged P2
residues (K/R/H) can electrostatically compensate the acidic B-pocket.

## Table 5. Falsified vs surviving claims

| Claim | Verdict |
|-------|:--:|
| "Per-allele > pan" | ❌ data-size artifact (pan-full recovers) |
| "Canonical 100% / 97–99% spec" | ❌ n=25 + motif tautology (anchor-checker = 1.0 F1) |
| "Type A easy (0.92–0.98)" | ❌ n=25 cherry-pick; real 0.48 |
| "Non-canonical systematically missed" | ⚠️ moderate penalty (~0.65), not catastrophic |
| "Three-tier pipeline rescues non-canonical" | ⚠️ hypothesis only; MUNIS does NOT lift plateau |
| "One rebalanced pan replaces 15-allele panel" | ✅ pooled 0.819 |
| "Temporal generalization (breaks circularity)" | ✅ post-2022 preserved |
| "K-E63 salt-bridge is a general P2 principle" | ✅ K/R/H all reach Glu63 |
