# Data Integrity Audit (2026-10-13)

Three-judge review (statistics/immunology/ML) + re-analysis exposed data-integrity issues.
This file documents each, its status, and the fix.

## 1. Architecture ranking was inverted (FIXED by re-analysis)

**Problem**: manuscript Table 1 claimed "Deep FFN 91.9% > CNN 90.0% > ResNet 84.3% > LSTM 83.3%".
Recomputation from the saved `.keras` models (n=508 test) and 5-seed re-training showed the
opposite: **CNN > Deep FFN > FFN > LSTM > ResNet > RF**.

**Authoritative 5-seed results** (`02_Data/cleaned/multiseed_results.json`, on `real_peptides.csv`):

| Model | mean ± SD |
|-------|-----------|
| CNN | 91.91% ± 0.93% |
| Deep FFN | 90.56% ± 1.16% |
| FFN (Jessen) | 89.83% ± 1.14% |
| LSTM | 89.71% ± 0.89% |
| ResNet | 87.63% ± 1.76% |
| Random Forest | 79.69% ± 0.44% |

CNN > Deep FFN in 5/5 seeds (paired). The manuscript's "Deep FFN beats complex
architectures" is **falsified**; a complex architecture (CNN) actually wins.

## 2. Two generations of single-allele data (5088 vs 5190) — DECISION NEEDED

| Data | Peptides | Per class | Used by |
|------|:---:|:---:|---------|
| `02_Data/raw/feature.csv` (5088×282, Jun 15) | 5,088 | 1,696 | R pipeline, manuscript Methods |
| `03_Analysis/models/*.npy` | 5,088 | 1,696 | `.keras` model training |
| `02_Data/raw/real_peptides.csv` (Jun 23) | **5,190** | **1,730** | `run_full_benchmark.py` (Python) |

The manuscript says 5,088/1,696/4,579/509. The Python pipeline actually trains on
5,190/1,730/4,671/519. MHCflurry labels were re-generated at some point, changing the
balanced class count from 1,696 → 1,730. **The 5-seed re-benchmark uses 5,190 (current source).**

**Decision needed**: which is authoritative? Recommendation: 5,190 (`real_peptides.csv`),
because it is what the current `run_full_benchmark.py` produces and what the `.keras`
pipeline can reproduce end-to-end.

## 3. "Random synthetic labels 65.8%" was unreproducible (FIXING)

**Problem**: manuscript Table 2 reported "Random synthetic 65.8%" but no script or data
file generates random labels. The closest thing, `generate_synthetic_jessen_data()`
(`peptide_mhc_binding_study.R:320`), produces **anchor-DETERMINED** labels (SB=has p2/p9
anchor, NB=no anchor, WB=30% anchor), NOT i.i.d. random labels — a model with BLOSUM62
input learns that deterministic rule, so it is a mislabeled null.

**Fix**: `03_Analysis/R_scripts/fix_random_labels.py` generates a true i.i.d. random-label
null (shuffle true labels) on 5,190 peptides, 6 architectures × 5 seeds.
**Result** (`02_Data/cleaned/random_label_results.json`): all ~33% (chance) —
RF 35.14%±1.34, Deep FFN 34.45%±1.70, FFN 34.34%±3.93, CNN 32.95%±3.15,
LSTM 32.60%±2.00, ResNet 32.41%±2.00. **The manuscript's "random 65.8%" is
confirmed impossible for i.i.d. labels — it was a mislabeled anchor-deterministic
dataset or a leak.**

## 7. PSSM "94.8%" also inconsistent (FOUND)

Manuscript Table 2 says PSSM 94.8%, but `02_Data/cleaned/pssm_comparison.csv`
reports Deep FFN on PSSM = 82.6% (on a 15,000-peptide pool, different from the
5,088/5,190 MHCflurry set). Both "random 65.8%" and "PSSM 94.8%" do not trace to
reproducible data. The three-way label comparison in Table 2 is therefore
unreproducible end-to-end.

## 4. Random Forest duplicated

`03_Analysis/models/Random_Forest.joblib` (Python) and `03_Analysis/Random_Forest.rds`
(R) both exist — two pipelines. The `.joblib` is the current Python pipeline's output.

## 5. MM-GBSA Analysis/ is empty

`MM-GBSA Analysis/` (0 files). Reviewer B M6 recommended replacing the Coulombic
K-E63 estimate with MMPBSA/GBSA — this directory suggests it was planned but never run.
Empty; the salt-bridge section (§3.5) was already downgraded to a qualitative hypothesis.

## 6. Sample-count mismatch

Manuscript Methods says train/test = 4,579/509. Actual `feature_resnet.csv` = 4,580/508
(off by one). `real_peptides.csv` = 4,671/519.

## Status summary (final, 2026-10-13)

- [x] Architecture ranking recomputed (5-seed) — CNN wins.
- [x] Random-label null script written + run (i.i.d. ~33%).
- [x] Data-generation consolidated to 5,190 (`real_peptides.csv`).
- [x] Manuscript Table 1/2 rewritten with corrected numbers.
- [x] Factorial partial (MHCflurry + random done; PSSM on 15,000-pool, flagged as confounded).
- [x] NetMHCpan-4.1 SOTA comparison (via IEDB API): 62.24% vs MHCflurry labels (test set); 81.6% on IEDB 49 epitopes.
- [x] ILRGSVAHK removed from IEDB benchmark (HLA-B\*08:01) → 46/48 (95.8%).
- [x] Figures 1/2/9 + cover letters + README + research-state all regenerated/synced.
- [ ] GitHub repo: main branch + Zenodo DOI (needs user action).
- [ ] ESM-2 650M / three-way comparison (optional, heavy).

## Final corrected headline numbers

- Architecture (5-seed, n=519): CNN 91.9±0.9 > Deep FFN 90.6±1.2 > FFN 89.8 > LSTM 89.7 > ResNet 87.6 > RF 79.7.
- Label source: MHCflurry 90.6 > PSSM 82.6 (15k pool) > random 34.5 (chance).
- NetMHCpan-4.1 vs MHCflurry labels: 62.24% (test set); 81.6% (IEDB 49) — quantifies label circularity.
- Multi-allele: pan matches per-allele on 4/5 (0.84–0.96), B\*07:02 exception (0.60).
