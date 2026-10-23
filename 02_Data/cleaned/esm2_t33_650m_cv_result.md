# ESM-2 t33 (650M) encoding — 5-fold CV (traceability re-run)

Re-run 2026-10-23 to back the Table 3 claim `92.4 ± 0.5% (CV)`.

- **Data**: `02_Data/raw/real_peptides.csv` — n=5,190, HLA-A\*02:01, 1,730 each SB/WB/NB
- **Model**: `facebook/esm2_t33_650M_UR50D` (fp16, RTX 2000 Ada 8GB)
- **Encoding**: per-position embeddings, 9 × 1280 = 11,520 dim
- **Head**: FFN 256→128→64→3, Adam lr 1e-3, 150 epochs/fold
- **Seed**: 42 (`StratifiedKFold(5, shuffle=True, random_state=42)`)
- **Script**: `esm2_650m_benchmark_local.py`

| Fold | Accuracy (%) |
|------|-------------:|
| 1    | 92.00 |
| 2    | 91.91 |
| 3    | 93.35 |
| 4    | 92.39 |
| 5    | 92.20 |

**Mean ± SD: 92.37% ± 0.52%** → reported as **92.4 ± 0.5 (CV)** in Table 3 / Figure 7.
