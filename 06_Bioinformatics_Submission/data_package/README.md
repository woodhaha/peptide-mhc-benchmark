# Peptide-MHC Binding Prediction — Data Package

Accompanying data for "Benchmark Hygiene Reveals the Non-Canonical Frontier of Peptide-MHC-I Binding Prediction" (submitted to *Bioinformatics*).

## Contents

- `data/` — curated multi-allele dataset and predictor scores
  - `type_ab_dataset.json` — 14,488 IEDB positives (5 alleles; canonical 11,225 / Type A 151 / Type B 3,112)
  - `type_ab_deepffn_scores.json` — Deep FFN binding scores per peptide
  - `type_ab_mhcflurry_scores.json` — MHCflurry affinity scores per peptide
  - `multiseed_results.json` — rebalanced pan model 5-seed non-canonical sensitivities
  - `munis_typeab_scores.csv` — MUNIS (sequence-based) scores
- `models/` — trained Keras models
  - `deep_ffn_pan_seed{0,1,42,123,2024}.keras` — rebalanced pan model (5 seeds)
  - `deep_ffn_pan16.keras` — 15-allele panel model
  - `deep_ffn_pan_b7x4.keras` — B*07:02 upsampled 4× variant
- `code/` — analysis scripts
  - `collect_type_ab.py` — multi-allele data curation
  - `eval_pan_multiseed.py` — 5-seed pan evaluation
  - `eval_pan16.py` — 15-allele panel evaluation
  - `score_deepffn_type_ab.py` — Deep FFN scoring
  - `temporal_holdout.py` — temporal holdout validation
- `figures/` — publication-quality figure PDFs

## Reproducibility

Full pipeline available at https://github.com/woodhaha/peptide-mhc-benchmark.
Random seed for train/test split: 42.
