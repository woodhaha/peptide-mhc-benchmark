# Multi-Allele Peptide-MHC Benchmark — Companion Data Package

Accompanying data for "A Rebalanced Pan-Specific Model and the Non-Canonical Frontier: A Multi-Allele Benchmark of Peptide-MHC-I Binding Prediction" (submitted to *PLOS Computational Biology*).

## Contents

- `data/` — curated multi-allele dataset and predictor scores
  - `type_ab_dataset.json` — 14,488 IEDB positives (5 alleles; canonical 11,225 / Type A 151 / Type B 3,112)
  - `type_ab_deepffn_scores.json` — Deep FFN binding scores per peptide
  - `type_ab_mhcflurry_scores.json` — MHCflurry affinity scores per peptide
  - `multiseed_results.json` — rebalanced pan model 5-seed non-canonical sensitivities
  - `munis_typeab_scores.csv` — MUNIS (structure-based) scores
- `models/` — trained Keras models
  - `deep_ffn_pan_seed{0,1,42,123,2024}.keras` — rebalanced pan model (5 seeds)
  - `deep_ffn_pan16.keras` — 15-allele panel model
- `code/` — analysis scripts
  - `collect_type_ab.py` — multi-allele data curation
  - `eval_pan_multiseed.py` — 5-seed pan evaluation
  - `eval_pan16.py` — 15-allele panel evaluation
  - `temporal_holdout.py` — temporal holdout validation
  - `score_deepffn_type_ab.py` — Deep FFN scoring
- `figures/` — manuscript figures (Fig 1-4, vector PDF)

## Reproducibility

Full pipeline available at https://github.com/woodhaha/peptide-mhc-benchmark.
Fixed thresholds: 0.95 (pan model); 500 nM (MHCflurry); 5th percentile of canonical scores (MUNIS). Wilson 95% CIs for pooled estimates.
