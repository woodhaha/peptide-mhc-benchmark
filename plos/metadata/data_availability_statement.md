# Data and Code Availability Statement

All code, trained models, feature matrices, prediction results, and the experimental-label benchmark (`experimental_label_scores.json`: 6,524 peptides — 3,438 binders with measured IC50 <500 nM and 3,086 non-binders with IC50 >5,000 nM) are available at https://github.com/woodhaha/peptide-mhc-benchmark and archived at Zenodo (DOI: 10.5281/zenodo.23615309).

The Deep FFN is distributed as the `pmhcbenchmark` Python package, with a `pmhc-predict` command-line tool and a standalone Windows GUI (`pmhc_predict.exe`); exported weights (`pmhc_models.npz`) enable tensorflow-free inference. The single-allele benchmark-hygiene analysis constitutes the complete study.

A manuscript-consistency checker (`peptide_audit.py` with `audit_config.json`) automates the reference-order, factual-tripwire, and figure-inventory checks used in this audit.

