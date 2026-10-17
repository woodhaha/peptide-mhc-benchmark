"""Peptide-MHC-I binding prediction with a per-allele Deep FFN (BLOSUM62).

Scores 9-mer peptides for five common HLA class I alleles. The score is
P(weak binder) + P(strong binder); >= 0.5 is classified as a binder.

Reference: "The Non-Canonical Frontier in Peptide-MHC-I Binding Prediction:
A Benchmark-Hygiene Case Study" (Zhou et al., 2026).
"""
from .predictor import ALLELES, load_predictor, predict

__version__ = "1.0.0"
__all__ = ["predict", "load_predictor", "ALLELES", "__version__"]
