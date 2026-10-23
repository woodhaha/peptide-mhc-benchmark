#!/usr/bin/env python3
"""Anchor-rule baseline: sensitivity/specificity on the corrected IEDB set.

Rule (manuscript A*02:01 motif): classify a 9-mer as a binder iff
p2 in {L, M, I, V} AND p9 in {V, L}. Computed on the 42 A*02:01-verified
positives + 20 homopolymer negatives (iedb_epitope_verification.csv).
"""
import csv
from collections import Counter

P2 = {"L", "M", "I", "V"}
P9 = {"V", "L"}

ver = {r["peptide"]: r["status"] for r in csv.DictReader(open("02_Data/cleaned/iedb_epitope_verification.csv", encoding="utf-8"))}
bench = {r["peptide"]: r for r in csv.DictReader(open("02_Data/cleaned/iedb_benchmark_results_v2.csv", encoding="utf-8"))}

def anchor(p):
    return bench[p]["p2"] in P2 and bench[p]["p9"] in P9

pos = [p for p in bench if bench[p]["true_label"] == "POS" and ver[p] == "A2-POS"]
neg = [p for p in bench if bench[p]["true_label"] == "NEG"]

tp = sum(1 for p in pos if anchor(p))
fn = len(pos) - tp
tn = sum(1 for p in neg if not anchor(p))
fp = len(neg) - tn

print(f"anchor rule  sens = {tp}/{len(pos)} = {tp/len(pos):.3f}")
print(f"anchor rule  spec = {tn}/{len(neg)} = {tn/len(neg):.3f}")
print(f"missed non-canonical binders (n={fn}):")
for p in pos:
    if not anchor(p):
        print(f"  {p}  p2={bench[p]['p2']} p9={bench[p]['p9']}")
