#!/usr/bin/env python3
"""Falsification Experiment A (Judge 1): anchor-checker baseline.

A zero-parameter rule predicts "binder iff P2+P9 (or P3/P9 for A*01:01) match the
allele's canonical motif". If this trivial rule reproduces the Deep FFN's canonical
F1, the "100% sensitivity" claim is near-tautological (learned anchors only).

Also computes the model's performance on the SAME negative set as the checker so
they're directly comparable. Uses protein-background-matched negatives where the
judges' protein-likeness concern applies.

Usage: python src/falsify_anchor_checker.py   (runs under venv; no training)
"""
import csv, os, sys
import numpy as np

AA_ORDER = "ARNDCQEGHILKMFPSTWYV"
AA_PROT = np.array([0.086,0.027,0.013,0.023,0.059,0.020,0.041,0.067,0.033,0.060,
                    0.077,0.062,0.032,0.012,0.033,0.060,0.058,0.046,0.033,0.016],
                   dtype=np.float64)  # roughly human proteome frequencies
AA_LIST = list(AA_ORDER)
SAFE = lambda a: a.replace("*", "_").replace(":", "")

# canonical anchor positions: {allele: {pos1indexed: allowed_residues}}
ANCHORS = {
    "HLA-A*02:01": {2: "LMIVAT", 9: "LVIM"},
    "HLA-A*01:01": {3: "DE", 9: "YF"},
    "HLA-A*03:01": {2: "LVIM", 9: "KR"},
    "HLA-B*07:02": {2: "P", 9: "LIMF"},
    "HLA-B*44:03": {2: "ED", 9: "FWYL"},
}
ALLELES = list(ANCHORS)

def anchor_checker(pep, allele):
    """Trivial rule: binder iff all anchor positions match the canonical motif."""
    if len(pep) != 9: return 0
    for p1, allowed in ANCHORS[allele].items():
        if pep[p1-1] not in allowed: return 0
    return 1

def load_benchmark(allele):
    with open(f"02_Data/cleaned/iedb_benchmark_{SAFE(allele)}.csv", encoding="utf-8") as f:
        return [r["peptide"] for r in csv.DictReader(f)]

def protein_bg_negatives(n, rng, excl):
    """9-mers sampled at proteome frequencies (protein-likeness matched)."""
    out = []
    while len(out) < n:
        s = ''.join(rng.choice(AA_LIST, size=9, p=AA_PROT/AA_PROT.sum()))
        if s not in excl: out.append(s); excl.add(s)
    return out

def uniform_negatives(n, rng, excl):
    out = []
    while len(out) < n:
        s = ''.join(rng.choice(AA_LIST, size=9))
        if s not in excl: out.append(s); excl.add(s)
    return out

def main():
    rng = np.random.RandomState(42)
    print(f"{'allele':<12}{'model':<10}{'neg':<9}{'sens':>6}{'spec':>6}{'F1':>6}  note")
    print("-"*70)
    for allele in ALLELES:
        pos = load_benchmark(allele)
        excl = set(pos)
        # same two negative types the judges flagged
        for negname, negfn in [("uniform", uniform_negatives), ("protein", protein_bg_negatives)]:
            neg = negfn(100, rng, set(pos))
            all_p = pos + neg
            true = [1]*len(pos)+[0]*len(neg)
            pred = [anchor_checker(p, allele) for p in all_p]
            tp=sum(1 for t,p in zip(true,pred) if t and p); fp=sum(1 for t,p in zip(true,pred) if not t and p)
            fn=sum(1 for t,p in zip(true,pred) if t and not p); tn=sum(1 for t,p in zip(true,pred) if not t and not p)
            sens=tp/(tp+fn) if tp+fn else 0; spec=tn/(tn+fp) if tn+fp else 0
            prec=tp/(tp+fp) if tp+fp else 0; f1=2*prec*sens/(prec+sens) if prec+sens else 0
            print(f"{allele:<12}{'checker':<10}{negname:<9}{sens:>5.3f}{spec:>6.3f}{f1:>6.3f}  rule-only")

if __name__ == "__main__":
    main()
