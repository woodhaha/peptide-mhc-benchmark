#!/usr/bin/env python3
"""Score the Type A/B dataset with MHCflurry (independent predictor). Base python.

Reads data/type_ab_dataset.json, predicts affinity per peptide (per-allele),
writes data/type_ab_mhcflurry_scores.json: {allele: {peptide: aff_nM}}.

Usage (base python, mhcflurry works): python src/score_mhcflurry_type_ab.py
"""
import json, time
import numpy as np

def main():
    from mhcflurry import Class1AffinityPredictor
    predictor = Class1AffinityPredictor.load()
    print("predictor loaded")
    with open("data/type_ab_dataset.json", encoding="utf-8") as f:
        ds = json.load(f)
    out = {}
    for allele, info in ds["alleles"].items():
        entries = info["entries"]
        peps = [e["peptide"] for e in entries]
        # batch predict in chunks to avoid huge single calls
        affs = []
        for i in range(0, len(peps), 500):
            chunk = peps[i:i+500]
            affs.extend(predictor.predict(peptides=chunk, alleles=[allele]*len(chunk)).tolist())
        out[allele] = dict(zip(peps, [round(a, 1) for a in affs]))
        print(f"{allele}: scored {len(peps)} peptides")
    with open("data/type_ab_mhcflurry_scores.json", "w") as f:
        json.dump(out, f, indent=1)
    print("saved data/type_ab_mhcflurry_scores.json")

if __name__ == "__main__":
    main()
