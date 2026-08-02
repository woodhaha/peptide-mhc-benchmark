#!/usr/bin/env python3
"""Aggregate Type A/B large-sample results across both predictors.

Deep FFN binder = binding_score (P_SB+P_WB) >= 0.95.
MHCflurry binder = affinity < 500 nM.
Reports per-allele, per-type sensitivity + Wilson 95% CI, pooled across alleles.

Usage (any python): python src/summarize_type_ab.py
"""
import json
import numpy as np

ALLELES = ["HLA-A*02:01", "HLA-A*01:01", "HLA-A*03:01", "HLA-B*07:02", "HLA-B*44:03"]
TYPES = ["canonical", "type_a", "type_b"]
DFF_TH = 0.95
MHC_TH = 500.0

def wilson(k, n):
    if n == 0: return (0, 0)
    p = k/n; z = 1.96
    denom = 1 + z*z/n
    centre = (p + z*z/(2*n))/denom
    half = z*np.sqrt(p*(1-p)/n + z*z/(4*n*n))/denom
    return (max(0, centre-half), min(1, centre+half))

def main():
    with open("data/type_ab_dataset.json", encoding="utf-8") as f:
        ds = json.load(f)
    with open("data/type_ab_deepffn_scores.json", encoding="utf-8") as f:
        dff = json.load(f)
    with open("data/type_ab_mhcflurry_scores.json", encoding="utf-8") as f:
        mhc = json.load(f)

    # pooled counters per type
    pool = {t: {"dff_n": 0, "dff_hit": 0, "mhc_n": 0, "mhc_hit": 0} for t in TYPES}
    print(f"{'allele':<14}{'type':<10}{'n':>5} | {'DFF sens':>9}{'CI':>22} | {'MHC sens':>9}{'CI':>22}")
    print("-"*95)
    per_allele = {}
    for allele in ALLELES:
        info = ds["alleles"][allele]
        dff_s = dff.get(allele, {}); mhc_s = mhc.get(allele, {})
        by_type = {t: [] for t in TYPES}
        for e in info["entries"]:
            by_type[e["type"]].append(e["peptide"])
        row = {}
        for t in TYPES:
            peps = by_type[t]
            dff_hits = sum(1 for p in peps if dff_s.get(p, 0) >= DFF_TH)
            mhc_hits = sum(1 for p in peps if mhc_s.get(p, 9e9) < MHC_TH)
            n = len(peps)
            dff_sens = dff_hits/n if n else 0
            mhc_sens = mhc_hits/n if n else 0
            dff_ci = wilson(dff_hits, n); mhc_ci = wilson(mhc_hits, n)
            pool[t]["dff_n"] += n; pool[t]["dff_hit"] += dff_hits
            pool[t]["mhc_n"] += n; pool[t]["mhc_hit"] += mhc_hits
            row[t] = {"n": n, "dff_sens": round(dff_sens,3), "mhc_sens": round(mhc_sens,3)}
            print(f"{allele:<14}{t:<10}{n:>5} | {dff_sens:>8.3f}  [{dff_ci[0]:.3f},{dff_ci[1]:.3f}] | {mhc_sens:>8.3f}  [{mhc_ci[0]:.3f},{mhc_ci[1]:.3f}]")
        per_allele[allele] = row
    print("\n=== POOLED across 5 alleles ===")
    for t in TYPES:
        dff_s = pool[t]["dff_hit"]/pool[t]["dff_n"] if pool[t]["dff_n"] else 0
        mhc_s = pool[t]["mhc_hit"]/pool[t]["mhc_n"] if pool[t]["mhc_n"] else 0
        dff_ci = wilson(pool[t]["dff_hit"], pool[t]["dff_n"])
        mhc_ci = wilson(pool[t]["mhc_hit"], pool[t]["mhc_n"])
        print(f"{t:<10} n={pool[t]['dff_n']:>5} | DFF {dff_s:.3f} [{dff_ci[0]:.3f},{dff_ci[1]:.3f}] | MHC {mhc_s:.3f} [{mhc_ci[0]:.3f},{mhc_ci[1]:.3f}]")

if __name__ == "__main__":
    main()
