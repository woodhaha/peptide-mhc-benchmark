#!/usr/bin/env python3
"""Re-measure the surviving kernel under corrected methodology.

Fixes demanded by the three judges:
  (1) protein-background negatives (not uniform-random) so the model can't win on
      protein-likeness;
  (2) broadened motifs (validated supertype anchors) so "non-canonical" isn't a
      dict-coverage artifact.

Re-runs per-allele vs pan-full on the broadened non-canonical positives, reporting
sens/spec/F1 at thresholds [0.5, 0.8, 0.9, 0.95]. Honest n + Wilson CI upper bound
for small n.

Usage (venv TF): python src/reeval_noncannonical_corrected.py
"""
import csv, json, sys
import numpy as np

AA_ORDER = "ARNDCQEGHILKMFPSTWYV"
AA_LIST = list(AA_ORDER)
AA_PROT = np.array([0.086,0.027,0.013,0.023,0.059,0.020,0.041,0.067,0.033,0.060,
                    0.077,0.062,0.032,0.012,0.033,0.060,0.058,0.046,0.033,0.016],
                   dtype=np.float64)
ALLELES = ["HLA-A*02:01", "HLA-A*01:01", "HLA-A*03:01", "HLA-B*07:02", "HLA-B*44:03"]
SAFE = lambda a: a.replace("*", "_").replace(":", "")
BROAD = {
    "HLA-A*02:01": {2: "LMIVATQ", 9: "LVIM"},
    "HLA-A*01:01": {3: "DE", 9: "YF"},
    "HLA-A*03:01": {2: "LVIM", 9: "KRYH"},
    "HLA-B*07:02": {2: "P", 9: "LIMFV"},
    "HLA-B*44:03": {2: "ED", 9: "FWYL"},
}
BLOSUM62 = np.array([
    [ 4,-1,-2,-2, 0,-1,-1, 0,-2,-1,-1,-1,-1,-2,-1, 1, 0,-3,-2, 0],
    [-1, 5, 0,-2,-3, 1, 0,-2, 0,-3,-2, 2,-1,-3,-2,-1,-1,-3,-2,-3],
    [-2, 0, 6, 1,-3, 0, 0, 0, 1,-3,-3, 0,-2,-3,-2, 1, 0,-4,-2,-3],
    [-2,-2, 1, 6,-3, 0, 2,-1,-1,-3,-4,-1,-3,-3,-1, 0,-1,-4,-3,-3],
    [ 0,-3,-3,-3, 9,-3,-4,-3,-3,-1,-1,-3,-1,-2,-3,-1,-1,-2,-2,-1],
    [-1, 1, 0, 0,-3, 5, 2,-2, 0,-3,-2, 1, 0,-3,-1, 0,-1,-2,-1,-2],
    [-1, 0, 0, 2,-4, 2, 5,-2, 0,-3,-3, 1,-2,-3,-1, 0,-1,-3,-2,-2],
    [ 0,-2, 0,-1,-3,-2,-2, 6,-2,-4,-4,-2,-3,-3,-2, 0,-2,-2,-3,-3],
    [-2, 0, 1,-1,-3, 0, 0,-2, 8,-3,-3,-1,-2,-1,-2,-1,-2,-2, 2,-3],
    [-1,-3,-3,-3,-1,-3,-3,-4,-3, 4, 2,-3, 1, 0,-3,-2,-1,-3,-1, 3],
    [-1,-2,-3,-4,-1,-2,-3,-4,-3, 2, 4,-2, 2, 0,-3,-2,-1,-2,-1, 1],
    [-1, 2, 0,-1,-3, 1, 1,-2,-1,-3,-2, 5,-1,-3,-1, 0,-1,-3,-2,-2],
    [-1,-1,-2,-3,-1, 0,-2,-3,-2, 1, 2,-1, 5, 0,-2,-1,-1,-1,-1, 1],
    [-2,-3,-3,-3,-2,-3,-3,-3,-1, 0, 0,-3, 0, 6,-4,-2,-2, 1, 3,-1],
    [-1,-2,-2,-1,-3,-1,-1,-2,-2,-3,-3,-1,-2,-4, 7,-1,-1,-4,-3,-2],
    [ 1,-1, 1, 0,-1, 0, 0, 0,-1,-2,-2, 0,-1,-2,-1, 4, 1,-3,-2,-2],
    [ 0,-1, 0,-1,-1,-1,-1,-2,-2,-1,-1,-1,-1,-2,-1, 1, 5,-2,-2, 0],
    [-3,-3,-4,-4,-2,-2,-3,-2, 2,-3,-2,-3,-1, 1,-4,-3,-2,11, 2,-3],
    [-2,-2,-2,-3,-2,-1,-2,-3, 2,-1,-1,-2,-1, 3,-3,-2,-2, 2, 7,-1],
    [ 0,-3,-3,-3,-1,-2,-2,-3,-3, 3, 1,-2, 1,-1,-2,-2, 0,-3,-1, 4],
], dtype=np.float32)
BLOSUM62_NORM = np.array([(row-row.min())/(row.max()-row.min()+1e-8) for row in BLOSUM62])

def encode(peptide):
    v = []
    for aa in peptide:
        if aa not in AA_ORDER: aa = 'A'
        v.extend(BLOSUM62_NORM[AA_ORDER.index(aa)].tolist())
    return np.array(v, dtype=np.float32)

def is_canon(pep, motif):
    return len(pep)==9 and all(pep[int(p)-1] in r for p, r in motif.items())

def load_benchmark(allele):
    with open(f"02_Data/cleaned/iedb_benchmark_{SAFE(allele)}.csv", encoding="utf-8") as f:
        return [r["peptide"] for r in csv.DictReader(f)]

def fetch_noncannonical(allele, motif, top=25):
    """Top-evidence IEDB positives that FAIL the (broad) canonical motif."""
    import urllib.request, urllib.parse, collections
    q = urllib.parse.urlencode({
        "select": "epitope__name,reference__pmid",
        "mhc_restriction__name": f"eq.{allele}",
        "epitope__object_type": "eq.Linear peptide",
        "assay__qualitative_measurement": "eq.Positive",
    }, safe="*")
    req = urllib.request.Request("https://query-api.iedb.org/mhc_export?"+q, headers={"Accept":"application/json"})
    with urllib.request.urlopen(req, timeout=120) as r:
        rows = json.load(r)
    by = collections.defaultdict(list)
    for x in rows:
        n = x.get("epitope__name")
        if not n or len(n) != 9 or is_canon(n, motif):
            continue
        if max(n.count(a) for a in set(n)) >= 4:
            continue
        by[n].append(x)
    ranked = sorted(by.items(), key=lambda kv: -len({z.get("reference__pmid") for z in kv[1] if z.get("reference__pmid")}))
    return [n for n, _ in ranked[:top]]

def protein_negatives(n, rng, excl):
    out = []
    while len(out) < n:
        s = ''.join(rng.choice(AA_LIST, size=9, p=AA_PROT/AA_PROT.sum()))
        if s not in excl: out.append(s); excl.add(s)
    return out

def wilson(k, n, z=1.96):
    if n == 0: return (0, 0)
    p = k/n
    denom = 1 + z*z/n
    centre = (p + z*z/(2*n))/denom
    half = z*np.sqrt(p*(1-p)/n + z*z/(4*n*n))/denom
    return (max(0, centre-half), min(1, centre+half))

def main():
    import keras
    rng = np.random.RandomState(7)
    models = {}
    out = {}
    for allele in ALLELES:
        pos = fetch_noncannonical(allele, BROAD[allele], top=25)
        if len(pos) < 5:
            out[allele] = {"n": len(pos), "skip": True}; continue
        # build shared negatives once
        excl = set(pos)
        neg = protein_negatives(100, rng, excl)
        all_p = pos + neg
        true = np.array([1]*len(pos)+[0]*len(neg))
        results = {}
        for name, mpath, is_pan in [
            ("per_allele", f"models/deep_ffn_{SAFE(allele)}.keras", False),
            ("pan_full", "models/deep_ffn_pan_full.keras", True),
        ]:
            m = models.get(mpath) or keras.models.load_model(mpath, compile=False)
            models[mpath] = m
            X = np.array([encode(p) for p in all_p], dtype=np.float32)
            if is_pan:
                o = np.zeros(len(ALLELES), dtype=np.float32); o[ALLELES.index(allele)] = 1.0
                X = np.concatenate([X, np.repeat(o[None,:], len(X), axis=0)], axis=1)
            probs = m.predict(X, verbose=0)
            scores = probs[:, 1] + probs[:, 2]
            per_th = {}
            for th in [0.5, 0.8, 0.9, 0.95]:
                pred = (scores >= th).astype(int)
                tp = int(((true==1)&(pred==1)).sum()); fp = int(((true==0)&(pred==1)).sum())
                fn = int(((true==1)&(pred==0)).sum()); tn = int(((true==0)&(pred==0)).sum())
                sens = tp/(tp+fn) if tp+fn else 0
                spec = tn/(tn+fp) if tn+fp else 0
                prec = tp/(tp+fp) if tp+fp else 0
                f1 = 2*prec*sens/(prec+sens) if prec+sens else 0
                lo, hi = wilson(tp, tp+fn)
                per_th[str(th)] = {"sens": round(sens,3), "spec": round(spec,3),
                                   "f1": round(f1,3), "tp": tp, "fn": fn,
                                   "sens_ci95": [round(lo,3), round(hi,3)]}
            results[name] = per_th
        out[allele] = {"n_noncanon_broad": len(pos), "results": results}
    print(json.dumps(out, indent=1))

if __name__ == "__main__":
    main()
