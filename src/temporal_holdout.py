#!/usr/bin/env python3
"""Temporal holdout validation — break the label circularity.

MHCflurry 2.x trained on IEDB data up to ~2021. Score epitopes published 2022+
(unseen by the labeler) vs pre-2021 (in-distribution). If the model's sensitivity
holds on post-2022 epitopes, the performance is NOT an artifact of training on the
same IEDB era.

Usage (venv TF): python src/temporal_holdout.py
"""
import json, urllib.request, urllib.parse, collections
import numpy as np
import sys
sys.path.insert(0, "src")

AA_ORDER = "ARNDCQEGHILKMFPSTWYV"
BASE = "https://query-api.iedb.org/mhc_export"
SAFE = lambda a: a.replace("*", "_").replace(":", "")
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

def encode(pep):
    v = []
    for aa in pep:
        if aa not in AA_ORDER: aa = 'A'
        v.extend(BLOSUM62_NORM[AA_ORDER.index(aa)].tolist())
    return np.array(v, dtype=np.float32)

def fetch(allele):
    q = urllib.parse.urlencode({
        "select": "epitope__name,reference__date",
        "mhc_restriction__name": f"eq.{allele}",
        "epitope__object_type": "eq.Linear peptide",
        "assay__qualitative_measurement": "eq.Positive",
    }, safe="*")
    req = urllib.request.Request(BASE+"?"+q, headers={"Accept":"application/json"})
    with urllib.request.urlopen(req, timeout=120) as r:
        return json.load(r)

def main():
    import keras
    for allele in ["HLA-A*02:01", "HLA-A*03:01", "HLA-B*07:02"]:
        rows = fetch(allele)
        by_bucket = collections.defaultdict(set)
        for x in rows:
            n = x.get("epitope__name"); yr = str(x.get("reference__date") or "")[:4]
            if not n or len(n) != 9: continue
            if max(n.count(c) for c in set(n)) >= 4: continue
            bucket = "pre2021" if yr and yr < "2021" else ("post2022" if yr and yr >= "2022" else "mid")
            by_bucket[bucket].add(n)
        m = keras.models.load_model(f"models/deep_ffn_{SAFE(allele)}.keras", compile=False)
        print(f"=== {allele} ===")
        for bucket in ["pre2021", "post2022"]:
            peps = list(by_bucket[bucket])
            if not peps: continue
            X = np.array([encode(p) for p in peps], dtype=np.float32)
            probs = m.predict(X, verbose=0)
            scores = probs[:,1] + probs[:,2]
            s05 = sum(1 for s in scores if s >= 0.5)/len(scores)
            s09 = sum(1 for s in scores if s >= 0.9)/len(scores)
            s095 = sum(1 for s in scores if s >= 0.95)/len(scores)
            print(f"  {bucket:<8} n={len(peps):>4}  sens@0.5={s05:.3f}  sens@0.9={s09:.3f}  sens@0.95={s095:.3f}")

if __name__ == "__main__":
    main()
