#!/usr/bin/env python3
"""Test H4: does the pan model dilute the B*07:02 P2-proline anchor?

B*07:02's motif is the most distinctive (P2=P) and its training set smallest (606).
Hypothesis: pan-full misses non-canonical B*07:02 binders because the shared pathway
under-weights the P2-proline anchor. If the misses are P2=P-present (they keep P2 but
fail P9), the pan model is under-using proline; if P2=P-absent, it's a different cause.

Usage (venv TF): python src/test_b0702_h4.py
"""
import json, urllib.request, urllib.parse, collections
import numpy as np

AA_ORDER = "ARNDCQEGHILKMFPSTWYV"
AA_LIST = list(AA_ORDER)
AA_PROT = np.array([0.086,0.027,0.013,0.023,0.059,0.020,0.041,0.067,0.033,0.060,
                    0.077,0.062,0.032,0.012,0.033,0.060,0.058,0.046,0.033,0.016],
                   dtype=np.float64)
ALLELES = ["HLA-A*02:01", "HLA-A*01:01", "HLA-A*03:01", "HLA-B*07:02", "HLA-B*44:03"]
BROAD_B7 = {2: "P", 9: "LIMFV"}

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

def fetch_noncannonical(allele, motif, top=25):
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

def main():
    import keras
    allele = "HLA-B*07:02"
    pos = fetch_noncannonical(allele, BROAD_B7, top=25)
    m = keras.models.load_model("models/deep_ffn_pan_full.keras", compile=False)
    o = np.zeros(len(ALLELES), dtype=np.float32); o[ALLELES.index(allele)] = 1.0
    X = np.array([np.concatenate([encode(p), o]) for p in pos], dtype=np.float32)
    probs = m.predict(X, verbose=0)
    scores = probs[:,1] + probs[:,2]
    print(f"B*07:02 non-canonical positives (broad): {len(pos)}")
    print(f"{'peptide':<12}{'P2':>3}{'P9':>4}{'score':>7}  {'P2=P?'}")
    print("-"*45)
    for p, s in sorted(zip(pos, scores), key=lambda x: x[1]):
        p2 = p[1]; p9 = p[8]
        flag = "  <-- P2=P miss" if p2=='P' and s<0.9 else ("  <-- P2-abs miss" if s<0.9 else "")
        print(f"{p:<12}{p2:>3}{p9:>4}{s:>7.3f}  {flag}")
    misses = [(p,s) for p,s in zip(pos,scores) if s<0.9]
    p2p_miss = sum(1 for p,s in misses if p[1]=='P')
    print(f"\nMisses @0.9: {len(misses)}/25; P2=P-present among misses: {p2p_miss}/{len(misses)}")

if __name__ == "__main__":
    main()
