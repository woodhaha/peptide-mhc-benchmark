#!/usr/bin/env python3
"""Evaluate the pan16 model (15-allele panel) on non-canonical Type B binders.

The pan16 model input = 180 BLOSUM62 + 15-dim allele one-hot. This script must use
the SAME allele order as training (PANEL_ANCHORS key order in scale_panel.py).

Usage (venv TF): python src/eval_pan16.py
"""
import json, urllib.request, urllib.parse, collections
import numpy as np

AA_ORDER = "ARNDCQEGHILKMFPSTWYV"
BASE = "https://query-api.iedb.org/mhc_export"
# MUST match scale_panel.PANEL_ANCHORS insertion order
PANEL = ["HLA-A*02:01","HLA-A*01:01","HLA-A*03:01","HLA-B*07:02","HLA-B*44:03",
         "HLA-A*11:01","HLA-A*24:02","HLA-B*08:01","HLA-B*15:01","HLA-B*27:05",
         "HLA-B*35:01","HLA-B*40:01","HLA-B*49:01","HLA-B*51:01","HLA-B*57:01"]
BROAD = {
    "HLA-A*02:01": {2: "LMIVATQ", 9: "LVIM"},
    "HLA-A*01:01": {3: "DE", 9: "YF"},
    "HLA-A*03:01": {2: "LVIM", 9: "KRYH"},
    "HLA-B*07:02": {2: "P", 9: "LIMFV"},
    "HLA-B*44:03": {2: "ED", 9: "FWYL"},
    "HLA-A*11:01": {2: "LVIM", 9: "KR"},
    "HLA-A*24:02": {2: "YF", 9: "FYLIM"},
    "HLA-B*08:01": {3: "RK", 9: "L"},
    "HLA-B*15:01": {2: "Q", 9: "FYWL"},
    "HLA-B*27:05": {2: "RK", 9: "FYWLIM"},
    "HLA-B*35:01": {2: "P", 9: "YFMWL"},
    "HLA-B*40:01": {2: "E", 9: "L"},
    "HLA-B*49:01": {2: "E", 9: "LIM"},
    "HLA-B*51:01": {2: "P", 9: "ILVW"},
    "HLA-B*57:01": {2: "ATS", 9: "FWYL"},
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

def encode(pep):
    v = []
    for aa in pep:
        if aa not in AA_ORDER: aa = 'A'
        v.extend(BLOSUM62_NORM[AA_ORDER.index(aa)].tolist())
    return np.array(v, dtype=np.float32)

def is_canon(pep, motif):
    if len(pep) != 9: return False
    return all(pep[int(p)-1] in r for p, r in motif.items())

def fetch_noncanon(allele, motif, top=25):
    q = urllib.parse.urlencode({
        "select": "epitope__name,reference__pmid",
        "mhc_restriction__name": f"eq.{allele}",
        "epitope__object_type": "eq.Linear peptide",
        "assay__qualitative_measurement": "eq.Positive",
    }, safe="*")
    req = urllib.request.Request(BASE+"?"+q, headers={"Accept":"application/json"})
    with urllib.request.urlopen(req, timeout=120) as r:
        rows = json.load(r)
    by = collections.defaultdict(list)
    for x in rows:
        n = x.get("epitope__name")
        if not n or len(n) != 9 or is_canon(n, motif): continue
        if max(n.count(a) for a in set(n)) >= 4: continue
        by[n].append(x)
    ranked = sorted(by.items(), key=lambda kv: -len({z.get("reference__pmid") for z in kv[1] if z.get("reference__pmid")}))
    return [n for n, _ in ranked[:top]]

def main():
    import keras
    m = keras.models.load_model("models/deep_ffn_pan16.keras", compile=False)
    print(f"{'allele':<14}{'n':>4}  sens@0.5  sens@0.9  sens@0.95")
    print("-"*50)
    total_n = total_hit = 0
    for i, allele in enumerate(PANEL):
        motif = BROAD.get(allele)
        if not motif: continue
        pos = fetch_noncanon(allele, motif, top=25)
        if not pos: continue
        o = np.zeros(len(PANEL), dtype=np.float32); o[i] = 1.0
        X = np.array([np.concatenate([encode(p), o]) for p in pos], dtype=np.float32)
        pr = m.predict(X, verbose=0)
        sc = pr[:,1] + pr[:,2]
        s05 = sum(1 for s in sc if s>=0.5)/len(pos)
        s09 = sum(1 for s in sc if s>=0.9)/len(pos)
        s095 = sum(1 for s in sc if s>=0.95)/len(pos)
        total_n += len(pos); total_hit += sum(1 for s in sc if s>=0.95)
        print(f"{allele:<14}{len(pos):>4}  {s05:.3f}   {s09:.3f}   {s095:.3f}")
    print(f"\nPOOLED: {total_hit}/{total_n} sens@0.95={total_hit/total_n:.3f}")

if __name__ == "__main__":
    main()
