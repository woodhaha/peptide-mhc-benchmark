#!/usr/bin/env python3
"""Evaluate 5-seed rebalanced pan models on non-canonical binders, aggregate mean+/-SD.

Usage (venv TF): python src/eval_pan_multiseed.py --seeds 0,1,42,123,2024
Outputs per-allele mean+/-SD non-canonical sens@0.95 across seeds.
"""
import argparse, json, urllib.request, urllib.parse, collections
import numpy as np

AA_ORDER = "ARNDCQEGHILKMFPSTWYV"
AA_LIST = list(AA_ORDER)
AA_PROT = np.array([0.086,0.027,0.013,0.023,0.059,0.020,0.041,0.067,0.033,0.060,
                    0.077,0.062,0.032,0.012,0.033,0.060,0.058,0.046,0.033,0.016],
                   dtype=np.float64)
ALLELES = ["HLA-A*02:01", "HLA-A*01:01", "HLA-A*03:01", "HLA-B*07:02", "HLA-B*44:03"]
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

def encode(p):
    v = []
    for aa in p:
        if aa not in AA_ORDER: aa = 'A'
        v.extend(BLOSUM62_NORM[AA_ORDER.index(aa)].tolist())
    return np.array(v, dtype=np.float32)

def is_canon(p, motif):
    return len(p)==9 and all(p[int(k)-1] in r for k, r in motif.items())

def fetch_noncanon(allele, motif, top=25):
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
        if not n or len(n)!=9 or is_canon(n, motif): continue
        if max(n.count(a) for a in set(n)) >= 4: continue
        by[n].append(x)
    ranked = sorted(by.items(), key=lambda kv: -len({z.get("reference__pmid") for z in kv[1] if z.get("reference__pmid")}))
    return [n for n,_ in ranked[:top]]

def protein_neg(n, rng, excl):
    out = []
    while len(out) < n:
        s = ''.join(rng.choice(AA_LIST, size=9, p=AA_PROT/AA_PROT.sum()))
        if s not in excl: out.append(s); excl.add(s)
    return out

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seeds", default="0,1,42,123,2024")
    ap.add_argument("--modeldir", default="models")
    args = ap.parse_args()
    import keras
    seeds = [int(s) for s in args.seeds.split(",")]
    rng = np.random.RandomState(7)
    out = {}
    for allele in ALLELES:
        pos = fetch_noncanon(allele, BROAD[allele], top=25)
        excl = set(pos); neg = protein_neg(100, rng, excl)
        all_p = pos + neg
        true = np.array([1]*len(pos)+[0]*len(neg))
        o = np.zeros(len(ALLELES), dtype=np.float32); o[ALLELES.index(allele)] = 1.0
        X = np.array([encode(p) for p in all_p], dtype=np.float32)
        Xp = np.concatenate([X, np.repeat(o[None,:], len(X), axis=0)], axis=1)
        sens_per_seed = []
        for seed in seeds:
            m = keras.models.load_model(os.path.join(args.modeldir, f"deep_ffn_pan_seed{seed}.keras"), compile=False)
            sc = m.predict(Xp, verbose=0)[:, 1] + m.predict(Xp, verbose=0)[:, 2]
            pred = (sc >= 0.95).astype(int)
            tp = int(((true==1)&(pred==1)).sum()); fn = int(((true==1)&(pred==0)).sum())
            sens_per_seed.append(tp/(tp+fn) if tp+fn else 0)
        out[allele] = {"sens_per_seed": [round(s,3) for s in sens_per_seed],
                       "mean": round(float(np.mean(sens_per_seed)),3),
                       "sd": round(float(np.std(sens_per_seed)),3)}
        print(f"{allele:<14} mean={out[allele]['mean']:.3f} ± {out[allele]['sd']:.3f}  {out[allele]['sens_per_seed']}")
    with open("06_Logs/multiseed_results.json", "w") as f:
        json.dump(out, f, indent=1)

if __name__ == "__main__":
    main()
