#!/usr/bin/env python3
"""Evaluate per-allele models against NON-canonical validated binders.

Tests the manuscript's core structural-blind-spot claim across alleles: models
trained partly on canonical-motif data may miss non-canonical binders (e.g. KRAS
G12V's P2-K). Uses IEDB positives that fail the allele's canonical anchor motif.
"""
import json, sys, urllib.request, urllib.parse, collections
import numpy as np

AA_ORDER = "ARNDCQEGHILKMFPSTWYV"
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

MOTIFS = {
    "HLA-A*02:01": {"p2": "LMIVAT", "p9": "LVIM"},
    "HLA-A*01:01": {"p3": "DE", "p9": "YF"},
    "HLA-A*03:01": {"p2": "LVIM", "p9": "KR"},
    "HLA-B*07:02": {"p2": "P", "p9": "LIMF"},
    "HLA-B*44:03": {"p2": "ED", "p9": "FWYL"},
}
ALLELES = ["HLA-A*02:01", "HLA-A*01:01", "HLA-A*03:01", "HLA-B*07:02", "HLA-B*44:03"]
BASE = "https://query-api.iedb.org/mhc_export"

def is_canonical(pep, motif):
    if len(pep) != 9: return None
    return all(pep[int(p[1:])-1] in a for p, a in motif.items())

def encode(peptide):
    vec = []
    for aa in peptide:
        if aa not in AA_ORDER: aa = 'A'
        vec.extend(BLOSUM62_NORM[AA_ORDER.index(aa)].tolist())
    return np.array(vec, dtype=np.float32)

def fetch_noncannonical(allele, top=25):
    motif = MOTIFS[allele]
    q = urllib.parse.urlencode({
        "select": "epitope__name,reference__pmid,assay__method",
        "mhc_restriction__name": f"eq.{allele}",
        "epitope__object_type": "eq.Linear peptide",
        "assay__qualitative_measurement": "eq.Positive",
    }, safe="*")
    with urllib.request.urlopen(urllib.request.Request(BASE+"?"+q, headers={"Accept":"application/json"}), timeout=120) as r:
        rows = json.load(r)
    by = collections.defaultdict(list)
    for row in rows:
        n = row.get("epitope__name")
        if not n:
            continue
        if not len(n) == 9:  # skip non-9-mers
            continue
        if is_canonical(n, motif):  # keep only NON-canonical
            continue
        if max(n.count(a) for a in set(n)) >= 4:  # skip homopolymer-ish
            continue
        by[n].append(row)
    ranked = sorted(by.items(), key=lambda kv: -len({x.get("reference__pmid") for x in kv[1] if x.get("reference__pmid")}))
    return [n for n, _ in ranked[:top]]

def main():
    import keras
    model_path, allele = sys.argv[1], sys.argv[2]
    pos = fetch_noncannonical(allele)
    rng = np.random.RandomState(7)
    aa = list("ACDEFGHIKLMNPQRSTVWY")
    excl = set(pos)
    neg = []
    while len(neg) < 40:
        s = ''.join(rng.choice(aa) for _ in range(9))
        if s not in excl: neg.append(s); excl.add(s)
    neg += [a*9 for a in aa]
    all_p = pos + neg
    true = [1]*len(pos) + [0]*len(neg)
    X = np.array([encode(p) for p in all_p], dtype=np.float32)
    if "--pan" in sys.argv:
        o = np.zeros(len(ALLELES), dtype=np.float32)
        o[ALLELES.index(allele)] = 1.0
        X = np.concatenate([X, np.repeat(o[None, :], len(X), axis=0)], axis=1)
    if "--pssm" in sys.argv:
        allowed = set()
        for _anpos, res in {"HLA-A*02:01": {2: "LMIVAT", 9: "LVIM"}, "HLA-A*01:01": {3: "DE", 9: "YF"},
                            "HLA-A*03:01": {2: "LVIM", 9: "KR"}, "HLA-B*07:02": {2: "P", 9: "LIMF"},
                            "HLA-B*44:03": {2: "ED", 9: "FWYL"}}[allele].items():
            allowed.update(res)
        v = np.array([1.0 if a in allowed else 0.05 for a in "ACDEFGHIKLMNPQRSTVWY"], dtype=np.float32)
        v /= v.sum()
        X = np.concatenate([X, np.repeat(v[None, :], len(X), axis=0)], axis=1)
    m = keras.models.load_model(model_path, compile=False)
    probs = m.predict(X, verbose=0)
    scores = probs[:, 1] + probs[:, 2]
    res = {}
    for th in [0.5, 0.7, 0.8, 0.9, 0.95, 0.97, 0.99]:
        pred = (scores >= th).astype(int)
        tp = sum(1 for t,p in zip(true,pred) if t and p); fp = sum(1 for t,p in zip(true,pred) if not t and p)
        fn = sum(1 for t,p in zip(true,pred) if t and not p)
        sens = tp/(tp+fn) if tp+fn else 0
        prec = tp/(tp+fp) if tp+fp else 0
        f1 = 2*prec*sens/(prec+sens) if prec+sens else 0
        res[f"t{th}"] = {"sens": round(sens,3), "prec": round(prec,3), "f1": round(f1,3), "tp": tp, "fn": fn}
    out = {"allele": allele, "n_noncannonical_pos": len(pos), "results": res}
    print(json.dumps(out))
    sys.exit(0)

if __name__ == "__main__":
    main()
