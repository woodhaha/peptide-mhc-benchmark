#!/usr/bin/env python3
"""Experimental-label benchmark — break the label circularity.

Query IEDB for HLA-A*02:01 peptides with EXPERIMENTAL quantitative binding
(IC50/EC50, not MS-elution / predictor-derived). Binder = IC50 < 500 nM;
non-binder = IC50 > 5000 nM. Score with the per-allele Deep FFN (BLOSUM62)
and report sensitivity / specificity / AUC against experimental ground truth.

Usage (venv TF): python experimental_label_benchmark.py
"""
import json, urllib.request, urllib.parse
import numpy as np

AA_ORDER = "ARNDCQEGHILKMFPSTWYV"
BASE = "https://query-api.iedb.org/mhc_export"
ALLELE = "HLA-A*02:01"
SAFE = lambda a: a.replace("*", "_").replace(":", "")
MODEL = f"models/deep_ffn_{SAFE(ALLELE)}.keras"

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

def fetch(quant_filter):
    q = urllib.parse.urlencode({
        "select": "epitope__name,assay__method,assay__quantitative_measurement,reference__pmid,reference__date",
        "mhc_restriction__name": f"eq.{ALLELE}",
        "epitope__object_type": "eq.Linear peptide",
        "assay__quantitative_measurement": quant_filter,
    }, safe="*")
    req = urllib.request.Request(BASE + "?" + q, headers={"Accept": "application/json"})
    with urllib.request.urlopen(req, timeout=300) as r:
        return json.load(r)

def collect():
    binders, nonbinders = {}, {}
    for row in fetch("lte.500"):
        p, ic = row.get("epitope__name"), row.get("assay__quantitative_measurement")
        if not p or len(p) != 9 or ic is None: continue
        if max(p.count(a) for a in set(p)) >= 4: continue
        binders[p] = min(binders.get(p, 1e9), float(ic))
    for row in fetch("gte.5000"):
        p, ic = row.get("epitope__name"), row.get("assay__quantitative_measurement")
        if not p or len(p) != 9 or ic is None: continue
        if max(p.count(a) for a in set(p)) >= 4: continue
        nonbinders[p] = min(nonbinders.get(p, 1e9), float(ic))
    for p in set(binders) & set(nonbinders):
        del binders[p], nonbinders[p]
    return binders, nonbinders

def main():
    import keras
    binders, nonbinders = collect()
    print(f"binders (IC50<500): {len(binders)} | non-binders (IC50>5000): {len(nonbinders)}")
    m = keras.models.load_model(MODEL, compile=False)
    def score(peps):
        X = np.array([encode(p) for p in peps], dtype=np.float32)
        probs = m.predict(X, verbose=0)
        return probs[:, 1] + probs[:, 2]
    sb, sn = score(list(binders)), score(list(nonbinders))
    for th in (0.5, 0.9, 0.95):
        sens = np.mean(sb >= th); spec = np.mean(sn < th)
        print(f"  thr={th}: sensitivity={sens:.3f} (n={len(sb)})  specificity={spec:.3f} (n={len(sn)})")
    y = np.concatenate([np.ones(len(sb)), np.zeros(len(sn))])
    s = np.concatenate([sb, sn])
    try:
        from sklearn.metrics import roc_auc_score
        print(f"  ROC AUC = {roc_auc_score(y, s):.3f}")
    except ImportError:
        order = np.argsort(s)[::-1]; y = y[order]
        tps = np.cumsum(y); auc = np.trapz(tps / tps[-1], np.arange(len(y)) / (len(y)-1)) if tps[-1] else float('nan')
        print(f"  ROC AUC (manual) = {auc:.3f}")
    out = {"binders": {p: round(float(v), 4) for p, v in zip(binders, sb)},
           "nonbinders": {p: round(float(v), 4) for p, v in zip(nonbinders, sn)}}
    with open("experimental_label_scores.json", "w", encoding="utf-8") as f:
        json.dump(out, f, indent=1)
    print("saved experimental_label_scores.json")

if __name__ == "__main__":
    main()
