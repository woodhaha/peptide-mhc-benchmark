#!/usr/bin/env python3
"""Evaluate a per-allele Deep FFN against the curated IEDB benchmark.

Usage: python src/evaluate_allele.py <model.keras> <allele> [--out result.json]
Runs under the r-tensorflow venv (TF/keras), reads 02_Data/cleaned/iedb_benchmark_{SAFE}.csv.
Prints JSON metrics and exits 0 on pass (sens>=0.90, spec_raw>=0.80, F1>=0.85).
"""
import argparse, csv, json, os, sys
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

SAFE = lambda a: a.replace("*", "_").replace(":", "")

def encode(peptide):
    vec = []
    for aa in peptide:
        if aa not in AA_ORDER: aa = 'A'
        vec.extend(BLOSUM62_NORM[AA_ORDER.index(aa)].tolist())
    return np.array(vec, dtype=np.float32)

def load_benchmark(allele, basedir="02_Data/cleaned"):
    path = os.path.join(basedir, f"iedb_benchmark_{SAFE(allele)}.csv")
    with open(path, encoding="utf-8") as f:
        return [r["peptide"] for r in csv.DictReader(f)]

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("model")
    ap.add_argument("allele")
    ap.add_argument("--out", default=None)
    ap.add_argument("--basedir", default="02_Data/cleaned")
    args = ap.parse_args()

    import keras
    pos = load_benchmark(args.allele, args.basedir)
    if len(pos) < 15:
        print(json.dumps({"pass": False, "error": f"benchmark too small: {len(pos)} positives"}))
        sys.exit(1)

    # Negative set: diverse (random 9-mers + protein windows + homopolymers),
    # same construction as the original evaluator for comparability.
    rng = np.random.RandomState(42)
    aa = list("ACDEFGHIKLMNPQRSTVWY")
    excl = set(pos)
    neg = []
    while len(neg) < 40:
        s = ''.join(rng.choice(aa) for _ in range(9))
        if s not in excl:
            neg.append(s); excl.add(s)
    # albumin windows
    alb = ("MKWVTFISLLFLFSSAYSRGVFRRDAHKSEVAHRFKDLGEENFKALVLIAFAQYLQQCPFEDHVK"
           "LVNEVTEFAKTCVADESAENCDKSLHTLFGDKLCTVATLRETYGEMADCCAKQEPERNECFLQHK"
           "DDNPNLPRLVRPEVDVMCTAFHDNEETFLKKYLYEIARRHPYFYAPELLFFAKRYKAAFTECCQA"
           "ADKAACLLPKLDELRDEGKASSAKQRLKCASLQKFGERAFKAWAVARLSQRFPKAEFAEVSKLVT")
    for i in range(len(alb)-8):
        w = alb[i:i+9]
        if all(c in AA_ORDER for c in w) and w not in excl:
            neg.append(w); excl.add(w)
        if len(neg) >= 80: break
    neg += [a*9 for a in aa]
    neg = list(dict.fromkeys(neg))

    all_p = pos + neg
    true = [1]*len(pos) + [0]*len(neg)
    X = np.array([encode(p) for p in all_p], dtype=np.float32)

    model = keras.models.load_model(args.model, compile=False)
    probs = model.predict(X, verbose=0)
    # binding score = P(SB) + P(WB); threshold 0.95 (from real-negative mission)
    scores = probs[:, 1] + probs[:, 2]

    best = None
    for th in [0.90, 0.95, 0.97, 0.99]:
        pred = (scores >= th).astype(int)
        tp = sum(1 for t,p in zip(true,pred) if t==1 and p==1)
        fp = sum(1 for t,p in zip(true,pred) if t==0 and p==1)
        fn = sum(1 for t,p in zip(true,pred) if t==1 and p==0)
        tn = sum(1 for t,p in zip(true,pred) if t==0 and p==0)
        sens = tp/(tp+fn) if tp+fn else 0
        spec = tn/(tn+fp) if tn+fp else 0
        prec = tp/(tp+fp) if tp+fp else 0
        f1 = 2*prec*sens/(prec+sens) if prec+sens else 0
        m = {"threshold": th, "sensitivity": round(sens,4), "specificity_raw": round(spec,4),
             "f1": round(f1,4), "tp": tp, "fp": fp, "tn": tn, "fn": fn,
             "n_pos": len(pos), "n_neg": len(neg)}
        if best is None or (m["f1"], m["sensitivity"]) > (best["f1"], best["sensitivity"]):
            best = m
    best["pass"] = best["sensitivity"] >= 0.90 and best["specificity_raw"] >= 0.80 and best["f1"] >= 0.85
    best.update({"allele": args.allele, "model": args.model})
    print(json.dumps(best))
    if args.out:
        with open(args.out, "w") as f:
            json.dump(best, f, indent=2)
    sys.exit(0 if best["pass"] else 1)

if __name__ == "__main__":
    main()
