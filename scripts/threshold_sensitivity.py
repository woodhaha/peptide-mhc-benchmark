#!/usr/bin/env python3
"""Threshold-sensitivity of the ~0.65 Type B plateau across 4 predictors.

Checks whether the reported Type B sensitivities (DFFN 0.607, MHCflurry 0.668,
MUNIS 0.647, BigMHC 0.660) are stable under threshold variation, or an artifact
of the three different operation points (0.95 / 500 nM / 5th-percentile).
"""
import csv, json
import numpy as np

# ---- load labels ----
d = json.load(open("06_Logs/type_ab_dataset.json"))
labels = {}  # (allele, peptide) -> type
for a, v in d["alleles"].items():
    for e in v["entries"]:
        labels[(a, e["peptide"])] = e["type"]

# ---- load scores ----
def _flatten(nested):
    out = {}
    for a, inner in nested.items():
        for p, s in inner.items():
            out[(a, p)] = s
    return out

dffn = _flatten(json.load(open("06_Logs/type_ab_deepffn_scores.json")))
mhcflurry = _flatten(json.load(open("06_Logs/type_ab_mhcflurry_scores.json")))

munis = {}
_lines = open("05_Submission/munis_results/munis_typeab_scores.csv").read().splitlines()
# first line is a stray filename; real header is the second line
for r in csv.DictReader(_lines[1:]):
    munis[(r["mhc"], r["pep"])] = float(r["score"])

bigmhc = {}
for f in ["A0101", "A0201", "A0301", "B0702", "B4403"]:
    for r in csv.DictReader(open(f"06_Bioinformatics_Submission/data_package/bigmhc_results/{f}.prd")):
        bigmhc[(r["mhc"], r["pep"])] = float(r["BigMHC_EL"])


def get_scores(pred):
    """Return (canonical_scores, typeb_scores) pooled across alleles."""
    can, tb = [], []
    for (a, p), t in labels.items():
        if t not in ("canonical", "type_b"):
            continue
        s = pred.get((a, p))
        if s is None:
            continue
        (can if t == "canonical" else tb).append(s)
    return np.array(can), np.array(tb)


preds = {
    "DFFN": dffn,
    "MHCflurry(nM)": mhcflurry,
    "MUNIS": munis,
    "BigMHC_EL": bigmhc,
}

# "binder" = score > threshold for DFFN/MUNIS/BigMHC; < threshold for MHCflurry affinity
HIGH = {"DFFN": True, "MHCflurry(nM)": False, "MUNIS": True, "BigMHC_EL": True}
# manuscript operation points
FIXED_THR = {"DFFN": 0.95, "MHCflurry(nM)": 500.0}

for name, pred in preds.items():
    can, tb = get_scores(pred)
    high = HIGH[name]
    can_sorted = np.sort(can)
    print(f"\n=== {name}: canonical n={len(can)}, typeB n={len(tb)} ===")

    # 1. reproduce manuscript threshold
    if name in FIXED_THR:
        thr = FIXED_THR[name]
        cs = (can > thr).mean() if high else (can < thr).mean()
        ts = (tb > thr).mean() if high else (tb < thr).mean()
        print(f"  MANUSCRIPT thr={thr:g}: canon-sens={cs:.3f}, TypeB sens={ts:.3f}")
    else:  # 5th percentile of canonical
        k = int(round(0.05 * len(can)))
        thr = can_sorted[k] if high else can_sorted[-1 - k]
        cs = (can > thr).mean() if high else (can < thr).mean()
        ts = (tb > thr).mean() if high else (tb < thr).mean()
        print(f"  MANUSCRIPT 5th-pctile thr={thr:.4g}: canon-sens={cs:.3f}, TypeB sens={ts:.3f}")

    # 2. aligned-threshold sweep (the key check)
    print("  aligned sweep (same canonical sensitivity):")
    for q in [0.90, 0.95, 0.98]:
        k = int(round((1 - q) * len(can)))
        k = max(0, min(len(can) - 1, k))
        thr = can_sorted[k] if high else can_sorted[-1 - k]
        ts = (tb > thr).mean() if high else (tb < thr).mean()
        print(f"    canon-sens={q:.2f}  TypeB sens={ts:.3f}")
