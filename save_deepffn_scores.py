#!/usr/bin/env python3
"""Score the experimental IC50 set with the single-allele MHCflurry-trained
Deep_FFN.keras and save alongside the per-allele scores."""
import json, numpy as np, keras
from experimental_label_benchmark import encode

d = json.load(open("experimental_label_scores.json", encoding="utf-8"))
binders = list(d["binders"]); non = list(d["nonbinders"])

m = keras.models.load_model("03_Analysis/models/Deep_FFN.keras", compile=False)
def score(peps):
    X = np.array([encode(p) for p in peps], dtype=np.float32)
    probs = m.predict(X, verbose=0)
    return probs[:, 1] + probs[:, 2]

sb = score(binders); sn = score(non)
out = {"binders": {p: round(float(v), 4) for p, v in zip(binders, sb)},
       "nonbinders": {p: round(float(v), 4) for p, v in zip(non, sn)}}
json.dump(out, open("experimental_label_scores_deepffn.json", "w", encoding="utf-8"), indent=1)
print(f"saved experimental_label_scores_deepffn.json  binders={len(binders)} nonbinders={len(non)}")
