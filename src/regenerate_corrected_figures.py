#!/usr/bin/env python3
"""Regenerate the headline figures to match the corrected results.

Figure 1 (architecture): CNN modestly best, 5-seed mean ± SD.
Figure 2 (label source): MHCflurry > PSSM > random (chance ~33%).
Companion Figure 1 (pan vs per-allele): 0.84-0.96 with B*07:02 exception.

Style: Nature palette + Arial font (matches the manuscript's figure spec).
"""
import json, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import rcParams

# Nature palette
BLUE, RED, ORANGE, GREEN, GRAY = "#377EB8", "#E41A1C", "#FF7F00", "#4DAF4A", "#999999"
for f in ("Arial", "Helvetica"):
    try:
        rcParams["font.family"] = f
        break
    except Exception:
        pass
rcParams["axes.linewidth"] = 0.8
rcParams["font.size"] = 10

OUT = "03_Analysis/figures"

def load_json(p):
    return json.load(open(p, encoding="utf-8"))

def save(fig, name):
    for ext in ("png", "pdf"):
        fig.savefig(f"{OUT}/{name}.{ext}", dpi=300, bbox_inches="tight")
    print(f"  saved {name}.png/.pdf")

# ---- Figure 1: architecture comparison ----
arch = load_json("02_Data/cleaned/multiseed_results.json")
order = ["CNN", "Deep FFN", "FFN (Jessen)", "LSTM", "ResNet", "Random Forest"]
means = [arch[m]["mean"] * 100 for m in order]
sds = [arch[m]["sd"] * 100 for m in order]
labels = [m if m != "FFN (Jessen)" else "FFN" for m in order]

fig, ax = plt.subplots(figsize=(5.5, 3.6))
y = np.arange(len(order))[::-1]
ax.barh(y, means, xerr=sds, color=[BLUE] + [GRAY] * 5, capsize=3, height=0.6)
ax.set_yticks(y); ax.set_yticklabels(labels)
ax.set_xlabel("Test accuracy (%)"); ax.set_xlim(70, 100)
ax.set_title("Architecture comparison (5-seed mean ± SD, n=519)")
for yi, m, s in zip(y, means, sds):
    ax.text(m + s + 0.4, yi, f"{m:.1f}", va="center", fontsize=9)
save(fig, "Figure1_model_comparison_corrected"); plt.close(fig)

# ---- Figure 2: label source ----
rand = load_json("02_Data/cleaned/random_label_results.json")
deep = rand["Deep FFN"]
lbl = [("MHCflurry 2.2.1", 90.56, 1.16, BLUE),
       ("PSSM", 82.6, None, ORANGE),
       ("Random (i.i.d.)", deep["mean"] * 100, deep["sd"] * 100, RED)]
fig, ax = plt.subplots(figsize=(5.0, 3.6))
x = np.arange(3)
for i, (name, m, s, c) in enumerate(lbl):
    ax.bar(i, m, yerr=s if s else None, color=c, capsize=4, width=0.55, label=name)
ax.axhline(33.3, ls="--", color=GRAY, lw=1, label="chance (3-class)")
ax.set_xticks(x); ax.set_xticklabels([l[0] for l in lbl])
ax.set_ylabel("Deep FFN accuracy (%)"); ax.set_ylim(0, 100)
ax.set_title("Label-source controls (mean ± SD)")
ax.legend(frameon=False, fontsize=8)
for i, (name, m, s, c) in enumerate(lbl):
    ax.text(i, m + (s or 0) + 2, f"{m:.1f}", ha="center", fontsize=9)
save(fig, "Figure2_label_comparison_corrected"); plt.close(fig)

# ---- Companion Figure 1: pan vs per-allele ----
reeval = load_json("data/reeval_corrected.json")
alleles = ["HLA-A*02:01", "HLA-A*01:01", "HLA-A*03:01", "HLA-B*07:02", "HLA-B*44:03"]
pa = [reeval[a]["results"]["per_allele"]["0.95"]["sens"] for a in alleles]
pan = [reeval[a]["results"]["pan_full"]["0.95"]["sens"] for a in alleles]
fig, ax = plt.subplots(figsize=(5.5, 3.4))
x = np.arange(len(alleles)); w = 0.38
ax.bar(x - w/2, pa, w, color=BLUE, label="per-allele")
ax.bar(x + w/2, pan, w, color=RED, label="pan")
ax.set_xticks(x); ax.set_xticklabels([a.replace("HLA-", "") for a in alleles], fontsize=9)
i = alleles.index("HLA-B*07:02")
ax.set_ylabel("Non-canonical sensitivity"); ax.set_ylim(0, 1.15)
ax.set_title("Pan vs per-allele (protein-background negatives, n=25)")
ax.legend(frameon=False, fontsize=8)
# highlight B*07:02 exception
ax.annotate("B*07:02\nexception", xy=(i + w/2, pan[i]), xytext=(i + w/2 + 0.6, pan[i] + 0.1),
            fontsize=8, color=RED, arrowprops=dict(arrowstyle="->", color=RED, lw=0.8))
save(fig, "Companion_Figure1_pan_vs_perallele_corrected"); plt.close(fig)

print("done")
