#!/usr/bin/env python3
"""Companion manuscript figures (Nature palette + Arial, consistent with main paper).

Fig 2: Predictor plateau — canonical vs Type B sensitivity (DFFN/MHCflurry/MUNIS)
Fig 3: Pan16 per-allele non-canonical sensitivity (15-allele panel)
Fig 4: Temporal holdout — pre-2021 vs post-2022 sensitivity

Usage: python src/plot_companion_figures.py  (base anaconda, matplotlib)
"""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

plt.rcParams["font.family"] = "Arial"
os.makedirs("06_Bioinformatics_Submission/companion_figures", exist_ok=True)

# ---- Fig 2: predictor plateau ----
fig, ax = plt.subplots(figsize=(6, 4.5))
types = ["Canonical", "Type B\n(non-canonical)"]
dff = [0.979, 0.607]
mhc = [0.986, 0.668]
munis = [0.950, 0.647]
x = np.arange(2); w = 0.25
ax.bar(x - w, dff, w, label="Deep FFN", color="#377EB8")
ax.bar(x, mhc, w, label="MHCflurry", color="#E41A1C")
ax.bar(x + w, munis, w, label="MUNIS (structure)", color="#FF7F00")
ax.set_xticks(x); ax.set_xticklabels(types)
ax.set_ylabel("Sensitivity")
ax.set_ylim(0, 1.1); ax.legend(fontsize=8)
ax.axhline(0.65, color="gray", ls="--", lw=0.8)
ax.text(1.45, 0.67, "plateau ~0.65", fontsize=8, color="gray", ha="right")
ax.set_title("Non-canonical plateau across predictor classes", fontsize=10)
fig.tight_layout()
fig.savefig("06_Bioinformatics_Submission/companion_figures/Fig2_plateau.png", dpi=300)
fig.savefig("06_Bioinformatics_Submission/companion_figures/Fig2_plateau.pdf")
plt.close(fig)
print("Fig2 saved")

# ---- Fig 3: pan16 per-allele ----
fig, ax = plt.subplots(figsize=(7, 5))
pan16 = {"A*02:01":0.88,"A*01:01":0.80,"A*03:01":1.00,"B*07:02":0.64,"B*44:03":1.00,
         "A*11:01":1.00,"A*24:02":0.76,"B*08:01":0.88,"B*15:01":0.96,"B*27:05":0.08,
         "B*35:01":0.84,"B*40:01":1.00,"B*49:01":1.00,"B*51:01":0.52,"B*57:01":0.92}
ks = list(pan16.keys()); vs = list(pan16.values())
order = np.argsort(vs)
ks = [ks[i] for i in order]; vs = [vs[i] for i in order]
colors = ["#E41A1C" if v < 0.6 else "#377EB8" for v in vs]
ax.barh(np.arange(len(ks)), vs, color=colors)
ax.set_yticks(np.arange(len(ks))); ax.set_yticklabels(ks, fontsize=8)
ax.set_xlabel("Non-canonical sensitivity @0.95")
ax.set_xlim(0, 1.1)
ax.axvline(0.819, color="gray", ls="--", lw=1.5)
ax.text(0.83, 14.3, "pooled 0.819", fontsize=9, color="gray")
ax.set_title("Pan16 — one model, 15-allele panel", fontsize=10)
fig.tight_layout()
fig.savefig("06_Bioinformatics_Submission/companion_figures/Fig3_pan16.png", dpi=300)
fig.savefig("06_Bioinformatics_Submission/companion_figures/Fig3_pan16.pdf")
plt.close(fig)
print("Fig3 saved")

# ---- Fig 4: temporal holdout ----
fig, ax = plt.subplots(figsize=(6, 4.5))
alleles = ["HLA-A*02:01", "HLA-A*03:01", "HLA-B*07:02"]
labels = [a.replace("HLA-", "").replace("*", "") for a in alleles]
pre = [0.879, 0.905, 0.898]
post = [0.824, 0.936, 0.973]
x = np.arange(len(alleles)); w = 0.3
ax.bar(x - w/2, pre, w, label="pre-2021 (in-distribution)", color="#999999")
ax.bar(x + w/2, post, w, label="post-2022 (unseen by labeler)", color="#377EB8")
ax.set_xticks(x); ax.set_xticklabels(labels)
ax.set_ylabel("Sensitivity @0.95")
ax.set_ylim(0, 1.1); ax.legend(fontsize=8)
ax.axhline(0.9, color="gray", ls="--", lw=0.8)
ax.set_title("Temporal holdout — generalizes to new era", fontsize=10)
fig.tight_layout()
fig.savefig("06_Bioinformatics_Submission/companion_figures/Fig4_temporal.png", dpi=300)
fig.savefig("06_Bioinformatics_Submission/companion_figures/Fig4_temporal.pdf")
plt.close(fig)
print("Fig4 saved")
