#!/usr/bin/env python3
"""Comprehensive findings panel — 4 key results in one figure.

Panels:
  A. Temporal holdout: pre-2021 vs post-2022 sensitivity (3 alleles)
  B. Pan16: per-allele non-canonical sensitivity (15-allele panel)
  C. P2 virtual mutagenesis: sidechain reachability to Glu63
  D. Predictor plateau: canonical vs Type B sensitivity (DeepFFN/MHCflurry/MUNIS)

Usage: python src/plot_findings_panel.py  (base anaconda, matplotlib)
"""
import json, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

os.makedirs("to_human", exist_ok=True)

fig, axes = plt.subplots(2, 2, figsize=(13, 9))

# ---- Panel A: temporal holdout ----
ax = axes[0, 0]
alleles = ["HLA-A*02:01", "HLA-A*03:01", "HLA-B*07:02"]
labels = [a.replace("HLA-", "").replace("*", "") for a in alleles]
pre = [0.879, 0.905, 0.898]
post = [0.824, 0.936, 0.973]
x = np.arange(len(alleles)); w = 0.3
ax.bar(x - w/2, pre, w, label="pre-2021 (in-distribution)", color="#718096")
ax.bar(x + w/2, post, w, label="post-2022 (unseen by labeler)", color="#2b6cb0")
ax.set_xticks(x); ax.set_xticklabels(labels)
ax.set_ylabel("Sensitivity @0.95")
ax.set_title("A. Temporal holdout — generalizes to new era", fontsize=11)
ax.set_ylim(0, 1.1); ax.legend(fontsize=8)
ax.axhline(0.9, color="gray", ls="--", lw=0.8)

# ---- Panel B: pan16 per-allele ----
ax = axes[0, 1]
# from eval_pan16.py output
pan16 = {"A*02:01":0.880,"A*01:01":0.800,"A*03:01":1.000,"B*07:02":0.640,"B*44:03":1.000,
         "A*11:01":1.000,"A*24:02":0.760,"B*08:01":0.880,"B*15:01":0.960,"B*27:05":0.080,
         "B*35:01":0.840,"B*40:01":1.000,"B*49:01":1.000,"B*51:01":0.520,"B*57:01":0.920}
ks = list(pan16.keys()); vs = list(pan16.values())
order = np.argsort(vs)
ks = [ks[i] for i in order]; vs = [vs[i] for i in order]
colors = ["#c53030" if v < 0.6 else "#2b6cb0" for v in vs]
ax.barh(np.arange(len(ks)), vs, color=colors)
ax.set_yticks(np.arange(len(ks))); ax.set_yticklabels(ks, fontsize=8)
ax.set_xlabel("Non-canonical sensitivity @0.95")
ax.set_title("B. Pan16 — one model, 15-allele panel (pooled 0.819)", fontsize=11)
ax.set_xlim(0, 1.1)
ax.axvline(0.819, color="gray", ls="--", lw=1.5)
ax.text(0.83, 14.5, "pooled 0.819", fontsize=9, color="gray")

# ---- Panel C: P2 virtual mutagenesis ----
ax = axes[1, 0]
scan = json.load(open("06_Logs/p2_virtual_scan.json"))
results = scan["results"]
# reachability vs reach length
ress = list(results.keys()); reaches = [results[r]["reach"] for r in ress]
reachable = [results[r]["reachable"] for r in ress]
colors = ["#38a169" if results[r]["reachable"] and r in "KRHDE" else
          ("#2b6cb0" if results[r]["reachable"] else "#e2e8f0") for r in ress]
order = np.argsort(reaches)
ress2 = [ress[i] for i in order]; reaches2 = [reaches[i] for i in order]; colors2 = [colors[i] for i in order]
ax.bar(np.arange(len(ress2)), reaches2, color=colors2)
ax.set_xticks(np.arange(len(ress2))); ax.set_xticklabels(ress2, fontsize=8)
ax.set_ylabel("Sidechain reach (Å)")
ax.axhline(3.83, color="red", ls="--", lw=1)
ax.text(0.5, 4.0, "P2CA→Glu63 = 3.83Å", fontsize=9, color="red")
ax.set_title("C. P2 residues that can reach Glu63 (green=charged)", fontsize=11)
ax.set_ylim(0, 8.5)

# ---- Panel D: predictor plateau ----
ax = axes[1, 1]
types = ["Canonical", "Type B\n(non-canonical)"]
dff = [0.979, 0.607]; mhc = [0.986, 0.668]; munis = [0.950, 0.647]
x = np.arange(2); w = 0.25
ax.bar(x - w, dff, w, label="Deep FFN", color="#2b6cb0")
ax.bar(x, mhc, w, label="MHCflurry", color="#e53e3e")
ax.bar(x + w, munis, w, label="MUNIS (structure)", color="#d69e2e")
ax.set_xticks(x); ax.set_xticklabels(types)
ax.set_ylabel("Sensitivity")
ax.set_title("D. Non-canonical plateau across predictor classes", fontsize=11)
ax.set_ylim(0, 1.1); ax.legend(fontsize=8)
ax.axhline(0.65, color="gray", ls="--", lw=0.8)
ax.text(1.25, 0.67, "plateau ~0.65", fontsize=9, color="gray")

fig.suptitle("Peptide-MHC Multi-Allele Research — Validated Findings", fontsize=14, y=0.99)
fig.tight_layout(rect=[0, 0, 1, 0.97])
fig.savefig("to_human/findings_panel.png", dpi=150)
fig.savefig("04_Manuscript/companion_multi_allele/fig1_findings_panel.pdf")
print("saved to_human/findings_panel.png")
