#!/usr/bin/env python3
"""Figure: rebalanced pan non-canonical sensitivity across 5 alleles (5-seed).

Headline positive result of the companion note. Shows only the seed-robust
rebalanced pan model (multiseed_results.json) — per-allele bars are NOT
overlaid because the per-allele benchmark (n=25 broad-motif) differs from the
large non-canonical set, and mixing them would be an apples-to-oranges
comparison. The per-allele parity claim is stated in text with its
same-benchmark source.

Usage: python src/plot_pan_vs_perallele.py  (base anaconda, matplotlib)
"""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

os.makedirs("04_Manuscript/companion_multi_allele", exist_ok=True)

# Verified from multiseed_results.json (rebalanced pan, B*07:02 x4)
pan = {
    "HLA-A*02:01": (0.984, 0.020),
    "HLA-A*01:01": (0.920, 0.044),
    "HLA-A*03:01": (0.976, 0.020),
    "HLA-B*07:02": (0.920, 0.051),
    "HLA-B*44:03": (1.000, 0.000),
}

alleles = list(pan.keys())
labels = [a.replace("HLA-", "").replace("*", "") for a in alleles]
means = [pan[a][0] for a in alleles]
sds = [pan[a][1] for a in alleles]

fig, ax = plt.subplots(figsize=(8, 5))
x = np.arange(len(alleles))
bars = ax.bar(x, means, 0.55, yerr=sds, capsize=5,
              color="#2b6cb0", edgecolor="black", linewidth=0.5)
for xi, m, sd in zip(x, means, sds):
    ax.text(xi, m + sd + 0.02, f"{m:.2f}±{sd:.2f}", ha="center", fontsize=9)

ax.set_xticks(x); ax.set_xticklabels(labels)
ax.set_ylabel("Non-canonical sensitivity @0.95")
ax.set_ylim(0, 1.15)
ax.axhline(0.9, color="gray", ls="--", lw=0.8)
ax.text(4.55, 0.915, "0.9 target", fontsize=8, color="gray", ha="right")
ax.set_title("Rebalanced pan model — non-canonical sensitivity, 5-seed mean ± SD",
             fontsize=11)
fig.tight_layout()
out = "04_Manuscript/companion_multi_allele/fig2_pan_vs_perallele.png"
fig.savefig(out, dpi=300)
fig.savefig("04_Manuscript/companion_multi_allele/fig2_pan_vs_perallele.pdf")
print(f"saved {out}")
