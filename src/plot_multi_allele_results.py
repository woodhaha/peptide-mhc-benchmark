#!/usr/bin/env python3
"""Generate the 5-allele multi-allele results figure (runs under base anaconda)."""
import json, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# Canonical benchmark results (from 06_Logs/multi_allele_results.md)
CANON = {
    "HLA-A*02:01": (1.0, 0.97, 0.943),
    "HLA-A*01:01": (1.0, 0.99, 0.980),
    "HLA-A*03:01": (1.0, 0.97, 0.943),
    "HLA-B*07:02": (1.0, 0.97, 0.943),
    "HLA-B*44:03": (1.0, 0.99, 0.980),
}
# Non-canonical sens @0.95 (from results log)
NONCAN = {
    "HLA-A*02:01": (1.0, 0.806, 0.893),
    "HLA-A*01:01": (0.88, 0.957, 0.917),
    "HLA-A*03:01": (1.0, 0.962, 0.980),
    "HLA-B*07:02": (1.0, 0.862, 0.926),
    "HLA-B*44:03": (0.96, 0.960, 0.960),
}
alleles = list(CANON.keys())
labels = [a.replace("HLA-", "").replace("*", "") for a in alleles]
x = np.arange(len(alleles)); w = 0.28

fig, axes = plt.subplots(1, 3, figsize=(13, 4.2))
metrics = ["sensitivity", "specificity_raw", "f1"]
for ax, m, title in zip(axes, range(3), ["Sensitivity", "Specificity", "F1"]):
    c = [CANON[a][m] for a in alleles]
    n = [NONCAN[a][m] for a in alleles]
    ax.bar(x - w/2, c, w, label="Canonical", color="#2b6cb0")
    ax.bar(x + w/2, n, w, label="Non-canonical", color="#e53e3e")
    ax.axhline(0.85, color="gray", ls="--", lw=0.8)
    ax.set_title(title, fontsize=11)
    ax.set_xticks(x); ax.set_xticklabels(labels, rotation=20, ha="right", fontsize=8)
    ax.set_ylim(0, 1.05)
    ax.legend(fontsize=7)
fig.suptitle("Multi-Allele Deep FFN: 5/5 alleles pass (canonical + non-canonical robustness)", fontsize=12)
fig.tight_layout()
os.makedirs("to_human", exist_ok=True)
fig.savefig("to_human/multi_allele_results.png", dpi=150)
print("saved to_human/multi_allele_results.png")
