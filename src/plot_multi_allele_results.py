#!/usr/bin/env python3
"""Generate 5-allele results figure incl. pan vs per-allele non-canonical comparison."""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

CANON = {
    "HLA-A*02:01": (1.0, 0.97, 0.943),
    "HLA-A*01:01": (1.0, 0.99, 0.980),
    "HLA-A*03:01": (1.0, 0.97, 0.943),
    "HLA-B*07:02": (1.0, 0.97, 0.943),
    "HLA-B*44:03": (1.0, 0.99, 0.980),
}
NONCAN_PA = {  # per-allele non-canonical @0.95
    "HLA-A*02:01": (1.0, 0.806, 0.893),
    "HLA-A*01:01": (0.88, 0.957, 0.917),
    "HLA-A*03:01": (1.0, 0.962, 0.980),
    "HLA-B*07:02": (1.0, 0.862, 0.926),
    "HLA-B*44:03": (0.96, 0.960, 0.960),
}
NONCAN_PAN = {  # pan model non-canonical @0.95
    "HLA-A*02:01": (0.80, 0.833, 0.816),
    "HLA-A*01:01": (0.76, 0.950, 0.844),
    "HLA-A*03:01": (0.48, 1.000, 0.649),
    "HLA-B*07:02": (0.96, 0.923, 0.941),
    "HLA-B*44:03": (1.00, 0.962, 0.980),
}
alleles = list(CANON.keys())
labels = [a.replace("HLA-", "").replace("*", "") for a in alleles]
x = np.arange(len(alleles)); w = 0.24

fig, axes = plt.subplots(1, 2, figsize=(12, 4.2))
# Panel 1: canonical sensitivity + specificity
for ax, m, title in zip([axes[0]], range(1), ["Canonical benchmark metrics"]):
    c_sens = [CANON[a][0] for a in alleles]
    c_spec = [CANON[a][1] for a in alleles]
    ax.bar(x - w/2, c_sens, w, label="Sensitivity", color="#2b6cb0")
    ax.bar(x + w/2, c_spec, w, label="Specificity", color="#63b3ed")
    ax.axhline(0.90, color="gray", ls="--", lw=0.8)
    ax.set_ylim(0, 1.05); ax.set_title(title, fontsize=11)
    ax.set_xticks(x); ax.set_xticklabels(labels, rotation=20, ha="right", fontsize=8)
    ax.legend(fontsize=7)

# Panel 2: non-canonical sensitivity — per-allele vs pan
pa = [NONCAN_PA[a][0] for a in alleles]
pn = [NONCAN_PAN[a][0] for a in alleles]
ax = axes[1]
ax.bar(x - w/2, pa, w, label="Per-allele models", color="#38a169")
ax.bar(x + w/2, pn, w, label="Pan-specific model", color="#d69e2e")
ax.axhline(0.90, color="gray", ls="--", lw=0.8)
ax.set_ylim(0, 1.05); ax.set_title("Non-canonical binder sensitivity (IEDB, @0.95)", fontsize=11)
ax.set_xticks(x); ax.set_xticklabels(labels, rotation=20, ha="right", fontsize=8)
ax.legend(fontsize=7)

fig.suptitle("Multi-Allele Deep FFN: canonical pass 5/5; per-allele > pan for non-canonical", fontsize=12)
fig.tight_layout()
os.makedirs("to_human", exist_ok=True)
fig.savefig("to_human/multi_allele_results.png", dpi=150)
print("saved to_human/multi_allele_results.png")
