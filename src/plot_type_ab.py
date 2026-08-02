#!/usr/bin/env python3
"""Type A/B large-sample summary figure (base python)."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# Match R figure font (redraw_all_nature.R uses family='sans' = Arial on Windows)
plt.rcParams["font.family"] = "Arial"

types = ["canonical", "Type A\n(tolerated)", "Type B\n(non-canonical)"]
dff = [0.979, 0.477, 0.607]
mhc = [0.986, 0.669, 0.668]
dff_lo = [0.976, 0.399, 0.590]; dff_hi = [0.982, 0.556, 0.624]
mhc_lo = [0.984, 0.590, 0.652]; mhc_hi = [0.988, 0.739, 0.685]

x = np.arange(3); w = 0.3
fig, ax = plt.subplots(figsize=(8, 5))
ax.bar(x-w/2, dff, w, yerr=[np.array(dff)-np.array(dff_lo), np.array(dff_hi)-np.array(dff)],
       capsize=4, label="Deep FFN (ours)", color="#377EB8")
ax.bar(x+w/2, mhc, w, yerr=[np.array(mhc)-np.array(mhc_lo), np.array(mhc_hi)-np.array(mhc)],
       capsize=4, label="MHCflurry", color="#E41A1C")
ax.set_xticks(x); ax.set_xticklabels(types)
ax.set_ylabel("Sensitivity on IEDB positives")
ax.set_ylim(0, 1.1)
ax.set_title("Non-canonical binding penalty is moderate, not catastrophic (n=14,488)")
ax.legend()
ax.axhline(0.6, color="gray", ls="--", lw=0.8)
fig.tight_layout()
fig.savefig("to_human/type_ab_largesample.png", dpi=150)
fig.savefig("to_human/type_ab_largesample.pdf")
print("saved to_human/type_ab_largesample.png")
