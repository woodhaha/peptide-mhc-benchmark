#!/usr/bin/env python3
"""Regenerate Figure 9 — non-canonical (Type B) sensitivity at matched canonical
operating points (0.90 / 0.95 / 0.98).

Shows that the apparent ~0.65 "plateau" is a threshold artifact: aligned to the
same canonical sensitivity, the four predictors diverge (0.37–0.77). Data from
02_Data/cleaned/threshold_sensitivity_results.txt. Publication-style matplotlib.
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

plt.rcParams.update({
    "font.family": "Arial",
    "font.size": 9,
    "axes.linewidth": 0.8,
    "axes.labelsize": 9,
    "axes.titlesize": 9,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "legend.fontsize": 8,
    "figure.dpi": 350,
})

# matched canonical operating points -> Type B sensitivity (from threshold_sensitivity_results.txt)
op_points = [0.90, 0.95, 0.98]
data = {
    "Deep FFN":  [0.369, 0.454, 0.610],
    "MHCflurry": [0.385, 0.497, 0.628],
    "MUNIS":     [0.521, 0.647, 0.773],
    "BigMHC":    [0.555, 0.653, 0.740],
}

# semantic palette (gray-safe)
colors = {
    "Deep FFN":  "#4C72B0",
    "MHCflurry": "#DD8452",
    "MUNIS":     "#55A868",
    "BigMHC":    "#C44E52",
}
markers = {"Deep FFN": "o", "MHCflurry": "s", "MUNIS": "^", "BigMHC": "D"}

fig, ax = plt.subplots(figsize=(4.8, 3.4))

for name, vals in data.items():
    ax.plot(op_points, vals, marker=markers[name], ms=4.5, lw=1.5,
            color=colors[name], label=name, zorder=3)

ax.set_xlabel("Matched canonical sensitivity (operating point)")
ax.set_ylabel("Type B (non-canonical) sensitivity")
ax.set_xlim(0.89, 0.99)
ax.set_ylim(0.30, 0.85)
ax.set_xticks(op_points)
ax.set_yticks(np.arange(0.3, 0.9, 0.1))
ax.legend(frameon=False, loc="upper left")

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

fig.tight_layout()
for out in ["06_Bioinformatics_Submission/figures/Figure9.pdf",
            "06_Bioinformatics_Submission/latex_template/figures/Figure9.pdf"]:
    fig.savefig(out)
    print("saved", out)
