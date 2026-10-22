#!/usr/bin/env python3
"""Regenerate Figure 9 (predictor plateau) with BigMHC added.

Canonical vs Type B sensitivity across 4 predictors; dashed line marks
the ~0.65 Type B plateau. Publication-style matplotlib (figures4papers).
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# figures4papers-style rcParams
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

predictors = ["Deep FFN", "MHCflurry", "MUNIS", "BigMHC"]
canonical = [0.979, 0.986, 0.950, 0.950]
type_b = [0.607, 0.668, 0.647, 0.660]

x = np.arange(len(predictors))
w = 0.35

fig, ax = plt.subplots(figsize=(4.5, 3.2))

# semantic palette, gray-safe
c_canon = "#4C72B0"   # blue
c_typeb = "#DD8452"   # orange

b1 = ax.bar(x - w/2, canonical, w, label="Canonical", color=c_canon)
b2 = ax.bar(x + w/2, type_b, w, label="Type B (non-canonical)", color=c_typeb)

ax.axhline(0.65, ls="--", lw=0.9, color="gray", zorder=3)
ax.text(len(predictors) - 0.5, 0.665, "~0.65 plateau", fontsize=7, color="gray", ha="right")

ax.set_ylabel("Sensitivity")
ax.set_xticks(x)
ax.set_xticklabels(predictors)
ax.set_ylim(0, 1.05)
ax.set_yticks(np.arange(0, 1.05, 0.2))
ax.legend(frameon=False, loc="lower left")

# value labels
for bars in (b1, b2):
    for r in bars:
        ax.text(r.get_x() + r.get_width()/2, r.get_height() + 0.02,
                f"{r.get_height():.2f}", ha="center", va="bottom", fontsize=7)

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

fig.tight_layout()
for out in ["06_Bioinformatics_Submission/figures/Figure9.pdf",
            "06_Bioinformatics_Submission/latex_template/figures/Figure9.pdf"]:
    fig.savefig(out)
    print("saved", out)
