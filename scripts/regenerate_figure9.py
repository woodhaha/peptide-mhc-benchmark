#!/usr/bin/env python3
"""Regenerate Figure 9 — cross-predictor rank-correlation matrix.

Shows that predictors sharing an objective correlate highly (binding-binding
0.96, presentation-presentation 0.90) while cross-objective correlation is
moderate (0.72-0.78). |Spearman rho| is shown; cross-objective signs are
inverted because presentation score and binding rank run in opposite directions.
"""
import csv
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import spearmanr

plt.rcParams.update({
    "font.family": "Arial", "font.size": 9, "axes.linewidth": 0.8,
    "figure.dpi": 350,
})

mhc = {r['peptide']: float(r['pct_rank']) for r in csv.DictReader(open('02_Data/raw/real_peptides.csv')) if r['data_type']=='test'}
nmc = {r['peptide']: float(r['netmhcpan_rank']) for r in csv.DictReader(open('02_Data/cleaned/netmhcpan_results.csv'))}
munis = {r['pep']: float(r['score']) for r in csv.DictReader(open('.munis_dl/predictions/pep519_munis_predictions.csv'))}
bigmhc = {}
for line in open('bigmhc_519.prd'):
    p = line.strip().split(',')
    if len(p) < 5 or p[0]=='mhc': continue
    bigmhc[p[1]] = float(p[4])

common = [p for p in mhc if p in nmc and p in munis and p in bigmhc]
names = ["MHCflurry", "NetMHCpan", "MUNIS", "BigMHC"]
data = {
    "MHCflurry": np.array([mhc[p] for p in common]),
    "NetMHCpan": np.array([nmc[p] for p in common]),
    "MUNIS": np.array([munis[p] for p in common]),
    "BigMHC": np.array([bigmhc[p] for p in common]),
}
n = len(names)
mat = np.zeros((n, n))
for i in range(n):
    for j in range(n):
        mat[i, j] = abs(spearmanr(data[names[i]], data[names[j]])[0])

fig, ax = plt.subplots(figsize=(3.8, 3.2))
im = ax.imshow(mat, cmap="viridis", vmin=0.70, vmax=1.0, origin="lower")
ax.set_xticks(range(n)); ax.set_xticklabels(names, rotation=45, ha="right")
ax.set_yticks(range(n)); ax.set_yticklabels(names)
for i in range(n):
    for j in range(n):
        ax.text(j, i, f"{mat[i,j]:.2f}", ha="center", va="center",
                color="white" if mat[i,j] < 0.85 else "black", fontsize=8)
# outline the two within-objective blocks
for rect in [(0.5, 0.5), (2.5, 2.5)]:
    ax.add_patch(plt.Rectangle((rect[0]-1, rect[1]-1), 2, 2,
                 fill=False, edgecolor="white", lw=1.5))
ax.set_title("|Spearman ρ| across predictors", fontsize=9, pad=8)
cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
cbar.set_label("|Spearman ρ|", fontsize=8)
ax.spines[:].set_visible(False)
fig.tight_layout()
for out in ["06_Bioinformatics_Submission/figures/Figure9.pdf",
            "06_Bioinformatics_Submission/latex_template/figures/Figure9.pdf"]:
    fig.savefig(out)
    print("saved", out)
