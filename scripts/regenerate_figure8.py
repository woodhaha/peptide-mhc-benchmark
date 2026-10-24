#!/usr/bin/env python3
"""Regenerate Figure 8 — cutpoint sensitivity (the "free knob").

Panel (a): binary binder agreement between MHCflurry and NetMHCpan as a
function of a shared percentile-rank threshold (the agreement dips at the
standard 2% cutpoint). Panel (b): agreement with MHCflurry's 2% binder call
as NetMHCpan's threshold is varied, showing the ~4x calibration offset.
"""
import csv
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

plt.rcParams.update({
    "font.family": "Arial", "font.size": 9, "axes.linewidth": 0.8,
    "axes.labelsize": 9, "axes.titlesize": 9, "xtick.labelsize": 8,
    "ytick.labelsize": 8, "figure.dpi": 350,
})

# --- recompute from data ---
mhc = {}
for r in csv.DictReader(open('02_Data/raw/real_peptides.csv')):
    if r['data_type'] == 'test':
        mhc[r['peptide']] = float(r['pct_rank'])
nmc = {}
for r in csv.DictReader(open('02_Data/cleaned/netmhcpan_results.csv')):
    nmc[r['peptide']] = float(r['netmhcpan_rank'])
common = [p for p in mhc if p in nmc]
a = [mhc[p] for p in common]; b = [nmc[p] for p in common]
n = len(common)

shared_ts = [0.5, 1.0, 2.0, 5.0, 10.0, 20.0]
shared_agree = [sum((x < t) == (y < t) for x, y in zip(a, b)) / n for t in shared_ts]
offset_ts = [0.5, 1.0, 2.0, 3.0, 5.0, 8.0, 12.0]
offset_agree = [sum((mhc[p] < 2.0) == (nmc[p] < t) for p in common) / n for t in offset_ts]

fig, axes = plt.subplots(1, 2, figsize=(7.0, 3.0))

# (a) agreement vs shared threshold
ax = axes[0]
ax.plot(shared_ts, shared_agree, marker='o', ms=4, lw=1.5, color='#4C72B0')
ax.axvline(2.0, ls='--', lw=0.8, color='gray')
ax.text(2.1, 0.72, 'standard 2%', fontsize=7, color='gray')
ax.set_xlabel('Shared binder threshold (% rank)')
ax.set_ylabel('Binary binder agreement')
ax.set_ylim(0.70, 0.95)
ax.set_xscale('log')
ax.set_xticks(shared_ts); ax.set_xticklabels([str(t) for t in shared_ts])
ax.set_title('(a)', loc='left', fontsize=9)

# (b) calibration offset
ax = axes[1]
ax.plot(offset_ts, offset_agree, marker='s', ms=4, lw=1.5, color='#DD8452')
ax.axhline(0.961, ls='--', lw=0.8, color='gray')
ax.axvline(8.0, ls=':', lw=0.8, color='#DD8452')
ax.text(8.1, 0.77, 'MHCflurry 2% ≈ NetMHCpan 8%', fontsize=7, color='#DD8452')
ax.set_xlabel('NetMHCpan threshold (% rank)')
ax.set_ylabel('Agreement with MHCflurry 2% binder call')
ax.set_ylim(0.50, 1.0)
ax.set_xscale('log')
ax.set_xticks([0.5, 1, 2, 3, 5, 8, 12]); ax.set_xticklabels(['0.5','1','2','3','5','8','12'])
ax.set_title('(b)', loc='left', fontsize=9)

for ax in axes:
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

fig.tight_layout()
for out in ["06_Bioinformatics_Submission/figures/Figure8.pdf",
            "06_Bioinformatics_Submission/latex_template/figures/Figure8.pdf"]:
    fig.savefig(out)
    print("saved", out)
