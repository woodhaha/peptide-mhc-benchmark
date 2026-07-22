#!/usr/bin/env python3
"""K-E63 Salt Bridge (1DUZ) — RMSF analysis plot.
Reference: D:/Researching/SMVT/MM-GBSA Analysis/_plots_123.py (flexibility map)
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path

plt.rcParams.update({'font.size': 9, 'figure.dpi': 150, 'savefig.bbox': 'tight'})

DATA = Path("D:/Researching/Peptide epitope/md results/salt_bridge")
OUT  = Path("D:/Researching/Peptide epitope/md pic")

d = np.load(str(DATA / "rmsf.npz"))
rmsf = d['rmsf']; chains = d['chains']; resnames = d['resnames']; resnums = d['resnums']
n_res = len(rmsf)

# Build x-axis: absolute residue index across chains
abs_idx = np.arange(n_res)

# Chain boundaries
ch_bounds = {}
for i in range(n_res):
    ch = int(chains[i])
    if ch not in ch_bounds:
        ch_bounds[ch] = [i, i]
    ch_bounds[ch][1] = i

# ============================================================
# FIGURE 1: Per-residue RMSF bar chart
# ============================================================
fig, axes = plt.subplots(3, 1, figsize=(16, 8), sharex=True,
    gridspec_kw={'height_ratios': [2.5, 1, 0.8]})
fig.suptitle('1DUZ K-E63 Salt Bridge — Calpha RMSF (100ns MD)',
             fontsize=14, fontweight='bold', y=0.98)

colors = {0: '#3498db', 1: '#2ecc71', 2: '#e74c3c'}
labels = {0: 'Chain A (MHC)', 1: 'Chain B (b2m)', 2: 'Chain C (peptide)'}

# 1a. Per-residue RMSF
ax = axes[0]
for ch in sorted(colors):
    mask = chains == ch
    ax.bar(abs_idx[mask], rmsf[mask], width=0.8, color=colors[ch], alpha=0.7, label=labels[ch])
ax.axhline(y=2.0, color='#95a5a6', ls='--', lw=0.8, alpha=0.5, label='Flexibility threshold')
ax.set_ylabel('Calpha RMSF (A)', fontsize=10)
ax.legend(fontsize=8, ncol=3)
ax.grid(axis='y', alpha=0.15)

# Annotate E63 and K2
for label, ch, seq_id, note in [
    ("E63 (salt bridge)", 0, 63, 'active site'),
    ("LYS2 (P2 ligand)", 2, 2, 'binding partner')]:
    for i in range(n_res):
        if int(chains[i]) == ch and int(resnums[i]) == seq_id:
            ax.annotate(f'{label}\n{rmsf[i]:.2f}A',
                        xy=(i, rmsf[i]), xytext=(i, rmsf[i]+1.5),
                        fontsize=7, ha='center', fontweight='bold',
                        arrowprops=dict(arrowstyle='->', color='#2c3e50', lw=1),
                        bbox=dict(boxstyle='round,pad=0.3', fc='white', alpha=0.8))
            break

# 1b. Chain average bar
ax = axes[1]
ch_means = []
ch_stds = []
for ch in sorted(colors):
    m = rmsf[chains == ch]
    ch_means.append(m.mean())
    ch_stds.append(m.std())
x_pos = np.arange(len(ch_means))
bars = ax.bar(x_pos, ch_means, yerr=ch_stds, color=[colors[i] for i in sorted(colors)],
              width=0.5, edgecolor='white', linewidth=1, capsize=3)
ax.set_ylabel('Mean RMSF (A)', fontsize=10)
ax.set_xticks(x_pos)
ax.set_xticklabels([f'{labels[c]}\n({int((chains==c).sum())} residues)' for c in sorted(colors)],
    fontsize=8)
ax.grid(axis='y', alpha=0.15)

for bar, v in zip(bars, ch_means):
    ax.text(bar.get_x()+bar.get_width()/2, bar.get_height()+0.1,
            f'{v:.2f}A', ha='center', fontsize=9, fontweight='bold')

# 1c. Chain boundary + annotation
ax = axes[2]; ax.axis('off')
y = 0.5
for ch in sorted(ch_bounds):
    s, e = ch_bounds[ch]
    ax.fill_between([s-0.5, e+0.5], y-0.3, y+0.3, color=colors[ch], alpha=0.3)
    ax.text((s+e)/2, y, f'{labels[ch]}\n{s}-{e}', ha='center', va='center', fontsize=8)
    ax.axvline(x=s-0.5, color='#7f8c8d', ls=':', lw=0.5)

# Mark residue range
ax.set_xlim(-1, n_res+1)
ax.set_ylim(0, 1)
ax.set_xlabel('Residue Index (sequential)', fontsize=9)

plt.savefig(OUT / 'sb_rmsf.png')
plt.close()
print('1. sb_rmsf.png done')

# ============================================================
# FIGURE 2: RMSF vs residue with chain shading + flexible regions
# ============================================================
fig, ax = plt.subplots(figsize=(16, 5))
fig.suptitle('1DUZ Calpha RMSF — Chain-Colored with Flexible Regions Highlighted',
             fontsize=14, fontweight='bold')

x = np.arange(n_res)
for ch in sorted(colors):
    mask = chains == ch
    ax.scatter(abs_idx[mask], rmsf[mask], c=colors[ch], s=12, alpha=0.7, label=labels[ch],
               zorder=3)
    # running average within chain
    idx_ch = abs_idx[mask]
    rmsf_ch = rmsf[mask]
    win = min(10, len(rmsf_ch))
    if win > 1:
        mu = np.convolve(rmsf_ch, np.ones(win)/win, mode='valid')
        ax.plot(idx_ch[win//2:win//2+len(mu)], mu, color=colors[ch], lw=1.5, alpha=0.6)

ax.axhline(y=2.0, color='#95a5a6', ls='--', lw=0.8, alpha=0.5)
ax.set_ylabel('Calpha RMSF (A)', fontsize=10)
ax.set_xlabel('Residue Index (sequential)', fontsize=10)
ax.legend(fontsize=8)
ax.grid(alpha=0.15)

# Flexible regions annotation (RMSF > 4.0)
high = rmsf > 4.0
if high.any():
    regions = []
    start = np.where(high)[0][0]
    for i in range(1, len(high)):
        if high[i] and not high[i-1]: start = i
        if not high[i] and high[i-1]: regions.append((start, i-1))
    if high[-1]: regions.append((start, len(high)-1))
    for s, e in regions:
        ax.axvspan(s, e, alpha=0.08, color='#e74c3c')
        ax.text((s+e)/2, ax.get_ylim()[1]*0.97, f'flex\n{s}-{e}',
                ha='center', fontsize=5.5, color='#e74c3c')

# Annotate key residues
for label, ch, seq_id in [
    ("E63", 0, 63), ("K2", 2, 2), ("peptide ends", None, None)]:
    if ch is None:
        continue
    for i in range(n_res):
        if int(chains[i]) == ch and int(resnums[i]) == seq_id:
            ax.annotate(f'  {label}={rmsf[i]:.2f}A',
                       xy=(i, rmsf[i]), fontsize=8, fontweight='bold',
                       color='#2c3e50',
                       bbox=dict(boxstyle='round,pad=0.2', fc='white', alpha=0.7))
            break

plt.savefig(OUT / 'sb_rmsf_detailed.png')
plt.close()
print('2. sb_rmsf_detailed.png done')

# ============================================================
# Summary
# ============================================================
print(f'\nRMSF Summary:')
for ch in sorted(colors):
    m = rmsf[chains == ch]
    print(f'  {labels[ch]}: mean={m.mean():.3f}A  max={m.max():.3f}A  min={m.min():.3f}A')
# E63
for i in range(n_res):
    if int(chains[i]) == 0 and int(resnums[i]) == 63:
        print(f'  A:E63 (salt bridge anchor): {rmsf[i]:.3f}A')
# K2 (LYS at position 2 in peptide chain)
for i in range(n_res):
    if int(chains[i]) == 2 and int(resnums[i]) == 2:
        print(f'  C:K2 (P2 ligand): {rmsf[i]:.3f}A')
print(f'\nOutput: {OUT}')
