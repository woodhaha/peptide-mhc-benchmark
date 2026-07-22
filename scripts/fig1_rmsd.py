#!/usr/bin/env python3
"""Fig1: RMSD — verify system equilibrium (all-protein backbone + per-chain + pocket)."""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path

DATA = Path("D:/Researching/Peptide epitope/md results/salt_bridge")
OUT  = Path("D:/Researching/Peptide epitope/md pic")

d = np.load(str(DATA / "rmsd.npz"))
t_ns = d['t_ns']

fig, axes = plt.subplots(3, 1, figsize=(14, 8), sharex=True,
    gridspec_kw={'height_ratios': [2, 1, 1]})
fig.suptitle('1DUZ K-E63 — Backbone RMSD (100ns MD on RTX 5090)',
             fontsize=13, fontweight='bold', y=0.98)

# 1a. All-protein backbone RMSD
ax = axes[0]
ax.plot(t_ns, d['backbone'], color='#2c3e50', lw=0.4, alpha=0.6, label='All protein backbone')
win = 50
mu = np.convolve(d['backbone'], np.ones(win)/win, mode='valid')
ax.plot(t_ns[win//2:win//2+len(mu)], mu, color='#e74c3c', lw=2, label=f'{win*0.1:.0f}ns running avg')
ax.axhline(d['backbone'].mean(), color='#3498db', ls='--', lw=1, alpha=0.7)
ax.set_ylabel('RMSD (Å)', fontsize=10)
ax.legend(fontsize=8)
ax.grid(alpha=0.15)
ax.text(0.01, 0.95, f'Mean={d["backbone"].mean():.2f}±{d["backbone"].std():.2f}Å  '
                     f'[{d["backbone"].min():.2f}–{d["backbone"].max():.2f}]',
        transform=ax.transAxes, fontsize=9, va='top',
        bbox=dict(boxstyle='round', fc='white', alpha=0.7))
# Stability region annotation
stable = t_ns > 30
if stable.sum() > 0:
    mu_stable = d['backbone'][stable].mean()
    ax.axvline(30, color='#27ae60', ls=':', lw=1.5)
    ax.text(30, ax.get_ylim()[1]*0.95, '30ns\n(stable)', fontsize=7, ha='center', color='#27ae60')

# 1b. Per-chain comparison
ax = axes[1]
for label, data, color in [('MHC α', d['MHC'], '#3498db'),
                            ('β2M', d['B2M'], '#2ecc71'),
                            ('Peptide', d['Peptide'], '#e74c3c')]:
    ax.plot(t_ns, data, color=color, lw=0.5, alpha=0.5)
    mu_ch = np.convolve(data, np.ones(win)/win, mode='valid')
    ax.plot(t_ns[win//2:win//2+len(mu_ch)], mu_ch, color=color, lw=2, label=f'{label} avg')
ax.set_ylabel('RMSD (Å)', fontsize=10)
ax.legend(fontsize=8)
ax.grid(alpha=0.15)

# 1c. Binding pocket RMSD
ax = axes[2]
ax.plot(t_ns, d['pocket'], color='#9b59b6', lw=0.5, alpha=0.6, label='Binding pocket (MHC 55-70 + peptide)')
mu_poc = np.convolve(d['pocket'], np.ones(win)/win, mode='valid')
ax.plot(t_ns[win//2:win//2+len(mu_poc)], mu_poc, color='#e74c3c', lw=2, label='running avg')
ax.axhline(d['pocket'].mean(), color='#9b59b6', ls='--', lw=1, alpha=0.5)
ax.set_xlabel('Time (ns)', fontsize=10)
ax.set_ylabel('RMSD (Å)', fontsize=10)
ax.legend(fontsize=8)
ax.grid(alpha=0.15)
ax.text(0.01, 0.95, f'Mean={d["pocket"].mean():.2f}±{d["pocket"].std():.2f}Å',
        transform=ax.transAxes, fontsize=9, va='top',
        bbox=dict(boxstyle='round', fc='white', alpha=0.7))

plt.tight_layout()
plt.savefig(OUT / 'Fig1_RMSD.png', dpi=300, bbox_inches='tight')
plt.savefig(OUT / 'Fig1_RMSD.pdf', bbox_inches='tight')
plt.close()
print('Fig1_RMSD.png/pdf done')
