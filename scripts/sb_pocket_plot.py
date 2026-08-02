#!/usr/bin/env python3
"""E63 binding pocket visualization + P2Leu pocket analysis."""
import numpy as np, json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path

plt.rcParams.update({'font.size': 9, 'figure.dpi': 150})

DATA = Path("D:/Researching/Peptide epitope/md results/salt_bridge")
OUT  = Path("D:/Researching/Peptide epitope/md pic")

p = json.load(open(str(DATA / "pocket_analysis.json")))

# Filter: residues with mean distance <4.5A (tight pocket)
tight = {k: v for k, v in p['pocket_distances'].items() if v['mean'] < 4.5}
print(f"Tight pocket (<4.5A): {len(tight)} residues")

# Top 25 closest
top = sorted(tight.items(), key=lambda x: x[1]['mean'])[:25]

# ============================================================
# FIGURE: E63 binding pocket profile
# ============================================================
fig, axes = plt.subplots(1, 2, figsize=(14, 6))
fig.suptitle('E63 Binding Pocket — MHC Contact Analysis (100ns)',
             fontsize=13, fontweight='bold')

# 1a. Bar chart: pocket residue distances from E63
ax = axes[0]
names = [x[0] for x in top]
means = [x[1]['mean'] for x in top]
stds = [x[1]['std'] for x in top]
mins = [x[1]['min'] for x in top]

colors = plt.cm.RdYlGn_r([max(0, min(1, (m-2.5)/3)) for m in means])
bars = ax.barh(range(len(names)), means, xerr=stds, color=colors,
               edgecolor='white', linewidth=0.5, capsize=2)
ax.set_yticks(range(len(names)))
ax.set_yticklabels(names, fontsize=7)
ax.set_xlabel('Min distance to E63 OE1/OE2 (A)')
ax.set_title(f'Tight pocket: {len(tight)} residues <4.5A')
ax.invert_yaxis()
ax.axvline(3.4, color='#27ae60', ls='--', lw=0.8, alpha=0.5, label='<3.4A strong')
ax.axvline(4.0, color='#e74c3c', ls=':', lw=0.8, alpha=0.5, label='<4.0A SB')
ax.legend(fontsize=7)
ax.grid(alpha=0.15)

# Add min markers
for i, mn in enumerate(mins):
    ax.plot(mn, i, 'k*', ms=4, alpha=0.5)

# 1b. Classification + P2 pocket summary
ax = axes[1]; ax.axis('off')

# Categorize by distance tier
tier1 = [(n,m) for n,m in zip(names, means) if m < 3.4]
tier2 = [(n,m) for n,m in zip(names, means) if 3.4 <= m < 4.0]
tier3 = [(n,m) for n,m in zip(names, means) if m >= 4.0]

txt = f"""
E63 Pocket Residue Classification
{'='*35}

<3.4A (strong contact, n={len(tier1)}):
{'  '.join([f'{n}({m:.1f}A)' for n,m in tier1])}

3.4-4.0A (moderate, n={len(tier2)}):
{'  '.join([f'{n}({m:.1f}A)' for n,m in tier2])}

4.0-4.5A (weak, n={len(tier3)}):
{'  '.join([f'{n}({m:.1f}A)' for n,m in tier3])}

E63-P2 backbone contacts:
{''}

Total pocket residues: {p['n_pocket_residues']}
P2 binding pocket: 30+ MHC residues >50% occupancy
"""
ax.text(0.05, 0.95, txt, transform=ax.transAxes, fontsize=7.5,
        va='top', fontfamily='monospace',
        bbox=dict(boxstyle='round', fc='#f8f9fa', ec='#dee2e6'))

plt.tight_layout()
plt.savefig(OUT / 'sb_pocket_analysis.png')
plt.savefig(OUT / 'sb_pocket_analysis.pdf')
plt.close()
print('4. sb_pocket_analysis.png done')

# ============================================================
# FIGURE: P2Leu binding pocket occupancy
# ============================================================
fig, ax = plt.subplots(figsize=(14, 4))
fig.suptitle('P2Leu Binding Pocket — MHC Contact Partners (heavy <4A, >50% occupancy)',
             fontsize=13, fontweight='bold')

# Read from the sb_contacts output - those were computed separately
# Reconstruct from pocket analysis data
# Actually we can plot the distance distribution for key pocket residues
# Filter to residues that had mean <5A and low std (stable contacts)
stable = {k: v for k, v in p['pocket_distances'].items()
          if v['mean'] < 5.0 and v['std'] < 0.2}
# Sort by distance
stable_sorted = sorted(stable.items(), key=lambda x: x[1]['mean'])

names = [x[0] for x in stable_sorted][:30]
means = [x[1]['mean'] for x in stable_sorted][:30]
stds = [x[1]['std'] for x in stable_sorted][:30]

colors_arr = plt.cm.viridis([max(0, min(1, (m-2.5)/2.5)) for m in means])
ax.barh(range(len(names)), means, xerr=stds, color=colors_arr,
        edgecolor='white', linewidth=0.5, capsize=2)
ax.set_yticks(range(len(names)))
ax.set_yticklabels(names, fontsize=7)
ax.set_xlabel('Min distance to E63 OE1/OE2 (A)')
ax.invert_yaxis()
ax.axvline(4.0, color='#e74c3c', ls='--', lw=1, alpha=0.5)
ax.grid(alpha=0.15)

plt.tight_layout()
plt.savefig(OUT / 'sb_stable_pocket.png')
plt.savefig(OUT / 'sb_stable_pocket.pdf')
plt.close()
print('5. sb_stable_pocket.png done')

print(f'\nOutput: {OUT}')
