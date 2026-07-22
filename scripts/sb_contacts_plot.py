#!/usr/bin/env python3
"""E63-P2Leu interface contact analysis plots.
Reference: D:/Researching/SMVT/MM-GBSA Analysis/_plots_123.py
"""
import numpy as np, json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path

plt.rcParams.update({'font.size': 9, 'figure.dpi': 150, 'savefig.bbox': 'tight'})

DATA = Path("D:/Researching/Peptide epitope/md results/salt_bridge")
OUT  = Path("D:/Researching/Peptide epitope/md pic")

d = np.load(str(DATA / "md_dist.npz"))
t_ns = d['t_ns']

try:
    cd = np.load(str(DATA / "sb_contacts.npz"))
    cs = json.load(open(str(DATA / "sb_contacts_summary.json")))
except:
    cd = None; cs = None

# ============================================================
# FIGURE: E63-P2Leu interaction modes
# ============================================================
fig, axes = plt.subplots(2, 2, figsize=(14, 8))
fig.suptitle('1DUZ K-E63 Salt Bridge — Interface Interaction Analysis (100ns)',
             fontsize=13, fontweight='bold')

# Panel 1: P2:CA - E63:OE (backbone contact, primary SB metric)
ax = axes[0, 0]
mo = d['mo']
ax.scatter(t_ns, mo, c='#3498db', s=4, alpha=0.3, edgecolors='none')
mu = np.convolve(mo, np.ones(100)/100, mode='valid')
ax.plot(t_ns[50:50+len(mu)], mu, color='#2c3e50', lw=1.5, label=f'Running avg')
ax.axhline(4.0, color='#27ae60', ls='--', lw=1, alpha=0.6, label='SB cutoff 4.0A')
ax.set_ylabel('P2:CA - E63:OE min (A)', fontsize=9)
ax.set_title(f"Backbone contact: {mo.mean():.2f}+-{mo.std():.2f}A  SB={(d['sb'].mean()*100):.1f}%")
ax.legend(fontsize=7)
ax.grid(alpha=0.15)

# Panel 2: E63 sidechain - P2Leu all heavy (sidechain contact)
ax = axes[0, 1]
if cd is not None:
    d_s = cd['d_min']
    ax.scatter(t_ns, d_s, c='#e74c3c', s=4, alpha=0.3, edgecolors='none')
    mu_s = np.convolve(d_s, np.ones(100)/100, mode='valid')
    ax.plot(t_ns[50:50+len(mu_s)], mu_s, color='#2c3e50', lw=1.5)
    ax.set_ylabel('E63 heavy - P2Leu heavy min (A)', fontsize=9)
    ax.set_title(f"Sidechain contact: {cs['mean']:.2f}+-{cs['std']:.2f}A  p34={cs['p34']:.0f}%")
else:
    ax.text(0.5, 0.5, 'no data', ha='center')
ax.grid(alpha=0.15)

# Panel 3: E63-binding mode schematic
ax = axes[1, 0]; ax.axis('off')
txt = f"""
E63-P2Leu Interaction Mode
{'='*30}

Primary contact (backbone):
  P2:CA - E63:OE  = {mo.mean():.2f}+-{mo.std():.2f}A
  Salt bridge:     {d['sb'].mean()*100:.1f}% occupancy

Secondary contact (sidechain):
  E63 heavy - P2Leu heavy = {cs['mean']:.1f}+-{cs['std']:.1f}A
  (<3.5A: {cs['p34']:.0f}%)

Interpretation:
- E63 carboxyl interacts with P2
  peptide backbone (CA/OE contact~4A)
- P2 leucine sidechain projects into
  MHC B-pocket (not contacting E63)
- This is typical for MHC-I: conserved
  E63 anchors peptide backbone at P2
"""
ax.text(0.05, 0.95, txt, transform=ax.transAxes, fontsize=8, va='top',
        fontfamily='monospace',
        bbox=dict(boxstyle='round', fc='#f8f9fa', ec='#dee2e6'))

# Panel 4: Structure schematic P2 in MHC pocket
ax = axes[1, 1]; ax.axis('off')
txt2 = f"""
Interface H-bonds (A-C): {cs['n_hbonds']} unique

Stability Summary:
  CA-E63: {'STABLE' if d['sb'].mean()>0.9 else 'VARIABLE'}
  Temp:   {d['T'].mean():.0f}+-{d['T'].std():.0f}K
  Energy: {d['Ep'].mean():.0f}+-{d['Ep'].std():.0f}

MM-GBSA: {json.load(open(str(DATA/'mmgbsa_summary.json')))['dG_mean']:.1f} kcal/mol
RMSF:    MHC={1.73}A  b2m={3.60}A  peptide={4.48}A

Replicas running: R2(seed=42), R3(seed=99)
"""
ax.text(0.05, 0.95, txt2, transform=ax.transAxes, fontsize=8, va='top',
        fontfamily='monospace',
        bbox=dict(boxstyle='round', fc='#f8f9fa', ec='#dee2e6'))

plt.tight_layout()
plt.savefig(OUT / 'sb_interface_analysis.png')
plt.close()
print('1. sb_interface_analysis.png done')
print(f'\nOutput: {OUT}')
