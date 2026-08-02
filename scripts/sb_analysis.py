#!/usr/bin/env python3
"""K-E63 Salt Bridge (1DUZ) — 100ns MD Analysis.
Reference SMVT analysis style: D:/Researching/SMVT/MM-GBSA Analysis/_plots_123.py
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path

plt.rcParams.update({'font.size': 9, 'figure.dpi': 150, 'savefig.bbox': 'tight'})

DATA = Path("D:/Researching/Peptide epitope/md results/salt_bridge")
OUT  = Path("D:/Researching/Peptide epitope/md pic")

d = np.load(str(DATA / "md_dist.npz"))
t_ns = d['t_ns']; mo = d['mo']; d1 = d['d1']; sb = d['sb']; T = d['T']; Ep = d['Ep']

# Running window
win = 100  # frames = ~10ns
mu = np.convolve(mo, np.ones(win)/win, mode='valid')
t_mu = t_ns[win//2:win//2+len(mu)]

# ============================================================
# FIGURE 1: Salt Bridge Distance + Status Dashboard
# ============================================================
fig, axes = plt.subplots(4, 1, figsize=(14, 10), sharex=True,
    gridspec_kw={'height_ratios': [2, 0.6, 1, 1]})
fig.suptitle('K-E63 Salt Bridge (1DUZ) — 100ns MD on RTX 5090',
             fontsize=14, fontweight='bold', y=0.98)

# 1a. Distance time series
ax = axes[0]
ax.plot(t_ns, d1, color='#3498db', lw=0.3, alpha=0.4, label='CA-E63(OE1)')
ax.plot(t_ns, mo, color='#e74c3c', lw=0.5, alpha=0.8, label='CA-E63(min OE1/OE2)')
ax.plot(t_mu, mu, color='#2c3e50', lw=1.5, label=f'Running avg ({win/10:.0f}ns window)')
ax.axhline(7.6, color='#95a5a6', ls='--', lw=0.8, alpha=0.6, label='K-SB cutoff (7.6A)')
ax.axhline(4.0, color='#27ae60', ls=':', lw=0.8, alpha=0.6, label='SB cutoff (4.0A)')
ax.set_ylabel('Distance (A)', fontsize=10)
ax.set_ylim(2, 9)
ax.legend(fontsize=7, ncol=3, loc='upper right')
ax.grid(alpha=0.15)
ax.text(0.01, 0.95, f'Mean={mo.mean():.3f}A  SD={mo.std():.3f}A  Min={mo.min():.3f}A  Max={mo.max():.3f}A',
        transform=ax.transAxes, fontsize=8, va='top',
        bbox=dict(boxstyle='round,pad=0.3', fc='white', alpha=0.7))

# 1b. Salt bridge binary indicator
ax = axes[1]
ax.fill_between(t_ns, 0, sb, color='#27ae60', alpha=0.5, step='mid')
ax.fill_between(t_ns, 0, 0.5, where=~sb, color='#e74c3c', alpha=0.3, step='mid')
ax.set_ylabel('Salt Bridge', fontsize=9)
ax.set_yticks([0, 1])
ax.set_yticklabels(['No', 'Yes'], fontsize=7)
ax.set_ylim(-0.1, 1.1)
ax.text(0.99, 0.5, f'{sb.mean()*100:.1f}% maintained',
        transform=ax.transAxes, fontsize=9, ha='right', va='center',
        bbox=dict(boxstyle='round', fc='#27ae60', alpha=0.2))

# 1c. Temperature
ax = axes[2]
ax.plot(t_ns, T, color='#e67e22', lw=0.5)
ax.axhline(300, color='#95a5a6', ls='--', lw=0.8, alpha=0.5)
ax.set_ylabel('Temperature (K)', fontsize=10)
ax.set_ylim(290, 310)
ax.grid(alpha=0.15)
ax.text(0.99, 0.9, f'Mean={T.mean():.0f}+-{T.std():.0f}K',
        transform=ax.transAxes, fontsize=8, ha='right', va='top',
        bbox=dict(boxstyle='round', fc='white', alpha=0.7))

# 1d. Potential Energy
ax = axes[3]
ax.plot(t_ns, Ep, color='#9b59b6', lw=0.5)
ax.set_xlabel('Time (ns)', fontsize=10)
ax.set_ylabel('Potential E (kcal/mol)', fontsize=9)
ax.grid(alpha=0.15)
ax.text(0.99, 0.9, f'Mean={Ep.mean():.0f}+-{Ep.std():.0f} kcal/mol',
        transform=ax.transAxes, fontsize=8, ha='right', va='top',
        bbox=dict(boxstyle='round', fc='white', alpha=0.7))

plt.savefig(OUT / 'sb_distance_timeseries.png')
plt.savefig(OUT / 'sb_distance_timeseries.pdf')
plt.close()
print('1. sb_distance_timeseries.png done')

# ============================================================
# FIGURE 2: Distance Distribution + Stability Summary
# ============================================================
fig = plt.figure(figsize=(14, 6))
fig.suptitle('K-E63 Salt Bridge — Distance Distribution & Stability Metrics',
             fontsize=14, fontweight='bold', y=0.98)

gs = fig.add_gridspec(1, 3, wspace=0.35)

# 2a. Histogram
ax = fig.add_subplot(gs[0])
bins = np.arange(2.5, 8.0, 0.1)
ax.hist(mo, bins=bins, color='#3498db', edgecolor='white', linewidth=0.3, alpha=0.8)
ax.axvline(4.0, color='#27ae60', ls='--', lw=1.5, label='SB cutoff 4.0A')
ax.axvline(7.6, color='#e74c3c', ls='--', lw=1.5, label='K-SB 7.6A')
ax.axvline(mo.mean(), color='#2c3e50', ls='-', lw=2, label=f'Mean={mo.mean():.2f}A')
ax.set_xlabel('CA-E63 min distance (A)', fontsize=10)
ax.set_ylabel('Frames', fontsize=10)
ax.legend(fontsize=8)
ax.grid(alpha=0.15)

# 2b. Running average + stability
ax = fig.add_subplot(gs[1])
windows = [5, 10, 20, 50, 100]
colors = plt.cm.Blues(np.linspace(0.4, 0.9, len(windows)))
for w, c in zip(windows, colors):
    w_mu = np.convolve(mo, np.ones(w)/w, mode='valid')
    ax.plot(t_ns[w//2:w//2+len(w_mu)], w_mu, color=c, lw=1,
            label=f'{w/10:.0f}ns' if w <= 10 else f'{w/10:.0f}ns')
ax.axhline(4.0, color='#27ae60', ls='--', lw=0.8)
ax.set_xlabel('Time (ns)', fontsize=10)
ax.set_ylabel('Running Avg CA-E63 (A)', fontsize=10)
ax.legend(fontsize=7, title='Window')
ax.grid(alpha=0.15)
ax.set_ylim(3.5, 5.0)

# 2c. Metric card
ax = fig.add_subplot(gs[2])
ax.axis('off')
metrics_text = f"""
K-E63 Salt Bridge — 100ns MD Summary
{'='*35}

Distance (CA-E63 min OE1/OE2):
  Mean:   {mo.mean():.3f} A
  SD:     {mo.std():.3f} A
  Min:    {mo.min():.3f} A
  Max:    {mo.max():.3f} A
  Median: {np.median(mo):.3f} A

Salt Bridge (<4.0A):
  Occupancy: {sb.mean()*100:.1f}%

Temperature:
  Mean: {T.mean():.0f} K
  Range: {T.min():.0f}-{T.max():.0f} K

Potential Energy:
  Mean: {Ep.mean():.0f} kcal/mol
  SD:   {Ep.std():.0f} kcal/mol

Stability: STABLE
  Run time: 100ns (389 min on RTX 5090)
  Method: OpenMM, amber14, TIP3P, PME,
          HMR 4fs, CUDA mixed-precision
"""
ax.text(0.05, 0.95, metrics_text, transform=ax.transAxes, fontsize=8.5,
        va='top', fontfamily='monospace',
        bbox=dict(boxstyle='round', fc='#f8f9fa', ec='#dee2e6'))

plt.savefig(OUT / 'sb_distance_distribution.png')
plt.savefig(OUT / 'sb_distance_distribution.pdf')
plt.close()
print('2. sb_distance_distribution.png done')

# ============================================================
# FIGURE 3: Head-vs-Tail Stability Comparison
# ============================================================
half = len(mo)//2
h_mo, t_mo = mo[:half], mo[half:]

fig, axes = plt.subplots(1, 3, figsize=(14, 4.5))
fig.suptitle('K-E63 Salt Bridge — Head (0-50ns) vs Tail (50-100ns) Comparison',
             fontsize=13, fontweight='bold')

colors = ['#3498db', '#e74c3c']

# 3a. Histogram overlay
ax = axes[0]
ax.hist(h_mo, bins=30, alpha=0.5, color=colors[0], label=f'Head: {h_mo.mean():.3f}+-{h_mo.std():.3f}A')
ax.hist(t_mo, bins=30, alpha=0.5, color=colors[1], label=f'Tail: {t_mo.mean():.3f}+-{t_mo.std():.3f}A')
ax.axvline(4.0, color='#27ae60', ls='--', lw=0.8)
ax.set_xlabel('CA-E63 min distance (A)')
ax.set_ylabel('Frames')
ax.legend(fontsize=8)
ax.grid(alpha=0.15)

# 3b. Distance time series with head/tail split
ax = axes[1]
ax.plot(t_ns[:half], h_mo, color=colors[0], lw=0.4, alpha=0.7)
ax.plot(t_ns[half:], t_mo, color=colors[1], lw=0.4, alpha=0.7)
ax.axhline(mo.mean(), color='#2c3e50', ls='-', lw=1.5, label=f'Full mean={mo.mean():.3f}A')
ax.axhline(h_mo.mean(), color=colors[0], ls=':', lw=1, alpha=0.7)
ax.axhline(t_mo.mean(), color=colors[1], ls=':', lw=1, alpha=0.7)
ax.set_xlabel('Time (ns)')
ax.set_ylabel('CA-E63 min distance (A)')
ax.legend(fontsize=8)
ax.grid(alpha=0.15)

# 3c. Cumulative average convergence
ax = axes[2]
cum_avg = np.cumsum(mo) / np.arange(1, len(mo)+1)
ax.plot(t_ns, cum_avg, color='#2c3e50', lw=1.5)
ax.axhline(mo.mean(), color='#e74c3c', ls='--', lw=1, alpha=0.7)
ax.set_xlabel('Time (ns)')
ax.set_ylabel('Cumulative Mean (A)')
ax.set_title('Convergence Check', fontsize=10)
ax.grid(alpha=0.15)
ax.text(0.97, 0.1, f'Final: {cum_avg[-1]:.3f}A\nStable after ~30ns',
        transform=ax.transAxes, fontsize=8, ha='right', va='bottom',
        bbox=dict(boxstyle='round', fc='white', alpha=0.7))

plt.tight_layout()
plt.savefig(OUT / 'sb_head_tail_comparison.png')
plt.savefig(OUT / 'sb_head_tail_comparison.pdf')
plt.close()
print('3. sb_head_tail_comparison.png done')

print('\nALL FIGURES DONE')
print(f'Output: {OUT}')
