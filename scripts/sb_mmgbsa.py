#!/usr/bin/env python3
"""K-E63 Salt Bridge (1DUZ) — MM-GBSA Binding Energy Analysis.
Reference: D:/Researching/SMVT/MM-GBSA Analysis/_compare_docking.py
"""
import numpy as np, json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path
from scipy import stats

plt.rcParams.update({'font.size': 9, 'figure.dpi': 150, 'savefig.bbox': 'tight'})

DATA = Path("D:/Researching/Peptide epitope/md results/salt_bridge")
OUT  = Path("D:/Researching/Peptide epitope/md pic")

d = np.load(str(DATA / "mmgbsa.npz"))
with open(str(DATA / "mmgbsa_summary.json")) as f:
    s = json.load(f)
with open(str(DATA / "sb_summary.json")) as f:
    sb = json.load(f)

Gc, Gr, Gl = d['E_c'], d['E_r'], d['E_l']
dG = d['dG']; dE = d['dE']; dSA = d['dSA']
frame = d['frame'] * 0.5  # each frame = 0.5ns (since N_SKIP=5, 1000/200=5 frames, times 0.1ns per frame...)
# Actually each frame index corresponds to 0.5ns: 1000 frames / 200 samples * 0.1ns per frame
# N_SKIP=5 means every 5th DCD frame. DCD was written every RINT=25000 steps, TSTEP=4fs
# time per DCD frame = 25000 * 4fs = 100ps = 0.1ns
# sampled frame time = 5 * 0.1ns = 0.5ns
t_ns = np.arange(len(dG)) * 0.5

# ============================================================
# FIGURE 1: Binding Free Energy Distribution + Time Series
# ============================================================
fig, axes = plt.subplots(2, 2, figsize=(14, 9))
fig.suptitle('1DUZ K-E63 Salt Bridge — MM-GBSA Binding Free Energy (200 frames/100ns)',
             fontsize=13, fontweight='bold')

# 1a. dG histogram
ax = axes[0, 0]
n_bins = 30
ax.hist(dG, bins=n_bins, color='#3498db', edgecolor='white', linewidth=0.5, alpha=0.8)
ax.axvline(s['dG_mean'], color='#e74c3c', ls='-', lw=2, label=f"Mean={s['dG_mean']:+.1f}")
ax.axvline(0, color='#2c3e50', ls='--', lw=1, alpha=0.5)
ax.set_xlabel('dG (kcal/mol)')
ax.set_ylabel('Frames')
ax.set_title('Binding Energy Distribution', fontsize=10)
ax.legend(fontsize=8)
ax.grid(alpha=0.15)
# Annotate favorable fraction
n_fav = int((dG < 0).sum())
ax.text(0.97, 0.95, f'Favorable: {n_fav}/{len(dG)} ({n_fav/len(dG)*100:.0f}%)',
        transform=ax.transAxes, fontsize=9, ha='right', va='top',
        bbox=dict(boxstyle='round', fc='white', alpha=0.7))

# 1b. dG time series
ax = axes[0, 1]
colors = np.where(dG < 0, '#27ae60', '#e74c3c')
ax.scatter(t_ns, dG, c=colors, s=8, alpha=0.5, edgecolors='none')
ax.axhline(s['dG_mean'], color='#2c3e50', ls='-', lw=1.5, label=f"Mean={s['dG_mean']:+.1f}")
ax.axhline(0, color='#7f8c8d', ls='--', lw=1, alpha=0.5)
# Running average
win = 20
mu_dG = np.convolve(dG, np.ones(win)/win, mode='valid')
ax.plot(t_ns[win//2:win//2+len(mu_dG)], mu_dG, color='#2c3e50', lw=2, label=f'{win*0.5:.0f}ns avg')
ax.set_xlabel('Time (ns)')
ax.set_ylabel('dG (kcal/mol)')
ax.set_title('Per-Frame Binding Free Energy', fontsize=10)
ax.legend(fontsize=8)
ax.grid(alpha=0.15)

# 1c. Energy components stack
ax = axes[1, 0]
ax.plot(t_ns, dE, color='#e67e22', lw=0.6, alpha=0.7, label=f"dE={s['dE_mean']:+.1f}")
ax.plot(t_ns, dSA, color='#9b59b6', lw=0.6, alpha=0.7, label=f"dSA={s['dSA_mean']:+.3f}")
ax.plot(t_ns, dG, color='#2c3e50', lw=0.6, alpha=0.5, label=f"dG={s['dG_mean']:+.1f}")
ax.axhline(0, color='#7f8c8d', ls='--', lw=0.8)
ax.set_xlabel('Time (ns)')
ax.set_ylabel('Energy (kcal/mol)')
ax.set_title('Energy Components', fontsize=10)
ax.legend(fontsize=7)
ax.grid(alpha=0.15)

# 1d. Metric summary card
ax = axes[1, 1]
ax.axis('off')
n_all = len(dG)
n_strong = int((dG < -20).sum())
n_mod = int(((dG >= -20) & (dG < -5)).sum())
n_neutral = int(((dG >= -5) & (dG < 5)).sum())
n_bad = int((dG >= 5).sum())

summary = f"""
MM-GBSA Results Summary
{'='*30}

Method: amber14/GBSAOBC/SA
         single-trajectory
         {s['n_frames']} frames / 100ns

dG (binding):  {s['dG_mean']:+7.1f} +/- {s['dG_std']:.1f} kcal/mol
dE (gas+GB):  {s['dE_mean']:+7.1f} +/- {s['dG_std']:.1f} kcal/mol
dG_SA:         {s['dSA_mean']:+.3f} kcal/mol

Classification:
  Favorable:    {n_fav:3d}/{n_all} ({n_fav/n_all*100:.0f}%)
    Strong:     {n_strong:3d}/{n_all} ({n_strong/n_all*100:.0f}%)
    Moderate:   {n_mod:3d}/{n_all} ({n_mod/n_all*100:.0f}%)
  Near-zero:    {n_neutral:3d}/{n_all} ({n_neutral/n_all*100:.0f}%)
  Unfavorable:  {n_bad:3d}/{n_all} ({n_bad/n_all*100:.0f}%)

Limitation: No entropy term (-TdS).
            Use for ranking, not absolute affinity.

Salt Bridge Result:
  Distance: {sb['ca_mean']:.2f}A, {sb['sb_pct']:.0f}% occupancy
  Stability: {'STABLE' if sb['stable'] else 'UNSTABLE'}
"""

ax.text(0.05, 0.95, summary, transform=ax.transAxes, fontsize=8.5,
        va='top', fontfamily='monospace',
        bbox=dict(boxstyle='round', fc='#f8f9fa', ec='#dee2e6'))

plt.savefig(OUT / 'sb_mmgbsa_analysis.png')
plt.close()
print('1. sb_mmgbsa_analysis.png done')

# ============================================================
# FIGURE 2: MM-GBSA vs Salt Bridge Distance Correlation
# ============================================================
# Load distance data
md = np.load(str(DATA / "md_dist.npz"))
sb_dist = md['mo']  # sampled to match mmgbsa frames
sb_sb = md['sb']

# Sample distance to match MM-GBSA frames
sample_idx = np.linspace(0, len(sb_dist)-1, len(dG), dtype=int)
sb_dist_sampled = sb_dist[sample_idx]
sb_sb_sampled = sb_sb[sample_idx]

fig, axes = plt.subplots(1, 3, figsize=(14, 5))
fig.suptitle('MM-GBSA vs Salt Bridge Distance — Cross-Correlation',
             fontsize=13, fontweight='bold')

# 2a. dG vs distance scatter
ax = axes[0]
sc = ax.scatter(sb_dist_sampled, dG, c=t_ns, cmap='viridis', s=15, alpha=0.6, edgecolors='none')
ax.axhline(0, color='#7f8c8d', ls='--', lw=0.8)
ax.axvline(4.0, color='#27ae60', ls='--', lw=0.8, alpha=0.5)
ax.set_xlabel('CA-E63 min distance (A)')
ax.set_ylabel('MM-GBSA dG (kcal/mol)')
ax.grid(alpha=0.15)
cbar = plt.colorbar(sc, ax=ax)
cbar.set_label('Time (ns)', fontsize=8)
# Add correlation
r, p = stats.pearsonr(sb_dist_sampled, dG)
ax.text(0.97, 0.05, f'Pearson r={r:.3f} (p={p:.3g})',
        transform=ax.transAxes, fontsize=8, ha='right', va='bottom',
        bbox=dict(boxstyle='round', fc='white', alpha=0.7))

# 2b. Joint histogram
ax = axes[1]
hb = ax.hexbin(sb_dist_sampled, dG, gridsize=20, cmap='Blues', mincnt=1, alpha=0.8)
ax.axhline(0, color='#7f8c8d', ls='--', lw=0.8)
ax.axvline(4.0, color='#27ae60', ls='--', lw=0.8, alpha=0.5)
ax.set_xlabel('CA-E63 min distance (A)')
ax.set_ylabel('dG (kcal/mol)')
cbar2 = plt.colorbar(hb, ax=ax)
cbar2.set_label('Frames', fontsize=8)

# Linear fit
slope, intercept, r_val, p_val, std_err = stats.linregress(sb_dist_sampled, dG)
x_fit = np.array([sb_dist_sampled.min(), sb_dist_sampled.max()])
ax.plot(x_fit, slope*x_fit + intercept, color='#e74c3c', ls='-', lw=1.5,
        label=f'Slope={slope:.1f}, r^2={r_val**2:.3f}')
ax.legend(fontsize=7)
ax.grid(alpha=0.15)

# 2c. Cumulative dG vs time with distance overlay
ax = axes[2]
t_overlay = np.arange(len(dG)) * 0.5
cum_dG = np.cumsum(dG) / np.arange(1, len(dG)+1)
ax.plot(t_overlay, cum_dG, color='#3498db', lw=1.5, label='Cumulative mean dG')
ax2 = ax.twinx()
ax2.plot(t_overlay, sb_dist_sampled, color='#e74c3c', lw=0.5, alpha=0.4, label='CA-E63 distance')
ax2.axhline(4.0, color='#27ae60', ls='--', lw=0.8, alpha=0.5)
ax.set_xlabel('Time (ns)')
ax.set_ylabel('Cumulative dG (kcal/mol)')
ax2.set_ylabel('Distance (A)')
lines1, labels1 = ax.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
ax.legend(lines1+lines2, labels1+labels2, fontsize=7, loc='upper right')
ax.grid(alpha=0.15)

plt.tight_layout()
plt.savefig(OUT / 'sb_mmgbsa_correlation.png')
plt.close()
print('2. sb_mmgbsa_correlation.png done')

# ============================================================
# FIGURE 3: Comprehensive 4-Panel Dashboard
# ============================================================
fig = plt.figure(figsize=(16, 8))
fig.suptitle('1DUZ K-E63 Salt Bridge — Comprehensive Analysis Dashboard',
             fontsize=14, fontweight='bold')

gs = fig.add_gridspec(2, 3, hspace=0.3, wspace=0.3)

# 3a. SB distance (top left, span 2 cols)
sb_t = md['t_ns']
sb_d1 = md['d1']; sb_mo = md['mo']; sb_T = md['T']; sb_Ep = md['Ep']
mu_sb = np.convolve(sb_mo, np.ones(100)/100, mode='valid')
ax = fig.add_subplot(gs[0, :2])
ax.plot(sb_t, sb_d1, color='#3498db', lw=0.3, alpha=0.4)
ax.plot(sb_t, sb_mo, color='#e74c3c', lw=0.5, alpha=0.8)
ax.plot(sb_t[50:50+len(mu_sb)], mu_sb, color='#2c3e50', lw=1.5)
ax.axhline(7.6, color='#95a5a6', ls='--', lw=0.8, alpha=0.5)
ax.axhline(4.0, color='#27ae60', ls=':', lw=0.8, alpha=0.5)
ax.set_ylabel('Distance (A)')
ax.set_title('CA-E63 Salt Bridge Distance (100ns)')
ax.legend(['OE1', 'min(OE1,OE2)', f'Running avg ({win/10:.0f}ns)'], fontsize=7)
ax.grid(alpha=0.15)
ax.text(0.01, 0.95, f'{sb_mo.mean():.3f}+-{sb_mo.std():.3f}A  SB={md["sb"].mean()*100:.1f}%',
        transform=ax.transAxes, fontsize=9, va='top',
        bbox=dict(boxstyle='round', fc='white', alpha=0.7))

# 3b. MM-GBSA dG (top right)
ax = fig.add_subplot(gs[0, 2])
ax.hist(dG, bins=25, color='#9b59b6', edgecolor='white', linewidth=0.5, alpha=0.8)
ax.axvline(s['dG_mean'], color='#e74c3c', ls='-', lw=2)
ax.axvline(0, color='#7f8c8d', ls='--', lw=1)
ax.set_xlabel('dG (kcal/mol)')
ax.set_ylabel('Frames')
ax.set_title(f"MM-GBSA: {s['dG_mean']:+.1f}+-{s['dG_std']:.1f}")
ax.grid(alpha=0.15)

# 3c. Stability (bottom left)
ax = fig.add_subplot(gs[1, 0])
cum_avg = np.cumsum(sb_mo) / np.arange(1, len(sb_mo)+1)
ax.plot(sb_t, cum_avg, color='#2c3e50', lw=1.5)
ax.axhline(sb_mo.mean(), color='#e74c3c', ls='--', lw=1)
ax.set_xlabel('Time (ns)')
ax.set_ylabel('Cumulative Mean (A)')
ax.set_title(f'Convergence: {sb_mo.mean():.3f}A final')
ax.grid(alpha=0.15)

# 3d. Temperature (bottom middle)
ax = fig.add_subplot(gs[1, 1])
ax.plot(sb_t, sb_T, color='#e67e22', lw=0.5)
ax.axhline(300, color='#95a5a6', ls='--', lw=0.8)
ax.set_xlabel('Time (ns)')
ax.set_ylabel('Temperature (K)')
ax.set_title(f'Mean {sb_T.mean():.0f}+-{sb_T.std():.0f}K')
ax.grid(alpha=0.15)

# 3e. Potential energy (bottom right)
ax = fig.add_subplot(gs[1, 2])
ax.plot(sb_t, sb_Ep, color='#9b59b6', lw=0.5)
ax.set_xlabel('Time (ns)')
ax.set_ylabel('PE (kcal/mol)')
ax.set_title(f'Mean {sb_Ep.mean():.0f}+-{sb_Ep.std():.0f}')
ax.grid(alpha=0.15)

plt.savefig(OUT / 'sb_comprehensive_dashboard.png')
plt.close()
print('3. sb_comprehensive_dashboard.png done')

print(f'\nALL FIGURES DONE -> {OUT}')
