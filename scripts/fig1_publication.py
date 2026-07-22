#!/usr/bin/env python3
"""Publication-ready composite figure: K-E63 Salt Bridge 100ns."""
import numpy as np, json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter

plt.rcParams.update({'font.size': 8, 'figure.dpi': 300,
    'font.family': 'sans-serif', 'axes.linewidth': 0.8})

DATA = "D:/Researching/Peptide epitope/md results/salt_bridge"
OUT  = "D:/Researching/Peptide epitope/md pic/Fig1_salt_bridge.png"

d = np.load(DATA + "/md_dist.npz")
t, mo, sb, T = d['t_ns'], d['mo'], d['sb'], d['T']
dg = np.load(DATA + "/mmgbsa.npz")['dG']
rf = np.load(DATA + "/rmsf.npz")
info = json.load(open(DATA + "/mmgbsa_summary.json"))
sb_i = json.load(open(DATA + "/sb_summary.json"))

R, KT = 0.001987, 300
bins = np.linspace(3.5, 6.0, 40)
H, xe, ye = np.histogram2d(mo, d['Ep'], bins=[bins, 40])
FEL = -R*KT*np.log(gaussian_filter(H.T+1e-10, 1)/H.max())

fig = plt.figure(figsize=(7.08, 5.5))
gs = fig.add_gridspec(2, 2, hspace=0.35, wspace=0.35, left=0.09, right=0.97, bottom=0.1, top=0.92)

# A
ax = fig.add_subplot(gs[0, 0])
ax.plot(t, d['d1'], '#3498db', lw=0.3, alpha=0.3)
ax.plot(t, mo, '#e74c3c', lw=0.5, alpha=0.7)
mu = np.convolve(mo, np.ones(100)/100, mode='valid')
ax.plot(t[50:50+len(mu)], mu, '#2c3e50', lw=1.2)
ax.axhline(4.0, color='#27ae60', ls=':', lw=0.8, alpha=0.5)
ax.set_xlim(0, 100); ax.set_ylim(3.2, 8)
ax.set_ylabel(r'CA$_{\rm P2}$-E63(OE) ($\AA$)')
ax.set_xlabel('Time (ns)')
ax.text(0.03, 0.95, 'Mean=%.2f +/- %.2f A' % (mo.mean(), mo.std()),
        transform=ax.transAxes, fontsize=7, va='top',
        bbox=dict(boxstyle='round', fc='white', alpha=0.7, ec='none'))
ax.text(-0.22, 1.02, 'A', transform=ax.transAxes, fontsize=11, fontweight='bold')
ax.grid(alpha=0.1, lw=0.3)

# B
ax = fig.add_subplot(gs[0, 1])
ax.hist(dg, bins=20, color='#9b59b6', edgecolor='white', lw=0.5, alpha=0.8)
ax.axvline(x=info['dG_mean'], color='#e74c3c', ls='-', lw=1.5)
ax.axvline(x=0, color='#7f8c8d', ls='--', lw=0.8)
ax.set_xlabel(r'DG$_{\rm bind}$ (kcal/mol)')
ax.set_ylabel('Frames')
ax.text(0.03, 0.95, 'DG=%+.1f +/- %.1f' % (info['dG_mean'], info['dG_std']),
        transform=ax.transAxes, fontsize=7, va='top',
        bbox=dict(boxstyle='round', fc='white', alpha=0.7, ec='none'))
ax.text(-0.22, 1.02, 'B', transform=ax.transAxes, fontsize=11, fontweight='bold')
ax.grid(alpha=0.1, lw=0.3)

# C
ax = fig.add_subplot(gs[1, 0])
extent = [xe[0], xe[-1], ye[0], ye[-1]]
FEL_c = np.clip(FEL, 0, 5)
im = ax.imshow(FEL_c, extent=extent, origin='lower', cmap='hot_r', aspect='auto', vmin=0, vmax=5)
ax.contour(FEL_c, levels=np.arange(0.5, 5.5, 0.5), colors='k', lw=0.3,
           extent=extent, origin='lower', alpha=0.3)
ax.set_xlabel(r'CA$_{\rm P2}$-E63(OE) ($\AA$)')
ax.set_ylabel('Potential E (kcal/mol)')
cb = plt.colorbar(im, ax=ax, shrink=0.7, pad=0.02)
cb.set_label('kT', fontsize=7)
ax.text(-0.22, 1.02, 'C', transform=ax.transAxes, fontsize=11, fontweight='bold')
ax.grid(alpha=0.1, lw=0.3)

# D
ax = fig.add_subplot(gs[1, 1])
x = np.arange(len(rf['rmsf']))
colors_ch = {0: '#3498db', 1: '#2ecc71', 2: '#e74c3c'}
lbls = {0: 'MHC', 1: 'b2m', 2: 'Peptide'}
for ch in [0, 1, 2]:
    m = rf['chains'] == ch
    ax.plot(x[m], rf['rmsf'][m], '.', color=colors_ch[ch], ms=2, alpha=0.5, label=lbls[ch])
for ch in [0, 1, 2]:
    m = rf['chains'] == ch
    if m.sum() > 10:
        mu_c = np.convolve(rf['rmsf'][m], np.ones(5)/5, mode='valid')
        ax.plot(x[m][:len(mu_c)], mu_c, color=colors_ch[ch], lw=1.2)
ax.axhline(2.0, color='#95a5a6', ls='--', lw=0.5, alpha=0.4)
ax.set_xlabel('Residue index (sequential)')
ax.set_ylabel(r'C$_\alpha$ RMSF ($\AA$)')
ax.legend(fontsize=6, markerscale=3, ncol=3, loc='upper right')
ax.set_xlim(0, len(rf['rmsf']))
ax.text(-0.22, 1.02, 'D', transform=ax.transAxes, fontsize=11, fontweight='bold')
ax.grid(alpha=0.1, lw=0.3)

# Footer
sb_pct = sb.mean() * 100
fig.text(0.5, 0.01,
    '1DUZ K-E63 Salt Bridge  |  100ns OpenMM (amber14, TIP3P, HMR 4fs)  |  '
    'SB occupancy: %.0f%%  |  CA-E63: %.2f A  |  DG_b: %+.1f kcal/mol  |  MHC RMSF: %.2f A' % (
    sb_pct, mo.mean(), info['dG_mean'], rf['rmsf'][rf['chains']==0].mean()),
    ha='center', va='center', fontsize=5.5, color='#555')

plt.savefig(OUT, dpi=300)
plt.close()
print('Figure saved:', OUT)
