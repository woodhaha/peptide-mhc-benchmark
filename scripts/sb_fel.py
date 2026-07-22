#!/usr/bin/env python3
"""FEL: Free Energy Landscape from salt bridge distance + potential energy.
DeltaG = -RT * ln(P) using 2D histogram.
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.ndimage import gaussian_filter

plt.rcParams.update({'font.size': 9, 'figure.dpi': 150, 'savefig.bbox': 'tight'})

DATA = Path("D:/Researching/Peptide epitope/md results/salt_bridge")
OUT  = Path("D:/Researching/Peptide epitope/md pic")
R = 0.001987  # kcal/mol/K
T = 300

d = np.load(str(DATA / "md_dist.npz"))
mo = d['mo']; Ep = d['Ep']; t_ns = d['t_ns']

# ============================================================
# FEL: 2D histogram over (distance, PE)
# ============================================================
fig, axes = plt.subplots(1, 3, figsize=(16, 5.5))
fig.suptitle('K-E63 Salt Bridge — Free Energy Landscape (100ns)',
             fontsize=13, fontweight='bold')

# 1a. 2D FEL
ax = axes[0]
bins_d = np.linspace(3.5, 6.0, 50)
bins_e = np.linspace(-292800, -291600, 50)
H, xedges, yedges = np.histogram2d(mo, Ep, bins=[bins_d, bins_e])
H = H.T + 1e-10  # avoid log(0)
H_smooth = gaussian_filter(H, sigma=1.0)
FEL = R * T * (-np.log(H_smooth) + np.log(H_smooth.max()))

# Clip to reasonable range
FEL = np.clip(FEL, 0, 10)

extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
im = ax.imshow(FEL, extent=extent, origin='lower', cmap='hot_r',
               aspect='auto', vmin=0, vmax=5)
ax.contour(FEL, levels=np.arange(0, 6, 0.5), colors='k', linewidths=0.4,
           extent=extent, origin='lower', alpha=0.3)
ax.set_xlabel('CA-E63 min distance (A)', fontsize=10)
ax.set_ylabel('Potential Energy (kcal/mol)', fontsize=9)
ax.set_title(f'FEL (kT units, 0-5 kT shown)', fontsize=10)
cbar = plt.colorbar(im, ax=ax)
cbar.set_label('k$_B$T', fontsize=8)

# Mark minimum
min_idx = np.unravel_index(FEL.argmin(), FEL.shape)
d_min = xedges[0] + (min_idx[1] + 0.5) * (xedges[1]-xedges[0])
e_min = yedges[0] + (min_idx[0] + 0.5) * (yedges[1]-yedges[0])
ax.plot(d_min, e_min, 'k*', ms=12, zorder=5)
ax.annotate(f'Minimum\n({d_min:.2f}A, {e_min:.0f})', (d_min, e_min),
            xytext=(d_min+0.3, e_min+100), fontsize=7, ha='left',
            arrowprops=dict(arrowstyle='->', color='k'))

# 1b. 1D free energy profile along distance
ax = axes[1]
H_d, _ = np.histogram(mo, bins=50, range=(3.5, 6.0))
H_d = H_d + 1e-10
FE_d = R * T * (-np.log(H_d) + np.log(H_d.max()))
bin_centers = np.linspace(3.525, 5.975, 50)
# Smooth
FE_d_s = gaussian_filter(FE_d, sigma=1.0)
ax.plot(bin_centers, FE_d_s, color='#3498db', lw=2)
ax.fill_between(bin_centers, 0, FE_d_s, color='#3498db', alpha=0.2)
ax.axvline(mo.mean(), color='#e74c3c', ls='--', lw=1, label=f"Mean={mo.mean():.2f}A")
ax.axvline(4.0, color='#27ae60', ls=':', lw=1, alpha=0.5, label='4.0A')
ax.set_xlabel('CA-E63 min distance (A)', fontsize=10)
ax.set_ylabel('Free Energy (kT)', fontsize=10)
ax.set_title('1D Profile along Distance', fontsize=10)
ax.legend(fontsize=8)
ax.grid(alpha=0.15)

# 1c. 1D free energy along PE
ax = axes[2]
H_e, _ = np.histogram(Ep, bins=50, range=(-292800, -291600))
H_e = H_e + 1e-10
FE_e = R * T * (-np.log(H_e) + np.log(H_e.max()))
FE_e_s = gaussian_filter(FE_e, sigma=1.0)
bc_e = np.linspace(-292775, -291625, 50)
ax.plot(bc_e, FE_e_s, color='#9b59b6', lw=2)
ax.fill_between(bc_e, 0, FE_e_s, color='#9b59b6', alpha=0.2)
ax.axvline(Ep.mean(), color='#e74c3c', ls='--', lw=1, label=f"Mean={Ep.mean():.0f}")
ax.set_xlabel('Potential Energy (kcal/mol)', fontsize=9)
ax.set_ylabel('Free Energy (kT)', fontsize=10)
ax.set_title('1D Profile along PE', fontsize=10)
ax.legend(fontsize=8)
ax.grid(alpha=0.15)

plt.tight_layout()
plt.savefig(OUT / 'sb_fel.png')
plt.close()
print('1. sb_fel.png done')

# ============================================================
# Additional: 1D FEL comparison (head vs tail 50ns)
# ============================================================
half = len(mo)//2
fig, ax = plt.subplots(figsize=(10, 4.5))
fig.suptitle('Free Energy Profile: Head (0-50ns) vs Tail (50-100ns)',
             fontsize=13, fontweight='bold')

for label, arr, color in [('0-50ns', mo[:half], '#3498db'), ('50-100ns', mo[half:], '#e74c3c')]:
    H_d, _ = np.histogram(arr, bins=50, range=(3.5, 6.0))
    H_d = H_d + 1e-10
    FE_d = R * T * (-np.log(H_d) + np.log(H_d.max()))
    FE_d_s = gaussian_filter(FE_d, sigma=1.0)
    ax.plot(bin_centers, FE_d_s, color=color, lw=2, alpha=0.8, label=f'{label}  mean={arr.mean():.2f}A')
    ax.fill_between(bin_centers, 0, FE_d_s, color=color, alpha=0.1)

ax.axvline(4.0, color='#27ae60', ls=':', lw=1, alpha=0.5)
ax.set_xlabel('CA-E63 min distance (A)', fontsize=10)
ax.set_ylabel('Free Energy (kT)', fontsize=10)
ax.legend(fontsize=9)
ax.grid(alpha=0.15)

plt.tight_layout()
plt.savefig(OUT / 'sb_fel_head_tail.png')
plt.close()
print('2. sb_fel_head_tail.png done')

print('\nFEL analysis done ->', OUT)
