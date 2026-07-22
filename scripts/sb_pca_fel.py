"""PCA + FEL + DCCM for main K-E63 trajectory.  Run on AutoDL (needs mdtraj)."""
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter
import mdtraj as md
plt.rcParams.update({'font.size': 9, 'figure.dpi': 150, 'font.family': 'sans-serif'})

D = '/root/autodl-tmp/docking'
FIG = '/root/autodl-tmp/docking'
KT = 0.596  # kcal/mol at 300K

# ── Load trajectory ──
print('Loading trajectory...', flush=True)
top = md.load(f'{D}/final.pdb')
traj = md.load(f'{D}/traj.dcd', top=f'{D}/final.pdb')
print(f'  {traj.n_frames} frames, {traj.n_atoms} atoms', flush=True)

# Select Cα
ca_idx = top.topology.select('name CA')
ca_coords = traj.xyz[:, ca_idx, :]
n_frames, n_ca = ca_coords.shape[:2]
ca_flat = ca_coords.reshape(n_frames, n_ca * 3)
print(f'  Cα: {n_ca} atoms, {ca_flat.shape[1]} dims', flush=True)

# ── PCA ──
mean = ca_flat.mean(axis=0)
centered = ca_flat - mean
U, S, Vt = np.linalg.svd(centered, full_matrices=False)
pc1 = centered @ Vt[0]; pc2 = centered @ Vt[1]
var_exp = S[:2]**2 / (S**2).sum()
print(f'  PC1={var_exp[0]*100:.1f}%, PC2={var_exp[1]*100:.1f}%', flush=True)
np.savez(f'{D}/pca.npz', pc1=pc1, pc2=pc2, var_exp=var_exp, mean=mean, components=Vt[:2])

# ── FEL ──
fig, axes = plt.subplots(1, 2, figsize=(10, 4.5))
fig.suptitle('K-E63 Salt Bridge — PCA Free Energy Landscape', fontsize=12, fontweight='bold')

ax = axes[0]
time = np.arange(n_frames) * 100.0 / n_frames
sc = ax.scatter(pc1, pc2, c=time, cmap='gist_heat_r', s=6, alpha=0.6, edgecolors='none')
plt.colorbar(sc, ax=ax, label='Time (ns)', shrink=0.8)
ax.plot(pc1[0], pc2[0], 'g*', ms=10, mew=1, label='Start')
ax.plot(pc1[-1], pc2[-1], 'rD', ms=6, mew=1, label='End')
ax.set_xlabel(f'PC1 ({var_exp[0]*100:.1f}%)')
ax.set_ylabel(f'PC2 ({var_exp[1]*100:.1f}%)')
ax.set_title('Trajectory Projection')
ax.legend(fontsize=7)

ax = axes[1]
H, xe, ye = np.histogram2d(pc1, pc2, bins=50, density=True)
H = gaussian_filter(H, sigma=1.2)
H[H < 1e-12] = 1e-12
FE = -KT * np.log(H)
FE -= FE.min()
X, Y = np.meshgrid((xe[:-1]+xe[1:])/2, (ye[:-1]+ye[1:])/2, indexing='ij')
cnt = ax.contourf(X, Y, FE.T, levels=np.arange(0, 7.1, 0.25), cmap='viridis', extend='max')
plt.colorbar(cnt, ax=ax, label='kcal/mol', shrink=0.8)
ax.set_xlabel(f'PC1 ({var_exp[0]*100:.1f}%)')
ax.set_ylabel(f'PC2 ({var_exp[1]*100:.1f}%)')
ax.set_title('Free Energy Landscape')

plt.tight_layout()
plt.savefig(f'{FIG}/Fig_pca_fel.png', dpi=300, bbox_inches='tight')
plt.savefig(f'{FIG}/Fig_pca_fel.pdf', bbox_inches='tight')
plt.close()
print('1. Fig_pca_fel.png/pdf', flush=True)

# ── DCCM ──
print('Computing DCCM...', flush=True)
ca_traj = traj.atom_slice(ca_idx)
frames = ca_traj.xyz  # n_frames × n_ca × 3
mean_xyz = frames.mean(axis=0)  # n_ca × 3
disp = frames - mean_xyz  # n_frames × n_ca × 3
# Loop over coords to avoid 1.2GB intermediate: 3 × (n_ca×n_frames @ n_frames×n_ca)
dccm_sum = np.zeros((n_ca, n_ca))
rmsf2 = np.zeros(n_ca)
for c in range(3):
    dc = disp[:,:,c]                    # n_frames × n_ca
    dccm_sum += dc.T @ dc               # n_ca × n_ca  (sum over frames)
    rmsf2 += (dc**2).mean(axis=0)       # n_ca
denom = np.sqrt(np.outer(rmsf2, rmsf2))
dccm_mean = np.clip(dccm_sum / n_frames / np.maximum(denom, 1e-12), -1, 1)

fig, axes = plt.subplots(1, 2, figsize=(12, 5.5))
fig.suptitle('K-E63 — Dynamic Cross-Correlation (Cα)', fontsize=12, fontweight='bold')

# Full
ds = max(1, n_ca // 100)
reduced = dccm_mean[::ds, ::ds]
ax = axes[0]
im = ax.imshow(reduced, cmap='RdBu_r', vmin=-1, vmax=1, aspect='auto')
plt.colorbar(im, ax=ax, label='Pearson R', shrink=0.85)
ax.set_xlabel('Cα residue index'); ax.set_ylabel('Cα residue index')
ax.set_title(f'Full ({n_ca} Cα, downsampled {ds}×)')

# Domain-blocked
block = 20
n_blocks = n_ca // block
blocked = np.zeros((n_blocks, n_blocks))
for i in range(n_blocks):
    for j in range(n_blocks):
        blocked[i,j] = dccm_mean[i*block:(i+1)*block, j*block:(j+1)*block].mean()
ax = axes[1]
im = ax.imshow(blocked, cmap='RdBu_r', vmin=-0.5, vmax=0.5, aspect='auto')
plt.colorbar(im, ax=ax, label='Mean Pearson R', shrink=0.85)
ax.set_xlabel(f'Block ({block} Cα)'); ax.set_ylabel(f'Block ({block} Cα)')
ax.set_title(f'Block-averaged ({block}-residue)')

plt.tight_layout()
plt.savefig(f'{FIG}/Fig_dccm.png', dpi=300, bbox_inches='tight')
plt.savefig(f'{FIG}/Fig_dccm.pdf', bbox_inches='tight')
plt.close()
np.save(f'{D}/dccm.npy', dccm_mean)
print('2. Fig_dccm.png/pdf', flush=True)
print('Done.', flush=True)
