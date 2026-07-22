#!/usr/bin/env python3
"""Block error analysis for MM-GBSA dG. Estimate statistical uncertainty."""
import numpy as np, json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path

plt.rcParams.update({'font.size': 9, 'figure.dpi': 150})

DATA = Path("D:/Researching/Peptide epitope/md results/salt_bridge")
OUT  = Path("D:/Researching/Peptide epitope/md pic")

d = np.load(str(DATA / "mmgbsa.npz"))
with open(str(DATA / "sb_summary.json")) as f:
    sb = json.load(f)

dG = d['dG']
n = len(dG)

# Block error estimation
max_block = n // 4
block_sizes = np.arange(1, max_block + 1)
block_errs = np.zeros(len(block_sizes))

for i, bs in enumerate(block_sizes):
    n_blocks = n // bs
    block_means = dG[:n_blocks*bs].reshape(n_blocks, bs).mean(axis=1)
    block_errs[i] = block_means.std() / np.sqrt(n_blocks - 1) if n_blocks > 1 else 0

# Automatic knee detection: where derivative stabilizes
d_err = np.gradient(block_errs)
# Find the point where relative change < 5%
knee = None
for i in range(len(block_errs)):
    if block_errs[i] <= 0:
        continue
    rel_change = abs(d_err[i]) / (block_errs[i] + 1e-10)
    if rel_change < 0.02 and block_errs[i] > 0:
        # Check that it stays flat
        if all(abs(d_err[i:]) / (block_errs[i:] + 1e-10) < 0.05):
            knee = block_sizes[i]
            break

# If knee not found, take where d(err)/d(block) first goes below threshold
if knee is None:
    for i in range(len(block_sizes)):
        if i > 0 and block_errs[i] > 0 and abs(d_err[i]) / block_errs[i] < 0.03:
            knee = block_sizes[i]
            break
if knee is None:
    knee = block_sizes[-1]

# Get the error at the knee
knee_idx = np.where(block_sizes == knee)[0][0] if knee in block_sizes else -1
estimated_error = block_errs[knee_idx] if knee_idx >= 0 else block_errs[-1]

# Also compute naive SEM for comparison
sem = dG.std() / np.sqrt(n)

print(f"Block error analysis:")
print(f"  Frames: {n}")
print(f"  Mean dG: {dG.mean():+.2f} kcal/mol")
print(f"  Naive SEM: {sem:.2f} kcal/mol")
print(f"  Block estimate: {estimated_error:.2f} kcal/mol")
print(f"  Knee block size: {knee}")
print(f"  True uncertainty (block-corrected): {estimated_error:.2f} kcal/mol")

# ============================================================
# FIGURE: Block error analysis
# ============================================================
fig, axes = plt.subplots(1, 3, figsize=(15, 5))
fig.suptitle('MM-GBSA dG — Block Error Estimation',
             fontsize=13, fontweight='bold')

ax = axes[0]
ax.plot(block_sizes, block_errs, color='#3498db', lw=1.5)
ax.axvline(knee, color='#e74c3c', ls='--', lw=1.5,
           label=f'Knee={knee} (err={estimated_error:.2f})')
ax.axhline(sem, color='#95a5a6', ls=':', lw=1.5, alpha=0.7,
           label=f'Naive SEM={sem:.2f}')
ax.set_xlabel('Block size (frames)')
ax.set_ylabel('Standard error (kcal/mol)')
ax.set_title('Block Error Convergence')
ax.legend(fontsize=8)
ax.grid(alpha=0.15)

ax = axes[1]
ax.plot(block_sizes, block_errs * np.sqrt(block_sizes), color='#e67e22', lw=1.5)
ax.axvline(knee, color='#e74c3c', ls='--', lw=1)
ax.set_xlabel('Block size (frames)')
ax.set_ylabel(r'$\sigma \times \sqrt{N_{block}}$')
ax.set_title('Scaled Error (plateau = convergence)')
ax.grid(alpha=0.15)

ax = axes[2]
ax.axis('off')
txt = f"""
Block Error Summary
{'='*30}

Mean dG:   {dG.mean():+7.2f} kcal/mol
Block est: +/- {estimated_error:.2f} kcal/mol
Naive SEM: +/- {sem:.2f} kcal/mol

Block method accounts for auto-
correlation. The true uncertainty
is {estimated_error/sem:.1f}x larger
than naive SEM.

Knee at block size={knee}
({knee * 0.5:.1f}ns per block,
{n//knee} blocks total)

Salt Bridge: {sb['sb_pct']:.1f}%
CA-E63: {sb['ca_mean']:.2f}A
Stable: {'YES' if sb['stable'] else 'NO'}
"""
ax.text(0.05, 0.95, txt, transform=ax.transAxes, fontsize=8.5,
        va='top', fontfamily='monospace',
        bbox=dict(boxstyle='round', fc='#f8f9fa', ec='#dee2e6'))

plt.tight_layout()
plt.savefig(OUT / 'sb_block_error.png')
plt.close()
print('3. sb_block_error.png done')

# Save estimate
json.dump({
    'mean_dG': float(dG.mean()), 'block_error': float(estimated_error),
    'naive_sem': float(sem), 'knee_blocks': int(knee),
    'n_frames': int(n), 'n_blocks': int(n//knee)
}, open(DATA/'block_error_estimate.json', 'w'), indent=2)
print(f'\nOutput: {OUT}')
print(f'Estimate saved: block_error_estimate.json')
