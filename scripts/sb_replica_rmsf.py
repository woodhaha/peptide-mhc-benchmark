"""RMSF for replica trajectories."""
import numpy as np
import mdtraj as md
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 10, 'figure.dpi': 150})

R2 = '/root/autodl-tmp/replica2'
R3 = '/root/autodl-tmp/replica3'
D  = '/root/autodl-tmp/docking'

colors = {'Main':'#e74c3c', 'Rep2':'#3498db', 'Rep3':'#2ecc71'}

print('Main RMSF...', flush=True)
m = md.load(f'{D}/traj.dcd', top=f'{D}/final.pdb')
r_main = md.rmsf(m, m, 0)  # per-atom RMSF

print('Rep2 RMSF...', flush=True)
r2 = md.load(f'{R2}/traj_s42.dcd', top=f'{R2}/final_s42.pdb')
r_r2 = md.rmsf(r2, r2, 0)

print('Rep3 RMSF...', flush=True)
r3 = md.load(f'{R3}/traj_s99.dcd', top=f'{R3}/final_s99.pdb')
r_r3 = md.rmsf(r3, r3, 0)

# Save
np.savez(f'{D}/rmsf_all.npz', main=r_main, rep2=r_r2, rep3=r_r3)

# Per-chain mean
for name, r in [('Main', r_main), ('Rep2', r_r2), ('Rep3', r_r3)]:
    print(f'{name}: mean={r.mean():.3f}Å')

# Plot
fig, ax = plt.subplots(figsize=(12, 3.5))
ax.plot(r_main, color=colors['Main'], lw=0.7, label='Main (GBSA)')
ax.plot(r_r2, color=colors['Rep2'], lw=0.7, alpha=0.7, label='Rep2 (PME)')
ax.plot(r_r3, color=colors['Rep3'], lw=0.7, alpha=0.7, label='Rep3 (PME)')
ax.set_xlabel('Atom index'); ax.set_ylabel('RMSF (Å)')
ax.set_title('Per-Atom RMSF — 3-Run Comparison')
ax.legend(fontsize=8)
plt.tight_layout()
plt.savefig(f'{D}/Fig_rmsf_3run.png', dpi=300, bbox_inches='tight')
plt.savefig(f'{D}/Fig_rmsf_3run.pdf', bbox_inches='tight')
plt.close()
print('Fig_rmsf_3run.png/pdf saved')
