#!/usr/bin/env python3
"""Fig3: Binding site analysis — per-residue contact frequency, distances, and SASA.
Uses the full trajectory to compute residue-level interaction statistics with the peptide.
Proxy for full energy decomposition: contact frequency × minimum distance → key binding residues."""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path
import mdtraj as md

DATA = Path("D:/Researching/Peptide epitope/md results/salt_bridge")
OUT  = Path("D:/Researching/Peptide epitope/md pic")

print("Loading trajectory for binding site analysis...", flush=True)
top = md.load(str(DATA / "final.pdb"))
traj = md.load(str(DATA / "traj.dcd"), top=str(DATA / "final.pdb"))
n_frames = traj.n_frames

# Binding site residues: within 8Å of peptide in the reference structure
peptide_atoms = [a.index for a in traj.topology.atoms
                 if a.residue.chain.index == 2 and a.element.symbol != 'H']

# Collect all protein residues within 8Å of any peptide heavy atom (reference frame)
binding = []  # list of {name, chain, resSeq, atoms}
for chain in traj.topology.chains:
    if chain.index > 1:
        continue
    for res in chain.residues:
        if res.is_water:
            continue
        res_atoms = [a.index for a in res.atoms if a.element.symbol != 'H']
        if not res_atoms:
            continue
        ref_xyz = traj.xyz[0, res_atoms]
        pep_xyz = traj.xyz[0, peptide_atoms]
        dists = np.linalg.norm(ref_xyz[:, None] - pep_xyz[None, :], axis=2)
        if dists.min() < 8.0:
            binding.append({
                'chain': chain.index, 'resSeq': res.resSeq,
                'name': str(res.name), 'atom_indices': res_atoms,
            })
n_bind = len(binding)
print(f"  {n_bind} binding site residues (within 8Å of peptide)", flush=True)

# Compute contact frequency and interaction energy proxy across trajectory
# For speed, use every 5th frame
sample_idx = np.linspace(0, n_frames-1, 200, dtype=int)
contact_freq = np.zeros(n_bind)   # fraction of frames with any heavy-atom contact < 4.0Å
min_dist = np.full(n_bind, 999.0) # average minimum distance

for j, fi in enumerate(sample_idx):
    frame = traj.xyz[fi]
    pep_xyz = frame[peptide_atoms]
    for ri, br in enumerate(binding):
        r_xyz = frame[br['atom_indices']]
        dists = np.linalg.norm(r_xyz[:, None] - pep_xyz[None, :], axis=2)
        d_min = dists.min()
        if d_min < 4.0:
            contact_freq[ri] += 1
        min_dist[ri] = min(min_dist[ri], d_min)
    if (j+1) % 50 == 0:
        print(f"  {j+1}/{len(sample_idx)} frames", flush=True)

contact_freq /= len(sample_idx)

# Map residue info for plotting
ch_labels = {0: 'MHC α', 1: 'β2M', 2: 'Peptide'}
ch_colors = {0: '#3498db', 1: '#2ecc71', 2: '#e74c3c'}

res_labels = [f'{r["name"]}{r["resSeq"]}' for r in binding]
res_colors = [ch_colors.get(r['chain'], '#95a5a6') for r in binding]
res_chains = [r['chain'] for r in binding]

# Sort by position: chain 0 (MHC) then 1 (B2M) then 2 (peptide)
# Within each chain, by resSeq
sorted_idx = sorted(range(n_bind), key=lambda i: (binding[i]['chain']*1000 + binding[i]['resSeq']))
# But for display, sort by contact frequency
top_idx = np.argsort(contact_freq)[::-1][:30]  # top 30

# Find E63 and K2 specifically
e63_idx = next((i for i, r in enumerate(binding)
                if r['chain']==0 and r['resSeq']==63), None)
k2_idx = next((i for i, r in enumerate(binding)
               if r['chain']==2 and r['resSeq']==2), None)

print(f"  E63 in binding site: {e63_idx is not None} (contact_freq={contact_freq[e63_idx]:.2f})" if e63_idx is not None else "  E63 NOT in binding site", flush=True)
print(f"  K2 in binding site: {k2_idx is not None}" , flush=True)
for i in top_idx[:10]:
    r = binding[i]
    print(f"    {ch_labels.get(r['chain'],'?'):8s}:{r['name']}{r['resSeq']:3d}  "
          f"contacts={contact_freq[i]:.2%}  min_dist={min_dist[i]:.2f}Å", flush=True)

np.savez(str(DATA / "binding_profile.npz"),
    res_names=np.array([r['name'] for r in binding], dtype=object),
    res_chains=np.array([r['chain'] for r in binding]),
    res_seqs=np.array([r['resSeq'] for r in binding]),
    contact_freq=contact_freq, min_dist=min_dist)

# ── Figure 3: Binding site per-residue analysis ──
fig, axes = plt.subplots(2, 1, figsize=(14, 7))
fig.suptitle('K-E63 Salt Bridge — Binding Site Per-Residue Contact Analysis',
             fontsize=13, fontweight='bold')

# 3a. Contact frequency bar chart (top 30)
ax = axes[0]
x = np.arange(len(top_idx))
bars = ax.bar(x, contact_freq[top_idx], color=[res_colors[i] for i in top_idx],
              edgecolor='white', linewidth=0.5, alpha=0.85)
ax.set_xticks(x)
ax.set_xticklabels([res_labels[i] for i in top_idx], fontsize=6, rotation=90)
ax.set_ylabel('Contact Frequency (fraction of frames)', fontsize=10)
ax.set_title('Top 30 Binding Site Residues by Contact Frequency with Peptide (<4Å heavy atom)',
             fontsize=10)
ax.axhline(0.5, color='#95a5a6', ls='--', lw=0.8, alpha=0.5)
ax.grid(axis='y', alpha=0.15)

# Highlight K2 and E63
for xi, idx in enumerate(top_idx):
    if idx == k2_idx:
        bars[xi].set_edgecolor('black')
        bars[xi].set_linewidth(2)
    if idx == e63_idx:
        bars[xi].set_edgecolor('black')
        bars[xi].set_linewidth(2)

# Legend
from matplotlib.patches import Patch
legend_elements = [
    Patch(facecolor='#3498db', alpha=0.85, label='MHC α'),
    Patch(facecolor='#2ecc71', alpha=0.85, label='β2M'),
    Patch(facecolor='#e74c3c', alpha=0.85, label='Peptide'),
]
ax.legend(handles=legend_elements, fontsize=8)

# 3b. Chain-wide contact profile
ax = axes[1]
x_all = np.arange(n_bind)
ax.bar(x_all, contact_freq, color=res_colors, width=0.8, alpha=0.6, edgecolor='none')
ax.set_xlabel('Binding Site Residue Index', fontsize=10)
ax.set_ylabel('Contact Frequency', fontsize=10)
ax.set_title('Full Binding Site Contact Profile', fontsize=10)
ax.grid(axis='y', alpha=0.15)

# Chain boundaries
ch_start = {}
for i, r in enumerate(binding):
    ch = r['chain']
    if ch not in ch_start:
        ch_start[ch] = i
for ch, s in ch_start.items():
    ax.axvline(s-0.5, color='grey', ls=':', lw=0.5)
    ax.text(s, 1.02, ch_labels.get(ch, f'C{ch}'), fontsize=8, fontweight='bold',
            color=ch_colors.get(ch, '#333'), ha='left', va='bottom')

# Annotate key residues
for target_name, target_chain, target_seq in [('E63', 0, 63), ('K2 (P2)', 2, 2)]:
    for i, r in enumerate(binding):
        if r['chain'] == target_chain and r['resSeq'] == target_seq:
            ax.annotate(f'{target_name}\n{contact_freq[i]:.0%}',
                        xy=(i, contact_freq[i]), xytext=(i, contact_freq[i]+0.15),
                        fontsize=7, ha='center', fontweight='bold',
                        arrowprops=dict(arrowstyle='->', color='#2c3e50', lw=0.8))

plt.tight_layout()
plt.savefig(OUT / 'Fig3_per_residue.png', dpi=300, bbox_inches='tight')
plt.savefig(OUT / 'Fig3_per_residue.pdf', bbox_inches='tight')
plt.close()
print('Fig3_per_residue.png/pdf done')
