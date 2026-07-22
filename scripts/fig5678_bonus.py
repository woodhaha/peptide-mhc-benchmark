#!/usr/bin/env python3
"""Bonus figures: H-bonds, Replica RMSF, PMF, DCCM.
Runs locally on md results/salt_bridge.
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path
from scipy import stats
import json

DATA = Path("D:/Researching/Peptide epitope/md results/salt_bridge")
REP2 = Path("D:/Researching/Peptide epitope/md results/salt_bridge_replica2")
REP3 = Path("D:/Researching/Peptide epitope/md results/salt_bridge_replica3")
OUT  = Path("D:/Researching/Peptide epitope/md pic")
KT = 0.596  # kcal/mol at 300K

plt.rcParams.update({'font.size': 9, 'figure.dpi': 200, 'savefig.bbox': 'tight'})

# ======================================================================
# 1. PMF from K-E63 distance distribution (Fig5)
# ======================================================================
print("[1/4] PMF from distance distribution...")
md_dist = np.load(str(DATA / "md_dist.npz"))
sb_dist = md_dist['mo']
t_ns = md_dist['t_ns']

# Boltzmann inversion: F(r) = -kT * ln(P(r))
bins = np.arange(3.0, 7.5, 0.05)
bin_centers = (bins[:-1] + bins[1:]) / 2
P, _ = np.histogram(sb_dist, bins=bins, density=True)
P[P < 1e-12] = 1e-12
FE = -KT * np.log(P)
FE -= FE.min()

fig, axes = plt.subplots(1, 2, figsize=(12, 5))
fig.suptitle('K-E63 Salt Bridge — Potential of Mean Force', fontsize=13, fontweight='bold')

ax = axes[0]
ax.plot(bin_centers, FE, 'o-', color='#2c3e50', ms=3, lw=1.5)
ax.axvline(4.0, color='#27ae60', ls='--', lw=1, alpha=0.7, label='Salt bridge cutoff (4.0Å)')
ax.axvline(7.6, color='#e74c3c', ls=':', lw=1, alpha=0.7, label='K-sidechain reach (7.6Å)')
# Mark minima
min_idx = np.argmin(FE)
ax.plot(bin_centers[min_idx], FE[min_idx], 'rD', ms=8, label=f'Minimum at {bin_centers[min_idx]:.2f}Å')
ax.set_xlabel('CA-K2 → E63 min distance (Å)')
ax.set_ylabel('PMF (kcal/mol)')
ax.set_title('1D PMF: F(r) = −kT ln[P(r)]')
ax.legend(fontsize=7)
ax.grid(alpha=0.15)

# Energy barrier
barrier_cutoff = 4.5
barrier_mask = bin_centers > barrier_cutoff
if barrier_mask.any():
    barrier_top = np.max(FE[~barrier_mask][-5:]) if (~barrier_mask).sum() > 5 else FE[-1]
    ax.annotate(f'Energy barrier\n~{barrier_top:.1f} kcal/mol',
                xy=(barrier_cutoff, barrier_top), fontsize=8, ha='center',
                bbox=dict(boxstyle='round', fc='#fee0d6', ec='#e74c3c'))

ax = axes[1]
# 2D view: distance vs time colored by energy
sc = ax.scatter(t_ns, sb_dist, c=-KT * np.log(np.maximum(np.interp(sb_dist, bin_centers, P), 1e-12)),
                cmap='viridis', s=5, alpha=0.6, edgecolors='none')
cbar = plt.colorbar(sc, ax=ax, label='Free Energy (kcal/mol)')
ax.axhline(4.0, color='#27ae60', ls='--', lw=0.8, alpha=0.5)
ax.set_xlabel('Time (ns)'); ax.set_ylabel('CA-K2 → E63 min distance (Å)')
ax.set_title('Trajectory Projected on PMF')
ax.grid(alpha=0.15)
ax.text(0.02, 0.95, f'ΔG‡ ≈ {FE[bin_centers > sb_dist.mean()].max() - FE.min():.1f} kcal/mol',
        transform=ax.transAxes, fontsize=9, va='top',
        bbox=dict(boxstyle='round', fc='white', alpha=0.7))

plt.tight_layout()
plt.savefig(OUT / 'Fig5_PMF.png', dpi=300, bbox_inches='tight')
plt.savefig(OUT / 'Fig5_PMF.pdf', bbox_inches='tight')
plt.close()
print("  Fig5_PMF.png/pdf done")

# ======================================================================
# 2. H-bond occupancy (Fig6)
# ======================================================================
print("\n[2/4] Hydrogen bond analysis...")
import mdtraj as md

top = md.load(str(DATA / "final.pdb"))
traj = md.load(str(DATA / "traj.dcd"), top=str(DATA / "final.pdb"))
# Extract only protein chains 0,1,2 for speed
prot_idx = [a.index for a in traj.topology.atoms
            if a.residue.chain.index in (0,1,2)]
traj_prot = traj.atom_slice(prot_idx)

# Get chain mapping for the sliced topology
# After slicing, we need map back to original chain IDs
chain_map = {}
for a_orig, a_new in zip([traj.topology.atom(i) for i in prot_idx],
                          traj_prot.topology.atoms):
    pass  # just verifying

# Compute H-bonds between peptide (chain 2) and MHC (chain 0)
# using the ORIGINAL topology to find the correct atom indices
print("  Computing H-bonds (Baker-Hubbard, 0.1 freq cutoff)...", flush=True)
hbonds = md.baker_hubbard(traj, freq=0.1, exclude_water=True)
print(f"  Total H-bonds found: {len(hbonds)}", flush=True)

# Filter: donor and acceptor must be in different chains
ac_hbonds = []
unique_pairs = {}
for hbond in hbonds:
    d, h, a = hbond[:3]
    d_res = traj.topology.atom(d).residue
    a_res = traj.topology.atom(a).residue
    d_chain = d_res.chain.index
    a_chain = a_res.chain.index
    # Peptide (chain 2) ↔ MHC (chain 0) only
    if set([d_chain, a_chain]) == {0, 2}:
        key = (d_chain, d_res.resSeq, a_chain, a_res.resSeq)
        if key not in unique_pairs:
            unique_pairs[key] = (d_res.name, a_res.name, str(d_res), str(a_res))

# Count occupancy per residue pair
hb_occupancy = {}
for hbond in hbonds:
    d, h, a = hbond[:3]
    d_res = traj.topology.atom(d).residue
    a_res = traj.topology.atom(a).residue
    if set([d_res.chain.index, a_res.chain.index]) == {0, 2}:
        pair = (d_res.chain.index, d_res.resSeq, d_res.name,
                a_res.chain.index, a_res.resSeq, a_res.name)
        hb_occupancy[pair] = hb_occupancy.get(pair, 0) + 1

# Compute framewise occupancy using the full trajectory H-bonds per frame
# mdtraj's baker_hubbard on full traj gives all hbonds that occur in ≥10% of frames
# For per-frame occupancy, we need a different approach
# Use the existing output format from sb_hbonds.py

# Use the full-frames approach with per-frame analysis (subsample 200 frames)
print("  Per-frame H-bond occupancy (200 frames)...", flush=True)
frame_idx = np.linspace(0, traj.n_frames-1, 200, dtype=int)
hbond_frames = md.baker_hubbard(traj[frame_idx], freq=0.0, exclude_water=True)
# freq=0 means show ALL hbonds found, not filtering by frequency
# But this still merges across frames... we need per-frame
# Let's compute pairwise minimum heavy-atom distance as H-bond proxy

# Simpler: compute per-residue pair minimum distances across trajectory
# to assess H-bond probability
peptide_chain2 = 2
mhc_chain = 0

# Identify backbone N and O atoms for MHC (chain 0) and peptide (chain 2)
mhc_bb_donors = []  # backbone N-H
mhc_bb_acceptors = []  # backbone O
pep_bb_donors = []
pep_bb_acceptors = []

for atom in traj.topology.atoms:
    if atom.residue.chain.index == 2:
        if atom.name == 'N':
            pep_bb_donors.append(atom.index)
        elif atom.name == 'O':
            pep_bb_acceptors.append(atom.index)
    elif atom.residue.chain.index == 0:
        if atom.name == 'N':
            mhc_bb_donors.append(atom.index)
        elif atom.name == 'O':
            mhc_bb_acceptors.append(atom.index)

# Compute distances for donor-acceptor pairs
# MHC donor → peptide acceptor, peptide donor → MHC acceptor
donor_sets = [('MHC→Pep', mhc_bb_donors, pep_bb_acceptors),
              ('Pep→MHC', pep_bb_donors, mhc_bb_acceptors)]

hb_data = {}  # (donor_resName_resSeq, acceptor_resName_resSeq) -> occupancy
# Sample 100 frames for speed
n_occ = len(frame_idx)
for label, donors, acceptors in donor_sets:
    for di in donors:
        d_atom = traj.topology.atom(di)
        d_res = d_atom.residue
        d_key = f"{d_res.name}{d_res.resSeq}"
        for ai in acceptors:
            a_atom = traj.topology.atom(ai)
            a_res = a_atom.residue
            a_key = f"{a_res.name}{a_res.resSeq}"
            key = (label, d_key, a_key, d_res.chain.index, a_res.chain.index,
                   d_res.resSeq, a_res.resSeq, d_atom.name, a_atom.name)

            count = 0
            for fi in frame_idx:
                d_xyz = traj.xyz[fi, di]
                a_xyz = traj.xyz[fi, ai]
                dist = np.linalg.norm(d_xyz - a_xyz)
                # H-bond: N-O distance < 3.5Å, N-H-O angle ~ 180°
                # Simple distance proxy: < 3.5Å
                if dist < 3.5:
                    count += 1

            occ = count / len(frame_idx)
            if occ > 0.05:  # >5% occupancy
                hb_data[key] = occ

# Sort by occupancy
sorted_hb = sorted(hb_data.items(), key=lambda x: x[1], reverse=True)
print(f"  Found {len(sorted_hb)} H-bond pairs (>5% occupancy):")
for (label, d_key, a_key, d_ch, a_ch, d_seq, a_seq, d_name, a_name), occ in sorted_hb[:15]:
    dir_label = "MHC→Pep" if d_ch == 0 else "Pep→MHC"
    print(f"    {dir_label} {d_key}({d_name}) → {a_key}({a_name}): {occ:.1%}")

# Save H-bond data
np.savez(str(DATA / "hbonds.npz"),
    n_pairs=len(sorted_hb),
    occ=np.array([o for _, o in sorted_hb]))

# Figure
fig, ax = plt.subplots(figsize=(12, max(5, len(sorted_hb[:20])*0.35)))
fig.suptitle('Peptide-MHC Hydrogen Bond Occupancy', fontsize=13, fontweight='bold')

labels_plot = []
occs = []
colors_hb = []
for (label, d_key, a_key, d_ch, a_ch, d_seq, a_seq, d_name, a_name), occ in sorted_hb[:20]:
    arrow = '→' if d_ch == 2 else '→'  # same
    labels_plot.append(f"{d_key}({d_name}) → {a_key}({a_name})")
    occs.append(occ)
    colors_hb.append('#e74c3c' if occ > 0.5 else '#f39c12' if occ > 0.2 else '#3498db')

y_pos = np.arange(len(labels_plot))
ax.barh(y_pos, occs, color=colors_hb, edgecolor='white', linewidth=0.5, height=0.6)
ax.set_yticks(y_pos)
ax.set_yticklabels(labels_plot, fontsize=7)
ax.invert_yaxis()
ax.set_xlabel('Occupancy (fraction of frames)')
ax.axvline(0.5, color='#27ae60', ls='--', lw=0.8, alpha=0.5)
ax.axvline(0.2, color='#f39c12', ls=':', lw=0.8, alpha=0.5)
ax.grid(axis='x', alpha=0.15)

plt.tight_layout()
plt.savefig(OUT / 'Fig6_hbonds.png', dpi=300, bbox_inches='tight')
plt.savefig(OUT / 'Fig6_hbonds.pdf', bbox_inches='tight')
plt.close()
print("  Fig6_hbonds.png/pdf done")

# ======================================================================
# 3. Replica RMSF comparison (Fig7)
# ======================================================================
print("\n[3/4] Replica RMSF comparison...")

rmsf_all = np.load(str(DATA / "rmsf_all.npz"))
r_main = rmsf_all['main']
r_rep2 = rmsf_all['rep2']
r_rep3 = rmsf_all['rep3']

# Get residue-level RMSF from per-atom data
# Map each atom to its residue
residue_names = []
residue_chains = []
residue_seqs = []

for atom in traj.topology.atoms:
    if atom.residue.chain.index <= 2:  # only protein chains
        residue_names.append(str(atom.residue.name))
        residue_chains.append(atom.residue.chain.index)
        residue_seqs.append(atom.residue.resSeq)

residue_names = np.array(residue_names)
residue_chains = np.array(residue_chains)
residue_seqs = np.array(residue_seqs)

# Main and rep2 have different atom counts (74345 vs 74195 vs 75238)
# But they should have similar protein regions. Align by chain/residue mapping.
# Use only the intersection of atoms that exist in all 3
# Since we can't align atom-by-atom across different topologies,
# compute per-residue mean RMSF for each trajectory using its own topology

# Actually the rmsf_all.npz is per-atom RMSF already computed on AutoDL
# The atom counts differ because of different water counts in GBSA vs PME
# We need per-residue averages for each

# Recompute per-residue RMSF from the saved per-atom data
# For the main traj:
def per_residue_rmsf(traj_obj, rmsf_array):
    """Compute per-residue mean RMSF from per-atom RMSF."""
    res_rmsf = {}
    for atom, rmsf_val in zip(traj_obj.topology.atoms, rmsf_array):
        res = atom.residue
        if res.chain.index > 2:
            continue
        key = (res.chain.index, res.resSeq, str(res.name))
        if key not in res_rmsf:
            res_rmsf[key] = []
        res_rmsf[key].append(rmsf_val)
    return res_rmsf

# Load rep trajectories
top_rep2 = md.load(str(REP2 / "final_s42.pdb"))
top_rep3 = md.load(str(REP3 / "final_s99.pdb"))

res_main = per_residue_rmsf(traj, r_main)

# For rep2/rep3 we don't have the trajectory objects loaded for their topologies
# Load them
rep2_atoms = sum(1 for _ in top_rep2.topology.atoms)
rep3_atoms = sum(1 for _ in top_rep3.topology.atoms)
print(f"  RMSF shapes: main={len(r_main)}, rep2={len(r_rep2)} ({rep2_atoms} atoms), rep3={len(r_rep3)} ({rep3_atoms} atoms)")

# Map rep2 RMSF to residues
res_rep2 = {}
for atom, rmsf_val in zip(top_rep2.topology.atoms, r_rep2):
    res = atom.residue
    if res.chain.index > 2:
        continue
    key = (res.chain.index, res.resSeq, str(res.name))
    res_rep2.setdefault(key, []).append(rmsf_val)

res_rep3 = {}
for atom, rmsf_val in zip(top_rep3.topology.atoms, r_rep3):
    res = atom.residue
    if res.chain.index > 2:
        continue
    key = (res.chain.index, res.resSeq, str(res.name))
    res_rep3.setdefault(key, []).append(rmsf_val)

# Find common residues across all 3 runs
common_keys = set(res_main.keys()) & set(res_rep2.keys()) & set(res_rep3.keys())
print(f"  Common residues across 3 runs: {len(common_keys)}")

# Build arrays
sorted_keys = sorted(common_keys, key=lambda k: (k[0]*1000 + k[1]))
r_main_res = np.array([np.mean(res_main[k]) for k in sorted_keys])
r_r2_res = np.array([np.mean(res_rep2[k]) for k in sorted_keys])
r_r3_res = np.array([np.mean(res_rep3[k]) for k in sorted_keys])
seq_labels = [f"{k[2]}{k[1]}" for k in sorted_keys]
chain_ids = [k[0] for k in sorted_keys]

# Per-chain averages
ch_labels = {0: 'MHC α', 1: 'β₂M', 2: 'Peptide'}
ch_colors_r = {0: '#3498db', 1: '#2ecc71', 2: '#e74c3c'}
for ch in [0, 1, 2]:
    mask_c = np.array(chain_ids) == ch
    print(f"  {ch_labels[ch]}: main={r_main_res[mask_c].mean():.2f}Å  rep2={r_r2_res[mask_c].mean():.2f}Å  rep3={r_r3_res[mask_c].mean():.2f}Å")

# Site-specific for E63 (chain 0, 63) and K2 (chain 2, 2)
for ch, seq, name in [(0, 63, 'A:E63'), (2, 2, 'C:K2')]:
    try:
        idx = sorted_keys.index((ch, seq, name.split(':')[1] if ':' in name else name))
        print(f"  {name}: main={r_main_res[idx]:.3f}Å  rep2={r_r2_res[idx]:.3f}Å  rep3={r_r3_res[idx]:.3f}Å")
    except ValueError:
        pass

# Figure: 3-panel RMSF comparison
fig, axes = plt.subplots(3, 1, figsize=(16, 8), sharex=True,
    gridspec_kw={'height_ratios': [2, 1, 1]})
fig.suptitle('K-E63 — Per-Residue Calpha RMSF: 3-Run Comparison (100ns each)',
             fontsize=13, fontweight='bold')

# Top: all runs overlaid
ax = axes[0]
x = np.arange(len(sorted_keys))
for lst, (name, data, color) in enumerate([
    ('Main (GBSA)', r_main_res, '#e74c3c'),
    ('Rep2 (PME)', r_r2_res, '#3498db'),
    ('Rep3 (PME)', r_r3_res, '#2ecc71')]):
    offset = (lst - 1) * 0.3
    ax.bar(x + offset, data, width=0.25, color=color, alpha=0.7, label=name)

# Per-chain shading
for ch, ch_color in ch_colors_r.items():
    ch_x = np.where(np.array(chain_ids) == ch)[0]
    if len(ch_x) > 0:
        ax.axvspan(ch_x[0]-1, ch_x[-1]+1, alpha=0.05, color=ch_color)

ax.set_ylabel('RMSF (Å)')
ax.legend(fontsize=8, ncol=3)
ax.set_title('Per-Residue RMSF (backbone, Calpha)')
ax.grid(axis='y', alpha=0.15)

# Middle: per-chain averages
ax = axes[1]
ch_idx = np.arange(3)
ch_colors_list = ['#3498db', '#2ecc71', '#e74c3c']
for lst, (name, data_arr) in enumerate([
    ('Main', r_main_res), ('Rep2', r_r2_res), ('Rep3', r_r3_res)]):
    means = [np.mean(data_arr[np.array(chain_ids) == ch]) for ch in [0, 1, 2]]
    ax.bar(ch_idx + lst*0.25, means, width=0.22, color=['#e74c3c','#3498db','#2ecc71'][lst],
           alpha=0.7, label=name)
ax.set_xticks(ch_idx + 0.25)
ax.set_xticklabels([f'{ch_labels[ch]}\n({(np.array(chain_ids)==ch).sum()} res)' for ch in [0,1,2]])
ax.set_ylabel('Mean RMSF (Å)')
ax.legend(fontsize=7)
ax.grid(axis='y', alpha=0.15)

# Bottom: RMSF difference (Rep - Main)
ax = axes[2]
ax.plot(x, r_r2_res - r_main_res, color='#3498db', lw=0.8, label='Rep2-Main Δ')
ax.plot(x, r_r3_res - r_main_res, color='#2ecc71', lw=0.8, alpha=0.7, label='Rep3-Main Δ')
ax.axhline(0, color='#7f8c8d', ls='-', lw=0.8)
ax.axhline(1.0, color='#95a5a6', ls='--', lw=0.5, alpha=0.5)
ax.axhline(-1.0, color='#95a5a6', ls='--', lw=0.5, alpha=0.5)
ax.set_xlabel('Residue Index')
ax.set_ylabel('ΔRMSF (Å)')
ax.set_ylim(-3, 3)
ax.legend(fontsize=7)
ax.grid(alpha=0.15)

plt.tight_layout()
plt.savefig(OUT / 'Fig7_RMSF_replica.png', dpi=300, bbox_inches='tight')
plt.savefig(OUT / 'Fig7_RMSF_replica.pdf', bbox_inches='tight')
plt.close()
print("  Fig7_RMSF_replica.png/pdf done")

# ======================================================================
# 4. DCCM — Dynamic Cross-Correlation Matrix (Fig8)
# ======================================================================
print("\n[4/4] DCCM plot...")

dccm = np.load(str(DATA / "dccm.npy"))
n_ca = dccm.shape[0]

# Get Ca residue info
ca_res_names = []
ca_res_seqs = []
ca_chains = []
ca_res_indices = []
for atom in traj.topology.atoms:
    if atom.name == 'CA' and atom.residue.chain.index <= 2:
        ca_res_names.append(str(atom.residue.name))
        ca_res_seqs.append(atom.residue.resSeq)
        ca_chains.append(atom.residue.chain.index)
        ca_res_indices.append(atom.residue.index)

n_ca_protein = len(ca_res_names)
print(f"  DCCM: {n_ca}x{n_ca} for {n_ca_protein} protein Ca atoms")

# Chain boundaries for plotting
ch_bounds = {}
for i, ch in enumerate(ca_chains):
    ch_bounds.setdefault(ch, []).append(i)

fig, axes = plt.subplots(1, 2, figsize=(14, 6.5))
fig.suptitle('K-E63 Salt Bridge — Dynamic Cross-Correlation (Cα, 100ns)',
             fontsize=13, fontweight='bold')

# 4a. Full DCCM
ax = axes[0]
im = ax.imshow(dccm, cmap='RdBu_r', vmin=-0.5, vmax=0.5, aspect='auto', interpolation='nearest')
plt.colorbar(im, ax=ax, label='Pearson Correlation', shrink=0.85)
# Chain boundary lines
for ch, idxs in ch_bounds.items():
    last = idxs[-1]
    ax.axhline(last+0.5, color='black', lw=0.8)
    ax.axvline(last+0.5, color='black', lw=0.8)
    mid = (idxs[0] + last) / 2
    labels = {0: 'MHC α', 1: 'β₂M', 2: 'Peptide'}
    ax.text(mid, -3, labels.get(ch, f'C{ch}'), fontsize=7, ha='center', va='bottom', rotation=45)
ax.set_xlabel('Cα residue index'); ax.set_ylabel('Cα residue index')
ax.set_title(f'Full DCCM ({n_ca} Cα atoms)')

# 4b. Zoom around binding pocket (MHC 50-80 vs Peptide)
# Find Ca indices
pocket_mhc = [i for i, (ch, seq) in enumerate(zip(ca_chains, ca_res_seqs))
              if ch == 0 and 55 <= seq <= 80]
pocket_peptide = [i for i, (ch, seq) in enumerate(zip(ca_chains, ca_res_seqs))
                  if ch == 2]
pocket_all = pocket_mhc + pocket_peptide
pocket_labels = [f'{ca_res_names[i]}{ca_res_seqs[i]}' for i in pocket_all]

if len(pocket_all) > 1:
    dccm_pocket = dccm[np.ix_(pocket_all, pocket_all)]
    ax = axes[1]
    im = ax.imshow(dccm_pocket, cmap='RdBu_r', vmin=-0.5, vmax=0.5, aspect='auto', interpolation='nearest')
    plt.colorbar(im, ax=ax, label='Pearson Correlation', shrink=0.85)
    ax.set_xticks(range(len(pocket_all)))
    ax.set_yticks(range(len(pocket_all)))
    ax.set_xticklabels(pocket_labels, fontsize=5, rotation=90)
    ax.set_yticklabels(pocket_labels, fontsize=5)
    # Separator line between MHC and peptide
    sep = len(pocket_mhc) - 0.5
    ax.axhline(sep, color='black', lw=1)
    ax.axvline(sep, color='black', lw=1)
    ax.set_xlabel('Residue'); ax.set_ylabel('Residue')
    ax.set_title('Binding Pocket Zoom: MHC(55-80) × Peptide')

    # Print key cross-correlations
    print(f"  Pocket DCCM: {len(pocket_mhc)} MHC × {len(pocket_peptide)} peptide Ca")
    # K2 (peptide resSeq 2) - E63 (MHC resSeq 63)
    for i, pi in enumerate(pocket_all):
        for j, pj in enumerate(pocket_all):
            if (ca_chains[pi] == 0 and ca_res_seqs[pi] == 63 and
                ca_chains[pj] == 2 and ca_res_seqs[pj] == 2):
                print(f"  A:E63(CA) × C:K2(CA) correlation: {dccm_pocket[i,j]:.3f}")

plt.tight_layout()
plt.savefig(OUT / 'Fig8_DCCM.png', dpi=300, bbox_inches='tight')
plt.savefig(OUT / 'Fig8_DCCM.pdf', bbox_inches='tight')
plt.close()
print("  Fig8_DCCM.png/pdf done")

print(f"\nAll bonus figures done → {OUT}")
