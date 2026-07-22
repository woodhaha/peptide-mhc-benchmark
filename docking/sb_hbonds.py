#!/usr/bin/env python3
"""Hydrogen bond + contact analysis for 1DUZ K-E63 salt bridge 100ns."""
import numpy as np, time, json
from pathlib import Path
import mdtraj as md

D = Path("/root/autodl-tmp/docking")
t0 = time.perf_counter()

traj = md.load_dcd(str(D / "traj.dcd"), top=str(D / "final.pdb"))
solvent = [a.index for a in traj.topology.atoms if a.residue.is_water]
if solvent:
    keep = [i for i in range(traj.n_atoms) if i not in solvent]
    traj = traj.atom_slice(keep)
print(f"  {traj.n_frames} frames, {traj.n_atoms} atoms", flush=True)

# 1. Salt bridge pair: E63 (A:63) - K2 (C:2)
# Find atom indices
e63_oe = []; k2_nz = []
for a in traj.topology.atoms:
    r = a.residue
    if r.chain.index == 0 and r.resSeq == 63 and a.name in ("OE1","OE2"):
        e63_oe.append(a.index)
    if r.chain.index == 2 and r.resSeq == 2 and a.name == "NZ":
        k2_nz.append(a.index)

print(f"  E63 OE atoms: {e63_oe}", flush=True)
print(f"  K2 NZ atoms:  {k2_nz}", flush=True)

if e63_oe and k2_nz:
    # Compute minimum distance E63 OE - K2 NZ per frame
    d_sb = np.zeros(traj.n_frames)
    for f in range(traj.n_frames):
        oe_pos = traj.xyz[f, e63_oe, :]
        nz_pos = traj.xyz[f, k2_nz[0], :]
        d_sb[f] = min(np.linalg.norm(oe_pos[0] - nz_pos), np.linalg.norm(oe_pos[-1] - nz_pos))
    
    print(f"\n  E63(OE)-K2(NZ) distance: mean={d_sb.mean():.3f}+-{d_sb.std():.3f}A", flush=True)
    print(f"  Salt bridge (<4.0A): {(d_sb<4.0).mean()*100:.1f}%", flush=True)

# 2. All hydrogen bonds between chain A and chain C
print("\n  Computing H-bonds (chain A-C interface)...", flush=True)
hbonds = md.baker_hubbard(traj, freq=0.1, exclude_water=True)
# Filter to A-C contacts
a_idx = set(a.index for a in traj.topology.atoms if a.residue.chain.index == 0)
c_idx = set(a.index for a in traj.topology.atoms if a.residue.chain.index == 2)

ac_hbonds = []
for hbond in hbonds:
    d, h, a = hbond[:3]
    if (d in a_idx and a in c_idx) or (d in c_idx and a in a_idx):
        ac_hbonds.append(hbond)

print(f"  Total Hbonds (A-C): {len(ac_hbonds)}", flush=True)
if ac_hbonds:
    ac_hbonds = np.array(ac_hbonds)
    # Get unique residue pair labels
    pairs = set()
    for hb in ac_hbonds:
        d_idx, h_idx, a_idx_val = int(hb[0]), int(hb[1]), int(hb[2])
        d_res = traj.topology.atom(d_idx).residue
        a_res = traj.topology.atom(a_idx_val).residue
        pairs.add((d_res.chain.index, d_res.name, d_res.resSeq,
                   a_res.chain.index, a_res.name, a_res.resSeq))
    for p in sorted(pairs):
        print(f"    C{p[0]}:{p[1]}{p[2]} - C{p[3]}:{p[4]}{p[5]}")

# 3. Per-frame contact map: residues within 3.5A of each other (A-C interface)
print("\n  Computing A-C residue contacts (CA<8A)...", flush=True)
ca_idx_a = [a.index for a in traj.topology.atoms if a.name == "CA" and a.residue.chain.index == 0]
ca_idx_c = [a.index for a in traj.topology.atoms if a.name == "CA" and a.residue.chain.index == 2]

if ca_idx_a and ca_idx_c:
    ca_a = traj.atom_slice(ca_idx_a)
    ca_c = traj.atom_slice(ca_idx_c)
    # Compute pairwise distances between all A CA and C CA
    contact_matrix = np.zeros((len(ca_idx_c), traj.n_frames))
    for i, ci in enumerate(ca_idx_c):
        c_pos = traj.xyz[:, ci:ci+1, :]  # (n_frames, 1, 3)
        a_pos = traj.xyz[:, ca_idx_a, :]  # (n_frames, n_A, 3)
        dists = np.sqrt(((c_pos - a_pos)**2).sum(axis=2))  # (n_frames, n_A)
        contact_matrix[i] = (dists.min(axis=1) < 8.0)

    # Occupancy per peptide residue
    print(f"  Peptide residue contacts with MHC (CA<8A, 100ns):")
    for i, ci in enumerate(ca_idx_c):
        a = traj.topology.atom(ci)
        r = a.residue
        occ = contact_matrix[i].mean() * 100
        bar = "█" * int(occ / 5)
        print(f"    C{r.resSeq}:{r.name}  {occ:5.1f}%  {bar}")

# Save
results = {
    "e63_k2_dist_mean": float(d_sb.mean()) if e63_oe and k2_nz else None,
    "e63_k2_dist_std": float(d_sb.std()) if e63_oe and k2_nz else None,
    "sb_pct": float((d_sb<4.0).mean()*100) if e63_oe and k2_nz else None,
    "n_ac_hbonds": len(ac_hbonds),
}
json.dump(results, open(D/"sb_contacts.json","w"), indent=2)
print(f"\n  DONE [{(time.perf_counter()-t0)/60:.1f}min]", flush=True)
