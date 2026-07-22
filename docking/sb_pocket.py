#!/usr/bin/env python3
"""1DUZ E63 binding pocket analysis: residues within 5/8A of E63."""
import numpy as np, time, json
from pathlib import Path
import mdtraj as md

D = Path("/root/autodl-tmp/docking")
t0 = time.perf_counter()

traj = md.load_dcd(str(D / "traj.dcd"), top=str(D / "final.pdb"))
# Remove water + ions
keep = [a.index for a in traj.topology.atoms 
    if not a.residue.is_water and a.residue.name not in ("NA","CL","K","MG","CA")]
traj = traj.atom_slice(keep)
print(f"  {traj.n_frames}frames {traj.n_atoms}atoms", flush=True)

# E63: chain 0, resSeq 63, atoms CD(978), OE1(979), OE2(980)
e63_idx = [978, 979, 980]

# 1. Per-frame: which residues have <5A contact with E63 any atom
print("[1] E63 contact partners (<5A, >50% occupancy)...", flush=True)
contact_occ = {}
for f in range(0, traj.n_frames, 10):
    for ai in range(traj.n_atoms):
        a = traj.topology.atom(ai)
        if a.residue.chain.index == 0 and a.element.symbol != "H":
            for ei in e63_idx:
                d = np.linalg.norm(traj.xyz[f, ai] - traj.xyz[f, ei])
                if d < 5.0:
                    k = (a.residue.name, a.residue.resSeq, a.name)
                    contact_occ[k] = contact_occ.get(k, 0) + 1
                    break

nsamp = traj.n_frames // 10
print(f"\n  Pocket residues (heavy contact with E63 >50%):")
pocket_5A = []
for (rn, rs, an), c in sorted(contact_occ.items(), key=lambda x: -x[1]):
    occ = c / nsamp * 100
    if occ > 50:
        bar = "█" * int(occ/20)
        print(f"    {rn}{rs}({an})  {occ:5.1f}%  {bar}", flush=True)
        pocket_5A.append((rn, rs, an))

# 2. Unique residues (group by resSeq)
unique_res = {}
for rn, rs, an in pocket_5A:
    if rs not in unique_res:
        unique_res[rs] = rn
print(f"\n  Unique residues in pocket: {len(unique_res)}")
print(f"    ", sorted([f"{rn}{rs}" for rs, rn in unique_res.items()]))

# 3. E63 OE1/OE2 distances to each pocket residue (average)
print(f"\n[2] E63-min distance to each pocket residue heavy atoms:")
pocket_dist = {}
for rs in unique_res:
    r_atoms = [ai for ai in range(traj.n_atoms) 
               if traj.topology.atom(ai).residue.resSeq == rs
               and traj.topology.atom(ai).residue.chain.index == 0
               and traj.topology.atom(ai).element.symbol != "H"]
    if not r_atoms:
        continue
    d_res = np.zeros(traj.n_frames)
    for f in range(traj.n_frames):
        best = 999.0
        for ei in e63_idx:
            for ri in r_atoms:
                d = np.linalg.norm(traj.xyz[f, ei] - traj.xyz[f, ri])
                if d < best: best = d
        d_res[f] = best
    pocket_dist[f"{unique_res[rs]}{rs}"] = {
        "mean": float(d_res.mean()), "std": float(d_res.std()),
        "min": float(d_res.min())
    }
    print(f"    {unique_res[rs]}{rs}: {d_res.mean():.2f}+-{d_res.std():.2f}A (min={d_res.min():.2f}A)")

# 4. Peptide-MHC contact surface: residues on chain A within 4A of chain C
print(f"\n[3] P2Leu binding pocket (chain A residues contacting P2Leu >50% <4A):")
p2_res = [a.index for a in traj.topology.atoms 
          if a.residue.chain.index == 2 and a.residue.resSeq == 2]
p2_pocket = {}
for f in range(0, traj.n_frames, 10):
    for ai in range(traj.n_atoms):
        a = traj.topology.atom(ai)
        if a.residue.chain.index != 0: continue
        for pi in p2_res:
            if np.linalg.norm(traj.xyz[f, ai] - traj.xyz[f, pi]) < 4.0:
                k = (a.residue.name, a.residue.resSeq)
                p2_pocket[k] = p2_pocket.get(k, 0) + 1
                break

for (rn, rs), c in sorted(p2_pocket.items(), key=lambda x: -x[1]):
    occ = c/nsamp*100
    if occ > 50:
        print(f"    {rn}{rs}: {occ:.0f}%", flush=True)

# 5. Distance matrix: E63 OE to nearest P2 backbone atoms
print(f"\n[4] E63-P2 backbone distances (100ns avg):")
for pi_name, pi_idx in [("N",6004),("CA",6006),("C",6008),("O",6009)]:
    d_all = np.zeros(traj.n_frames)
    for f in range(traj.n_frames):
        best = min(np.linalg.norm(traj.xyz[f, ei] - traj.xyz[f, pi_idx]) for ei in e63_idx)
        d_all[f] = best
    print(f"    P2Leu({pi_name}) - E63: {d_all.mean():.2f}+-{d_all.std():.2f}A", flush=True)

# Save
results = {
    "n_pocket_residues": len(unique_res),
    "pocket_residues": [f"{rn}{rs}" for rs, rn in sorted(unique_res.items())],
    "pocket_distances": pocket_dist,
}
json.dump(results, open(D/"pocket_analysis.json","w"), indent=2)
print(f"\n  DONE [{(time.perf_counter()-t0)/60:.1f}min]", flush=True)
