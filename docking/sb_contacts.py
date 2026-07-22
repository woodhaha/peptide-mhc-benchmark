#!/usr/bin/env python3
"""1DUZ E63-P2Leu interface: contact distance, H-bonds, pocket occupancy."""
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
print(f"  {traj.n_frames}frames {traj.n_atoms}atoms", flush=True)

# E63 sidechain heavy atoms, P2-leu all heavy atoms
e63 = [978, 979, 980]  # CD, OE1, OE2
p2 =  [6004, 6006, 6008, 6009, 6010, 6013, 6015, 6019]  # N,CA,C,O,CB,CG,CD1,CD2

# 1. E63-P2Leu min distance per frame
print("[1] E63-P2Leu distance...", flush=True)
d_min = np.zeros(traj.n_frames)
for f in range(traj.n_frames):
    best = 999.0
    for ei in e63:
        for pi in p2:
            d = np.linalg.norm(traj.xyz[f, ei] - traj.xyz[f, pi])
            if d < best: best = d
    d_min[f] = best

print(f"  E63-P2Leu: mean={d_min.mean():.3f}+-{d_min.std():.3f}A", flush=True)
print(f"  <3.4A: {(d_min<3.4).mean()*100:.1f}%  <4.0A: {(d_min<4.0).mean()*100:.1f}%", flush=True)

# 2. A-C H-bonds
print("[2] A-C H-bonds...", flush=True)
hbonds = md.baker_hubbard(traj, freq=0.1, exclude_water=True)
a_set = set(a.index for a in traj.topology.atoms if a.residue.chain.index == 0)
c_set = set(a.index for a in traj.topology.atoms if a.residue.chain.index == 2)
ac_hb = []
for hb in hbonds:
    d, h, a_id = int(hb[0]), int(hb[1]), int(hb[2])
    if (d in a_set and a_id in c_set) or (d in c_set and a_id in a_set):
        ac_hb.append(hb)
print(f"  Found {len(ac_hb)} A-C H-bonds", flush=True)

# 3. P2 contact partners (>50% occupancy, <3.5A)
print("[3] P2 pocket partners...", flush=True)
contact_cnt = {}
for f in range(0, traj.n_frames, 10):
    for a_idx in range(traj.n_atoms):
        a = traj.topology.atom(a_idx)
        if a.residue.chain.index != 0 or a.element.symbol == "H":
            continue
        for pi in p2:
            if np.linalg.norm(traj.xyz[f, a_idx] - traj.xyz[f, pi]) < 3.5:
                k = (a.residue.name, a.residue.resSeq)
                contact_cnt[k] = contact_cnt.get(k, 0) + 1
                break

nsamp = traj.n_frames // 10
for (rn, rs), c in sorted(contact_cnt.items(), key=lambda x: -x[1]):
    occ = c/nsamp*100
    if occ > 40:
        bar = "█" * int(occ/10)
        print(f"    {rn}{rs}  {occ:5.1f}%  {bar}", flush=True)

np.savez(D/"sb_contacts.npz", d_min=d_min)
json.dump({"mean":float(d_min.mean()),"std":float(d_min.std()),
    "p34":float((d_min<3.4).mean()*100),"p40":float((d_min<4.0).mean()*100),
    "n_hbonds":len(ac_hb)}, open(D/"sb_contacts_summary.json","w"), indent=2)
print(f"  DONE [{(time.perf_counter()-t0)/60:.1f}min]", flush=True)
