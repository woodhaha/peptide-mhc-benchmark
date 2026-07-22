#!/usr/bin/env python3
"""RMSF for 1DUZ K-E63 salt bridge 100ns - Calpha only."""
import numpy as np, time
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

bb = traj.topology.select("backbone")
traj.superpose(traj, frame=0, atom_indices=bb)

# Calpha RMSF - iterate residues
rmsf_by_res = []
res_info = []
for ch in traj.topology.chains:
    for res in ch.residues:
        ca_atoms = [a.index for a in res.atoms if a.name == "CA"]
        if not ca_atoms:
            continue
        ca_traj = traj.atom_slice(ca_atoms)
        rmsf_val = float(np.sqrt(np.mean((ca_traj.xyz - ca_traj.xyz[0:1])**2)))
        rmsf_by_res.append(rmsf_val)
        res_info.append((ch.index, res.name, res.index))

rmsf = np.array(rmsf_by_res)
print(f"  {len(rmsf)} residues, RMSF: mean={rmsf.mean():.3f}A median={float(np.median(rmsf)):.3f}A", flush=True)

flex = np.argsort(rmsf)[::-1][:30]
print(f"\n  Top 30 flexible residues:")
for i in flex:
    ci, rn, ri = res_info[i]
    print(f"    {'ABC'[ci]}:{rn}{ri}  {rmsf[i]:.3f}A", flush=True)

# E63 specifically
for i, (ci, rn, ri) in enumerate(res_info):
    if rn == "GLU" and ri == 63 and ci == 0:
        print(f"\n  Target: A:GLU63 Calpha RMSF = {rmsf[i]:.3f}A", flush=True)
    if rn == "LYS" and ri == 2 and ci == 2:
        print(f"  Ligand: C:LYS2 Calpha RMSF = {rmsf[i]:.3f}A", flush=True)

# Per-chain
for ch in range(3):
    ch_vals = [rmsf[i] for i, (ci, rn, ri) in enumerate(res_info) if ci == ch]
    if ch_vals:
        print(f"  Chain {['A(MHC)','B(b2m)','C(peptide)'][ch]}: mean={float(np.mean(ch_vals)):.3f}A "
              f"n={len(ch_vals)} residues", flush=True)

np.savez(D/"rmsf.npz", rmsf=rmsf,
         chains=np.array([x[0] for x in res_info]),
         resnames=np.array([x[1] for x in res_info]),
         resnums=np.array([x[2] for x in res_info]))
print(f"\n  DONE [{(time.perf_counter()-t0)/60:.1f}min]", flush=True)
