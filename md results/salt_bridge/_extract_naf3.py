#!/usr/bin/env python3
import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis import rms

DATA = "/root/autodl-tmp/SMVT_AutoDL_Package/trajectories/NAFTAZONE"
TOP = DATA + "/NAFTAZONE_final.pdb"
DCD = DATA + "/NAFTAZONE_100ns.dcd"
OUT = DATA + "/NAFTAZONE_100ns_full.csv"

print("Loading...", flush=True)
u = mda.Universe(TOP, DCD)
n = len(u.trajectory)
print(f"  {n} frames, {u.atoms.n_atoms} atoms", flush=True)

# Use FIRST frame as reference
u.trajectory[0]
ref_pos = u.atoms.positions.copy()
bb = u.select_atoms("backbone")

time_ns = np.empty(n)
rmsd = np.empty(n)

for i, ts in enumerate(u.trajectory):
    # ponytail: 0.1ns per frame = 100ns / 1000 frames (2fs*25000steps = 50ps wait no...)
    # Actually dt in DCD is stored per step. Let me just use i * (100/n) 
    time_ns[i] = i * 100.0 / n
    bb_now = bb.positions.copy()
    # RMSD vs first frame
    rmsd[i] = rms.rmsd(bb_now, ref_pos[bb.indices], superposition=True)
    if (i+1) % 200 == 0:
        print(f"  {i+1}/{n} t={time_ns[i]:.0f}ns RMSD={rmsd[i]:.3f}", flush=True)

np.savetxt(OUT, np.column_stack([time_ns, rmsd]), delimiter=",",
           header="Time_ns,RMSD_A", comments="")
print(f"\nSaved {OUT}")
print(f"RMSD: {rmsd.mean():.3f}+/-{rmsd.std():.3f} A [{rmsd.min():.2f}-{rmsd.max():.2f}]")
print(f"Time: {time_ns[0]:.1f}-{time_ns[-1]:.0f}ns")
