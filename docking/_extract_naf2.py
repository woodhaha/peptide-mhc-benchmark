#!/usr/bin/env python3
import numpy as np, pandas as pd
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

bb = u.select_atoms("backbone")
ref = mda.Universe(TOP); ref_bb = ref.select_atoms("backbone")
rmsd_v = np.empty(n); tn = np.empty(n)

for i, ts in enumerate(u.trajectory):
    rmsd_v[i] = rms.rmsd(bb.positions, ref_bb.positions, superposition=True)
    tn[i] = ts.time / 1000
    if (i+1) % 500 == 0:
        print(f"  {i+1}/{n} ({tn[i]:.0f}ns): RMSD={rmsd_v[i]:.3f}", flush=True)

np.savetxt(OUT, np.column_stack([tn, rmsd_v]), delimiter=",",
           header="Time_ns,RMSD_A", comments="")
print(f"\nSaved {OUT}")
print(f"RMSD: {rmsd_v.mean():.3f}±{rmsd_v.std():.3f} A [{rmsd_v.min():.2f}-{rmsd_v.max():.2f}]")
print(f"Range: {tn[0]:.1f}-{tn[-1]:.0f}ns ({len(tn)} frames)")
