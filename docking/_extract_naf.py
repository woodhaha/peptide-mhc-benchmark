#!/usr/bin/env python3
"""Extract full 100ns data from NAFTAZONE DCD + topology, save as clean CSV."""
import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis import rms
import sys, os

DATA = "/root/autodl-tmp/SMVT_AutoDL_Package/trajectories/NAFTAZONE"
TOP  = f"{DATA}/protein.pdb"
DCD  = f"{DATA}/NAFTAZONE_100ns.dcd"
OUT  = f"{DATA}/NAFTAZONE_100ns_full.csv"

print(f"Loading {DCD}...", flush=True)
u = mda.Universe(TOP, DCD)
n_frames = len(u.trajectory)
print(f"  {n_frames} frames, {u.atoms.n_atoms} atoms", flush=True)

# Backbone RMSD
bb = u.select_atoms("backbone")
ref = mda.Universe(TOP)
ref_bb = ref.select_atoms("backbone")
rmsd_vals = np.empty(n_frames)
time_ns = np.empty(n_frames)

for i, ts in enumerate(u.trajectory):
    rmsd_vals[i] = rms.rmsd(bb.positions, ref_bb.positions, superposition=True)
    time_ns[i] = ts.time / 1000  # ps -> ns
    if (i+1) % 500 == 0:
        print(f"  {i+1}/{n_frames} frames ({time_ns[i]:.0f}ns)", flush=True)

# Simple T approximation from kinetic energy not available without velocities
# Write time + RMSD
out = np.column_stack([time_ns, rmsd_vals])
header = "Time_ns,RMSD_A"
np.savetxt(OUT, out, delimiter=",", header=header, comments="")
print(f"Saved {OUT} — {len(time_ns)} rows", flush=True)
print(f"RMSD: {rmsd_vals.mean():.3f}±{rmsd_vals.std():.3f} Å  [{rmsd_vals.min():.2f}–{rmsd_vals.max():.2f}]", flush=True)
