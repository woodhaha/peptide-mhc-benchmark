#!/usr/bin/env python3
"""Extract Cα PCA projections from DCD for all compounds."""
import numpy as np, os
import MDAnalysis as mda
from MDAnalysis.analysis import pca

COMPOUNDS = ["NAFTAZONE","BIOTIN","ESKETAMINE","FUROSEMIDE",
             "GABAPENTIN_ENACARBIL","HYDROMORPHONE","PHENOBARBITAL","RIBOFLAVIN"]
BASE = "/root/autodl-tmp/SMVT_AutoDL_Package/trajectories"

for c in COMPOUNDS:
    d = f"{BASE}/{c}"
    top = f"{d}/{c}_final.pdb"
    dcd = f"{d}/{c}_100ns.dcd"
    out = f"{d}/{c}_pca.npz"
    if os.path.exists(out):
        print(f"SKIP {c}: exists"); continue

    print(f"\n=== {c} ===", flush=True)
    u = mda.Universe(top, dcd)
    n = len(u.trajectory)

    # Use Cα for PCA
    ca = u.select_atoms("name CA")
    n_ca = ca.n_atoms
    n_dim = n_ca * 3
    print(f"  {n} frames, {n_ca} Cα atoms ({n_dim} dims)", flush=True)

    # Build trajectory matrix: frames × (N*3)
    XYZ = np.empty((n, n_dim))
    for i, ts in enumerate(u.trajectory):
        XYZ[i] = ca.positions.flatten()
        if (i+1) % 200 == 0:
            print(f"  {i+1}/{n} loaded", flush=True)

    # Centering
    mean_xyz = XYZ.mean(axis=0)
    XYZ_c = XYZ - mean_xyz

    # PCA via SVD (faster than eigendecomposition for n_frames < n_dims usually)
    # For n_frames << n_dims (999 << ~600), SVD is efficient
    U, S, Vt = np.linalg.svd(XYZ_c, full_matrices=False)
    # Projections: frames × PCs
    projections = U * S  # same as XYZ_c @ Vt.T
    pc1 = projections[:, 0]
    pc2 = projections[:, 1]
    var_exp = (S**2) / (S**2).sum()

    print(f"  PC1: {var_exp[0]*100:.1f}%  PC2: {var_exp[1]*100:.1f}%  "
          f"Total: {(var_exp[0]+var_exp[1])*100:.1f}%", flush=True)

    # Save
    np.savez_compressed(out, pc1=pc1, pc2=pc2, var_exp=var_exp[:10])
    print(f"  -> {out}", flush=True)

print("\nALL DONE", flush=True)
