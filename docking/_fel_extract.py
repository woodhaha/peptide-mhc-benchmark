"""Extract RMSD + Rg + SASA from DCD for all compounds for FEL."""
import numpy as np, os
import MDAnalysis as mda
from MDAnalysis.analysis import rms

COMPOUNDS = ["BIOTIN","ESKETAMINE","FUROSEMIDE","GABAPENTIN_ENACARBIL",
             "HYDROMORPHONE","PHENOBARBITAL","RIBOFLAVIN","NAFTAZONE"]
BASE = "/root/autodl-tmp/SMVT_AutoDL_Package/trajectories"

for c in COMPOUNDS:
    d = f"{BASE}/{c}"
    top = f"{d}/{c}_final.pdb"
    dcd = f"{d}/{c}_100ns.dcd"
    out = f"{d}/{c}_fel.npz"
    if os.path.exists(out):
        print(f"SKIP {c}: exists"); continue
    print(f"\n=== {c} ===", flush=True)
    u = mda.Universe(top, dcd)
    n = len(u.trajectory)
    u.trajectory[0]; ref = u.atoms.positions.copy()
    bb = u.select_atoms("backbone")
    ca = u.select_atoms("name CA")
    prot = u.select_atoms("protein")
    ref_bb = ref[bb.indices]
    rmsd_v = np.empty(n); rg_v = np.empty(n)
    for i, ts in enumerate(u.trajectory):
        rmsd_v[i] = rms.rmsd(bb.positions, ref_bb, superposition=True)
        rg_v[i] = prot.radius_of_gyration()
        if (i+1) % 400 == 0:
            print(f"  {i+1}/{n} t={i*100/n:.0f}ns RMSD={rmsd_v[i]:.2f} Rg={rg_v[i]:.2f}", flush=True)
    np.savez_compressed(out, rmsd=rmsd_v, rg=rg_v)
    print(f"  -> {out}  RMSD={rmsd_v.mean():.2f}+/-{rmsd_v.std():.2f}  Rg={rg_v.mean():.2f}+/-{rg_v.std():.2f}", flush=True)
print("\nALL DONE", flush=True)
