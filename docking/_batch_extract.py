"""Extract RMSD from DCD for all compounds."""
import numpy as np, sys, os
import MDAnalysis as mda
from MDAnalysis.analysis import rms

COMPOUNDS = ["BIOTIN","ESKETAMINE","FUROSEMIDE","GABAPENTIN_ENACARBIL",
             "HYDROMORPHONE","PHENOBARBITAL","RIBOFLAVIN"]
BASE = "/root/autodl-tmp/SMVT_AutoDL_Package/trajectories"

for c in COMPOUNDS:
    d = f"{BASE}/{c}"
    top = f"{d}/{c}_final.pdb"
    dcd = f"{d}/{c}_100ns.dcd"
    out = f"{d}/{c}_100ns_full.csv"
    
    if not os.path.exists(dcd):
        print(f"SKIP {c}: no DCD"); continue
    if os.path.exists(out):
        sz = os.path.getsize(out)
        if sz > 100:
            print(f"SKIP {c}: already extracted ({sz}b)"); continue
    
    print(f"\n=== {c} ===", flush=True)
    u = mda.Universe(top, dcd)
    n = len(u.trajectory)
    u.trajectory[0]
    ref = u.atoms.positions.copy()
    bb = u.select_atoms("backbone")
    ref_bb = ref[bb.indices]
    print(f"  {n} frames, {u.atoms.n_atoms} atoms", flush=True)
    
    times, rmsds = [], []
    for i, ts in enumerate(u.trajectory):
        times.append(i * 100.0 / n)
        rmsds.append(float(rms.rmsd(bb.positions, ref_bb, superposition=True)))
        if (i+1) % 400 == 0:
            print(f"  {i+1}/{n} t={times[-1]:.0f}ns RMSD={rmsds[-1]:.2f}", flush=True)
    
    out_a = np.column_stack([np.array(times), np.array(rmsds)])
    np.savetxt(out, out_a, delimiter=",", header="Time_ns,RMSD_A", comments="")
    print(f"  -> {out}  RMSD={np.mean(rmsds):.2f}+/-{np.std(rmsds):.2f}A [{np.min(rmsds):.2f}-{np.max(rmsds):.2f}]", flush=True)

print("\n=== ALL DONE ===", flush=True)
