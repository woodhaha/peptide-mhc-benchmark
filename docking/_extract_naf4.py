"""Extract RMSD from DCD — use lists, not pre-allocated arrays."""
import numpy as np, sys
import MDAnalysis as mda
from MDAnalysis.analysis import rms

DATA = "/root/autodl-tmp/SMVT_AutoDL_Package/trajectories/NAFTAZONE"
u = mda.Universe(DATA + "/NAFTAZONE_final.pdb", DATA + "/NAFTAZONE_100ns.dcd")
n = len(u.trajectory)
u.trajectory[0]
ref = u.atoms.positions.copy()
bb = u.select_atoms("backbone")
ref_bb = ref[bb.indices]
print(f"{n} frames, {u.atoms.n_atoms} atoms")

times, rmsds = [], []
for i, ts in enumerate(u.trajectory):
    times.append(i * 100.0 / n)
    rmsds.append(float(rms.rmsd(bb.positions, ref_bb, superposition=True)))
    if (i+1) % 200 == 0:
        print(f"  {i+1}/{n}  t={times[-1]:.0f}ns  RMSD={rmsds[-1]:.3f}")

out = np.column_stack([np.array(times), np.array(rmsds)])
np.savetxt(DATA + "/NAFTAZONE_100ns_full.csv", out, delimiter=",",
           header="Time_ns,RMSD_A", comments="")
print(f"\nSaved. RMSD={np.mean(rmsds):.3f}+/-{np.std(rmsds):.3f}A [{np.min(rmsds):.2f}-{np.max(rmsds):.2f}]")
print(f"Time: {times[0]:.1f}-{times[-1]:.0f}ns")
