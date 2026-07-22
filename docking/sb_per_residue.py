#!/usr/bin/env python3
"""1DUZ interaction energy decomposition."""
import numpy as np, time, json
from pathlib import Path
import mdtraj as md
import openmm as mm, openmm.app as app, openmm.unit as unit

D = Path("/root/autodl-tmp/docking")
N_SKIP = 20
t0 = time.perf_counter()

def make_gb_sys(top):
    ff = app.ForceField("amber14-all.xml", "amber14/tip3pfb.xml")
    sys = ff.createSystem(top, nonbondedMethod=app.NoCutoff, constraints=None)
    charges = []
    for f in sys.getForces():
        if isinstance(f, mm.NonbondedForce):
            for i in range(f.getNumParticles()):
                q, sig, eps = f.getParticleParameters(i)
                charges.append(q)
            break
    gb = mm.GBSAOBCForce()
    for i in range(top.getNumAtoms()):
        c = float(charges[i].value_in_unit(unit.elementary_charge)) if i < len(charges) else 0.0
        gb.addParticle(c, 0.15, 0.5)
    gb.setSolventDielectric(78.5); gb.setSoluteDielectric(1.0)
    sys.addForce(gb)
    return sys

def gb_energy(sim, pos):
    sim.context.setPositions(pos)
    return sim.context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilocalories_per_mole)

traj = md.load_dcd(str(D / "traj.dcd"), top=str(D / "final.pdb"))
# Remove water AND ions
keep_idx = [a.index for a in traj.topology.atoms 
    if not a.residue.is_water and a.residue.name not in ("NA","CL","K","MG","CA")]
traj = traj.atom_slice(keep_idx)
n_tot = traj.n_frames
n = max(1, n_tot // N_SKIP)
idx = np.linspace(0, n_tot-1, n, dtype=int)
print(f"  {n_tot}f -> {n} ({traj.n_atoms} atoms)", flush=True)

lig_idx = [a.index for a in traj.topology.atoms if a.residue.chain.index == 2]
rec_idx = [a.index for a in traj.topology.atoms if a.residue.chain.index != 2]

print(f"  Rec: {len(rec_idx)} atoms, Lig: {len(lig_idx)} atoms", flush=True)

# Build systems
full_top = traj.topology.to_openmm()
rec_top = traj.atom_slice(rec_idx).topology.to_openmm()
lig_top = traj.atom_slice(lig_idx).topology.to_openmm()

print("  Building systems...", flush=True)
full_sys = make_gb_sys(full_top)
rec_sys = make_gb_sys(rec_top)
lig_sys = make_gb_sys(lig_top)

full_sim = app.Simulation(full_top, full_sys, mm.LangevinIntegrator(300, 1, 0.002))
rec_sim = app.Simulation(rec_top, rec_sys, mm.LangevinIntegrator(300, 1, 0.002))
lig_sim = app.Simulation(lig_top, lig_sys, mm.LangevinIntegrator(300, 1, 0.002))

print(f"  Computing {n} frames...", flush=True)
E_c = np.zeros(n); E_r = np.zeros(n); E_l = np.zeros(n)
for j, fi in enumerate(idx):
    E_c[j] = gb_energy(full_sim, traj.xyz[fi] * unit.nanometer)
    E_r[j] = gb_energy(rec_sim, traj.xyz[fi, rec_idx] * unit.nanometer)
    E_l[j] = gb_energy(lig_sim, traj.xyz[fi, lig_idx] * unit.nanometer)
    if (j+1) % max(1, n//5) == 0:
        print(f"  {j+1}/{n} [{time.perf_counter()-t0:.0f}s]", flush=True)

dE = E_c - E_r - E_l
print(f"\n  Interaction: {dE.mean():+.1f}+-{dE.std():.1f} kcal/mol", flush=True)

np.savez(D/"per_residue.npz", E_c=E_c, E_r=E_r, E_l=E_l, dE=dE)
json.dump({"n":n,"dE_mean":float(dE.mean()),"dE_std":float(dE.std()),
    "E_c":float(E_c.mean()),"E_r":float(E_r.mean()),"E_l":float(E_l.mean())},
    open(D/"decomp_summary.json","w"), indent=2)
print(f"  DONE [{(time.perf_counter()-t0)/60:.1f}min]", flush=True)
