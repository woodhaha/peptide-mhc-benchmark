#!/usr/bin/env python3
"""MM-GBSA for 1DUZ K-E63 salt bridge. OpenMM 8.5 compatible."""
import numpy as np, time, json
from pathlib import Path
import mdtraj as md
import openmm as mm, openmm.app as app, openmm.unit as unit

D = Path("/root/autodl-tmp/docking")
GAMMA = 0.00542
N_SKIP = 5
t0 = time.perf_counter()

def make_sys_gb(top):
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
    gb.setSolventDielectric(78.5)
    gb.setSoluteDielectric(1.0)
    sys.addForce(gb)
    return sys

def gb_energy(sim, pos):
    sim.context.setPositions(pos)
    state = sim.context.getState(getEnergy=True)
    return state.getPotentialEnergy().value_in_unit(unit.kilocalories_per_mole)

print("[1] Load DCD...", flush=True)
traj = md.load_dcd(str(D / "traj.dcd"), top=str(D / "final.pdb"))
solvent = [a.index for a in traj.topology.atoms if a.residue.is_water or a.residue.name in ("NA","CL","K","MG","CA")]
if solvent:
    keep = [i for i in range(traj.n_atoms) if i not in solvent]
    traj = traj.atom_slice(keep)
n_tot = traj.n_frames
n = max(1, n_tot // N_SKIP)
idx = np.linspace(0, n_tot-1, n, dtype=int)
print(f"  {n_tot} frames -> {n}", flush=True)

print("[2] Split rec (A,B) / lig (C)...", flush=True)
rec_idx = [a.index for a in traj.topology.atoms if a.residue.chain.index in (0,1)]
lig_idx = [a.index for a in traj.topology.atoms if a.residue.chain.index == 2]
t_rec = traj.atom_slice(rec_idx); t_lig = traj.atom_slice(lig_idx)
print(f"  Rec {t_rec.n_atoms}  Lig {t_lig.n_atoms}", flush=True)
top_c = traj.topology.to_openmm()
top_r = t_rec.topology.to_openmm()
top_l = t_lig.topology.to_openmm()

print("[3] Build GB systems...", flush=True)
print("  Complex...", flush=True); sys_c = make_sys_gb(top_c)
print("  Receptor...", flush=True); sys_r = make_sys_gb(top_r)
print("  Ligand...", flush=True); sys_l = make_sys_gb(top_l)

intg_c = mm.LangevinIntegrator(300*unit.kelvin, 1/unit.picoseconds, 0.002*unit.picoseconds)
intg_r = mm.LangevinIntegrator(300*unit.kelvin, 1/unit.picoseconds, 0.002*unit.picoseconds)
intg_l = mm.LangevinIntegrator(300*unit.kelvin, 1/unit.picoseconds, 0.002*unit.picoseconds)
sim_c = app.Simulation(top_c, sys_c, intg_c)
sim_r = app.Simulation(top_r, sys_r, intg_r)
sim_l = app.Simulation(top_l, sys_l, intg_l)
print("  Done", flush=True)

print(f"[4] Compute {n} frames...", flush=True)
dt = np.dtype([("frame","i4"),("E_c","f4"),("E_r","f4"),("E_l","f4"),
    ("SASA_c","f4"),("SASA_r","f4"),("SASA_l","f4")])
out = np.zeros(n, dtype=dt)

for j, fi in enumerate(idx):
    pos_c = traj.xyz[fi] * unit.nanometer
    Ec = gb_energy(sim_c, pos_c)
    pos_r = t_rec.xyz[fi] * unit.nanometer
    Er = gb_energy(sim_r, pos_r)
    pos_l = t_lig.xyz[fi] * unit.nanometer
    El = gb_energy(sim_l, pos_l)
    sasa_c = md.shrake_rupley(traj[fi], probe_radius=0.14, n_sphere_points=100)[0].sum()
    sasa_r = md.shrake_rupley(t_rec[fi], probe_radius=0.14, n_sphere_points=100)[0].sum()
    sasa_l = md.shrake_rupley(t_lig[fi], probe_radius=0.14, n_sphere_points=100)[0].sum()
    out[j] = (fi, Ec, Er, El, sasa_c, sasa_r, sasa_l)
    if (j+1) % max(1, n//10) == 0:
        elapsed = time.perf_counter() - t0
        print(f"  {j+1}/{n} [{elapsed:.0f}s]", flush=True)

Gc, Gr, Gl = out["E_c"], out["E_r"], out["E_l"]
SAc, SAr, SAl = out["SASA_c"], out["SASA_r"], out["SASA_l"]
dE = Gc - Gr - Gl
dSA = GAMMA * (SAc - SAr - SAl)
dG = dE + dSA

print()
sep = "=" * 60
print(sep)
print(f"  MM-GBSA (amber14/GBSAOBC/SA) single-trajectory, {len(idx)} frames")
for label, arr in [("G_complex",Gc),("G_rec",Gr),("G_lig",Gl),
    ("dE (vacuum+GB)",dE),("dG_SA",dSA),("dG_binding",dG)]:
    mn = float(np.mean(arr)); sd = float(np.std(arr)); lo = float(np.min(arr)); hi = float(np.max(arr))
    print(f"  {label:>20}:  {mn:9.1f}  {sd:9.1f}  {lo:9.1f}  {hi:9.1f}")
print(sep)
print("  No entropy term. Use for ranking, not absolute affinity.")

np.savez(D/"mmgbsa.npz", frame=out["frame"], E_c=Gc, E_r=Gr, E_l=Gl,
         SASA_c=SAc, SASA_r=SAr, SASA_l=SAl, dG=dG, dE=dE, dSA=dSA)
json.dump({"n_frames":int(len(idx)),
    "dG_mean":float(dG.mean()),"dG_std":float(dG.std()),
    "dE_mean":float(dE.mean()),"dSA_mean":float(dSA.mean()),
    "Gc_mean":float(Gc.mean()),"Gr_mean":float(Gr.mean()),"Gl_mean":float(Gl.mean())},
    open(D/"mmgbsa_summary.json","w"), indent=2)
tot = (time.perf_counter() - t0) / 60
print(f"  DONE [{tot:.1f}min]", flush=True)
