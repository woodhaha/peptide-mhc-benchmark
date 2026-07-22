#!/usr/bin/env python3
import numpy as np, json, time, os, sys
from pathlib import Path
import pdbfixer
import openmm as mm, openmm.app as app, openmm.unit as unit

SEED = int(sys.argv[1]) if len(sys.argv) > 1 else 42
D = Path("/root/autodl-tmp/replica2")
D.mkdir(parents=True, exist_ok=True)
os.chdir(str(D))
TSTEP, RINT, NPROD = 4.0, 25000, int(50 * 1e6 / 4.0)
t0 = time.perf_counter()
print(f"SALT BRIDGE REPLICA seed={SEED} -- GBSA 4fs HMR x 50ns", flush=True)

print("[1] Build GBSA system...", flush=True)
fix = pdbfixer.PDBFixer(str(D / "1DUZ.pdb"))
fix.findMissingResidues(); fix.findMissingAtoms()
fix.addMissingAtoms(); fix.addMissingHydrogens(7.0)
mod = app.Modeller(fix.topology, fix.positions)
rem = []
for c in mod.topology.chains():
    if c.id not in ("A","B","C"): rem.extend(c.residues())
for c in mod.topology.chains():
    if c.id in ("A","B","C"): rem.extend(x for x in c.residues() if x.name == "HOH")
mod.delete(rem)

ff = app.ForceField("amber14-all.xml")
sys = ff.createSystem(mod.topology, nonbondedMethod=app.NoCutoff,
    constraints=app.HBonds, ignoreExternalBonds=False)
gb = mm.GBSAOBCForce()
for a in mod.topology.atoms():
    el = a.element; c = el.atomic_number if el else 0
    gb.addParticle(c, 0.15 if c == 1 else 0.17, 0.5)
gb.setSolventDielectric(78.5); gb.setSoluteDielectric(1.0)
sys.addForce(gb)

al = list(mod.topology.atoms())
ad = [[] for _ in range(mod.topology.getNumAtoms())]
for b in mod.topology.bonds(): ad[b[0].index].append(b[1].index); ad[b[1].index].append(b[0].index)
n_h = 0
for a in al:
    if not a.element or a.element.symbol != "H": continue
    hv = next((n for n in ad[a.index] if al[n].element and al[n].element.symbol != "H"), None)
    if hv is None: continue
    hm = sys.getParticleMass(hv)
    sys.setParticleMass(a.index, 3.0*unit.amu)
    sys.setParticleMass(hv, hm - (3.0*unit.amu - sys.getParticleMass(a.index)))
    n_h += 1
print(f"  {n_h}H repartitioned, {sys.getNumParticles()} particles", flush=True)

print("[2] Index...", flush=True)
p2 = e1 = e2 = None
for i, a in enumerate(mod.topology.atoms()):
    r = a.residue
    if r.chain.id == "C" and str(r.id) == "2" and a.name == "CA": p2 = i
    if r.chain.id == "A" and str(r.id) == "63":
        if a.name == "OE1": e1 = i
        elif a.name == "OE2": e2 = i
print(f"  P2:CA={p2} E63:OE1={e1} OE2={e2}", flush=True)
assert None not in (p2, e1, e2)

try:
    platform = mm.Platform.getPlatformByName("CUDA")
    plat_props = {"CudaPrecision": "mixed"}
    print("  Using CUDA", flush=True)
except:
    platform = mm.Platform.getPlatformByName("CPU")
    plat_props = {"Threads": "15"}
    print("  Using CPU", flush=True)

print("[3] EQ (0.5fs, short)...", flush=True)
intg = mm.LangevinMiddleIntegrator(300*unit.kelvin, 1.0/unit.picosecond, 0.5*unit.femtosecond)
sim = app.Simulation(mod.topology, sys, intg, platform, plat_props)
sim.context.setPositions(mod.positions)
sim.minimizeEnergy(maxIterations=5000)
pe = lambda: sim.context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilocalories_per_mole)
print(f"  min: PE={pe():.0f}", flush=True)
for tmp in [50, 100, 150, 200, 250, 300]:
    sim.integrator.setTemperature(tmp*unit.kelvin); sim.step(5000)
    pv = pe()
    if np.isnan(pv) or np.isinf(pv): raise RuntimeError(f"NaN at {tmp}K")
k0 = 20.0*unit.kilocalories_per_mole/(unit.angstroms**2)
bf = mm.CustomExternalForce("k*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
bf.addGlobalParameter("k", k0)
bf.addPerParticleParameter("x0"); bf.addPerParticleParameter("y0"); bf.addPerParticleParameter("z0")
for i, a in enumerate(mod.topology.atoms()):
    if a.name in ("N","CA","C","O") and a.residue.chain.id in ("A","B","C"):
        p = mod.positions[i]
        bf.addParticle(i, [p[0]._value, p[1]._value, p[2]._value])
sys.addForce(bf)
for kv in [10, 5, 2, 1, 0]:
    bf.setGlobalParameterDefaultValue(0, kv*unit.kilocalories_per_mole/(unit.angstroms**2))
    bf.updateParametersInContext(sim.context); sim.step(10000)
    pv = pe()
    if np.isnan(pv) or np.isinf(pv): raise RuntimeError(f"NaN at k={kv}")
print(f"  Eq ok [{(time.perf_counter()-t0)/60:.1f}min]", flush=True)

print(f"[4] Prod 50ns...", flush=True)
st = sim.context.getState(getPositions=True, getVelocities=True, getEnergy=True)
pint = mm.LangevinMiddleIntegrator(300*unit.kelvin, 1.0/unit.picosecond, TSTEP*unit.femtosecond)
pint.setRandomSeed(SEED)
prod = app.Simulation(mod.topology, sys, pint, platform, plat_props)
prod.context.setState(st)
prod.context.setVelocitiesToTemperature(300*unit.kelvin)
prod.reporters.append(app.DCDReporter(str(D/f"traj_s{SEED}.dcd"), RINT))

dt = np.dtype([("step","i4"),("t_ns","f4"),("d1","f4"),("d2","f4"),("minOE","f4"),
               ("nz","f4"),("sb","?"),("T","f4"),("Ep","f4")])
buf = np.zeros(NPROD//RINT, dtype=dt); ix = 0
for step in range(0, NPROD, RINT):
    prod.step(RINT)
    s = prod.context.getState(getPositions=True, getEnergies=True)
    posA = s.getPositions(asNumpy=True).value_in_unit(unit.nanometer)*10
    ke = s.getKineticEnergy()._value; ndof = sys.getNumParticles()*3-3; Tk = 2*ke/(ndof*0.008314462618)
    pv = posA[p2]; mo = min(float(np.linalg.norm(pv-posA[e1])), float(np.linalg.norm(pv-posA[e2])))
    nz = max(0.0, mo-4.0+1.5); sb = nz < 4.0
    buf[ix] = (step+RINT, (step+RINT)*TSTEP/1000, 0.0, 0.0, mo, nz, sb, Tk,
               s.getPotentialEnergy()._value)
    ix += 1
    if ix % 20 == 0:
        el = time.perf_counter()-t0
        print(f"  {step/NPROD*100:.0f}% t={buf[ix-1][t_ns]:.1f}ns CA-E63={mo:.1f}A "
              f"SB={Y if sb else N} [{el/60:.1f}min]", flush=True)

v = buf[:ix]; cd = v["minOE"]
print(f"\n  {v[t_ns][-1]:.1f}ns: mean={cd.mean():.1f}+/-{cd.std():.1f}A SB={v[sb].mean()*100:.1f}%", flush=True)
np.savez(D/f"md_dist_s{SEED}.npz", t_ns=v["t_ns"], minOE=v["minOE"], sb=v["sb"])
json.dump({"seed":SEED,"ns":float(v["t_ns"][-1]),"sb_pct":float(v["sb"].mean()*100),
    "ca_mean":float(cd.mean()),"ca_std":float(cd.std())}, open(D/f"sb_summary_s{SEED}.json","w"), indent=2)
print(f"  Total: {(time.perf_counter()-t0)/60:.1f}min", flush=True)
