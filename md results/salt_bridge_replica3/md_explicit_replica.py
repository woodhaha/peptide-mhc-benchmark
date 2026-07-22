#!/usr/bin/env python3
"""K-E63 explicit solvent replica — PME 4fs HMR, 20ns, all-CUDA."""
import numpy as np, json, time, os, sys
from pathlib import Path
import pdbfixer
import openmm as mm, openmm.app as app, openmm.unit as unit

SEED = int(sys.argv[1]) if len(sys.argv) > 1 else 42
D = Path("/root/autodl-tmp/replica3")
D.mkdir(parents=True, exist_ok=True)
os.chdir(str(D))
TSTEP, RINT, NPROD = 4.0, 25000, int(100 * 1e6 / 4.0)  # 20ns
t0 = time.perf_counter()
print(f"REPLICA seed={SEED} — explicit PME 4fs 20ns", flush=True)

print("[1] Build...", flush=True)
fix = pdbfixer.PDBFixer(str(D / "1DUZ.pdb"))
fix.findMissingResidues(); fix.findMissingAtoms()
fix.addMissingAtoms(); fix.addMissingHydrogens(7.0)
for r in fix.topology.residues():
    if r.name == "CYS": r.name = "CYM"
mod = app.Modeller(fix.topology, fix.positions)
# Keep chains A/B/C only, remove water
rem = []
for c in mod.topology.chains():
    if c.id not in ("A","B","C"): rem.extend(c.residues())
for c in mod.topology.chains():
    if c.id in ("A","B","C"): rem.extend(x for x in c.residues() if x.name == "HOH")
mod.delete(rem)
ff = app.ForceField("amber14-all.xml", "amber14/tip3pfb.xml")
mod.addSolvent(ff, padding=1.0*unit.nanometer, ionicStrength=0.15*unit.molar)
sys = ff.createSystem(mod.topology, nonbondedMethod=app.PME,
    nonbondedCutoff=1.0*unit.nanometer, constraints=app.HBonds, ignoreExternalBonds=False)
sys.addForce(mm.MonteCarloBarostat(1.0*unit.atmosphere, 300.0*unit.kelvin))

# HMR
al = list(mod.topology.atoms()); ad = [[] for _ in range(mod.topology.getNumAtoms())]
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
print(f"  {n_h}H, {sys.getNumParticles()} particles", flush=True)

print("[2] Index...", flush=True)
p2 = e1 = e2 = None
for i, a in enumerate(mod.topology.atoms()):
    r = a.residue
    if r.chain.id == "C" and str(r.id) == "2" and a.name == "CA": p2 = i
    if r.chain.id == "A" and str(r.id) == "63":
        if a.name == "OE1": e1 = i
        elif a.name == "OE2": e2 = i
print(f"  P2:CA={p2} E63:OE1={e1} OE2={e2}", flush=True)
assert None not in (p2, e1, e2), "Index failed"

# ── EQ (all on CUDA, quick) ──
print("[3] CUDA EQ...", flush=True)
plat = mm.Platform.getPlatformByName("CUDA")
intg = mm.LangevinMiddleIntegrator(300*unit.kelvin, 1.0/unit.picosecond, 2.0*unit.femtosecond)
sim = app.Simulation(mod.topology, sys, intg, plat, {"CudaPrecision": "mixed"})
sim.context.setPositions(mod.positions)
sim.minimizeEnergy(maxIterations=5000)
pe = lambda: sim.context.getState(getEnergy=True).getPotentialEnergy()._value
print(f"  min: PE={pe():.0f}", flush=True)
sim.step(50000); print(f"  100ps NPT: PE={pe():.0f}", flush=True)

# ── Prod ──
print(f"[4] Prod 100ns (4fs HMR)...", flush=True)
st = sim.context.getState(getPositions=True, getVelocities=True, getEnergy=True)
pint = mm.LangevinMiddleIntegrator(300*unit.kelvin, 1.0/unit.picosecond, TSTEP*unit.femtosecond)
prod = app.Simulation(mod.topology, sys, pint, plat, {"CudaPrecision": "mixed"})
prod.context.setState(st)
prod.context.setVelocitiesToTemperature(300*unit.kelvin)
prod.reporters.append(app.DCDReporter(str(D/f"traj_s{SEED}.dcd"), RINT))

dt = np.dtype([("step","i4"),("t_ns","f4"),("minOE","f4"),("sb","?"),("T","f4")])
buf = np.zeros(NPROD//RINT, dtype=dt); ix = 0
for step in range(0, NPROD, RINT):
    prod.step(RINT)
    s = prod.context.getState(getPositions=True, getEnergy=True)
    posA = s.getPositions(asNumpy=True)._value*10
    ke = s.getKineticEnergy()._value; ndof = sys.getNumParticles()*3-3; Tk = 2*ke/(ndof*0.008314462618)
    pv = posA[p2]
    mo = min(float(np.linalg.norm(pv-posA[e1])), float(np.linalg.norm(pv-posA[e2])))
    sb = mo < 4.0
    buf[ix] = (step+RINT, (step+RINT)*TSTEP/1000, mo, sb, Tk); ix += 1
    if ix % 20 == 0:
        el = time.perf_counter()-t0
        print(f"  {step/NPROD*100:.0f}% t={buf[ix-1]['t_ns']:.1f}ns "
              f"CA-E63={mo:.1f}A SB={'Y' if sb else 'N'} [{el/60:.1f}min]", flush=True)

v = buf[:ix]; cd = v['minOE']
print(f"\n  {v['t_ns'][-1]:.1f}ns: mean={cd.mean():.1f}+/-{cd.std():.1f}A SB={v['sb'].mean()*100:.1f}%", flush=True)
np.savez(D/f"md_dist_s{SEED}.npz", t_ns=v['t_ns'], minOE=v['minOE'], sb=v['sb'])
json.dump({"seed":SEED,"ns":float(v['t_ns'][-1]),"sb_pct":float(v['sb'].mean()*100),
    "ca_mean":float(cd.mean()),"ca_std":float(cd.std())}, open(D/f"sb_summary_s{SEED}.json","w"), indent=2)
with open(D/f"final_s{SEED}.pdb","w") as fh:
    app.PDBFile.writeFile(prod.topology, prod.context.getState(getPositions=True).getPositions(), fh)
print(f"  Total: {(time.perf_counter()-t0)/60:.1f}min", flush=True)
