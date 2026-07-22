#!/usr/bin/env python3
"""K-E63 salt bridge: CPU eq (proven) + CUDA prod 100ns 4fs HMR."""
import numpy as np, json, time
from pathlib import Path
import pdbfixer
import openmm as mm, openmm.app as app, openmm.unit as unit

D = Path("/root/autodl-tmp/docking")
TSTEP, K_SB, SB_MAX, RINT = 4.0, 7.6, 4.0, 25000
N_SAMP = int(100 * 1e6 / TSTEP / RINT)
t0 = time.perf_counter()
print("K-E63 SALT BRIDGE -- CPU eq + CUDA prod 100ns", flush=True)

# ── 1. Build ──
print("[1] Build...", flush=True)
fix = pdbfixer.PDBFixer(str(D / "1DUZ.pdb"))
fix.findMissingResidues(); fix.findMissingAtoms()
fix.addMissingAtoms(); fix.addMissingHydrogens(7.0)
for r in fix.topology.residues():
    if r.name == "CYS": r.name = "CYM"
mod = app.Modeller(fix.topology, fix.positions)
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
print(f"  {n_h}H repartitioned, {sys.getNumParticles()} particles", flush=True)

# ── 2. Index ──
print("[2] Index...", flush=True)
p2 = e1 = e2 = eCD = None; bb = []
for i, a in enumerate(mod.topology.atoms()):
    r = a.residue
    if r.chain.id == "C" and str(r.id) == "2" and a.name == "CA": p2 = i
    if r.chain.id == "A" and str(r.id) == "63":
        if a.name == "OE1": e1 = i
        elif a.name == "OE2": e2 = i
        elif a.name == "CD": eCD = i
    if r.chain.id in ("A","B","C") and a.name in ("N","CA","C","O"): bb.append(i)
print(f"  P2:CA={p2} E63:OE1={e1} OE2={e2}", flush=True)
assert p2 and eCD

k0 = 20.0 * unit.kilocalories_per_mole / (unit.angstroms**2)
bf = mm.CustomExternalForce("k*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
bf.addGlobalParameter("k", k0)
bf.addPerParticleParameter("x0"); bf.addPerParticleParameter("y0"); bf.addPerParticleParameter("z0")
for i in bb: p = mod.positions[i]; bf.addParticle(i, [p.x, p.y, p.z])
sys.addForce(bf)

def chk(s, lbl):
    e = s.context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilocalories_per_mole)
    ok = not (np.isnan(e) or np.isinf(e))
    print(f"  {lbl}: PE={e:.0f} {'OK' if ok else 'NaN!'}", flush=True)

# ── 3. CPU eq (proven stable) ──
print("[3] CPU eq...", flush=True)
cpu = mm.Platform.getPlatformByName("CPU")
intg = mm.LangevinMiddleIntegrator(300.0*unit.kelvin, 1.0/unit.picosecond, 2.0*unit.femtosecond)
sim = app.Simulation(mod.topology, sys, intg, cpu, {"Threads": "15"})
sim.context.setPositions(mod.positions)
sim.minimizeEnergy(maxIterations=5000); chk(sim, "minimize")
for t in [50, 100, 150, 200, 250, 300]:
    sim.integrator.setTemperature(t*unit.kelvin); sim.step(10000); chk(sim, f"NVT {t}K")
sim.integrator.setTemperature(300.0*unit.kelvin); sim.step(50000); chk(sim, "NVT eq")
# faster annealing than original
for kv, st in [(20,10000),(10,10000),(5,20000),(2,20000),(1,20000),(0,10000)]:
    bf.setGlobalParameterDefaultValue(0, kv*unit.kilocalories_per_mole/(unit.angstroms**2))
    bf.updateParametersInContext(sim.context); sim.step(st); chk(sim, f"rst k={kv}")
bf.setGlobalParameterDefaultValue(0, 0.0*unit.kilocalories_per_mole/(unit.angstroms**2))
bf.updateParametersInContext(sim.context); sim.step(100000); chk(sim, "NPT free")
print(f"  Eq done [{(time.perf_counter()-t0)/60:.0f}min]", flush=True)

# ── 4. CUDA production ──
print(f"[4] CUDA prod {100}ns 4fs...", flush=True)
st = sim.context.getState(getPositions=True, getVelocities=True, getEnergy=True, getParameters=False)
plat = mm.Platform.getPlatformByName("CUDA")
pint = mm.LangevinMiddleIntegrator(300.0*unit.kelvin, 1.0/unit.picosecond, TSTEP*unit.femtosecond)
prod = app.Simulation(mod.topology, sys, pint, plat, {"CudaPrecision": "mixed"})
prod.context.setState(st)
prod.reporters.append(app.StateDataReporter(str(D/"md_log.csv"), RINT,
    step=True, time=True, potentialEnergy=True, temperature=True, density=True))

PROD = int(100 * 1_000_000 / TSTEP)
dt = np.dtype([("step","i4"),("t_ns","f4"),("d1","f4"),("d2","f4"),("dCD","f4"),
    ("minOE","f4"),("k_reach","?"),("nz","f4"),("sb","?"),("T","f4"),("Ep","f4")])
buf = np.zeros(N_SAMP, dtype=dt); ix = 0

for step in range(0, PROD, RINT):
    prod.step(RINT)
    s = prod.context.getState(getPositions=True, getEnergy=True)
    posA = s.getPositions(asNumpy=True).value_in_unit(unit.nanometer)*10
    # T from KE (LangevinMiddle has no computeSystemTemperature)
    ke = s.getKineticEnergy().value_in_unit(unit.kilojoules_per_mole)
    ndof = sys.getNumParticles() * 3 - 3
    Tk = 2.0 * ke / (ndof * 0.008314462618)
    pv = posA[p2]; d1 = float(np.linalg.norm(pv-posA[e1])); d2 = float(np.linalg.norm(pv-posA[e2]))
    dCD = float(np.linalg.norm(pv-posA[eCD]))
    mo = min(d1,d2); kr = mo < K_SB; nz = max(0.0, mo-K_SB+1.5); sb = nz < SB_MAX
    buf[ix] = (step+RINT, (step+RINT)*TSTEP/1000, d1, d2, dCD, mo, kr, nz, sb, Tk,
               s.getPotentialEnergy().value_in_unit(unit.kilocalories_per_mole))
    ix += 1
    if ix % 25 == 0:
        r = buf[ix-1]; el = time.perf_counter()-t0
        print(f"  {(step+RINT)/PROD*100:.0f}% t={r['t_ns']:.1f}ns CA-E63={r['minOE']:.1f}A "
              f"K={'Y' if r['k_reach'] else 'N'} NZ={r['nz']:.1f}A T={r['T']:.0f}K [{el/60:.0f}min]", flush=True)

# ── 5. Results ──
v = buf[:ix]; cd = v["minOE"]
print(f"\n  {v['t_ns'][-1]:.1f}ns: CA-E63 {cd.mean():.1f}+/-{cd.std():.1f}A min={cd.min():.1f}A")
print(f"  Salt bridge: {v['sb'].mean()*100:.1f}%")
stable = v["T"].max() < 400 and v["Ep"].std()/abs(v["Ep"].mean()) < 0.1
print(f"  {'STABLE' if stable else 'UNSTABLE'}")
np.savez(D/"md_dist.npz", step=v["step"], t_ns=v["t_ns"], d1=v["d1"], mo=v["minOE"],
         sb=v["sb"], T=v["T"], Ep=v["Ep"])
json.dump({"ns":float(v["t_ns"][-1]),"sb_pct":float(v["sb"].mean()*100),
    "ca_mean":float(cd.mean()),"ca_min":float(cd.min()),"stable":bool(stable)},
    open(D/"sb_summary.json","w"), indent=2)
print(f"  Total: {(time.perf_counter()-t0)/60:.0f}min")
