#!/usr/bin/env python3
"""
MD: K-E63 salt bridge — full CUDA (no CPU eq, SMVT is done).
100ns, HMR + 4fs + mixed precision.
"""
import pickle, numpy as np, json, time, sys
from pathlib import Path
import pdbfixer
import openmm as mm
import openmm.app as app
import openmm.unit as unit

DATA_DIR = Path("/root/autodl-tmp/docking")
MD_LENGTH_NS = 100
TIMESTEP = 4.0
DT_EQ = 2.0
K_SIDECHAIN = 7.6
SALT_BRIDGE_MAX = 4.0
REPORT_INTERVAL = 25000
N_SAMPLES = int(MD_LENGTH_NS * 1e6 / TIMESTEP / REPORT_INTERVAL)

ckpt_file = DATA_DIR / "cuda_checkpoint.pkl"
resume = ckpt_file.exists()

t0 = time.perf_counter()
print(f"K-E63 SALT BRIDGE — full CUDA {MD_LENGTH_NS}ns", flush=True)

plat = mm.Platform.getPlatformByName("CUDA")
print(f"Platform: CUDA ({plat.getSpeed()})", flush=True)

# ── 1. Build system ──
print("[1/5] Building system...", flush=True)
fixer = pdbfixer.PDBFixer(str(DATA_DIR / "1DUZ.pdb"))
fixer.findMissingResidues(); fixer.findMissingAtoms()
fixer.addMissingAtoms(); fixer.addMissingHydrogens(7.0)
for r in fixer.topology.residues():
    if r.name == "CYS": r.name = "CYM"

modeller = app.Modeller(fixer.topology, fixer.positions)
to_delete = []
for chain in modeller.topology.chains():
    if chain.id not in ("A", "B", "C"): to_delete.extend(list(chain.residues()))
for chain in modeller.topology.chains():
    if chain.id in ("A", "B", "C"):
        to_delete.extend([r for r in chain.residues() if r.name == "HOH"])
modeller.delete(to_delete)

ff = app.ForceField("amber14-all.xml", "amber14/tip3pfb.xml")
modeller.addSolvent(ff, model="tip3p", padding=1.0 * unit.nanometer, ionicStrength=0.15 * unit.molar)
system = ff.createSystem(modeller.topology, nonbondedMethod=app.PME,
                         nonbondedCutoff=1.0 * unit.nanometer,
                         constraints=app.HBonds, ignoreExternalBonds=False)
system.addForce(mm.MonteCarloBarostat(1.0 * unit.atmosphere, 300.0 * unit.kelvin))

# HMR
atom_list = list(modeller.topology.atoms())
adj = [[] for _ in range(modeller.topology.getNumAtoms())]
for bond in modeller.topology.bonds():
    adj[bond[0].index].append(bond[1].index); adj[bond[1].index].append(bond[0].index)
n_hmr = 0
for atom in atom_list:
    if atom.element and atom.element.symbol == 'H':
        h_idx = atom.index
        heavy = next((nbr for nbr in adj[h_idx] if atom_list[nbr].element and atom_list[nbr].element.symbol != 'H'), None)
        if heavy is None: continue
        hm = system.getParticleMass(h_idx); hv = system.getParticleMass(heavy)
        system.setParticleMass(h_idx, 3.0 * unit.amu)
        system.setParticleMass(heavy, hv - (3.0 * unit.amu - hm))
        n_hmr += 1
print(f"  {n_hmr} H repartitioned, {system.getNumParticles()} particles", flush=True)

# ── 2. Index ──
print("[2/5] Indexing...", flush=True)
p2_ca = e63_oe1 = e63_oe2 = e63_cd = None
bb_idx = []
for i, a in enumerate(modeller.topology.atoms()):
    r = a.residue
    if r.chain.id == "C" and str(r.id) == "2" and a.name == "CA": p2_ca = i
    if r.chain.id == "A" and str(r.id) == "63":
        if a.name == "OE1": e63_oe1 = i
        elif a.name == "OE2": e63_oe2 = i
        elif a.name == "CD": e63_cd = i
    if r.chain.id in ("A","B","C") and a.name in ("N","CA","C","O"): bb_idx.append(i)
assert p2_ca and e63_cd, "Key atoms not found"

# Backbone restraints
k_rest = 20.0 * unit.kilocalories_per_mole / (unit.angstroms ** 2)
bbf = mm.CustomExternalForce("k*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
bbf.addGlobalParameter("k", k_rest)
bbf.addPerParticleParameter("x0"); bbf.addPerParticleParameter("y0"); bbf.addPerParticleParameter("z0")
for i in bb_idx:
    p = modeller.positions[i]
    bbf.addParticle(i, [p.x, p.y, p.z])
bbf_idx = system.addForce(bbf)

def check(label):
    pe = sim.context.getState(getEnergy=True).getPotentialEnergy()
    v = pe.value_in_unit(unit.kilocalories_per_mole)
    ok = not (np.isnan(v) or np.isinf(v))
    print(f"  {label}: PE={v:.0f} {'OK' if ok else 'NaN!'}", flush=True)
    return ok

# ── 3. Equilibrate on CUDA (SMVT done, GPU free) ──
intg = mm.LangevinMiddleIntegrator(300.0 * unit.kelvin, 1.0 / unit.picosecond, DT_EQ * unit.femtosecond)
sim = app.Simulation(modeller.topology, system, intg, plat, {"CudaPrecision": "mixed"})
sim.context.setPositions(modeller.positions)

sim.minimizeEnergy(maxIterations=5000); check("minimize")
for t in [50, 100, 150, 200, 250, 300]:
    intg.setTemperature(t * unit.kelvin); sim.step(10_000); check(f"NVT {t}K")
intg.setTemperature(300.0 * unit.kelvin); sim.step(50_000); check("NVT eq")

for k_val, steps in [(20,25000),(10,25000),(5,50000),(2,50000),(1,50000),(0,25000)]:
    bbf.setGlobalParameterDefaultValue(0, k_val * unit.kilocalories_per_mole / (unit.angstroms ** 2))
    bbf.updateParametersInContext(sim.context); sim.step(steps); check(f"k={k_val}")

bbf.setGlobalParameterDefaultValue(0, 0.0 * unit.kilocalories_per_mole / (unit.angstroms ** 2))
bbf.updateParametersInContext(sim.context); sim.step(100_000); check("NPT unrestrained")
t_eq = (time.perf_counter() - t0) / 60
print(f"  CUDA eq done [{t_eq:.0f}min]", flush=True)

ckpt = {"state": sim.context.getState(getPositions=True, getVelocities=True, getEnergy=True),
        "topology": modeller.topology, "step": "eq_done"}
with open(str(ckpt_file), "wb") as f: pickle.dump(ckpt, f)

# ── 4. Production ──
print(f"[4/5] Production {MD_LENGTH_NS}ns (4fs CUDA)...", flush=True)
intg.setStepSize(TIMESTEP * unit.femtosecond)
sim.reporters.append(app.StateDataReporter(
    str(DATA_DIR / "md_log.csv"), REPORT_INTERVAL,
    step=True, time=True, potentialEnergy=True, temperature=True, density=True))

prod_steps = int(MD_LENGTH_NS * 1_000_000 / TIMESTEP)
sf = REPORT_INTERVAL

dtype = np.dtype([("step","i4"),("time_ns","f4"),("d_ca_oe1","f4"),("d_ca_oe2","f4"),
                  ("d_ca_cd","f4"),("min_oe","f4"),("k_ca_reach","?"),("est_nz_dist","f4"),
                  ("salt_bridge","?"),("T","f4"),("E_pot","f4")])
dist = np.zeros(N_SAMPLES, dtype=dtype)
si = 0

for step in range(0, prod_steps, sf):
    sim.step(sf)
    s = sim.context.getState(getPositions=True, getEnergy=True)
    pa = s.getPositions(asNumpy=True).value_in_unit(unit.nanometer) * 10
    Tk = intg.computeSystemTemperature(sim.context)
    c = pa[p2_ca]
    d1, d2 = float(np.linalg.norm(c - pa[e63_oe1])), float(np.linalg.norm(c - pa[e63_oe2]))
    mo = min(d1, d2); kr = mo < K_SIDECHAIN; enz = max(0.0, mo - K_SIDECHAIN + 1.5)
    dist[si] = (step+sf, (step+sf)*TIMESTEP/1000, d1, d2,
                float(np.linalg.norm(c - pa[e63_cd])), mo, kr, enz, enz < SALT_BRIDGE_MAX,
                Tk, s.getPotentialEnergy().value_in_unit(unit.kilocalories_per_mole))
    si += 1
    if si % 25 == 0:
        pct = (step+sf)/prod_steps*100; el = time.perf_counter()-t0
        r = dist[si-1]
        print(f"  {pct:.0f}% t={r['time_ns']:.1f}ns P2CA-E63={r['min_oe']:.1f}A "
              f"{'SB' if r['salt_bridge'] else 'NO'} T={r['T']:.0f}K [{el/60:.0f}min]", flush=True)

# ── 5. Results ──
v = dist[:si]
cd = v["min_oe"]; sb = v["salt_bridge"]
print(f"\n{'='*60}")
print(f"  MD {v['time_ns'][-1]:.1f}ns RESULTS")
print(f"  P2:CA->E63: {cd.mean():.1f}+/-{cd.std():.1f}A  min={cd.min():.1f}A")
print(f"  Salt bridge: {sb.mean()*100:.1f}%  {v['T'].max():.0f}K max")

np.savez(DATA_DIR / "md_distances.npz",
         step=v["step"], time_ns=v["time_ns"], d_ca_oe1=v["d_ca_oe1"],
         d_ca_oe2=v["d_ca_oe2"], min_oe=v["min_oe"], k_ca_reach=v["k_ca_reach"],
         salt_bridge=v["salt_bridge"], T=v["T"], E_pot=v["E_pot"])
json.dump({"duration_ns": float(v["time_ns"][-1]), "salt_bridge_pct": float(sb.mean()*100),
           "ca_mean": float(cd.mean()), "ca_min": float(cd.min())},
          open(DATA_DIR / "salt_bridge_summary.json","w"), indent=2)
print(f"  DONE. Total: {(time.perf_counter()-t0)/60:.0f}min", flush=True)
