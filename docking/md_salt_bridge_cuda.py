#!/usr/bin/env python3
"""
MD simulation: K-E63 salt bridge — CUDA/RTX 5090 + HMR + 4fs.
Template: 1DUZ (LLFGYPVYV) — CA-based K reachability.
100ns production for statistically meaningful salt bridge occupancy.

ponytail: v13 SMVT protocol adapted for docking — HMR + mixed precision + barostat from start.
"""
import numpy as np
import json
import time
from pathlib import Path
import pdbfixer
import openmm as mm
import openmm.app as app
import openmm.unit as unit

DATA_DIR = Path("/root/autodl-tmp/docking")
MD_LENGTH_NS = 100
TIMESTEP = 4.0             # fs (HMR-safe)
DT_EQ = 2.0                # fs equilibration
K_SIDECHAIN = 7.6          # Å, CA→NZ for lysine
SALT_BRIDGE_MAX = 4.0      # Å
REPORT_INTERVAL = 25000    # steps (100 ps at 4fs)
N_SAMPLES = int(MD_LENGTH_NS * 1e6 / TIMESTEP / REPORT_INTERVAL)

t0 = time.perf_counter()
print(f"  K-E63 SALT BRIDGE MD — {MD_LENGTH_NS}ns CUDA + HMR")
print(f"  Template: 1DUZ (LLFGYPVYV)")

# 1. Build system
print("[1/5] Building system...")
fixer = pdbfixer.PDBFixer(str(DATA_DIR / "1DUZ.pdb"))
fixer.findMissingResidues(); fixer.findMissingAtoms()
fixer.addMissingAtoms(); fixer.addMissingHydrogens(7.0)

modeller = app.Modeller(fixer.topology, fixer.positions)
to_delete = []
for chain in modeller.topology.chains():
    if chain.id not in ("A", "B", "C"):
        to_delete.extend(list(chain.residues()))
for chain in modeller.topology.chains():
    if chain.id in ("A", "B", "C"):
        to_delete.extend([r for r in chain.residues() if r.name == "HOH"])
modeller.delete(to_delete)

ff = app.ForceField("amber14-all.xml", "amber14/tip3pfb.xml")
modeller.addSolvent(ff, model="tip3p", padding=1.0 * unit.nanometer,
                    ionicStrength=0.15 * unit.molar)

system = ff.createSystem(modeller.topology, nonbondedMethod=app.PME,
                         nonbondedCutoff=1.0 * unit.nanometer,
                         constraints=app.HBonds, ignoreExternalBonds=True)
# barostat from start (ponytail: v13 fix)
system.addForce(mm.MonteCarloBarostat(1.0 * unit.atmosphere, 300.0 * unit.kelvin))

# HMR
atom_list = list(modeller.topology.atoms())
adj = [[] for _ in range(modeller.topology.getNumAtoms())]
for bond in modeller.topology.bonds():
    adj[bond[0].index].append(bond[1].index)
    adj[bond[1].index].append(bond[0].index)
n_hmr = 0
for atom in atom_list:
    if atom.element is not None and atom.element.symbol == 'H':
        h_idx = atom.index
        heavy = next((nbr for nbr in adj[h_idx]
                      if atom_list[nbr].element is not None
                      and atom_list[nbr].element.symbol != 'H'), None)
        if heavy is None: continue
        h_mass = system.getParticleMass(h_idx)
        heavy_mass = system.getParticleMass(heavy)
        new_h_mass = 3.0 * unit.amu
        system.setParticleMass(h_idx, new_h_mass)
        system.setParticleMass(heavy, heavy_mass - (new_h_mass - h_mass))
        n_hmr += 1
print(f"  {n_hmr} H repartitioned, {system.getNumParticles()} particles")

# 2. Index atoms
print("[2/5] Indexing...")
p2_ca_idx = e63_oe1_idx = e63_oe2_idx = e63_cd_idx = None
bb_indices = []
for i, atom in enumerate(modeller.topology.atoms()):
    res = atom.residue
    if res.chain.id == "C" and str(res.id) == "2" and atom.name == "CA":
        p2_ca_idx = i
    if res.chain.id == "A" and str(res.id) == "63":
        if atom.name == "OE1": e63_oe1_idx = i
        elif atom.name == "OE2": e63_oe2_idx = i
        elif atom.name == "CD": e63_cd_idx = i
    if res.chain.id in ("A", "B", "C") and atom.name in ("N", "CA", "C", "O"):
        bb_indices.append(i)
print(f"  P2:CA={p2_ca_idx}  E63:OE1={e63_oe1_idx} OE2={e63_oe2_idx}")
assert p2_ca_idx and e63_cd_idx, "Key atoms not found"

# Backbone restraints
k_rest = 20.0 * unit.kilocalories_per_mole / (unit.angstroms ** 2)
bb_force = mm.CustomExternalForce("k*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
bb_force.addGlobalParameter("k", k_rest)
bb_force.addPerParticleParameter("x0"); bb_force.addPerParticleParameter("y0")
bb_force.addPerParticleParameter("z0")
for i in bb_indices:
    p = modeller.positions[i]
    bb_force.addParticle(i, [p.x, p.y, p.z])
bb_idx = system.addForce(bb_force)

# 3. Minimize + equilibrate
print("[3/5] Minimize + equilibrate...")
plat = mm.Platform.getPlatformByName("CUDA")
integrator = mm.LangevinMiddleIntegrator(300.0 * unit.kelvin,
                                         1.0 / unit.picosecond, DT_EQ * unit.femtosecond)
sim = app.Simulation(modeller.topology, system, integrator, plat,
                     {"CudaPrecision": "mixed"})
sim.context.setPositions(modeller.positions)

def check(label):
    pe = sim.context.getState(getEnergy=True).getPotentialEnergy()
    v = pe.value_in_unit(unit.kilocalories_per_mole)
    ok = not (np.isnan(v) or np.isinf(v))
    print(f"  {label}: PE={v:.0f} {'OK' if ok else 'NaN!'}")
    return ok

sim.minimizeEnergy(maxIterations=5000)
check("minimize")

# NVT heating
for t in [50, 100, 150, 200, 250, 300]:
    sim.integrator.setTemperature(t * unit.kelvin)
    sim.step(10_000)
    check(f"NVT {t}K")

sim.integrator.setTemperature(300.0 * unit.kelvin)
sim.step(50_000)
check("NVT eq")

# NPT equil with restraint annealing
anneal = [(20, 25000), (10, 25000), (5, 50000), (2, 50000), (1, 50000), (0, 25000)]
for k_val, steps in anneal:
    bb_force.setGlobalParameterDefaultValue(0, k_val * unit.kilocalories_per_mole / (unit.angstroms ** 2))
    bb_force.updateParametersInContext(sim.context)
    sim.step(steps)
    check(f"restraint k={k_val}")

system.removeForce(bb_idx)
sim.context.reinitialize(preserveState=True)
sim.step(100_000)
check("NPT unrestrained")

# 4. Production
print(f"[4/5] Production {MD_LENGTH_NS}ns (4fs)...")
sim.integrator.setStepSize(TIMESTEP * unit.femtosecond)
sim.reporters.append(app.StateDataReporter(
    str(DATA_DIR / "md_log.csv"), REPORT_INTERVAL,
    step=True, time=True, potentialEnergy=True, temperature=True, density=True))

prod_steps = int(MD_LENGTH_NS * 1_000_000 / TIMESTEP)
save_freq = REPORT_INTERVAL

dtype = np.dtype([
    ("step", "i4"), ("time_ns", "f4"),
    ("d_ca_oe1", "f4"), ("d_ca_oe2", "f4"), ("d_ca_cd", "f4"),
    ("min_oe", "f4"), ("k_ca_reach", "?"),
    ("est_nz_dist", "f4"), ("salt_bridge", "?"),
    ("T", "f4"), ("E_pot", "f4"),
])
dist_data = np.zeros(N_SAMPLES, dtype=dtype)
sample_idx = 0

for step in range(0, prod_steps, save_freq):
    sim.step(save_freq)
    state = sim.context.getState(getPositions=True, getEnergy=True)
    pos_a = state.getPositions(asNumpy=True).value_in_unit(unit.nanometer) * 10
    T_k = integrator.computeSystemTemperature(sim.context)

    p2_ca = pos_a[p2_ca_idx]
    d_oe1 = float(np.linalg.norm(p2_ca - pos_a[e63_oe1_idx]))
    d_oe2 = float(np.linalg.norm(p2_ca - pos_a[e63_oe2_idx]))
    d_cd = float(np.linalg.norm(p2_ca - pos_a[e63_cd_idx]))
    min_oe = min(d_oe1, d_oe2)
    k_ca_reach = min_oe < K_SIDECHAIN
    est_nz = max(0.0, min_oe - K_SIDECHAIN + 1.5)
    salt_bridge = est_nz < SALT_BRIDGE_MAX

    dist_data[sample_idx] = (
        step + save_freq, (step + save_freq) * TIMESTEP / 1000,
        d_oe1, d_oe2, d_cd, min_oe, k_ca_reach, est_nz, salt_bridge,
        T_k, state.getPotentialEnergy().value_in_unit(unit.kilocalories_per_mole),
    )
    sample_idx += 1

    if sample_idx % 25 == 0:
        pct = (step + save_freq) / prod_steps * 100
        elapsed = time.perf_counter() - t0
        r = dist_data[sample_idx - 1]
        print(f"  {pct:.0f}%  t={r['time_ns']:.1f}ns  "
              f"P2CA-E63={r['min_oe']:.1f}Å  "
              f"K_reach={'Y' if r['k_ca_reach'] else 'N'}  "
              f"est_NZ={r['est_nz_dist']:.1f}Å  T={r['T']:.0f}K  [{elapsed/60:.0f}min]")

# 5. Results
valid = dist_data[:sample_idx]
ca_dists = valid["min_oe"]
sb = valid["salt_bridge"]
T_max = valid["T"].max()
E_var = valid["E_pot"].std() / abs(valid["E_pot"].mean())

print(f"\n{'='*60}")
print(f"  MD {valid['time_ns'][-1]:.1f}ns RESULTS")
print(f"  P2:CA→E63: {ca_dists.mean():.1f}±{ca_dists.std():.1f}Å  min={ca_dists.min():.1f}Å")
print(f"  Salt bridge: {sb.mean()*100:.1f}%  Tmax={T_max:.0f}K  E_CV={E_var:.1%}")
print(f"  {'STABLE' if T_max < 400 and E_var < 0.1 else 'UNSTABLE — need rerun'}")

np.savez(DATA_DIR / "md_distances.npz",
         step=valid["step"], time_ns=valid["time_ns"],
         d_ca_oe1=valid["d_ca_oe1"], d_ca_oe2=valid["d_ca_oe2"],
         min_oe=valid["min_oe"], k_ca_reach=valid["k_ca_reach"],
         salt_bridge=valid["salt_bridge"], T=valid["T"], E_pot=valid["E_pot"])
json.dump({"duration_ns": float(valid["time_ns"][-1]),
           "salt_bridge_pct": float(sb.mean() * 100),
           "ca_mean": float(ca_dists.mean()), "ca_min": float(ca_dists.min()),
           "stable": T_max < 400 and E_var < 0.1},
          open(DATA_DIR / "salt_bridge_summary.json", "w"), indent=2)
print(f"  Saved. Total: {(time.perf_counter()-t0)/60:.0f}min")
