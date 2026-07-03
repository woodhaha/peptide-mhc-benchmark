#!/usr/bin/env python3
"""
MD simulation: K-E63 salt bridge feasibility for KRAS G12V neoepitope.
Uses 1DUZ template (LLFGYPVYV) — P2:CA position is identical for L and K.
Tracks P2:CA→E63 distance over trajectory; infers K reachability geometrically.

Key insight: LEU and LYS share the same backbone (CA, N, C, O).
If P2:CA stays within 7.6Å of E63 over MD, then a K sidechain CAN reach.

Merges: md_kras_k_at_p2.py + md_optimized.py (fixes crash: T→630K, frozen ρ)

Output: md_distances.npz, md_log.csv, salt_bridge_summary.json
"""
import numpy as np
import json
import time
from pathlib import Path
import pdbfixer
import openmm as mm
import openmm.app as app
import openmm.unit as unit

DOCK_DIR = Path(__file__).parent
K_SIDECHAIN = 7.6          # Å, CA→NZ for lysine
SALT_BRIDGE_MAX = 4.0      # Å, ideal N-O salt bridge
MD_LENGTH_NS = 1           # ns production (practical for CPU: ~8h)
TIMESTEP = 2.0             # fs
N_THREADS = 15             # CPU threads for OpenMM
GPU_OPENCL = False          # CPU-only (OpenCL NaN on RTX 2000 Ada)
REPORT_INTERVAL = 5000     # steps (10 ps)
N_SAMPLES = int(MD_LENGTH_NS * 1e6 / TIMESTEP / REPORT_INTERVAL)

print("=" * 70)
print(f"  KRAS G12V — K-E63 SALT BRIDGE MD ({MD_LENGTH_NS}ns)")
print(f"  Platform: {'OpenCL GPU (NVIDIA RTX 2000 Ada)' if GPU_OPENCL else f'CPU x {N_THREADS} threads'} | {TIMESTEP}fs | ~30K atoms")
print(f"  Template: 1DUZ (LLFGYPVYV) — CA-based K reachability")
print("=" * 70)

t0 = time.perf_counter()

# ═══════════════════════════════════════════════════════════════════════════════
#  1. BUILD SYSTEM (pdbfixer → OpenMM)
# ═══════════════════════════════════════════════════════════════════════════════
print("\n[1/5] Building system from 1DUZ (pdbfixer → OpenMM)...")

# 1a. Fix PDB with pdbfixer
fixer = pdbfixer.PDBFixer(str(DOCK_DIR / "1DUZ.pdb"))
fixer.findMissingResidues()
fixer.findMissingAtoms()
fixer.addMissingAtoms()
fixer.addMissingHydrogens(7.0)  # pH 7.0
print(f"  PDBFixer: {len(fixer.missingResidues)} missing residues, "
      f"{len(fixer.missingAtoms)} missing atoms fixed")

# 1b. Keep only chains A (heavy), B (B2M), C (peptide)
modeller = app.Modeller(fixer.topology, fixer.positions)
to_delete = []
for chain in modeller.topology.chains():
    if chain.id not in ("A", "B", "C"):
        to_delete.extend(list(chain.residues()))
for chain in modeller.topology.chains():
    if chain.id in ("A", "B", "C"):
        to_delete.extend([r for r in chain.residues() if r.name == "HOH"])
modeller.delete(to_delete)

# 1c. Forcefield + solvate
ff = app.ForceField("amber14-all.xml", "amber14/tip3pfb.xml")
modeller.addSolvent(
    ff, model="tip3p", padding=1.0 * unit.nanometer,
    ionicStrength=0.15 * unit.molar,
)

system = ff.createSystem(
    modeller.topology,
    nonbondedMethod=app.PME,
    nonbondedCutoff=1.0 * unit.nanometer,
    constraints=app.HBonds,
)
n_atoms = sum(1 for _ in modeller.topology.atoms())
print(f"  System: {n_atoms} atoms, {len(list(modeller.topology.chains()))} chains")

# ═══════════════════════════════════════════════════════════════════════════════
#  2. PRE-COMPUTE ATOM INDICES
# ═══════════════════════════════════════════════════════════════════════════════
print("[2/5] Indexing atoms...")

p2_ca_idx = None
e63_oe1_idx = e63_oe2_idx = e63_cd_idx = None
bb_indices = []

for i, atom in enumerate(modeller.topology.atoms()):
    res = atom.residue
    if res.chain.id == "C":
        if str(res.id) == "2" and atom.name == "CA":
            p2_ca_idx = i
    if res.chain.id == "A" and str(res.id) == "63":
        if atom.name == "OE1":
            e63_oe1_idx = i
        elif atom.name == "OE2":
            e63_oe2_idx = i
        elif atom.name == "CD":
            e63_cd_idx = i
    if res.chain.id in ("A", "B", "C") and atom.name in ("N", "CA", "C", "O"):
        bb_indices.append(i)

print(f"  P2:CA={p2_ca_idx}  E63:OE1={e63_oe1_idx} OE2={e63_oe2_idx} CD={e63_cd_idx}")
print(f"  Backbone atoms: {len(bb_indices)}")

if p2_ca_idx is None or e63_cd_idx is None:
    print("  ❌ FATAL: Could not find key atoms. Check PDB structure.")
    exit(1)

# Backbone restraints with gradual annealing (k=20→10→5→2→1→0)
RESTRAINT_K = 20.0  # start value (will be ramped down)
force = mm.CustomExternalForce("k*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
force.addGlobalParameter("k", RESTRAINT_K * unit.kilojoules_per_mole / unit.nanometer**2)
force.addPerParticleParameter("x0")
force.addPerParticleParameter("y0")
force.addPerParticleParameter("z0")
for i in bb_indices:
    pos = modeller.positions[i]
    force.addParticle(i, [pos.x, pos.y, pos.z])
restraint_idx = system.addForce(force)

# ═══════════════════════════════════════════════════════════════════════════════
#  3. MINIMIZE + EQUILIBRATE (gradual restraint annealing)
# ═══════════════════════════════════════════════════════════════════════════════
print("[3/5] Minimizing + equilibrating (gradual restraint release)...")

integrator = mm.LangevinMiddleIntegrator(
    300 * unit.kelvin, 1.0 / unit.picosecond, TIMESTEP * unit.femtosecond
)
if GPU_OPENCL:
    sim = app.Simulation(
        modeller.topology, system, integrator,
        mm.Platform.getPlatformByName('OpenCL'),
        {'OpenCLPlatformIndex': '0', 'OpenCLDeviceIndex': '0', 'Precision': 'single'},
    )
else:
    sim = app.Simulation(
        modeller.topology, system, integrator,
        mm.Platform.getPlatformByName('CPU'),
        {'Threads': str(N_THREADS)},
    )
sim.context.setPositions(modeller.positions)

# Energy minimization (with restraints at k=20)
e0 = sim.context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
sim.minimizeEnergy(maxIterations=5000)
e1 = sim.context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
print(f"  Minimized: {e0:.0f} → {e1:.0f} kJ/mol")

# NVT warmup (50 ps, k=20)
sim.step(25000)

# NPT equilibration with gradual restraint annealing
system.addForce(mm.MonteCarloBarostat(1 * unit.atmosphere, 300 * unit.kelvin, 25))
sim.context.reinitialize(preserveState=True)

def set_restraint_k(sim, system, k_val):
    """Gradually update restraint force constant."""
    for i in range(system.getNumForces()):
        f = system.getForce(i)
        if isinstance(f, mm.CustomExternalForce):
            f.setGlobalParameterDefaultValue(0, k_val * unit.kilojoules_per_mole / unit.nanometer**2)
            f.updateParametersInContext(sim.context)
            return

anneal_schedule = [
    (20.0, 25000,  "k=20"),    # 50 ps
    (10.0, 25000,  "k=10"),    # 50 ps
    (5.0,  50000,  "k=5"),     # 100 ps
    (2.0,  50000,  "k=2"),     # 100 ps
    (1.0,  50000,  "k=1"),     # 100 ps
    (0.0,  25000,  "k=0"),     # 50 ps — fully released for last bit
]
total_eq_ps = sum(steps for _, steps, _ in anneal_schedule) * TIMESTEP / 1000
print(f"  Annealing: {total_eq_ps:.0f}ps NPT ({sum(steps for _,steps,_ in anneal_schedule)} steps)")

for k_val, steps, label in anneal_schedule:
    set_restraint_k(sim, system, k_val)
    sim.step(steps)
    state = sim.context.getState(getEnergy=True)
    e = state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
    print(f"    {label}  E={e:.0f} kJ/mol")

print(f"  Equilibration done [{time.perf_counter() - t0:.0f}s]")

# ═══════════════════════════════════════════════════════════════════════════════
#  4. PRODUCTION MD
# ═══════════════════════════════════════════════════════════════════════════════
print(f"[4/5] Production MD ({MD_LENGTH_NS} ns, {TIMESTEP} fs timestep)...")
print(f"  Samples: {N_SAMPLES}  Interval: {REPORT_INTERVAL} steps (10 ps)")

sim.reporters.append(app.StateDataReporter(
    str(DOCK_DIR / "md_log.csv"), REPORT_INTERVAL,
    step=True, time=True, potentialEnergy=True, temperature=True, density=True,
))

dtype = np.dtype([
    ("step", "i4"), ("time_ns", "f4"),
    ("d_ca_oe1", "f4"), ("d_ca_oe2", "f4"), ("d_ca_cd", "f4"),
    ("min_oe", "f4"),
    ("k_ca_reach", "?"),      # CA within sidechain reach
    ("est_nz_dist", "f4"),    # estimated K:NZ → E63:OE
    ("salt_bridge", "?"),     # est_nz < SALT_BRIDGE_MAX
    ("T", "f4"), ("E_pot", "f4"),  # for crash detection
])
dist_data = np.zeros(N_SAMPLES, dtype=dtype)

N_STEPS = int(MD_LENGTH_NS * 1e6 / TIMESTEP)
sample_idx = 0

for step in range(0, N_STEPS, REPORT_INTERVAL):
    sim.step(REPORT_INTERVAL)

    state = sim.context.getState(getPositions=True, getEnergy=True)
    pos_a = state.getPositions(asNumpy=True).value_in_unit(unit.nanometer) * 10
    T = integrator.computeSystemTemperature(sim.context) if hasattr(integrator, 'computeSystemTemperature') else 0

    p2_ca = pos_a[p2_ca_idx]
    d_oe1 = float(np.linalg.norm(p2_ca - pos_a[e63_oe1_idx]))
    d_oe2 = float(np.linalg.norm(p2_ca - pos_a[e63_oe2_idx]))
    d_cd  = float(np.linalg.norm(p2_ca - pos_a[e63_cd_idx]))
    min_oe = min(d_oe1, d_oe2)

    k_ca_reach = min_oe < K_SIDECHAIN
    est_nz = max(0.0, min_oe - K_SIDECHAIN + 1.5)
    salt_bridge = est_nz < SALT_BRIDGE_MAX

    dist_data[sample_idx] = (
        step + REPORT_INTERVAL,
        (step + REPORT_INTERVAL) * TIMESTEP / 1000,
        d_oe1, d_oe2, d_cd,
        min_oe, k_ca_reach, est_nz, salt_bridge,
        T, state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole),
    )
    sample_idx += 1

    if (step // REPORT_INTERVAL) % 50 == 0:
        pct = (step + REPORT_INTERVAL) / N_STEPS * 100
        elapsed = time.perf_counter() - t0
        eta = elapsed / max(pct, 0.1) * (100 - pct) / 60 if pct > 0 else 0
        r = dist_data[sample_idx - 1]
        print(f"  {pct:.0f}%  t={r['time_ns']:.1f}ns  "
              f"P2CA-E63={r['min_oe']:.1f}Å  "
              f"K_reach={'✓' if r['k_ca_reach'] else '✗'}  "
              f"est_NZ={r['est_nz_dist']:.1f}Å  "
              f"T={r['T']:.0f}K  [{elapsed:.0f}s / ETA {eta:.0f}min]")

# ═══════════════════════════════════════════════════════════════════════════════
#  5. RESULTS
# ═══════════════════════════════════════════════════════════════════════════════
print(f"\n{'='*70}")
print(f"  MD RESULTS ({dist_data['time_ns'][sample_idx - 1]:.1f} ns)")
print(f"  {'='*66}")

valid = dist_data[:sample_idx]

# Quality checks
T_max = valid["T"].max()
T_min = valid["T"].min()
E_var = valid["E_pot"].std() / abs(valid["E_pot"].mean()) if valid["E_pot"].mean() != 0 else 99

print(f"\n  Quality:")
print(f"    T range: {T_min:.0f}–{T_max:.0f} K  "
      f"{'✅' if T_max < 400 else '⚠️ HIGH' if T_max < 600 else '❌ CRASHED'}")
print(f"    E_pot CV: {E_var:.1%}  "
      f"{'✅' if E_var < 0.1 else '⚠️ unstable' if E_var < 0.3 else '❌ CRASHED'}")

if T_max > 500 or E_var > 0.3:
    print(f"\n  ❌ Simulation unstable! Results are NOT reliable.")
    print(f"  Causes: insufficient equilibration, bad starting geometry, forcefield clash")
    print(f"  Fix: longer NPT, weaker restraints, check for steric clashes")

# CA-based distances
ca_dists = valid["min_oe"]
ca_reach = valid["k_ca_reach"]
est_nz = valid["est_nz_dist"]
sb = valid["salt_bridge"]

print(f"\n  P2:CA → E63:OE (template has L at P2, same backbone as K):")
print(f"    Mean: {ca_dists.mean():.1f} ± {ca_dists.std():.1f} Å")
print(f"    Min:  {ca_dists.min():.1f} Å  Max: {ca_dists.max():.1f} Å")
print(f"    K CA reach (< {K_SIDECHAIN}Å):  {ca_reach.mean() * 100:.1f}%")

print(f"\n  Estimated K:NZ → E63:OE (CA reach + sidechain geometry):")
print(f"    Mean: {est_nz.mean():.1f} ± {est_nz.std():.1f} Å")
print(f"    Min:  {est_nz.min():.1f} Å")
print(f"    Salt bridge (< {SALT_BRIDGE_MAX}Å): {sb.mean() * 100:.1f}% of frames")

if T_max < 400 and E_var < 0.1:
    if sb.mean() >= 0.5:
        print(f"\n  VERDICT: ✅ K-E63 salt bridge is STABLE ({sb.mean() * 100:.0f}%)")
    elif sb.mean() >= 0.2:
        print(f"\n  VERDICT: ⚠️ K-E63 salt bridge forms transiently ({sb.mean() * 100:.0f}%)")
    elif ca_reach.mean() >= 0.5:
        print(f"\n  VERDICT: ⚠️ K can reach E63 but salt bridge rarely forms")
        print(f"     ({ca_reach.mean() * 100:.0f}% reachable, {sb.mean() * 100:.0f}% salt bridge)")
    else:
        print(f"\n  VERDICT: ❌ K cannot reliably reach E63 ({ca_reach.mean() * 100:.0f}%)")
else:
    print(f"\n  VERDICT: ⚠️ SIMULATION UNSTABLE — results unreliable, need rerun")

# Save
np.savez(
    DOCK_DIR / "md_distances.npz",
    step=valid["step"], time_ns=valid["time_ns"],
    d_ca_oe1=valid["d_ca_oe1"], d_ca_oe2=valid["d_ca_oe2"],
    d_ca_cd=valid["d_ca_cd"], min_oe=valid["min_oe"],
    k_ca_reach=valid["k_ca_reach"], est_nz_dist=valid["est_nz_dist"],
    salt_bridge=valid["salt_bridge"],
    T=valid["T"], E_pot=valid["E_pot"],
    k_sidereach=K_SIDECHAIN, salt_bridge_cutoff=SALT_BRIDGE_MAX,
)

summary = {
    "duration_ns": float(valid["time_ns"][-1]),
    "template": "1DUZ (LLFGYPVYV)", "p2_residue": "LEU (same CA as LYS)",
    "k_sidereach_angstrom": K_SIDECHAIN, "salt_bridge_cutoff_angstrom": SALT_BRIDGE_MAX,
    "quality": {"T_max": float(T_max), "T_min": float(T_min), "E_CV": float(E_var), "stable": T_max < 400 and E_var < 0.1},
    "ca_based": {
        "mean": float(ca_dists.mean()), "std": float(ca_dists.std()),
        "min": float(ca_dists.min()), "max": float(ca_dists.max()),
        "k_ca_reachable_pct": float(ca_reach.mean() * 100),
    },
    "nz_estimated": {
        "mean": float(est_nz.mean()), "std": float(est_nz.std()),
        "min": float(est_nz.min()),
        "salt_bridge_pct": float(sb.mean() * 100),
    },
}
with open(DOCK_DIR / "salt_bridge_summary.json", "w") as f:
    json.dump(summary, f, indent=2)

print(f"\n  ✓ Saved: md_distances.npz, md_log.csv, salt_bridge_summary.json")
print(f"  Total time: {(time.perf_counter() - t0) / 60:.1f} min")
