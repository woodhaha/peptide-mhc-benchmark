"""Quick MD validation: implicit solvent GBSA, 1ns, ~15 min CPU.

ponytail: explicit solvent was taking hours to build. GBSA is the standard
first-pass for salt bridge stability — no water to equilibrate, 30× faster.
"""
import numpy as np, json, time
from pathlib import Path
import pdbfixer
import openmm as mm
import openmm.app as app
import openmm.unit as unit

DOCK_DIR = Path(__file__).parent
K_SIDECHAIN = 7.6          # Å
SALT_BRIDGE_MAX = 4.0      # Å
MD_LENGTH_NS = 1           # ns
TIMESTEP = 4.0             # fs (GBSA allows longer timestep)
REPORT_INTERVAL = 5000     # steps (20 ps)
N_SAMPLES = int(MD_LENGTH_NS * 1e6 / TIMESTEP / REPORT_INTERVAL)

t0 = time.perf_counter()
print("=" * 60)
print(f"  K-E63 SALT BRIDGE MD ({MD_LENGTH_NS}ns GBSA implicit solvent)")
print(f"  Template: 1DUZ (LLFGYPVYV) — CA-based K reachability")
print(f"  Timestep: {TIMESTEP}fs | Samples: {N_SAMPLES} | ~15 min CPU")
print("=" * 60)

# ── 1. Build system (pdbfixer, no solvent) ──
print("\n[1/4] Building system...")
fixer = pdbfixer.PDBFixer(str(DOCK_DIR / "1DUZ.pdb"))
fixer.findMissingResidues()
fixer.findMissingAtoms()
fixer.addMissingAtoms()
fixer.addMissingHydrogens(7.0)
print(f"  Missing: {len(fixer.missingResidues)} residues, {len(fixer.missingAtoms)} atoms fixed")

# Keep only chains A, B, C
modeller = app.Modeller(fixer.topology, fixer.positions)
to_delete = []
for chain in modeller.topology.chains():
    if chain.id not in ("A", "B", "C"):
        to_delete.extend(list(chain.residues()))
for chain in modeller.topology.chains():
    if chain.id in ("A", "B", "C"):
        to_delete.extend([r for r in chain.residues() if r.name == "HOH"])
modeller.delete(to_delete)

# GBSA implicit solvent — no water box needed
ff = app.ForceField("amber14-all.xml", "amber14/tip3pfb.xml")
system = ff.createSystem(
    modeller.topology,
    nonbondedMethod=app.NoCutoff,          # GBSA: no cutoff needed
    constraints=app.HBonds,
    implicitSolvent=app.GBn2,             # GBSA model
    soluteDielectric=1.0,
    solventDielectric=78.5,
)
n_atoms = sum(1 for _ in modeller.topology.atoms())
del system  # will recreate after indexing

# ── 2. Index atoms ──
print(f"[2/4] System: {n_atoms} atoms (no solvent). Indexing...")

p2_ca_idx = e63_oe1_idx = e63_oe2_idx = e63_cd_idx = None
for i, atom in enumerate(modeller.topology.atoms()):
    res = atom.residue
    if res.chain.id == "C" and str(res.id) == "2" and atom.name == "CA":
        p2_ca_idx = i
    if res.chain.id == "A" and str(res.id) == "63":
        if atom.name == "OE1": e63_oe1_idx = i
        elif atom.name == "OE2": e63_oe2_idx = i
        elif atom.name == "CD": e63_cd_idx = i

if p2_ca_idx is None or e63_cd_idx is None:
    print("FATAL: Key atoms not found")
    exit(1)
print(f"  P2:CA={p2_ca_idx}  E63:OE1={e63_oe1_idx} OE2={e63_oe2_idx}")

# Recreate system for real (GBSA)
system = ff.createSystem(
    modeller.topology,
    nonbondedMethod=app.NoCutoff,
    constraints=app.HBonds,
    implicitSolvent=app.GBn2,
    soluteDielectric=1.0,
    solventDielectric=78.5,
)

# ── 3. Minimize + equilibrate ──
print("[3/4] Minimizing + equilibrating...")
integrator = mm.LangevinMiddleIntegrator(
    300 * unit.kelvin, 1.0 / unit.picosecond, TIMESTEP * unit.femtosecond
)
sim = app.Simulation(modeller.topology, system, integrator,
                     mm.Platform.getPlatformByName('CPU'),
                     {'Threads': '15'})
sim.context.setPositions(modeller.positions)

# Minimize
e0 = sim.context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
sim.minimizeEnergy(maxIterations=5000)
e1 = sim.context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
print(f"  Energy: {e0:.0f} → {e1:.0f} kJ/mol")

# Short equilibration: 100ps NVT
print("  Equilibrating (100ps NVT)...")
sim.step(25000)   # 100ps at 4fs

# ── 4. Production ──
print(f"[4/4] Production MD ({MD_LENGTH_NS}ns, {TIMESTEP}fs)...")
sim.reporters.append(app.StateDataReporter(
    str(DOCK_DIR / "md_log.csv"), REPORT_INTERVAL,
    step=True, time=True, potentialEnergy=True, temperature=True,
))

dtype = np.dtype([
    ("step", "i4"), ("time_ns", "f4"),
    ("d_ca_oe1", "f4"), ("d_ca_oe2", "f4"), ("d_ca_cd", "f4"),
    ("min_oe", "f4"),
    ("k_ca_reach", "?"), ("est_nz_dist", "f4"), ("salt_bridge", "?"),
    ("T", "f4"), ("E_pot", "f4"),
])
dist_data = np.zeros(N_SAMPLES, dtype=dtype)

N_STEPS = int(MD_LENGTH_NS * 1e6 / TIMESTEP)
sample_idx = 0

for step in range(0, N_STEPS, REPORT_INTERVAL):
    sim.step(REPORT_INTERVAL)
    state = sim.context.getState(getPositions=True, getEnergy=True)
    pos_a = state.getPositions(asNumpy=True).value_in_unit(unit.nanometer) * 10

    p2_ca = pos_a[p2_ca_idx]
    d_oe1 = float(np.linalg.norm(p2_ca - pos_a[e63_oe1_idx]))
    d_oe2 = float(np.linalg.norm(p2_ca - pos_a[e63_oe2_idx]))
    d_cd  = float(np.linalg.norm(p2_ca - pos_a[e63_cd_idx]))
    min_oe = min(d_oe1, d_oe2)

    k_ca_reach = min_oe < K_SIDECHAIN
    est_nz = max(0.0, min_oe - K_SIDECHAIN + 1.5)
    salt_bridge = est_nz < SALT_BRIDGE_MAX

    T = (2 * state.getKineticEnergy() / (3 * n_atoms * unit.MOLAR_GAS_CONSTANT_R)).value_in_unit(unit.kelvin)

    dist_data[sample_idx] = (
        step + REPORT_INTERVAL,
        (step + REPORT_INTERVAL) * TIMESTEP / 1000,
        d_oe1, d_oe2, d_cd, min_oe, k_ca_reach, est_nz, salt_bridge,
        T, state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole),
    )
    sample_idx += 1

    if (step // REPORT_INTERVAL) % 25 == 0:
        pct = (step + REPORT_INTERVAL) / N_STEPS * 100
        r = dist_data[sample_idx - 1]
        elapsed = time.perf_counter() - t0
        eta = elapsed / max(pct, 0.1) * (100 - pct) / 60 if pct > 0 else 0
        print(f"  {pct:.0f}%  t={r['time_ns']:.1f}ns  "
              f"P2CA-E63={r['min_oe']:.1f}Å  K_reach={'✓' if r['k_ca_reach'] else '✗'}  "
              f"salt_bridge={'✓' if r['salt_bridge'] else '✗'}  "
              f"T={r['T']:.0f}K  [{elapsed:.0f}s ETA{eta:.0f}min]")

# ── 5. Results ──
valid = dist_data[:sample_idx]
ca_dists = valid["min_oe"]
sb = valid["salt_bridge"]
ca_reach = valid["k_ca_reach"]

print(f"\n{'='*60}")
print(f"  MD RESULTS ({valid['time_ns'][-1]:.1f} ns GBSA)")
print(f"  {'='*56}")
print(f"\n  P2:CA → E63:OE:")
print(f"    Mean: {ca_dists.mean():.1f} ± {ca_dists.std():.1f} Å")
print(f"    Min:  {ca_dists.min():.1f} Å  Max: {ca_dists.max():.1f} Å")
print(f"    K CA reach (<{K_SIDECHAIN}Å): {ca_reach.mean()*100:.1f}%")
print(f"    Salt bridge (<{SALT_BRIDGE_MAX}Å est_NZ): {sb.mean()*100:.1f}% of frames")

if sb.mean() >= 0.5:
    print(f"\n  VERDICT: ✅ STABLE ({sb.mean()*100:.0f}%)")
elif sb.mean() >= 0.2:
    print(f"\n  VERDICT: ⚠️ TRANSIENT ({sb.mean()*100:.0f}%)")
elif ca_reach.mean() >= 0.5:
    print(f"\n  VERDICT: ⚠️ Reachable but rarely forms ({sb.mean()*100:.0f}%)")
else:
    print(f"\n  VERDICT: ❌ Cannot reliably reach ({ca_reach.mean()*100:.0f}%)")

# Save
np.savez(DOCK_DIR / "md_distances.npz",
    step=valid["step"], time_ns=valid["time_ns"],
    d_ca_oe1=valid["d_ca_oe1"], d_ca_oe2=valid["d_ca_oe2"],
    d_ca_cd=valid["d_ca_cd"], min_oe=valid["min_oe"],
    k_ca_reach=valid["k_ca_reach"], est_nz_dist=valid["est_nz_dist"],
    salt_bridge=valid["salt_bridge"],
    T=valid["T"], E_pot=valid["E_pot"],
    method="GBSA implicit solvent",
)

json.dump({
    "method": "GBSA implicit solvent (GBn2)",
    "duration_ns": float(valid["time_ns"][-1]),
    "ca_mean": float(ca_dists.mean()), "ca_std": float(ca_dists.std()),
    "ca_min": float(ca_dists.min()), "ca_max": float(ca_dists.max()),
    "k_reachable_pct": float(ca_reach.mean() * 100),
    "salt_bridge_pct": float(sb.mean() * 100),
    "verdict": ("STABLE" if sb.mean()>=0.5 else "TRANSIENT" if sb.mean()>=0.2
                else "REACHABLE" if ca_reach.mean()>=0.5 else "UNSTABLE"),
}, open(DOCK_DIR / "salt_bridge_summary.json", "w"), indent=2)

print(f"\n  Total time: {time.perf_counter()-t0:.0f}s")
print(f"  Saved: md_distances.npz, salt_bridge_summary.json")
