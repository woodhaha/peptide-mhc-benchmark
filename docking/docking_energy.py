#!/usr/bin/env python3
"""
Unified peptide-MHC energy analysis.
Merges: binding_energy.py + energy_optimized.py + mmgbsa_neo_vs_wt.py

Single-pass OpenMM: loads each complex once, computes:
  1. Full complex potential energy (minimized)
  2. Receptor-only energy (peptide deleted)
  3. Peptide-only strain energy
  4. ΔG_binding = E(complex) − E(receptor) − E(peptide)
  5. K-E63 distance from minimized structure
"""
import numpy as np
import json
from pathlib import Path
import time
import openmm as mm
import openmm.app as app
import openmm.unit as unit

DOCK_DIR = Path(__file__).parent
HDOCK_DIR = DOCK_DIR / "hdock_models"

FORCEFIELD = app.ForceField("amber14-all.xml", "amber14/tip3pfb.xml")
INTEGRATOR = mm.LangevinMiddleIntegrator(300 * unit.kelvin, 1.0 / unit.picosecond, 2.0 * unit.femtosecond)
K_SIDECHAIN = 7.6


def clear_non_mhc(modeller):
    """Remove everything except chains A, B, C and delete waters."""
    to_delete = []
    for chain in modeller.topology.chains():
        if chain.id not in ("A", "B", "C"):
            to_delete.extend(list(chain.residues()))
        elif chain.residues:
            to_delete.extend([r for r in chain.residues() if r.name == "HOH"])
    modeller.delete(to_delete)


def make_system(topology):
    return FORCEFIELD.createSystem(
        topology,
        nonbondedMethod=app.NoCutoff,
        constraints=app.HBonds,
        hydrogenMass=1.5 * unit.amu,
    )


def minimize_energy(pdb_path):
    """Load PDB, add H, minimize, return (energy_kj, modeller, simulation)."""
    pdb = app.PDBFile(str(pdb_path))
    modeller = app.Modeller(pdb.topology, pdb.positions)
    clear_non_mhc(modeller)
    modeller.addHydrogens(forcefield=FORCEFIELD)

    system = make_system(modeller.topology)
    sim = app.Simulation(modeller.topology, system, INTEGRATOR)
    sim.context.setPositions(modeller.positions)
    sim.minimizeEnergy(maxIterations=2000)
    energy = sim.context.getState(getEnergy=True).getPotentialEnergy()
    return energy.value_in_unit(unit.kilojoules_per_mole), modeller, sim


def split_energy(pdb_path):
    """Compute ΔG = E(complex) − E(receptor) − E(peptide)."""
    # 1. Full complex
    e_complex, _, _ = minimize_energy(pdb_path)

    # 2. Receptor only (delete peptide chain C)
    pdb = app.PDBFile(str(pdb_path))
    modeller = app.Modeller(pdb.topology, pdb.positions)
    clear_non_mhc(modeller)
    # delete peptide
    pep_residues = []
    for chain in modeller.topology.chains():
        if chain.id == "C":
            pep_residues.extend(list(chain.residues()))
    modeller.delete(pep_residues)
    modeller.addHydrogens(forcefield=FORCEFIELD)
    system = make_system(modeller.topology)
    sim = app.Simulation(modeller.topology, system, INTEGRATOR)
    sim.context.setPositions(modeller.positions)
    sim.minimizeEnergy(maxIterations=2000)
    e_receptor = sim.context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)

    # 3. Peptide only (delete everything except chain C)
    pdb = app.PDBFile(str(pdb_path))
    modeller = app.Modeller(pdb.topology, pdb.positions)
    non_pep = []
    for chain in modeller.topology.chains():
        if chain.id != "C":
            non_pep.extend(list(chain.residues()))
        elif chain.id == "C":
            non_pep.extend([r for r in chain.residues() if r.name == "HOH"])
    modeller.delete(non_pep)

    if len(list(modeller.topology.atoms())) == 0:
        return None

    modeller.addHydrogens(forcefield=FORCEFIELD)
    system = make_system(modeller.topology)
    sim = app.Simulation(modeller.topology, system, INTEGRATOR)
    sim.context.setPositions(modeller.positions)
    sim.minimizeEnergy(maxIterations=2000)
    e_peptide = sim.context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)

    delta_g = e_complex - e_receptor - e_peptide
    return {
        "e_complex": round(e_complex, 1),
        "e_receptor": round(e_receptor, 1),
        "e_peptide": round(e_peptide, 1),
        "delta_g_kj": round(delta_g, 1),
        "delta_g_kcal": round(delta_g / 4.184, 1),
    }


def measure_k_e63(pdb_path):
    """Measure P2→E63 distances from minimized structure."""
    pdb = app.PDBFile(str(pdb_path))
    modeller = app.Modeller(pdb.topology, pdb.positions)
    clear_non_mhc(modeller)
    modeller.addHydrogens(forcefield=FORCEFIELD)
    system = make_system(modeller.topology)
    sim = app.Simulation(modeller.topology, system, INTEGRATOR)
    sim.context.setPositions(modeller.positions)
    sim.minimizeEnergy(maxIterations=2000)

    positions = sim.context.getState(getPositions=True).getPositions(asNumpy=True)
    pos_a = positions.value_in_unit(unit.nanometer) * 10  # → Angstrom

    p2_ca = e63_oe1 = e63_oe2 = e63_cd = None
    for i, atom in enumerate(modeller.topology.atoms()):
        res = atom.residue
        if res.chain.id == "C" and res.id == "2" and atom.name == "CA":
            p2_ca = np.array(pos_a[i])
        if res.chain.id == "A" and res.id == "63":
            if atom.name == "OE1":
                e63_oe1 = np.array(pos_a[i])
            elif atom.name == "OE2":
                e63_oe2 = np.array(pos_a[i])
            elif atom.name == "CD":
                e63_cd = np.array(pos_a[i])

    if p2_ca is None or e63_cd is None:
        return {}

    d1 = float(np.linalg.norm(p2_ca - e63_oe1)) if e63_oe1 is not None else None
    d2 = float(np.linalg.norm(p2_ca - e63_oe2)) if e63_oe2 is not None else None
    dc = float(np.linalg.norm(p2_ca - e63_cd))
    min_oe = min(d1, d2) if d1 and d2 else None
    k_reach = min_oe < K_SIDECHAIN if min_oe else False
    est_nz = max(0.0, min_oe - K_SIDECHAIN + 1.5) if min_oe else None

    return {
        "p2_ca_to_e63_oe1": round(d1, 1) if d1 else None,
        "p2_ca_to_e63_oe2": round(d2, 1) if d2 else None,
        "p2_ca_to_e63_cd":  round(dc, 1),
        "min_oe":            round(min_oe, 1) if min_oe else None,
        "k_reachable":       bool(k_reach),
        "est_nz_dist":       round(est_nz, 1) if est_nz else None,
        "salt_bridge_constraint": "2.5-4.0Å for N-O salt bridge",
    }


# ═══════════════════════════════════════════════════════════════════════════════
#  MAIN
# ═══════════════════════════════════════════════════════════════════════════════
print("=" * 70)
print("  PEPTIDE-MHC BINDING ENERGY + K-E63 GEOMETRY")
print("=" * 70)

t0 = time.perf_counter()

complexes = [
    ("KRAS_G12V_neo", HDOCK_DIR / "KRAS_G12V_neo_complex_top10.pdb"),
    ("KRAS_G12V_wt",  HDOCK_DIR / "KRAS_G12V_wt_complex_top10.pdb"),
    ("p53_R248W_neo", HDOCK_DIR / "p53_R248W_neo_complex_top10.pdb"),
    ("p53_R248W_wt",  HDOCK_DIR / "p53_R248W_wt_complex_top10.pdb"),
]

results = {}
for name, fpath in complexes:
    print(f"\n── {name} ──")
    if not fpath.exists():
        print(f"  ⚠ {fpath.name} not found, skipping")
        continue

    energy = split_energy(str(fpath))
    geom = measure_k_e63(str(fpath)) if energy else {}

    result = {**energy, **geom} if energy else geom
    results[name] = result

    if energy:
        print(f"  E_complex:  {energy['e_complex']:>10.1f} kJ/mol")
        print(f"  E_receptor: {energy['e_receptor']:>10.1f} kJ/mol")
        print(f"  E_peptide:  {energy['e_peptide']:>10.1f} kJ/mol")
        print(f"  ΔG_bind:    {energy['delta_g_kcal']:>10.1f} kcal/mol  ({energy['delta_g_kj']:.0f} kJ/mol)")
        sign = "favorable ✓" if energy['delta_g_kj'] < 0 else "unfavorable ✗"
        print(f"  Binding:    {sign}")
    if geom.get("min_oe"):
        print(f"  P2:CA→E63:  {geom['min_oe']:.1f} Å  K reachable: {'YES' if geom['k_reachable'] else 'NO'}  est NZ dist: {geom.get('est_nz_dist', '?')} Å")

# Neo vs WT comparison
print(f"\n{'='*70}")
print(f"  NEO vs WT COMPARISON")
print(f"  {'='*66}")

for pair_label, neo_key, wt_key in [
    ("KRAS G12V", "KRAS_G12V_neo", "KRAS_G12V_wt"),
    ("p53 R248W", "p53_R248W_neo", "p53_R248W_wt"),
]:
    if neo_key in results and wt_key in results:
        neo, wt = results[neo_key], results[wt_key]
        dg = neo.get("delta_g_kcal", 0) - wt.get("delta_g_kcal", 0)
        print(f"\n  {pair_label}:")
        print(f"    ΔΔG(neo−WT) = {dg:+.1f} kcal/mol  "
              f"({'neo binds STRONGER ✅' if dg < 0 else 'WT binds stronger' if dg > 0 else 'same'})")
        if neo.get("min_oe") and wt.get("min_oe"):
            print(f"    P2→E63:      {neo['min_oe']:.1f}Å (neo) vs {wt['min_oe']:.1f}Å (wt)")

# K-E63 salt bridge verdict
print(f"\n{'='*70}")
print(f"  K-E63 SALT BRIDGE ENERGETICS")
print(f"  {'='*66}")

if "KRAS_G12V_neo" in results:
    r = results["KRAS_G12V_neo"]
    if r.get("k_reachable"):
        if r["delta_g_kj"] < -20:
            print(f"  VERDICT: ✅ Salt bridge likely stable (ΔG={r['delta_g_kcal']} kcal/mol, reachable)")
        else:
            print(f"  VERDICT: ⚠️ Geometry ok but ΔG marginal ({r['delta_g_kcal']} kcal/mol)")
    else:
        print(f"  VERDICT: ❌ K cannot reach E63 in minimized structure")

# Save
with open(DOCK_DIR / "energy_binding.json", "w") as f:
    json.dump(results, f, indent=2)
print(f"\n  ✓ Saved: energy_binding.json")
print(f"  Total: {time.perf_counter() - t0:.0f}s")
