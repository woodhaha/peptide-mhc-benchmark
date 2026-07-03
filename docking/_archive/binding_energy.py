#!/usr/bin/env python3
"""
Compute peptide-MHC interaction energy via OpenMM single-point.
Aligns neo/WT peptides to template position, computes E(complex) - E(rec) - E(pep).
"""
import numpy as np
import openmm as mm
import openmm.app as app
import openmm.unit as unit

DOCK_DIR = r'D:\Researching\Peptide epitope\docking'

def energy_of_pdb(pdb_path):
    """Load, add H, minimize, return potential energy (kJ/mol)."""
    pdb = app.PDBFile(pdb_path)
    forcefield = app.ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')

    modeller = app.Modeller(pdb.topology, pdb.positions)
    to_del = []
    for chain in modeller.topology.chains():
        to_del.extend([r for r in chain.residues() if r.name == 'HOH'])
    modeller.delete(to_del)
    modeller.addHydrogens(forcefield=forcefield)

    system = forcefield.createSystem(modeller.topology, nonbondedMethod=app.NoCutoff,
                                     constraints=app.HBonds, hydrogenMass=1.5*unit.amu)
    integrator = mm.LangevinMiddleIntegrator(300*unit.kelvin, 1.0/unit.picosecond, 2.0*unit.femtosecond)
    sim = app.Simulation(modeller.topology, system, integrator)
    sim.context.setPositions(modeller.positions)
    sim.minimizeEnergy(maxIterations=2000)
    energy = sim.context.getState(getEnergy=True).getPotentialEnergy()
    return energy.value_in_unit(unit.kilojoules_per_mole)

# The HDOCK complex models have receptor+peptide together.
# Compute interaction energy = E(complex) - E(receptor_alone) - E(peptide_alone)
# But separating them from the PDB requires chain deletion.

# Simpler: compute per-peptide energies and compare.
# The peptide alone energy difference reflects intrinsic stability.

print("=" * 60)
print("  PEPTIDE STRAIN ENERGY: NEO vs WT")
print("=" * 60)

peptides = {
    'KRAS_G12V_neo': 'peptide_KRAS_G12V_neo.pdb',
    'KRAS_G12V_wt':  'peptide_KRAS_G12V_wt.pdb',
    'p53_R248W_neo': 'peptide_p53_R248W_neo.pdb',
    'p53_R248W_wt':  'peptide_p53_R248W_wt.pdb',
}

results = {}
for name, fname in peptides.items():
    e = energy_of_pdb(fname)
    results[name] = e
    print(f"  {name}: E = {e:.1f} kJ/mol")

# Compare
print(f"\n{'='*60}")
print("  COMPARISON")
print(f"{'='*60}")
if 'KRAS_G12V_neo' in results:
    d = results['KRAS_G12V_neo'] - results['KRAS_G12V_wt']
    print(f"  KRAS: neo-WT = {d:+.1f} kJ/mol ({'neo more stable' if d<0 else 'WT more stable'})")
if 'p53_R248W_neo' in results:
    d = results['p53_R248W_neo'] - results['p53_R248W_wt']
    print(f"  p53:  neo-WT = {d:+.1f} kJ/mol ({'neo more stable' if d<0 else 'WT more stable'})")

# Now try the full complex from HDOCK
print(f"\n{'='*60}")
print("  COMPLEX INTERACTION ENERGY (HDOCK top pose)")
print(f"{'='*60}")
for name in ['KRAS_G12V_neo', 'KRAS_G12V_wt', 'p53_R248W_neo', 'p53_R248W_wt']:
    fname = f'hdock_models/{name}_complex_top10.pdb'
    try:
        e = energy_of_pdb(fname)
        print(f"  {name}: E_complex = {e:.0f} kJ/mol")
    except Exception as ex:
        print(f"  {name}: ERROR - {ex}")

print("\nDone.")
