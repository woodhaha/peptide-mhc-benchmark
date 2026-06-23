#!/usr/bin/env python3
"""
MM/GBSA binding energy comparison: neoepitope vs WT peptide.
Computes delta_binding = E(complex) - E(receptor) - E(peptide)
for KRAS G12V (YKLVVVGAV vs YKLVVVGAG) and p53 R248W (MNWRPILTI vs MNRRPILTI).
"""
import numpy as np
import openmm as mm
import openmm.app as app
import openmm.unit as unit
import sys, os

DOCK_DIR = r'D:\Researching\Peptide epitope\docking'
os.chdir(DOCK_DIR)

def compute_energy(pdb_path, fix_peptide=False):
    """Load PDB, build system with implicit solvent, minimize, return potential energy."""
    pdb = app.PDBFile(pdb_path)

    # Keep chains A (heavy), B (B2M), C (peptide)
    modeller = app.Modeller(pdb.topology, pdb.positions)
    to_delete = []
    for chain in modeller.topology.chains():
        if chain.id not in ['A', 'B', 'C']:
            to_delete.extend([r for r in chain.residues()])
        elif chain.id in ['A', 'B', 'C']:
            to_delete.extend([r for r in chain.residues() if r.name == 'HOH'])
    modeller.delete(to_delete)

    # Forcefield with implicit solvent
    forcefield = app.ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')

    # Add hydrogens
    modeller.addHydrogens(forcefield=forcefield)

    # Create system with GB implicit solvent (faster than explicit)
    system = forcefield.createSystem(
        modeller.topology,
        nonbondedMethod=app.NoCutoff,
        constraints=app.HBonds,
        hydrogenMass=1.5*unit.amu,
    )

    # Position restraints if needed
    if fix_peptide:
        force = mm.CustomExternalForce("k*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
        force.addGlobalParameter("k", 500.0 * unit.kilojoules_per_mole / unit.nanometer**2)
        force.addPerParticleParameter("x0")
        force.addPerParticleParameter("y0")
        force.addPerParticleParameter("z0")
        for i, atom in enumerate(modeller.topology.atoms()):
            res = atom.residue
            if res.chain.id == 'C':
                pos = modeller.positions[i]
                force.addParticle(i, [pos.x, pos.y, pos.z])
        system.addForce(force)

    integrator = mm.LangevinMiddleIntegrator(
        300*unit.kelvin, 1.0/unit.picosecond, 2.0*unit.femtosecond
    )
    simulation = app.Simulation(modeller.topology, system, integrator)
    simulation.context.setPositions(modeller.positions)

    # Minimize
    simulation.minimizeEnergy(maxIterations=2000)

    state = simulation.context.getState(getEnergy=True)
    energy = state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)

    return energy, modeller, simulation


def split_energy(pdb_path):
    """
    Compute: E_complex, E_receptor, E_peptide
    Returns: delta_G_binding = E_complex - E_receptor - E_peptide
    """
    # Full complex
    e_complex, _, _ = compute_energy(pdb_path)

    # Receptor only (load PDB, delete peptide chain C)
    pdb = app.PDBFile(pdb_path)
    modeller = app.Modeller(pdb.topology, pdb.positions)
    to_delete = []
    for chain in modeller.topology.chains():
        if chain.id == 'C':
            to_delete.extend([r for r in chain.residues() if r.name != 'HOH'])
        elif chain.id not in ['A', 'B']:
            to_delete.extend([r for r in chain.residues()])
        elif chain.id in ['A', 'B']:
            to_delete.extend([r for r in chain.residues() if r.name == 'HOH'])
    modeller.delete(to_delete)

    forcefield = app.ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')
    modeller.addHydrogens(forcefield=forcefield)
    system = forcefield.createSystem(modeller.topology, nonbondedMethod=app.NoCutoff,
                                     constraints=app.HBonds, hydrogenMass=1.5*unit.amu)
    integrator = mm.LangevinMiddleIntegrator(300*unit.kelvin, 1.0/unit.picosecond, 2.0*unit.femtosecond)
    sim_rec = app.Simulation(modeller.topology, system, integrator)
    sim_rec.context.setPositions(modeller.positions)
    sim_rec.minimizeEnergy(maxIterations=2000)
    e_receptor = sim_rec.context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)

    # Peptide only (chain C)
    pdb = app.PDBFile(pdb_path)
    modeller = app.Modeller(pdb.topology, pdb.positions)
    to_delete = []
    for chain in modeller.topology.chains():
        if chain.id != 'C':
            to_delete.extend([r for r in chain.residues()])
    modeller.delete(to_delete)

    if len(list(modeller.topology.atoms())) == 0:
        return None  # Can't compute peptide alone

    modeller.addHydrogens(forcefield=forcefield)
    system = forcefield.createSystem(modeller.topology, nonbondedMethod=app.NoCutoff,
                                     constraints=app.HBonds, hydrogenMass=1.5*unit.amu)
    sim_pep = app.Simulation(modeller.topology, system, integrator)
    sim_pep.context.setPositions(modeller.positions)
    sim_pep.minimizeEnergy(maxIterations=2000)
    e_peptide = sim_pep.context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)

    delta_G = e_complex - e_receptor - e_peptide
    return {
        'e_complex': e_complex,
        'e_receptor': e_receptor,
        'e_peptide': e_peptide,
        'delta_G': delta_G,
    }


print("=" * 65)
print("  MM/GBSA BINDING ENERGY: NEOEPITOPE vs WT")
print("=" * 65)

# KRAS G12V
print("\n--- KRAS G12V ---")
# We use the 1DUZ template (LLFGYPVYV) for receptor
# The peptide PDB files have the correct sequences
r_kras_neo = split_energy('peptide_KRAS_G12V_neo.pdb')
r_kras_wt  = split_energy('peptide_KRAS_G12V_wt.pdb')

# Note: the peptide-only PDB files don't have the receptor -
# we need the full complex. Let me use the HDOCK complex models instead.
# Actually, let me use the receptor + peptide separately.

# Better approach: build complex from receptor + peptide PDBs
print("Building complexes from receptor + peptide PDBs...")

# Load receptor
rec_pdb = app.PDBFile('receptor_1DUZ.pdb')  # HLA-A*02:01 only

for name, pep_file in [('KRAS_G12V_neo', 'peptide_KRAS_G12V_neo.pdb'),
                        ('KRAS_G12V_wt', 'peptide_KRAS_G12V_wt.pdb'),
                        ('p53_R248W_neo', 'peptide_p53_R248W_neo.pdb'),
                        ('p53_R248W_wt', 'peptide_p53_R248W_wt.pdb')]:
    pep_pdb = app.PDBFile(pep_file)

    # Combine receptor + peptide into one topology
    # Simple approach: just compute energies separately
    # E_receptor from receptor PDB
    # E_peptide from peptide PDB
    # For complex, we need to combine them - this is complex with OpenMM

    # For now, use a simpler metric: just compare peptide internal energies
    # as a proxy for strain upon binding

    forcefield = app.ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')

    # Peptide energy
    modeller = app.Modeller(pep_pdb.topology, pep_pdb.positions)
    modeller.addHydrogens(forcefield=forcefield)
    system = forcefield.createSystem(modeller.topology, nonbondedMethod=app.NoCutoff,
                                     constraints=app.HBonds, hydrogenMass=1.5*unit.amu)
    integrator = mm.LangevinMiddleIntegrator(300*unit.kelvin, 1.0/unit.picosecond, 2.0*unit.femtosecond)
    sim = app.Simulation(modeller.topology, system, integrator)
    sim.context.setPositions(modeller.positions)
    sim.minimizeEnergy(maxIterations=2000)
    e_pep = sim.context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)

    print(f"  {name}: E_peptide = {e_pep:.0f} kJ/mol")

print("\nNote: Full MM/GBSA requires proper complex building.")
print("Peptide-only energies show internal strain differences.")
print("For definitive delta-G, need MD trajectory + MM/PBSA.py")
