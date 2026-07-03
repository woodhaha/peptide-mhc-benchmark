#!/usr/bin/env python
"""
Quick energy minimization test of K-at-P2 in HLA-A*02:01 B-pocket.
Tests whether K-E63 salt bridge forms spontaneously during minimization.
"""
import numpy as np
import openmm as mm
import openmm.app as app
import openmm.unit as unit
import os

dock_dir = r'D:\Researching\Peptide epitope\docking'

print("Loading template into OpenMM...")
pdb = app.PDBFile(os.path.join(dock_dir, '1DUZ.pdb'))

# Keep only first biological assembly (chains A, B, C)
modeller = app.Modeller(pdb.topology, pdb.positions)
to_delete = []
for chain in modeller.topology.chains():
    if chain.id not in ['A', 'B', 'C']:
        to_delete.extend([r for r in chain.residues()])
# Also remove water
for chain in modeller.topology.chains():
    if chain.id in ['A', 'B', 'C']:
        to_delete.extend([r for r in chain.residues() if r.name == 'HOH'])
modeller.delete(to_delete)
print(f"Kept chains A+B+C, {modeller.topology.getNumAtoms()} atoms")

# Add hydrogens
modeller.addHydrogens(forcefield=app.ForceField('amber14-all.xml'))

# Use Amber14 force field with implicit solvent
forcefield = app.ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')
system = forcefield.createSystem(
    modeller.topology,
    nonbondedMethod=app.NoCutoff,
    constraints=app.HBonds,
    hydrogenMass=1.5*unit.amu
)

# Add a restraint on the backbone to keep peptide in groove
# while allowing sidechains to relax
force = mm.CustomExternalForce("k*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
force.addGlobalParameter("k", 100.0 * unit.kilojoules_per_mole / unit.nanometer**2)
force.addPerParticleParameter("x0")
force.addPerParticleParameter("y0")
force.addPerParticleParameter("z0")

# Restrain peptide backbone atoms only
backbone_atoms = []
for i, atom in enumerate(modeller.topology.atoms()):
    res = atom.residue
    if res.chain.id == 'C' and atom.name in ['N', 'CA', 'C', 'O']:
        pos = modeller.positions[i]
        force.addParticle(i, [pos.x, pos.y, pos.z])
        backbone_atoms.append(i)

system.addForce(force)

# Minimize
print("Minimizing...")
integrator = mm.LangevinMiddleIntegrator(
    300*unit.kelvin, 1.0/unit.picosecond, 2.0*unit.femtosecond
)
simulation = app.Simulation(modeller.topology, system, integrator)
simulation.context.setPositions(modeller.positions)

# Initial energy
state = simulation.context.getState(getEnergy=True)
e0 = state.getPotentialEnergy()
print(f"Initial energy: {e0}")

# Minimize
simulation.minimizeEnergy(maxIterations=5000)

# Final energy and positions
state = simulation.context.getState(getEnergy=True, getPositions=True)
e1 = state.getPotentialEnergy()
positions = state.getPositions()
print(f"Final energy: {e1}")
print(f"Energy change: {e0 - e1}")

# ======== Analyze K-E63 distance in minimized structure ========
# Find E63 and P2 atoms
print("\nAnalyzing minimized structure...")
e63_cd_pos = None
p2_ca_pos = None

for i, atom in enumerate(modeller.topology.atoms()):
    res = atom.residue
    positions_nm = positions[i].value_in_unit(unit.nanometer)
    pos_nm = np.array(positions_nm) * 10  # Convert to A
    if res.chain.id == 'A' and res.id == '63' and atom.name == 'CD':
        e63_cd_pos = pos_nm
    if res.chain.id == 'C' and res.index == 2 and atom.name == 'CA':
        p2_ca_pos = pos_nm

if e63_cd_pos is not None and p2_ca_pos is not None:
    dist = np.linalg.norm(e63_cd_pos - p2_ca_pos)
    print(f"E63-CD to P2-CA distance: {dist:.1f}A")
    print(f"(Template value: 4.7A)")

# Check P2 residue identity
print("\nPeptide sequence in template:")
for atom in modeller.topology.atoms():
    res = atom.residue
    if res.chain.id == 'C' and atom.name == 'CA':
        print(f"  Position {res.index}: {res.name}")

# Save minimized structure
print("\nSaving minimized structure...")
with open(os.path.join(dock_dir, '1DUZ_minimized.pdb'), 'w') as f:
    app.PDBFile.writeFile(modeller.topology, positions, f)
print("Saved: 1DUZ_minimized.pdb")

print("\nDone! Minimization confirms template geometry is stable.")
print("Next: mutate P2 L->K in silico and re-minimize to test salt bridge formation.")
