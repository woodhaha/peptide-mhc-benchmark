#!/usr/bin/env python3
"""Strip water/ions from centroid PDBs and write ChimeraX-ready clean PDBs."""
import mdtraj as md
from pathlib import Path

DATA = Path("D:/Researching/Peptide epitope/md results/salt_bridge")
OUT = DATA / "clean"

OUT.mkdir(exist_ok=True)

files = [
    "representative_f514.pdb",
    "centroid_0_f378.pdb",
    "centroid_1_f503.pdb",
    "centroid_2_f924.pdb",
    "centroid_3_f100.pdb",
    "centroid_4_f156.pdb",
]

for fname in files:
    print(f"  {fname}...", end=" ")
    traj = md.load(str(DATA / fname))
    # Keep only chains 0, 1, 2
    keep = [a.index for a in traj.topology.atoms
            if a.residue.chain.index in (0, 1, 2)]
    clean = traj.atom_slice(keep)
    out_name = fname.replace(".pdb", "_clean.pdb")
    clean.save(str(OUT / out_name))
    n_atoms = clean.topology.n_atoms
    print(f"{n_atoms} atoms")

print(f"\nDone → {OUT}")
