#!/usr/bin/env python
"""Analyze p53 DBD structures for R248W docking"""

from Bio.PDB import PDBParser, NeighborSearch
import math

parser = PDBParser(QUIET=True)

# ===== 1. Analyze 6ZNC: p53 DBD + DNA + MQ (APR-246) =====
print("=" * 60)
print("6ZNC: p53 WT DBD + DNA + MQ (APR-246 covalent adduct)")
print("=" * 60)
struct = parser.get_structure('6ZNC', r'D:\Researching\Peptide epitope\docking\p53\6ZNC.pdb')

for chain in struct[0]:
    std_res = [r for r in chain if r.id[0] == ' ']
    het_res = set(r.resname for r in chain if r.id[0] != ' ' and r.resname != 'HOH')
    print(f"Chain {chain.id}: {len(std_res)} aa, het: {het_res}")

# Find MQ
for model in struct:
    for chain in model:
        for res in chain:
            if 'MQ' in res.resname:
                print(f"\nMQ found: chain {chain.id}, resname={res.resname}, id={res.id}")
                for atom in res:
                    if atom.name in ['C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7', 'N', 'O', 'SG', 'S1', 'S2'] + \
                       [a for a in [a.name for a in res] if 'C' in a.name]:
                        print(f"  {atom.name:4s} {atom.coord}")

# ===== 2. Analyze R248 position from 2XWR =====
print()
print("=" * 60)
print("2XWR: p53 WT DBD — R248 position")
print("=" * 60)
struct2 = parser.get_structure('2XWR', r'D:\Researching\Peptide epitope\docking\p53\2XWR.pdb')

for chain in struct2[0]:
    if chain.id == 'A':
        for res in chain:
            if res.id[1] == 248:
                print(f"R248: chain {chain.id}, {res.resname}{res.id[1]}")
                for atom in res:
                    if atom.name in ['CA', 'CB', 'CG', 'CD', 'NE', 'CZ', 'NH1', 'NH2']:
                        print(f"  {atom.name:4s} {atom.coord}")

                # Find nearby DNA-binding residues
                all_atoms = [a for r in chain for a in r if a.element != 'H']
                ns = NeighborSearch(all_atoms)

                print(f"\n  Nearby residues (within 5A of R248 CA):")
                ca = [a for a in res if a.name == 'CA'][0]
                nearby = ns.search(ca.coord, 5.0, level='R')
                for nr in sorted(nearby, key=lambda r: r.id[1]):
                    if nr.id[1] != 248:
                        print(f"    {nr.resname}{nr.id[1]}")

# ===== 3. Analyze C124 and C277 (APR-246 covalent targets) =====
print()
print("=" * 60)
print("C124 and C277 — APR-246 covalent targets")
print("=" * 60)

for pdb_id in ['6ZNC', '2XWR', '2OCJ']:
    struct = parser.get_structure(pdb_id, rf'D:\Researching\Peptide epitope\docking\p53\{pdb_id}.pdb')
    print(f"\n{pdb_id}:")
    for chain in struct[0]:
        if chain.id == 'A':
            for res in chain:
                if res.id[1] in [124, 277, 248] and res.resname in ['CYS', 'ARG']:
                    sg = [a for a in res if a.name == 'SG']
                    nz = [a for a in res if a.name in ['NH1', 'NH2', 'CZ']]
                    print(f"  {res.resname}{res.id[1]}: ", end='')
                    if sg:
                        print(f"SG at {sg[0].coord}")
                    elif nz:
                        print(f"NH1/NH2 at {[a.coord for a in nz]}")
                    else:
                        print("CA at", [a.coord for a in res if a.name == 'CA'][0] if [a for a in res if a.name == 'CA'] else "no CA")

# ===== 4. Key distances =====
print()
print("=" * 60)
print("Key Distances in 6ZNC")
print("=" * 60)
struct = parser.get_structure('6ZNC', r'D:\Researching\Peptide epitope\docking\p53\6ZNC.pdb')

for chain in struct[0]:
    if chain.id == 'A':
        # Find C124, C277, R248, and MQ
        c124_sg = None
        c277_sg = None
        r248_atoms = []
        mq_atoms = []
        zn_atom = None

        for res in chain:
            if res.id[1] == 124 and res.resname == 'CYS':
                c124_sg = [a for a in res if a.name == 'SG'][0]
            if res.id[1] == 277 and res.resname == 'CYS':
                c277_sg = [a for a in res if a.name == 'SG'][0]
            if res.id[1] == 248:
                r248_atoms = [a for a in res if a.name in ['CA', 'CZ', 'NH1', 'NH2']]
            if res.resname == 'ZN':
                zn_atom = [a for a in res][0]

        # Find MQ
        for res in chain:
            if 'MQ' in res.resname:
                mq_atoms = [a for a in res]

        if c124_sg:
            print(f"C124-C277 distance: {c124_sg - c277_sg:.1f}A" if c277_sg else "C277 not found")
            if mq_atoms:
                mq_dist_to_c124 = min([c124_sg - a for a in mq_atoms])
                mq_dist_to_c277 = min([c277_sg - a for a in mq_atoms])
                print(f"MQ closest to C124: {mq_dist_to_c124:.1f}A")
                print(f"MQ closest to C277: {mq_dist_to_c277:.1f}A")

        if zn_atom:
            for atom_name in ['C176_SG', 'H179_NE2', 'C238_SG', 'C242_SG']:
                pass
            print(f"Zn coordination center at {zn_atom.coord}")

# ===== 5. R248W mutation analysis =====
print()
print("=" * 60)
print("R248W Mutation Analysis")
print("=" * 60)
print("""
R248 (Arg, positively charged, long flexible side chain):
  - Normally contacts DNA minor groove via guanidinium group
  - Located in L3 loop (DNA-recognition loop)
  - Mutation to W: removes positive charge, adds bulky hydrophobic indole ring
  - Expected effect:
    1. Loss of DNA contact → DNA binding affinity drops
    2. Tryptophan indole ring may create hydrophobic patch on protein surface
    3. Potential druggable pocket: the empty space where Arg used to reach DNA

Docking strategy for R248W:
  1. In silico mutate R248→W in 2XWR DBD structure
  2. Remove DNA (use 2OCJ/2XWR apo structures)
  3. Define binding box around the R248 region (L3 loop)
  4. Dock small molecules to find:
     a) Stabilizers that fill the Arg-void and restore WT-like conformation
     b) Allosteric modulators at the R248W-induced surface patch
""")
