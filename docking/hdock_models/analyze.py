#!/usr/bin/env python3
from Bio.PDB import PDBParser
import numpy as np

parser = PDBParser(QUIET=True)
aa_map = {'ALA':'A','CYS':'C','ASP':'D','GLU':'E','PHE':'F','GLY':'G','HIS':'H',
          'ILE':'I','LYS':'K','LEU':'L','MET':'M','ASN':'N','PRO':'P','GLN':'Q',
          'ARG':'R','SER':'S','THR':'T','VAL':'V','TRP':'W','TYR':'Y'}

files = [
    ('KRAS_G12V_neo', 'KRAS_G12V_neo_complex_top10.pdb'),
    ('KRAS_G12V_wt',  'KRAS_G12V_wt_complex_top10.pdb'),
    ('p53_R248W_neo', 'p53_R248W_neo_complex_top10.pdb'),
    ('p53_R248W_wt',  'p53_R248W_wt_complex_top10.pdb'),
]

print("=" * 65)
print("  HDOCK DOCKING — TOP 10 POSE ANALYSIS")
print("=" * 65)

results = {}
for name, fname in files:
    struct = parser.get_structure('x', fname)

    # Analyze all 10 models
    k_e63_dists = []
    p7_contacts_list = []

    for model in struct:
        chains = list(model.get_chains())
        chain_a = chain_c = None
        for c in chains:
            residues = [r for r in c if r.id[0]==' ']
            if len(residues) > 100: chain_a = c
            elif len(residues) < 20: chain_c = c

        if not chain_a or not chain_c:
            continue

        pep_res = [r for r in chain_c if r.id[0]==' ']
        seq = ''.join(aa_map.get(r.resname.strip(),'X') for r in pep_res)

        # E63
        e63_oe1 = e63_oe2 = None
        for r in chain_a:
            if r.id[1] == 63 and isinstance(r.id[1], int):
                for a in r:
                    if a.name == 'OE1': e63_oe1 = a.coord
                    if a.name == 'OE2': e63_oe2 = a.coord

        # P2 CA, P7 CA
        p2_ca = p7_ca = None
        for r in pep_res:
            if r.id[1] == 2:
                for a in r:
                    if a.name == 'CA': p2_ca = a.coord
            if r.id[1] == 7:
                for a in r:
                    if a.name == 'CA': p7_ca = a.coord

        # K-E63 distance
        if p2_ca and e63_oe1 and e63_oe2:
            d = min(np.linalg.norm(p2_ca-e63_oe1), np.linalg.norm(p2_ca-e63_oe2))
            k_e63_dists.append(d)

        # P7 contacts
        if p7_ca:
            contacts = 0
            for r in chain_a:
                if r.id[0] != ' ': continue
                for a in r:
                    if a.element == 'H': continue
                    if np.linalg.norm(p7_ca - a.coord) < 5.0:
                        contacts += 1
            p7_contacts_list.append(contacts)

    if k_e63_dists:
        d_arr = np.array(k_e63_dists)
        results[name] = {
            'k_e63_mean': d_arr.mean(),
            'k_e63_std': d_arr.std(),
            'k_e63_min': d_arr.min(),
            'k_reachable': (d_arr < 7.6).mean() * 100,
            'p7_contacts_mean': np.array(p7_contacts_list).mean() if p7_contacts_list else 0,
            'seq': seq[:9] if 'seq' in dir() else '?',
        }

# Print results
for name, r in results.items():
    print(f"\n{name}:")
    print(f"  Peptide: {r.get('seq','?')}")
    print(f"  P2-E63: mean={r['k_e63_mean']:.1f}+/-{r['k_e63_std']:.1f}A  min={r['k_e63_min']:.1f}A")
    print(f"  K reachable (<7.6A): {r['k_reachable']:.0f}%")
    print(f"  P7 contacts <5A: {r['p7_contacts_mean']:.0f}")

# Neo vs WT comparison
print(f"\n{'='*65}")
print("  K-E63: NEO vs WT COMPARISON")
print(f"{'='*65}")
if 'KRAS_G12V_neo' in results and 'KRAS_G12V_wt' in results:
    neo = results['KRAS_G12V_neo']
    wt = results['KRAS_G12V_wt']
    print(f"  KRAS G12V NEO: P2-E63 = {neo['k_e63_mean']:.1f}A")
    print(f"  KRAS G12V WT:  P2-E63 = {wt['k_e63_mean']:.1f}A")
    print(f"  Delta: {neo['k_e63_mean']-wt['k_e63_mean']:+.1f}A")

print(f"\n{'='*65}")
print("  P7 CONTACTS: NEO vs WT COMPARISON")
print(f"{'='*65}")
if 'p53_R248W_neo' in results and 'p53_R248W_wt' in results:
    neo = results['p53_R248W_neo']
    wt = results['p53_R248W_wt']
    print(f"  p53 R248W NEO: P7 contacts = {neo['p7_contacts_mean']:.0f}")
    print(f"  p53 R248W WT:  P7 contacts = {wt['p7_contacts_mean']:.0f}")
    print(f"  Delta: {neo['p7_contacts_mean']-wt['p7_contacts_mean']:+.0f} contacts")
