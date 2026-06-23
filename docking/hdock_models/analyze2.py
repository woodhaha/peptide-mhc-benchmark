#!/usr/bin/env python3
"""Analyze HDOCKlite complex models — K-E63 distance and P7 contacts."""
from Bio.PDB import PDBParser
import numpy as np
parser = PDBParser(QUIET=True)
AA = {'ALA':'A','CYS':'C','ASP':'D','GLU':'E','PHE':'F','GLY':'G','HIS':'H',
      'ILE':'I','LYS':'K','LEU':'L','MET':'M','ASN':'N','PRO':'P','GLN':'Q',
      'ARG':'R','SER':'S','THR':'T','VAL':'V','TRP':'W','TYR':'Y'}

def analyze_pdb(fname):
    """Extract K-E63 distance and P7 contacts from top 10 poses."""
    struct = parser.get_structure('x', fname)
    k_e63_all, p7_contacts_all = [], []
    pep_seq = "?"

    for model in struct:
        chains = list(model.get_chains())
        chain_a = chain_c = None
        for c in chains:
            nres = sum(1 for r in c if r.id[0]==' ')
            if nres > 100: chain_a = c
            elif nres < 20: chain_c = c
        if chain_a is None or chain_c is None:
            continue

        # Peptide seq
        pep_res = [r for r in chain_c if r.id[0]==' ']
        pep_seq = ''.join(AA.get(r.resname.strip(),'X') for r in pep_res)

        # P2 CA from peptide
        p2_coord = None
        for r in pep_res:
            rid = r.id[1]
            if isinstance(rid, tuple): rid = rid[0]
            if rid == 2:
                for a in r:
                    if a.name == 'CA': p2_coord = np.array(a.coord)
                break
        if p2_coord is None: continue

        # P7 CA
        p7_coord = None
        for r in pep_res:
            rid = r.id[1]
            if isinstance(rid, tuple): rid = rid[0]
            if rid == 7:
                for a in r:
                    if a.name == 'CA': p7_coord = np.array(a.coord)
                break

        # E63 from chain A
        e63_oe = []
        for r in chain_a:
            rid = r.id[1]
            if isinstance(rid, tuple): rid = rid[0]
            if rid == 63:
                for a in r:
                    if a.name in ('OE1','OE2'):
                        e63_oe.append(np.array(a.coord))

        if not e63_oe or p2_coord is None: continue

        # K-E63 min distance
        d_min = min(np.linalg.norm(p2_coord - oe) for oe in e63_oe)
        k_e63_all.append(float(d_min))

        # P7 contacts
        if p7_coord is not None:
            contacts = 0
            for r in chain_a:
                if r.id[0] != ' ': continue
                for a in r:
                    if a.element == 'H': continue
                    if np.linalg.norm(p7_coord - np.array(a.coord)) < 5.0:
                        contacts += 1
            p7_contacts_all.append(contacts)

    return {
        'pep_seq': pep_seq,
        'k_e63': np.array(k_e63_all),
        'p7_contacts': np.array(p7_contacts_all),
        'n_models': len(k_e63_all),
    }

# Run
FILES = [
    ('KRAS_G12V_neo', 'KRAS_G12V_neo_complex_top10.pdb'),
    ('KRAS_G12V_wt',  'KRAS_G12V_wt_complex_top10.pdb'),
    ('p53_R248W_neo', 'p53_R248W_neo_complex_top10.pdb'),
    ('p53_R248W_wt',  'p53_R248W_wt_complex_top10.pdb'),
]

results = {}
for name, fname in FILES:
    r = analyze_pdb(fname)
    results[name] = r
    print(f"\n{name}: {r['pep_seq'][:9]} ({r['n_models']} poses)")
    print(f"  P2-E63: {r['k_e63'].mean():.1f}+/-{r['k_e63'].std():.1f}A  "
          f"min={r['k_e63'].min():.1f}A  K_reach={(r['k_e63']<7.6).mean()*100:.0f}%")
    if len(r['p7_contacts']) > 0:
        print(f"  P7 contacts <5A: {r['p7_contacts'].mean():.0f}")

# Comparison
print(f"\n{'='*60}")
print("  NEO vs WT COMPARISON")
print(f"{'='*60}")

if 'KRAS_G12V_neo' in results and 'KRAS_G12V_wt' in results:
    neo = results['KRAS_G12V_neo']['k_e63']
    wt  = results['KRAS_G12V_wt']['k_e63']
    print(f"  KRAS K-E63: NEO={neo.mean():.1f}A  WT={wt.mean():.1f}A  delta={neo.mean()-wt.mean():+.1f}A")

if 'p53_R248W_neo' in results and 'p53_R248W_wt' in results:
    neo = results['p53_R248W_neo']['p7_contacts']
    wt  = results['p53_R248W_wt']['p7_contacts']
    if len(neo) and len(wt):
        print(f"  p53 P7 contacts: NEO={neo.mean():.0f}  WT={wt.mean():.0f}  delta={neo.mean()-wt.mean():+.0f}")

print("\nDone.")
