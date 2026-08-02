#!/usr/bin/env python3
"""P2 virtual mutagenesis scan — which P2 residues can reach Glu63?

Uses the 1DUZ template (HLA-A*02:01 + LLFGYPVYV peptide). Measures the P2 C-alpha
to Glu63 carboxyl-OE distance from the crystal structure, then for each candidate
P2 residue computes whether its sidechain (C-alpha to terminal atom length) can
reach the carboxyl. Identifies how common energy-compensating non-canonical P2
residues are.

Sidechain reach (C-alpha to terminal atom, from amino-acid geometry):
  K: 7.6 A (N-zeta), R: 7.3 A (N-eta), H: 4.4 A (N-delta), D: 3.0 A, E: 4.4 A,
  Q: 4.5 A, N: 3.6 A, plus neutral small (L/M/V/I/A) ~2-3 A.

Usage: python src/p2_virtual_mutagenesis.py
"""
import json, os
import numpy as np

# P2 sidechain reach from C-alpha to terminal functional atom (A), from rotamer geometry
SIDECHAIN_REACH = {
    # charged / polar — potential electrostatics with Glu63
    "K": 7.6, "R": 7.3, "H": 4.4, "D": 3.0, "E": 4.4, "Q": 4.5, "N": 3.6,
    # canonical A*02:01 P2 (hydrophobic)
    "L": 2.6, "M": 3.6, "I": 2.5, "V": 2.0, "A": 1.5, "T": 2.5,
    # others
    "S": 2.0, "G": 0.0, "P": 2.5, "F": 4.4, "W": 5.3, "Y": 4.9, "C": 3.0,
}
SALT_BRIDGE_MAX = 4.0  # A, N-O contact threshold
CORE_CHARGED = {"K", "R", "H", "D", "E"}

def load_p2_glu63_distance(pdb_path):
    """Read PDB, find peptide P2 CA and HLA Glu63 carboxyl (OE1/OE2), return min distance."""
    residues = {}  # (chain, resSeq) -> (resname, atoms{name: (x,y,z)})
    # 1DUZ: chain A = HLA, chain P or C = peptide (LLFGYPVYV). Glu63 in HLA.
    # P2 = second peptide residue (L in LLFGYPVYV).
    atoms = []
    with open(pdb_path) as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")) and len(line) >= 54:
                name = line[12:16].strip()
                resname = line[17:20].strip()
                chain = line[21]
                resseq = int(line[22:26].strip())
                x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
                atoms.append((chain, resseq, resname, name, (x, y, z)))
    # 1DUZ: chain C (or F) = peptide (residues 1-9, LLFGYPVYV), chain A (or D) = HLA
    # Use the FIRST complex: peptide chain C, HLA chain A.
    p2_ca = None; glu_oe = []
    for chain, resseq, resname, name, coord in atoms:
        # peptide P2 (chain C, resSeq 2, CA atom)
        if chain == "C" and resseq == 2 and name == "CA":
            p2_ca = coord
        # HLA Glu63 (chain A, GLU resSeq 63, carboxyl OE1/OE2)
        if chain == "A" and resname == "GLU" and resseq == 63 and name in ("OE1", "OE2"):
            glu_oe.append(coord)
    if p2_ca is None or not glu_oe:
        return None
    d = min(np.linalg.norm(np.array(p2_ca) - np.array(oe)) for oe in glu_oe)
    return d

def main():
    pdb = "docking/1DUZ.pdb"
    if not os.path.exists(pdb):
        print(f"ERROR: {pdb} not found"); return
    d_p2ca_glu = load_p2_glu63_distance(pdb)
    if d_p2ca_glu is None:
        print("ERROR: could not find P2 CA / Glu63 OE in template"); return
    print(f"P2 C-alpha -> Glu63 carboxyl distance: {d_p2ca_glu:.2f} A\n")
    print(f"{'P2 res':<7}{'reach(A)':<10}{'CA-Glu dist':<13}{'NZ-OE est':<11}{'reachable?':<12}{'class'}")
    print("-"*65)
    results = {}
    for res, reach in sorted(SIDECHAIN_REACH.items(), key=lambda kv: -kv[1]):
        # est terminal-atom-to-Glu = |P2CA-Glu| - reach (sidechain extends toward pocket)
        est = d_p2ca_glu - reach  # if positive, terminal atom falls short; if <=0, reaches past
        # salt bridge requires terminal atom within 4A of Glu OE; P2CA->OE minus reach
        nz_oe = max(d_p2ca_glu - reach, 0.0)  # simplified: reach subtracts directly
        reachable = nz_oe <= SALT_BRIDGE_MAX
        klass = "charged/polar" if res in CORE_CHARGED else "hydrophobic/other"
        results[res] = {"reach": float(reach), "nz_oe_est": float(round(nz_oe, 2)), "reachable": bool(reachable)}
        flag = "  <== can reach Glu63" if reachable else ""
        print(f"{res:<7}{reach:<10.1f}{d_p2ca_glu:<13.2f}{nz_oe:<11.2f}{str(reachable):<12}{klass}{flag}")
    n_charged_reach = sum(1 for r in CORE_CHARGED if results.get(r, {}).get("reachable"))
    n_all_reach = sum(1 for r, v in results.items() if v["reachable"])
    print(f"\nSummary: {n_charged_reach}/{len(CORE_CHARGED)} charged residues can reach Glu63; {n_all_reach}/{len(results)} total P2 residues reach.")
    print("Interpretation: if multiple charged residues reach, energy-compensating non-canonical P2 is NOT unique to K.")
    with open("06_Logs/p2_virtual_scan.json", "w") as f:
        json.dump({"p2_ca_glu63": round(d_p2ca_glu, 2), "results": results}, f, indent=1)

if __name__ == "__main__":
    main()
