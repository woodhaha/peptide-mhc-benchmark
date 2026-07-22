#!/usr/bin/env python3
"""Autodock Vina: re-dock 1DUZ peptide to MHC to validate binding mode."""
import subprocess, os, json, tempfile, shutil, warnings
from pathlib import Path
import numpy as np

warnings.filterwarnings('ignore')

DATA = Path("D:/Researching/Peptide epitope/md results/salt_bridge")
OUT  = Path("D:/Researching/Peptide epitope/scripts")
PDBQT_DIR = OUT / "pdbqt"
PDBQT_DIR.mkdir(exist_ok=True)

# Use final.pdb (stripped) as receptor, extract peptide as ligand
FINAL_PDB = DATA / "final.pdb"

# 1. Create receptor PDB (chains A+B)
print("[1] Extract receptor (MHC + b2m)...", flush=True)
rec_pdb = PDBQT_DIR / "1DUZ_receptor.pdb"
rec_pdbqt = PDBQT_DIR / "1DUZ_receptor.pdbqt"

# Use grep to extract chains A and B (ATOM records)
with open(FINAL_PDB) as f:
    lines = f.readlines()

rec_lines = []
for l in lines:
    if l.startswith("ATOM") or l.startswith("HETATM"):
        chain = l[21:22]
        if chain in ("A", "B"):
            rec_lines.append(l)

with open(rec_pdb, "w") as f:
    f.writelines(rec_lines)
print(f"  Receptor: {sum(1 for l in rec_lines if l.startswith('ATOM'))} atoms", flush=True)

# 2. Create peptide PDB (chain C)
print("[2] Extract ligand (peptide chain C)...", flush=True)
lig_pdb = PDBQT_DIR / "1DUZ_peptide.pdb"
lig_pdbqt = PDBQT_DIR / "1DUZ_peptide.pdbqt"

lig_lines = []
for l in lines:
    if l.startswith("ATOM") or l.startswith("HETATM"):
        chain = l[21:22]
        if chain == "C":
            lig_lines.append(l)

with open(lig_pdb, "w") as f:
    f.writelines(lig_lines)
print(f"  Ligand: {sum(1 for l in lig_lines if l.startswith('ATOM'))} atoms", flush=True)

# 3. Convert to PDBQT using prepare_receptor/prepare_ligand (ADT) or obabel
print("[3] Converting to PDBQT...", flush=True)

def run_obabel(inp, outp):
    cmd = f'obabel "{inp}" -O "{outp}" -h --partialcharge gasteiger'
    r = subprocess.run(cmd, shell=True, capture_output=True, text=True, cwd=str(PDBQT_DIR))
    if r.returncode != 0:
        print(f"  obabel error: {r.stderr[:200]}", flush=True)
    return r.returncode == 0

rec_ok = run_obabel(str(rec_pdb), str(rec_pdbqt))
lig_ok = run_obabel(str(lig_pdb), str(lig_pdbqt))

if not os.path.exists(str(rec_pdbqt)) or not os.path.exists(str(lig_pdbqt)):
    print("Trying alternate obabel syntax...", flush=True)
    for inp, outp, lbl in [(rec_pdb, rec_pdbqt, "receptor"), (lig_pdb, lig_pdbqt, "ligand")]:
        cmd = f'obabel -ipdb "{inp}" -opdbqt -O "{outp}" -h --partialcharge gasteiger'
        r = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        if r.returncode == 0:
            print(f"  {lbl}: OK", flush=True)
        else:
            print(f"  {lbl}: {r.stderr[:200]}", flush=True)

if not os.path.exists(str(rec_pdbqt)) or not os.path.exists(str(lig_pdbqt)):
    print("PDBQT conversion failed.", flush=True)
    exit(1)

# 4. Define search box (centered on P2 binding pocket -> E63 region)
print("[4] Defining search box...", flush=True)
# Use last frame to find E63 CA position
# From the trajectory, E63:CA = atom 968
e63_ca = None
for l in lines:
    if l.startswith("ATOM") and "  CA " in l and l[22:27].strip() == "968":
        x, y, z = float(l[30:38]), float(l[38:46]), float(l[46:54])
        e63_ca = (x, y, z)
        break

if e63_ca is None:
    # Fallback: center of binding pocket
    xs, ys, zs = [], [], []
    for l in lig_lines:
        xs.append(float(l[30:38]))
        ys.append(float(l[38:46]))
        zs.append(float(l[46:54]))
    e63_ca = (np.mean(xs), np.mean(ys), np.mean(zs))

cx, cy, cz = e63_ca
box_size = 22  # Angstrom, large enough for full peptide
print(f"  Box center: {cx:.1f} {cy:.1f} {cz:.1f}", flush=True)
print(f"  Box size: {box_size}A", flush=True)

# 5. Run Vina
print("[5] Running Vina docking...", flush=True)
vina_out = PDBQT_DIR / "vina_output.pdbqt"
vina_cmd = (
    f'vina --receptor "{rec_pdbqt}" --ligand "{lig_pdbqt}" '
    f'--center_x {cx:.3f} --center_y {cy:.3f} --center_z {cz:.3f} '
    f'--size_x {box_size} --size_y {box_size} --size_z {box_size} '
    f'--exhaustiveness 16 --num_modes 9 '
    f'--out "{vina_out}"'
)

print(f"  Running: vina (exhaustiveness=16, 9 modes)", flush=True)
r = subprocess.run(vina_cmd, shell=True, capture_output=True, text=True, timeout=3600)
print(r.stdout[-500:] if len(r.stdout) > 500 else r.stdout, flush=True)

if r.returncode != 0:
    print(f"  Vina error: {r.stderr[:500]}", flush=True)
    exit(1)

# 6. Parse results
print("\n[6] Results:", flush=True)
modes = []
for l in r.stdout.split('\n'):
    if l.strip().startswith('1') or 'REMARK VINA RESULT:' in l:
        parts = l.strip().split()
        try:
            if 'RESULT' in l:
                mode_num = 1
                score = float(parts[3])
            else:
                mode_num = int(parts[0])
                score = float(parts[1])
            modes.append((mode_num, score))
        except:
            pass

print(f"  Docking results (kcal/mol):")
results = {}
for mn, sc in modes[:9]:
    print(f"    Mode {mn}: {sc:.2f}")
    results[f"mode_{mn}"] = sc

json.dump(results, open(PDBQT_DIR / "vina_results.json", "w"), indent=2)

print(f"\n  Output: {vina_out}")
print("  PDBQT files in:", PDBQT_DIR)
