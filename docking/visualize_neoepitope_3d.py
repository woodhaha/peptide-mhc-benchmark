"""Generate 3D py3Dmol view with neoepitope YKLVVVGAV + K-E63 salt bridge.

Mutates P2 Leu→Lys in 1DUZ, places K sidechain toward E63, and renders
with distance labels.
"""
import json, sys, os
import numpy as np
from pathlib import Path

DOCK_DIR = Path(__file__).parent
FIG_DIR = DOCK_DIR / "figures"
FIG_DIR.mkdir(exist_ok=True)

# ═══════════════════════════════════════════════════════════════════════════════
# 1. Mutate P2 (Leu → Lys) in 1DUZ — replace sidechain atoms, keep backbone
# ═══════════════════════════════════════════════════════════════════════════════
def mutate_p2_to_lys(pdb_path):
    """Replace P2 LEU sidechain with LYS, orient toward E63 OE2."""
    with open(pdb_path) as f:
        lines = f.readlines()

    # Find P2 (chain C, residue 2) backbone and sidechain atoms
    p2_atoms = []
    for i, l in enumerate(lines):
        if not l.startswith('ATOM'):
            continue
        chain = l[21]
        resi = l[22:26].strip()
        if chain == 'C' and resi == '2':
            p2_atoms.append((i, l))

    # Get backbone atoms (N, CA, C, O) — keep these
    backbone = {'N', 'CA', 'C', 'O'}
    backbone_lines = {i: l for i, l in p2_atoms if l[12:16].strip() in backbone}

    # Get CA position for placing NZ
    ca_line = None
    for i, l in p2_atoms:
        if l[12:16].strip() == 'CA':
            ca_line = l
            ca_pos = np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])])
            break

    if ca_line is None:
        print("ERROR: No CA found for P2")
        return None

    # Get E63 OE2 position (chain A, residue 63)
    e63_oe2 = None
    for l in lines:
        if (l.startswith('ATOM') and l[21] == 'A' and l[22:26].strip() == '63'
            and l[12:16].strip() == 'OE2'):
            e63_oe2 = np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])])
            break

    if e63_oe2 is None:
        print("ERROR: No E63 OE2 found")
        return None

    # Direction from CA toward OE2
    direction = e63_oe2 - ca_pos
    dist = np.linalg.norm(direction)
    unit_dir = direction / dist
    print(f"P2 CA → E63 OE2: {dist:.1f} Å")

    # Build K sidechain: CB, CG, CD, CE, NZ along the direction toward E63
    # Standard bond lengths ~1.5 Å per C-C bond
    cb = ca_pos + unit_dir * 1.53
    cg = cb + unit_dir * 1.52
    cd = cg + unit_dir * 1.52
    ce = cd + unit_dir * 1.52
    nz = ce + unit_dir * 1.49  # NZ at ~7.6 Å from CA

    print(f"K NZ → E63 OE2: {np.linalg.norm(nz - e63_oe2):.1f} Å")
    print(f"K sidechain length: {np.linalg.norm(nz - ca_pos):.1f} Å")

    # Now build the mutated PDB: keep backbone, add K sidechain
    new_lines = []
    atom_num = 0
    for i, l in enumerate(lines):
        if not l.startswith('ATOM'):
            new_lines.append(l)
            continue

        chain = l[21]
        resi = l[22:26].strip()

        if chain == 'C' and resi == '2':
            atom_name = l[12:16].strip()
            if atom_name in backbone:
                # Keep backbone, rename residue to LYS
                new_l = l[:17] + 'LYS' + l[20:]
                new_lines.append(new_l)
                atom_num = max(atom_num, int(l[6:11]))
            # Skip original LEU sidechain atoms
        else:
            new_lines.append(l)

    # Add K sidechain atoms
    def format_atom(atom_num, name, resn, chain, resi, x, y, z):
        return (f"ATOM  {atom_num:5d}  {name:<3s} {resn} {chain}{resi:4d}    "
                f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           C  ")

    k_sidechain = [
        ('CB', cb), ('CG', cg), ('CD', cd), ('CE', ce), ('NZ', nz)
    ]

    # Insert after P2 CA line — find its index in new_lines
    insert_idx = None
    for idx, l in enumerate(new_lines):
        if (l.startswith('ATOM') and l[21] == 'C' and l[22:26].strip() == '2'
            and l[12:16].strip() == 'CA'):
            insert_idx = idx + 1
            break

    if insert_idx is None:
        insert_idx = len(new_lines)

    for atom_name, pos in k_sidechain:
        atom_num += 1
        new_lines.insert(insert_idx, format_atom(atom_num, atom_name, 'LYS', 'C', 2,
                                                  pos[0], pos[1], pos[2]))
        insert_idx += 1

    # Update atom numbering
    final_lines = []
    atom_num = 0
    for l in new_lines:
        if l.startswith('ATOM'):
            atom_num += 1
            final_lines.append(l[:6] + f"{atom_num:5d}" + l[11:])
        else:
            final_lines.append(l)

    output_path = DOCK_DIR / "1DUZ_K2_mutant.pdb"
    with open(output_path, 'w') as f:
        f.write(''.join(final_lines))
    print(f"✓ Mutant PDB saved: {output_path}")
    return output_path


# ═══════════════════════════════════════════════════════════════════════════════
# 2. Generate py3Dmol views
# ═══════════════════════════════════════════════════════════════════════════════
def render_3d():
    try:
        import py3Dmol
    except ImportError:
        print("py3Dmol not found, skipping")
        return

    mutant_pdb = mutate_p2_to_lys(DOCK_DIR / "1DUZ_minimized.pdb")
    if mutant_pdb is None:
        return

    with open(mutant_pdb) as f:
        pdb_str = f.read()

    # ── View 1: Overview + zoom ──
    view = py3Dmol.view(width=1400, height=900)
    view.setBackgroundColor('white')
    view.addModel(pdb_str, 'pdb')

    # Receptor: light gray cartoon
    view.setStyle({'chain': 'A'}, {'cartoon': {'color': '#BDC3C7', 'opacity': 0.5}})
    view.setStyle({'chain': 'B'}, {'cartoon': {'color': '#D5D8DC', 'opacity': 0.4}})

    # MHC α1-α2 helices — highlight groove
    groove_helix = list(range(49, 85)) + list(range(138, 180))
    view.addStyle({'chain': 'A', 'resi': groove_helix},
                  {'cartoon': {'color': '#95A5A6', 'opacity': 0.5}})

    # Neoepitope peptide: sticks + cartoon
    view.setStyle({'chain': 'C'},
                  {'cartoon': {'color': '#F39C12', 'opacity': 0.6, 'thickness': 0.4},
                   'stick': {'color': '#F39C12', 'radius': 0.2}})

    # P2 K: highlighted in blue with thicker sticks
    view.addStyle({'chain': 'C', 'resi': 2},
                  {'stick': {'color': '#3498DB', 'radius': 0.4}})

    # P9 V: highlighted
    view.addStyle({'chain': 'C', 'resi': 9},
                  {'stick': {'color': '#2ECC71', 'radius': 0.3}})

    # E63: highlighted in red
    view.addStyle({'chain': 'A', 'resi': 63},
                  {'stick': {'color': '#E74C3C', 'radius': 0.4}})

    # B-pocket residues: surface
    pocket_res = [5, 7, 9, 45, 63, 66, 67, 99, 159]
    view.addStyle({'chain': 'A', 'resi': pocket_res},
                  {'stick': {'color': 'salmon', 'radius': 0.12, 'opacity': 0.5}})
    view.addSurface(py3Dmol.VDW, {'chain': 'A', 'resi': pocket_res,
                    'opacity': 0.25, 'color': 'salmon'})

    # Distance line: P2 K NZ → E63 OE2
    # Note: need to find the actual atom positions
    view.addLabel('K(P2)\nYKLVVVGAV',
                  {'position': {'chain': 'C', 'resi': 2, 'atom': 'CA'},
                   'fontColor': '#3498DB', 'fontSize': 12,
                   'backgroundColor': 'white', 'backgroundOpacity': 0.8})
    view.addLabel('E63',
                  {'position': {'chain': 'A', 'resi': 63, 'atom': 'OE2'},
                   'fontColor': '#E74C3C', 'fontSize': 12,
                   'backgroundColor': 'white', 'backgroundOpacity': 0.8})
    view.addLabel('P9 (Val)',
                  {'position': {'chain': 'C', 'resi': 9, 'atom': 'CA'},
                   'fontColor': '#2ECC71', 'fontSize': 10,
                   'backgroundColor': 'white', 'backgroundOpacity': 0.7})

    # Zoom to peptide + B-pocket
    view.zoomTo({'chain': 'C'})
    view.zoom(1.3)

    # Rotate for best view of pocket
    view.rotate(-20, 'x')
    view.rotate(30, 'y')

    # Save
    html_path = FIG_DIR / "K63_neoepitope_3D.html"
    with open(html_path, 'w') as f:
        f.write(view._make_html())
    print(f"\n✓ 3D view saved: {html_path}")
    print(f"  Open in browser: file:///{html_path}")

    # ── View 2: Surface representation ──
    view2 = py3Dmol.view(width=1200, height=800)
    view2.setBackgroundColor('white')
    view2.addModel(pdb_str, 'pdb')

    # Surface for receptor
    view2.addSurface(py3Dmol.VDW, {'chain': 'A', 'opacity': 0.4, 'color': 'white'})
    # Highlight B-pocket surface
    view2.addSurface(py3Dmol.VDW, {'chain': 'A', 'resi': pocket_res,
                     'opacity': 0.5, 'color': 'salmon'})
    # Peptide as sticks through surface
    view2.setStyle({'chain': 'C'},
                   {'stick': {'color': '#F39C12', 'radius': 0.25}})
    view2.addStyle({'chain': 'C', 'resi': 2},
                   {'stick': {'color': '#3498DB', 'radius': 0.4}})
    view2.addStyle({'chain': 'A', 'resi': 63},
                   {'stick': {'color': '#E74C3C', 'radius': 0.35}})

    view2.addLabel('K(P2)\nYKLVVVGAV',
                   {'position': {'chain': 'C', 'resi': 2, 'atom': 'CA'},
                    'fontColor': '#3498DB', 'fontSize': 11,
                    'backgroundColor': 'white', 'backgroundOpacity': 0.7})
    view2.addLabel('E63',
                   {'position': {'chain': 'A', 'resi': 63, 'atom': 'OE2'},
                    'fontColor': '#E74C3C', 'fontSize': 11,
                    'backgroundColor': 'white', 'backgroundOpacity': 0.7})

    view2.zoomTo({'chain': 'C'})
    view2.zoom(1.2)
    view2.rotate(-30, 'x')
    view2.rotate(20, 'y')

    html2_path = FIG_DIR / "K63_neoepitope_surface_3D.html"
    with open(html2_path, 'w') as f:
        f.write(view2._make_html())
    print(f"✓ Surface view saved: {html2_path}")

    print(f"\nDone. YKLVVVGAV (P2=K) modeled into 1DUZ HLA-A*02:01.")


if __name__ == "__main__":
    render_3d()
