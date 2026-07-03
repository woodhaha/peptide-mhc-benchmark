#!/usr/bin/env python
"""Generate publication figures using PyMOL for Phase A neoepitope analysis."""
import pymol2, os

dock_dir = r'D:\Researching\Peptide epitope\docking'
out_dir = os.path.join(dock_dir, 'figures')
os.makedirs(out_dir, exist_ok=True)

with pymol2.PyMOL() as p:
    # ===== Panel A: Binding groove overview =====
    p.cmd.load(os.path.join(dock_dir, '1DUZ.pdb'), 'complex')
    p.cmd.remove('resn HOH')
    p.cmd.remove('chain D+E+F')
    p.cmd.bg_color('white')
    p.cmd.set('ray_opaque_background', 0)
    p.cmd.set('antialias', 2)

    # Receptor cartoon
    p.cmd.show('cartoon', 'chain A')
    p.cmd.color('gray80', 'chain A')
    p.cmd.show('cartoon', 'chain B')
    p.cmd.color('lightblue', 'chain B')

    # Peptide as yellow sticks
    p.cmd.show('sticks', 'chain C')
    p.cmd.color('yellow', 'chain C')
    p.cmd.set('stick_radius', 0.2, 'chain C')

    # B-pocket (salmon) and F-pocket (orange) residues
    p.cmd.select('bp', 'resi 7+9+45+63+66+67+99+159 and chain A')
    p.cmd.show('sticks', 'bp')
    p.cmd.color('salmon', 'bp')
    p.cmd.select('fp', 'resi 77+80+81+84+116+143+146 and chain A')
    p.cmd.show('sticks', 'fp')
    p.cmd.color('lightorange', 'fp')

    # Labels
    p.cmd.label('resi 63 and chain A and name CA', '"E63"')
    p.cmd.label('resi 2 and chain C and name CA', '"P2(L)"')
    p.cmd.label('resi 9 and chain C and name CA', '"P9(V)"')
    p.cmd.set('label_size', 18)
    p.cmd.set('label_color', 'black')
    p.cmd.set('label_font_id', 5)

    # E63-P2 distance
    p.cmd.distance('d_e63p2', 'resi 63 and chain A and name CD', 'resi 2 and chain C and name CA')
    p.cmd.color('red', 'd_e63p2')

    # Orient and render
    p.cmd.orient('chain C')
    p.cmd.zoom('chain A or chain B or chain C', 1.2)
    p.cmd.ray(1600, 1200)
    p.cmd.png(os.path.join(out_dir, 'FigA_groove_overview.png'), dpi=300)
    print('FigA done')

    # ===== Panel B: B-pocket close-up with semitransparent surface =====
    p.cmd.hide('labels')
    p.cmd.hide('everything')
    p.cmd.show('cartoon', 'chain A and resi 1-180')
    p.cmd.color('gray90', 'chain A')
    p.cmd.show('sticks', 'bp')
    p.cmd.color('salmon', 'bp')
    p.cmd.show('sticks', 'chain C and resi 1+2+3')
    p.cmd.color('yellow', 'chain C and resi 1+2+3')

    # Semi-transparent surface on B-pocket
    p.cmd.show('surface', 'chain A within 10 of bp')
    p.cmd.set('transparency', 0.5)

    p.cmd.label('resi 63 and chain A and name CA', '"E63"')
    p.cmd.label('resi 2 and chain C and name CA', '"P2"')
    p.cmd.set('label_size', 22)

    p.cmd.zoom('bp', 3)
    p.cmd.ray(1600, 1200)
    p.cmd.png(os.path.join(out_dir, 'FigB_b_pocket.png'), dpi=300)
    print('FigB done')

    # ===== Panel C: Schematic model — K-at-P2 vs L-at-P2 =====
    p.cmd.hide('everything')
    p.cmd.show('cartoon', 'chain A and resi 55-72')
    p.cmd.color('gray80', 'chain A')
    p.cmd.show('sticks', 'resi 63 and chain A')
    p.cmd.color('red', 'resi 63 and chain A')  # E63 in red

    p.cmd.show('sticks', 'chain C and resi 2')
    p.cmd.color('yellow', 'chain C and resi 2')

    p.cmd.label('resi 63 and chain A and name OE1', '"E63"')
    p.cmd.set('label_size', 24)

    p.cmd.zoom('resi 63 or (chain C and resi 2)', 5)
    p.cmd.ray(1600, 1200)
    p.cmd.png(os.path.join(out_dir, 'FigC_e63_anchor.png'), dpi=300)
    print('FigC done')

    print(f'\nAll figures saved to {out_dir}')
    print('Files:')
    for f in sorted(os.listdir(out_dir)):
        sz = os.path.getsize(os.path.join(out_dir, f))
        print(f'  {f}: {sz:,} bytes')
