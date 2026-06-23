#!/usr/bin/env python
"""Generate p53 R248W mutant via PyMOL for Phase B."""
import pymol2, os

dock_dir = r'D:\Researching\Peptide epitope\docking\p53'
out_dir = os.path.join(dock_dir, 'figures')
os.makedirs(out_dir, exist_ok=True)

with pymol2.PyMOL() as p:
    # Load WT DBD
    p.cmd.load(os.path.join(dock_dir, '2XWR.pdb'), 'wt')
    p.cmd.remove('resn HOH')

    # Create R248W mutant
    p.cmd.wizard('mutagenesis')
    p.cmd.do('refresh_wizard')
    p.cmd.get_wizard().set_mode('TRP')
    p.cmd.get_wizard().do_select('resi 248 and chain A')
    p.cmd.get_wizard().apply()
    p.cmd.set_wizard()

    # Save mutant
    p.cmd.save(os.path.join(dock_dir, '2XWR_R248W.pdb'), 'wt')
    print('R248W mutant saved: 2XWR_R248W.pdb')

    # ===== Figure: WT vs R248W comparison =====
    p.cmd.bg_color('white')
    p.cmd.set('ray_opaque_background', 0)

    # Show both as cartoon
    p.cmd.show('cartoon')
    p.cmd.color('gray70', 'all')

    # Highlight position 248
    p.cmd.select('site248', 'resi 248')
    p.cmd.show('sticks', 'site248')
    p.cmd.color('red', 'site248')
    p.cmd.set('stick_radius', 0.3, 'site248')

    # Nearby residues
    p.cmd.select('nearby', 'byres (all within 5 of resi 248)')
    p.cmd.show('sticks', 'nearby')
    p.cmd.color('lightblue', 'nearby')

    p.cmd.label('resi 248 and name CA', '"248(W)"')
    p.cmd.set('label_size', 20)

    p.cmd.zoom('nearby', 3)
    p.cmd.ray(1600, 1200)
    p.cmd.png(os.path.join(out_dir, 'p53_R248W_mutant.png'), dpi=300)
    print('Figure saved: p53_R248W_mutant.png')

    # Also save surface representation of mutation site
    p.cmd.show('surface', 'nearby')
    p.cmd.set('transparency', 0.3)
    p.cmd.ray(1600, 1200)
    p.cmd.png(os.path.join(out_dir, 'p53_R248W_surface.png'), dpi=300)
    print('Figure saved: p53_R248W_surface.png')

    print(f'\nDone! Files in {dock_dir}:')
    for f in sorted(os.listdir(dock_dir)):
        sz = os.path.getsize(os.path.join(dock_dir, f))
        print(f'  {f}: {sz:,} bytes')
