#!/usr/bin/env python3
"""Generate interactive 3D HTML views of centroids + representative via py3Dmol.
Opens automatically in browser. Also saves static images.
"""
import py3Dmol, json, os, webbrowser
import mdtraj as md
import numpy as np
from pathlib import Path

DATA = Path("D:/Researching/Peptide epitope/md results/salt_bridge/clean")
OUT = Path("D:/Researching/Peptide epitope/md pic")

files = [
    ("representative_f783", "Representative (78.3ns, C3 main)"),
    ("centroid_0_f378", "Centroid 0 (37.8ns, C1)"),
    ("centroid_1_f503", "Centroid 1 (50.3ns, C2)"),
    ("centroid_2_f924", "Centroid 2 (92.4ns, C3 main)"),
    ("centroid_3_f100", "Centroid 3 (10.0ns, C4 early)"),
    ("centroid_4_f156", "Centroid 4 (15.6ns, C5)"),
]

def make_view(pdb_file, title, zoom_to_peptide=True):
    """Create a 3Dmol view with best automatic orientation."""
    view = py3Dmol.view(width=800, height=600)

    with open(pdb_file) as f:
        pdb = f.read()

    # Use py3Dmol.js built-in viewer
    view.addModel(pdb, "pdb")

    # Chain colors
    view.setStyle({'chain': 'A'}, {'cartoon': {'color': '#3498db', 'opacity': 0.7}})
    view.setStyle({'chain': 'B'}, {'cartoon': {'color': '#2ecc71', 'opacity': 0.7}})
    view.setStyle({'chain': 'C'}, {'stick': {'color': '#e74c3c', 'radius': 0.15}})

    # Highlight K2 (P2) and E63
    view.addStyle({'chain': 'C', 'resi': 2},
                  {'stick': {'color': '#ff7f0e', 'radius': 0.25}})
    view.addStyle({'chain': 'A', 'resi': 63},
                  {'stick': {'color': '#d62728', 'radius': 0.25}})

    # Labels
    view.addLabel("K2-P2", {"fontSize": 14, "fontColor": "#ff7f0e",
                            "showBackground": True, "backgroundColor": "white",
                            "opacity": 0.8},
                  {'chain': 'C', 'resi': 2})
    view.addLabel("E63", {"fontSize": 14, "fontColor": "#d62728",
                          "showBackground": True, "backgroundColor": "white",
                          "opacity": 0.8},
                  {'chain': 'A', 'resi': 63})

    if zoom_to_peptide:
        view.zoomTo({'chain': 'C'})

    # Best view: center on binding site
    view.center({'chain': 'C'})
    view.zoom(0.9)

    # Add title
    view.addLabel(f"<b>{title}</b>", {"fontSize": 16, "fontColor": "black",
                  "showBackground": True, "backgroundColor": "white",
                  "useOx": True, "x": 0.05, "y": 0.95})

    return view

# Generate individual HTML files
html_parts = []
for fstem, label in files:
    pdb_file = DATA / f"{fstem}_clean.pdb"
    if not pdb_file.exists():
        print(f"  SKIP {fstem} (not found)")
        continue
    out_html = OUT / f"view_{fstem}.html"
    out_png = "file:///" + str(out_html.resolve()).replace('\\', '/')

    view = make_view(str(pdb_file), label)
    view.zoomTo({'chain': 'C'})
    html = view._make_html()
    with open(out_html, 'w') as f:
        f.write(html)
    html_parts.append((label, out_html))
    print(f"  saved: {out_html.name}")

# ===== Superposition: all centroids =====
print("\nGenerating superposition...")
super_view = py3Dmol.view(width=1000, height=700)

colors = ['#7f7f7f', '#e74c3c', '#3498db', '#2ecc71', '#f39c12', '#9b59b6']
names = ['Rep (78.3ns)', 'C0 (37.8ns)', 'C1 (50.3ns)', 'C2 (92.4ns)', 'C3 (10.0ns)', 'C4 (15.6ns)']

for idx, (fstem, label) in enumerate(files):
    pdb_file = DATA / f"{fstem}_clean.pdb"
    if not pdb_file.exists():
        continue
    # Load with mdtraj and align to first
    traj = md.load(str(pdb_file))
    if idx == 0:
        ref = traj
        coords = ref.xyz[0]
    else:
        # Align CA to reference
        ref_ca = ref.topology.select('name CA and chainid 0')
        mov_ca = traj.topology.select('name CA and chainid 0')
        if len(ref_ca) == len(mov_ca):
            traj.superpose(ref, atom_indices=mov_ca, ref_atom_indices=ref_ca)
        coords = traj.xyz[0]

    # Write aligned PDB
    aligned_pdb = DATA / f"{fstem}_aligned.pdb"
    traj[0].save(str(aligned_pdb))

    with open(aligned_pdb) as f:
        pdb_content = f.read()

    super_view.addModel(pdb_content, "pdb" if idx == 0 else f"pdb")
    # Only show peptide
    super_view.setStyle({'chain': 'C'}, {'stick': {'color': colors[idx], 'radius': 0.12}})
    super_view.setStyle({'chain': ['A', 'B']}, {'cartoon': {'color': colors[idx], 'opacity': 0.05}})

super_view.center({'chain': 'C'})
super_view.zoom(1.1)
super_view.addLabel("<b>Centroid Superposition (aligned on MHC Ca)</b>",
                    {"fontSize": 16, "fontColor": "black",
                     "showBackground": True, "backgroundColor": "white",
                     "useOx": True, "x": 0.05, "y": 0.95})

out_super = OUT / "view_superposition.html"
with open(out_super, 'w') as f:
    f.write(super_view._make_html())
print(f"  saved: {out_super.name}")

# ===== Summary HTML =====
summary = """<!DOCTYPE html><html><head>
<meta charset="utf-8"><title>K-E63 Salt Bridge — Structural Snapshots</title>
<style>
body { font-family: sans-serif; max-width: 1200px; margin: auto; padding: 20px; background: #f5f5f5; }
h1 { color: #2c3e50; border-bottom: 2px solid #3498db; padding-bottom: 10px; }
.grid { display: grid; grid-template-columns: 1fr 1fr; gap: 20px; }
.card { background: white; border-radius: 8px; box-shadow: 0 2px 8px rgba(0,0,0,0.1); overflow: hidden; }
.card iframe { width: 100%; height: 500px; border: none; }
.card h3 { margin: 10px; color: #2c3e50; }
.card p { margin: 5px 10px 10px; color: #666; font-size: 0.9em; }
</style></head><body>
<h1>K-E63 Salt Bridge — Structural Snapshots (100ns MD)</h1>
<p>5 KMeans centroids + representative frame. Chain A (MHC α) = blue, Chain B (β₂M) = green, Chain C (peptide) = red. K2 (P2) = orange, E63 = dark red.</p>
<div class="grid">
"""

for label, out_html in html_parts:
    rel = out_html.name
    summary += f'<div class="card"><h3>{label}</h3>'
    summary += f'<iframe src="{rel}"></iframe></div>\n'

summary += '</div>'
summary += f'<h2>Superposition</h2>'
summary += f'<div class="card"><iframe src="view_superposition.html" style="height:600px"></iframe>'
summary += '<p>All centroids aligned on MHC Cα. Each peptide colored by centroid.</p></div>'
summary += '<p><i>Generated by py3Dmol — rotate/zoom with mouse</i></p>'
summary += '</body></html>'

with open(OUT / "index_structures.html", 'w') as f:
    f.write(summary)

print(f"\nAll done → {OUT}")
print(f"  Open: {OUT}/index_structures.html")
webbrowser.open(str(OUT / "index_structures.html"))
