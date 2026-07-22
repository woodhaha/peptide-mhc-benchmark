#!/usr/bin/env python3
"""Generate superposition HTML + summary index for centroid PDBs."""
import py3Dmol, webbrowser
from pathlib import Path

DATA = Path("D:/Researching/Peptide epitope/md results/salt_bridge/clean")
OUT = Path("D:/Researching/Peptide epitope/md pic")

files = [
    ("representative_f783", "Rep (78.3ns, main)"),
    ("centroid_0_f378", "C0 (37.8ns, C1)"),
    ("centroid_1_f503", "C1 (50.3ns, C2)"),
    ("centroid_2_f924", "C2 (92.4ns, C3 main)"),
    ("centroid_3_f100", "C3 (10.0ns, C4 early)"),
    ("centroid_4_f156", "C4 (15.6ns, C5)"),
]
colors = ["#7f7f7f", "#e74c3c", "#3498db", "#2ecc71", "#f39c12", "#9b59b6"]
color_names = ["Grey", "Red", "Blue", "Green", "Orange", "Purple"]

print("Building superposition...")
sv = py3Dmol.view(width=1000, height=700)
for idx, (fstem, _) in enumerate(files):
    p = DATA / f"{fstem}_aligned.pdb"
    if not p.exists():
        continue
    with open(p) as f:
        sv.addModel(f.read(), "pdb")
    sv.setStyle({"chain": "C"}, {"stick": {"color": colors[idx], "radius": 0.12}})
    sv.setStyle({"chain": ["A", "B"]}, {"cartoon": {"color": colors[idx], "opacity": 0.03}})

sv.center({"chain": "C"})
sv.zoom(1.1)
with open(OUT / "view_superposition.html", "w") as f:
    f.write(sv._make_html())
print("  view_superposition.html")

# Summary landing page
summary = f"""<!DOCTYPE html><html><head>
<meta charset="utf-8"><title>K-E63 Salt Bridge — Structural Snapshots</title>
<style>
body {{ font-family: sans-serif; max-width: 1200px; margin: auto; padding: 20px; background: #f5f5f5; }}
h1 {{ color: #2c3e50; border-bottom: 2px solid #e74c3c; padding-bottom: 10px; }}
.grid {{ display: grid; grid-template-columns: 1fr 1fr; gap: 20px; }}
.card {{ background: white; border-radius: 8px; box-shadow: 0 2px 8px rgba(0,0,0,0.1); overflow: hidden; }}
.card iframe {{ width: 100%; height: 500px; border: none; }}
.card h3 {{ margin: 10px; color: #2c3e50; }}
.card p {{ margin: 5px 10px 10px; color: #666; font-size: 0.9em; }}
</style></head><body>
<h1>K-E63 Salt Bridge — 6 Snapshots from 100ns MD</h1>
<p>5 KMeans centroids + representative frame (main cluster middle).<br>
Chain A (MHC &alpha;) = blue cartoon, Chain B (&beta;<sub>2</sub>M) = green cartoon, Chain C (peptide) = red stick.<br>
<b>K2</b> (P2, orange) and <b>E63</b> (dark red) highlighted. <i>Drag to rotate, scroll to zoom.</i></p>
<div class="grid">"""

for (fstem, label), cn in zip(files, color_names):
    summary += f'<div class="card"><h3>{label}</h3>'
    summary += f'<iframe src="view_{fstem}.html"></iframe></div>'

summary += """</div>
<h2>Superposition (aligned on MHC C&alpha;)</h2>
<div class="card" style="max-width:100%"><iframe src="view_superposition.html" style="height:600px"></iframe>
<p>Grey=rep &middot; Red=C0 &middot; Blue=C1 &middot; Green=C2 &middot; Orange=C3 &middot; Purple=C4</p></div>
<p><i>Generated with py3Dmol. Rotate with mouse, zoom with scroll wheel.</i></p>
</body></html>"""

with open(OUT / "index_structures.html", "w") as f:
    f.write(summary)
print("  index_structures.html")

print(f"\nOpen: {OUT / 'index_structures.html'}")
webbrowser.open(str(OUT / "index_structures.html"))
