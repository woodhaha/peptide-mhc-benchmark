#!/usr/bin/env python3
"""Superposition via py3Dmol direct - proven to work for individual views."""
import py3Dmol
from pathlib import Path

OUT = Path("D:/Researching/Peptide epitope/md pic")
TMP = OUT

colors = ["#7f7f7f", "#e74c3c", "#3498db", "#2ecc71", "#f39c12", "#9b59b6"]

v = py3Dmol.view(width=1000, height=700, linked=False)

for i in range(6):
    pdb_file = TMP / f"_pep_{i}.pdb"
    if not pdb_file.exists():
        print(f"  SKIP {pdb_file}")
        continue
    with open(pdb_file) as f:
        pdb_text = f.read()
    v.addModel(pdb_text, "pdb")
    v.setStyle({"model": i}, {"stick": {"color": colors[i], "radius": 0.14}})

v.center({"chain": "C"})
v.zoom(1.2)

html = v._make_html()
html = html.replace("3dmolviewer_", "3dmolsuper_")

with open(OUT / "view_superposition.html", "w") as f:
    f.write(html)

print(f"Done: view_superposition.html ({len(html)} bytes)")
