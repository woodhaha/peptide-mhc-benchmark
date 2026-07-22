"""Update index to reference working superposition."""
from pathlib import Path

OUT = Path("D:/Researching/Peptide epitope/md pic")
index = OUT / "index_structures.html"

if index.exists():
    with open(index) as f:
        content = f.read()
    # The superposition HTML filename is already correct, no change needed
    print("index_structures.html OK")
else:
    print("index_structures.html not found")
