"""Fix all py3Dmol generated HTML files: add DOCTYPE, fix CDN URL, remove Quirks Mode."""
from pathlib import Path

OUT = Path("D:/Researching/Peptide epitope/md pic")

files = list(OUT.glob("view_*.html"))
files += [OUT / "index_structures.html"]

LOCAL = "3Dmol-min.js"
cdn_remap = {
    "https://cdn.jsdelivr.net/npm/3dmol@2.5.5/build/3Dmol-min.js": LOCAL,
    "https://3Dmol.org/build/3Dmol-min.js": LOCAL,
    "https://cdnjs.cloudflare.com/ajax/libs/3Dmol/2.5.5/3Dmol-min.js": LOCAL,
}

for f in files:
    if not f.exists():
        continue
    content = f.read_text(encoding="utf-8")

    # Add DOCTYPE if missing
    if not content.startswith("<!DOCTYPE"):
        content = "<!DOCTYPE html>\n" + content

    # Fix CDN URL
    for old, new in cdn_remap.items():
        content = content.replace(old, new)

    # Remove the warning <p> element completely (just the tag, inner text stays inert)
    import re
    content = re.sub(
        r'<p\s+id="3dmolwarning[^"]*"[^>]*>.*?<br></p>',
        "",
        content,
    )

    f.write_text(content, encoding="utf-8")
    print(f"Fixed: {f.name}")

print("All done")
