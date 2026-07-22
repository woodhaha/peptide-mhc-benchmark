#!/usr/bin/env python3
"""Superposition HTML — proper CDN async loading via $3Dmolpromise pattern."""
import base64
from pathlib import Path

OUT = Path("D:/Researching/Peptide epitope/md pic")
pdb_file = OUT / "tmp_superposition.pdb"

with open(pdb_file) as f:
    raw = f.read()

models_raw = raw.split("ENDMDL")
models = []
for m in models_raw:
    m = m.strip()
    if not m:
        continue
    lines = m.split("\n")
    clean = [l for l in lines if not l.startswith("MODEL")]
    models.append("\n".join(clean) + "\nENDMDL\n")
models = models[:6]

colors = ["#7f7f7f", "#e74c3c", "#3498db", "#2ecc71", "#f39c12", "#9b59b6"]
labels = [
    "Rep (78.3ns, main)", "C0 (37.8ns)", "C1 (50.3ns)",
    "C2 (92.4ns, main)", "C3 (10.0ns)", "C4 (15.6ns)",
]

# Build JS array of base64-encoded models
models_js = "var pdbData = [\n"
for m in models:
    b64 = base64.b64encode(m.encode()).decode()
    models_js += f'  atob("{b64}"),\n'
models_js += "];\n"

html = f"""<!DOCTYPE html><html><head>
<meta charset="utf-8"><title>K-E63 Superposition</title>
<style>
body {{ margin: 0; font-family: sans-serif; background: #f0f0f0; overflow: hidden; }}
#container {{ width: 100vw; height: 100vh; position: relative; }}
#viewer {{ width: 100%; height: 100%; }}
#legend {{ position: absolute; top: 20px; right: 20px; z-index: 100;
  background: rgba(255,255,255,0.95); padding: 12px 16px;
  border-radius: 8px; box-shadow: 0 2px 12px rgba(0,0,0,0.15);
  font-size: 14px; line-height: 1.6; min-width: 180px; }}
#legend span {{ display: inline-block; width: 14px; height: 14px; border-radius: 3px;
  margin-right: 6px; vertical-align: middle; }}
#legend .row {{ display: flex; align-items: center; margin: 2px 0; }}
#title {{ position: absolute; top: 20px; left: 20px; z-index: 100;
  background: rgba(255,255,255,0.9); padding: 8px 16px;
  border-radius: 8px; font-weight: bold; font-size: 16px; }}
</style></head><body>
<div id="container">
<div id="title">K-E63 Salt Bridge &mdash; Peptide Superposition<br><span style="font-size:12px;font-weight:normal">6 centroids aligned on MHC C&alpha;</span></div>
<div id="viewer"></div>
<div id="legend"><div style="font-weight:bold;margin-bottom:4px;border-bottom:1px solid #ddd;padding-bottom:4px">Legend</div>"""

for i in range(6):
    html += f'<div class="row"><span style="background:{colors[i]}"></span>{labels[i]}</div>\n'

html += """</div></div>
<script>
var loadScriptAsync = function(uri){
  return new Promise((resolve, reject) => {
    var tag = document.createElement('script');
    tag.src = uri;
    tag.async = true;
    tag.onload = () => resolve();
    var first = document.getElementsByTagName('script')[0];
    first.parentNode.insertBefore(tag, first);
  });
};
if(typeof $3Dmolpromise === 'undefined') {
  $3Dmolpromise = null;
  $3Dmolpromise = loadScriptAsync('https://cdn.jsdelivr.net/npm/3dmol@2.5.5/build/3Dmol-min.js');
}
"""

html += models_js

html += """
var colors = ["#7f7f7f", "#e74c3c", "#3498db", "#2ecc71", "#f39c12", "#9b59b6"];
var labels = ["Rep (78.3ns, main)", "C0 (37.8ns)", "C1 (50.3ns)", "C2 (92.4ns, main)", "C3 (10.0ns)", "C4 (15.6ns)"];

$3Dmolpromise.then(function() {
  var viewer = $3Dmol.createViewer("viewer", {defaultColors: {backgroundColor: "white"}});
  for (var i = 0; i < pdbData.length; i++) {
    viewer.addModel(pdbData[i], "pdb");
    viewer.setStyle({model: i}, {stick: {color: colors[i], radius: 0.14}});
  }
  viewer.zoomTo();
  viewer.center();
  viewer.render();
});
</script>
</body></html>
"""

with open(OUT / "view_superposition.html", "w") as f:
    f.write(html)
print(f"Written: {OUT / 'view_superposition.html'} ({len(html)} bytes, {len(models)} models)")
