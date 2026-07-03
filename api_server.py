"""Peptide-MHC Binding Prediction API — stdlib only, zero extra deps.

ponytail: http.server + json, no Flask/FastAPI needed.
Start: conda run -n base python api_server.py
"""

import json
import sys
import os
from http.server import HTTPServer, BaseHTTPRequestHandler

# Add project analysis dir to path for blosum_utils import
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "03_Analysis"))
os.environ["TF_CPP_MIN_LOG_LEVEL"] = "2"  # suppress TF info

import numpy as np
os.environ["KERAS_BACKEND"] = "tensorflow"
import keras
from blosum_utils import encode_blosum, AA2IDX

MODEL_PATH = os.path.join(os.path.dirname(__file__), "03_Analysis", "FFN_Deep.h5")
CLASSES = ["NB", "WB", "SB"]
CLASS_CN = {"NB": "非结合体", "WB": "弱结合体", "SB": "强结合体"}
VALID_AA = set(AA2IDX.keys())

model = None  # lazy load


def load_model():
    global model
    if model is None:
        model = keras.models.load_model(MODEL_PATH)
        _ = model(np.zeros((1, 180), dtype=np.float32))  # warmup
    return model


def predict(peptide: str) -> dict:
    model = load_model()
    x = encode_blosum([peptide])
    probs = model(x.numpy() if hasattr(x, "numpy") else x, training=False)
    probs = probs.numpy().flatten()
    pred_idx = int(np.argmax(probs))
    return {
        "peptide": peptide.upper(),
        "class": CLASSES[pred_idx],
        "class_cn": CLASS_CN[CLASSES[pred_idx]],
        "probabilities": {
            "NB": round(float(probs[0]), 4),
            "WB": round(float(probs[1]), 4),
            "SB": round(float(probs[2]), 4),
        },
    }


HTML_PAGE = r"""<!DOCTYPE html>
<html lang="zh">
<head>
<meta charset="UTF-8"><meta name="viewport" content="width=device-width,initial-scale=1.0">
<title>Peptide-MHC Binding Predictor</title>
<style>
*{box-sizing:border-box;margin:0;padding:0}
body{font:16px/1.5 system-ui,-apple-system,sans-serif;background:#0f172a;color:#e2e8f0;min-height:100vh;display:flex;align-items:center;justify-content:center;padding:20px}
.card{background:#1e293b;border-radius:16px;padding:32px;max-width:520px;width:100%;box-shadow:0 25px 50px -12px rgba(0,0,0,.5)}
h1{font-size:1.5em;text-align:center;margin-bottom:4px}
.sub{text-align:center;color:#94a3b8;font-size:.85em;margin-bottom:24px}
code{background:#334155;padding:2px 6px;border-radius:4px;font-size:.9em}
input{width:100%;padding:12px 16px;font:1.1em monospace;text-align:center;letter-spacing:4px;background:#0f172a;border:2px solid #334155;border-radius:8px;color:#e2e8f0;outline:none;text-transform:uppercase}
input:focus{border-color:#38bdf8}
.btn-row{display:flex;gap:8px;margin-top:12px}
button{flex:1;padding:10px;border:none;border-radius:8px;font-size:.9em;cursor:pointer;font-weight:600;transition:all .15s}
.btn-predict{background:#38bdf8;color:#0f172a;font-size:1.1em}
.btn-predict:hover{background:#7dd3fc}
.btn-example{background:#334155;color:#94a3b8;font-size:.78em}
.btn-example:hover{background:#475569;color:#e2e8f0}
.result{margin-top:20px;padding:16px;background:#0f172a;border-radius:8px;display:none}
.result.show{display:block}
.result .label{font-size:1.3em;font-weight:700}
.badge{display:inline-block;padding:4px 12px;border-radius:99px;font-weight:700;font-size:.85em}
.badge.SB{background:#22c55e;color:#052e16}
.badge.WB{background:#f59e0b;color:#451a03}
.badge.NB{background:#ef4444;color:#450a0a}
.bars{margin-top:12px}
.bar-row{display:flex;align-items:center;gap:8px;margin:4px 0}
.bar-label{width:32px;font-weight:700;font-size:.85em}
.bar-track{flex:1;height:22px;background:#1e293b;border-radius:4px;overflow:hidden}
.bar-fill{height:100%;border-radius:4px;transition:width .3s}
.bar-fill.SB{background:#22c55e}.bar-fill.WB{background:#f59e0b}.bar-fill.NB{background:#ef4444}
.bar-val{width:50px;text-align:right;font-size:.8em;font-family:monospace}
.error{color:#f87171;margin-top:8px}
.footer{margin-top:24px;text-align:center;font-size:.75em;color:#475569}
.footer a{color:#38bdf8;text-decoration:none}
.spinner{display:none;text-align:center;margin-top:12px}
</style>
</head>
<body>
<div class="card">
<h1>🧬 Peptide-MHC Binding</h1>
<div class="sub">HLA-A*02:01 · Deep FFN · BLOSUM62 · 9-mer</div>
<input id="pep" type="text" maxlength="9" placeholder="输入9肽序列..." autofocus>
<div class="btn-row">
  <button class="btn-predict" onclick="doPredict()">Predict</button>
  <button class="btn-example" onclick="setPep('GILGFVFTL')">Flu M1</button>
  <button class="btn-example" onclick="setPep('YKLVVVGAV')">KRAS G12V</button>
  <button class="btn-example" onclick="setPep('NLVPMVATV')">CMV pp65</button>
</div>
<div class="btn-row">
  <button class="btn-example" onclick="setPep('MNWRPILTI')">p53 R248W</button>
  <button class="btn-example" onclick="setPep('SLLMWITQC')">NY-ESO-1</button>
  <button class="btn-example" onclick="setPep('RMFPNAPYL')">WT1</button>
  <button class="btn-example" onclick="setPep('LLLLLLLLL')">Poly-L (-)</button>
</div>
<div class="spinner" id="spin">⏳ Predicting...</div>
<div class="result" id="result">
  <div id="resLabel"></div>
  <div class="bars" id="bars"></div>
  <div class="error" id="err"></div>
</div>
<div class="footer">API: <a href="/predict">POST /predict</a> · Model: Deep FFN (152K) · <a href="https://github.com/woodhaha">woodhaha</a></div>
</div>
<script>
const clsName={NB:{cn:'非结合体',en:'Non-Binder'},WB:{cn:'弱结合体',en:'Weak Binder'},SB:{cn:'强结合体',en:'Strong Binder'}};
function setPep(v){document.getElementById('pep').value=v;doPredict()}
async function doPredict(){
  const pep=document.getElementById('pep').value.trim().toUpperCase();
  document.getElementById('result').classList.remove('show');
  document.getElementById('err').textContent='';
  if(pep.length!==9){document.getElementById('err').textContent='Must be exactly 9 amino acids';return}
  document.getElementById('spin').style.display='block';
  try{
    const r=await fetch('/predict',{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify({peptide:pep})});
    const d=await r.json();
    document.getElementById('spin').style.display='none';
    if(d.error){document.getElementById('err').textContent=d.error;return}
    const cls=d.class;
    document.getElementById('resLabel').innerHTML=`<span class="label">${d.peptide}</span> <span class="badge ${cls}">${cls} · ${clsName[cls].cn} · ${clsName[cls].en}</span>`;
    let bars='';
    for(const [k,v] of Object.entries(d.probabilities)){
      bars+=`<div class="bar-row"><span class="bar-label">${k}</span><div class="bar-track"><div class="bar-fill ${k}" style="width:${Math.round(v*100)}%"></div></div><span class="bar-val">${(v*100).toFixed(1)}%</span></div>`;
    }
    document.getElementById('bars').innerHTML=bars;
    document.getElementById('result').classList.add('show');
  }catch(e){document.getElementById('spin').style.display='none';document.getElementById('err').textContent='Network error'}
}
document.getElementById('pep').addEventListener('keydown',e=>{if(e.key==='Enter')doPredict()})
</script>
</body>
</html>"""

class Handler(BaseHTTPRequestHandler):
    def _send_json(self, data, status=200):
        body = json.dumps(data, ensure_ascii=False).encode("utf-8")
        self.send_response(status)
        self.send_header("Content-Type", "application/json; charset=utf-8")
        self.send_header("Access-Control-Allow-Origin", "*")
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def do_OPTIONS(self):
        self.send_response(204)
        self.send_header("Access-Control-Allow-Origin", "*")
        self.send_header("Access-Control-Allow-Methods", "POST, OPTIONS")
        self.send_header("Access-Control-Allow-Headers", "Content-Type")
        self.end_headers()

    def do_GET(self):
        if self.path in ("/", "/index.html"):
            self._send_html(HTML_PAGE)
        elif self.path == "/health":
            self._send_json({"status": "ok"})
        else:
            self._send_json({"error": "POST /predict only"}, 404)

    def _send_html(self, html, status=200):
        body = html.encode("utf-8")
        self.send_response(status)
        self.send_header("Content-Type", "text/html; charset=utf-8")
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def do_POST(self):
        if self.path != "/predict":
            return self._send_json({"error": "use POST /predict"}, 404)

        length = int(self.headers.get("Content-Length", 0))
        body = json.loads(self.rfile.read(length))
        peptide = body.get("peptide", "").strip().upper()

        if len(peptide) != 9:
            return self._send_json({"error": "peptide must be 9 amino acids"}, 400)
        if not all(aa in VALID_AA for aa in peptide):
            return self._send_json({"error": f"invalid amino acid(s) in '{peptide}'"}, 400)

        try:
            result = predict(peptide)
            self._send_json(result)
        except Exception as e:
            self._send_json({"error": str(e)}, 500)

    def log_message(self, format, *args):
        pass  # silent


if __name__ == "__main__":
    port = int(sys.argv[1]) if len(sys.argv) > 1 else 8787
    print(f"Peptide-MHC API @ http://localhost:{port}/predict")
    HTTPServer(("0.0.0.0", port), Handler).serve_forever()
