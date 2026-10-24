#!/usr/bin/env python3
"""One-command build + reproducibility release.

Rebuilds the manuscript (authoryear -> numbered -> docx -> pdf), runs the
consistency audit, then snapshots the key data/code (SHA256 manifest),
pins the environment, and writes a Zenodo-ready metadata file.

Usage:
    python scripts/build_release.py            # build + package (no upload)
    ZENODO_TOKEN=... python scripts/build_release.py --upload  # also upload
"""
import argparse, hashlib, json, os, subprocess, sys
from pathlib import Path
from datetime import date

ROOT = Path(__file__).resolve().parent.parent
RELEASE = ROOT / "release"

# Files whose integrity matters for reproducibility (data + code + manuscript).
SNAPSHOT_FILES = [
    "02_Data/raw/real_peptides.csv",
    "02_Data/cleaned/iedb_epitope_verification.csv",
    "02_Data/cleaned/iedb_benchmark_results_v2.csv",
    "02_Data/cleaned/netmhcpan_results.csv",
    "02_Data/cleaned/multiseed_results.json",
    "experimental_label_scores.json",
    "netmhcpan_dl/netmhcpan_42_standalone.xls",
    "netmhcpan_dl/netmhcpan_6524.xls",
    "bigmhc_519.prd",
    ".munis_dl/predictions/pep519_munis_predictions.csv",
    "scripts/convert_refs_numbered.py",
    "scripts/regenerate_figure8.py",
    "scripts/build_release.py",
]


def run(cmd, cwd=None):
    print("  $", " ".join(cmd) if isinstance(cmd, list) else cmd)
    r = subprocess.run(cmd, cwd=cwd or ROOT, shell=isinstance(cmd, str))
    if r.returncode != 0:
        print(f"  !! failed (exit {r.returncode})")
        sys.exit(r.returncode)


def build_manuscript():
    print("[1/4] convert authoryear -> numbered")
    run([sys.executable, "scripts/convert_refs_numbered.py"])
    print("[2/4] consistency audit")
    run([sys.executable, "peptide_audit.py"])
    print("[3/4] docx + pdf")
    sub = ROOT / "06_Bioinformatics_Submission"
    run(["pandoc", "manuscript_plos_numbered.md", "-o", "manuscript_plos.docx"], cwd=sub)
    run(["pandoc", "manuscript_plos_numbered.md", "-f", "markdown", "-t", "html5",
         "-s", "--metadata", "title=", "-o", "manuscript_plos.html"], cwd=sub)
    run([sys.executable, "-c",
         "from weasyprint import HTML; HTML('manuscript_plos.html').write_pdf('manuscript_plos.pdf')"],
         cwd=sub)
    (sub / "manuscript_plos.html").unlink(missing_ok=True)


def snapshot():
    print("[4/4] snapshot + package")
    RELEASE.mkdir(exist_ok=True)
    manifest = {}
    for rel in SNAPSHOT_FILES:
        p = ROOT / rel
        if p.exists():
            h = hashlib.sha256(p.read_bytes()).hexdigest()
            manifest[rel] = {"sha256": h, "bytes": p.stat().st_size}
    (RELEASE / "manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")

    # Pin environment
    r = subprocess.run([sys.executable, "-m", "pip", "freeze"], capture_output=True, text=True)
    (RELEASE / "requirements.txt").write_text(r.stdout, encoding="utf-8")

    # Package the reproducibility-critical files into a tarball for the Zenodo deposit
    import tarfile
    with tarfile.open(RELEASE / "source_data.tar.gz", "w:gz") as tf:
        for rel in SNAPSHOT_FILES:
            p = ROOT / rel
            if p.exists():
                tf.add(p, arcname=rel)

    # Zenodo metadata from CITATION.cff (best-effort)
    try:
        import yaml
        cff = yaml.safe_load((ROOT / "CITATION.cff").read_text(encoding="utf-8"))
        creators = []
        for a in cff.get("authors", []):
            fam = a.get("family-names", ""); giv = a.get("given-names", "")
            creators.append({"person_or_org": {
                "type": "personal",
                "name": f"{fam}, {giv}".strip(", "),
                "given_name": giv, "family_name": fam,
            }})
        zenodo = {
            "metadata": {
                "title": cff.get("title", ""),
                "description": cff.get("abstract", ""),
                "publication_date": date.today().isoformat(),
                "creators": creators,
                "publisher": "Zenodo",
                "resource_type": {"id": "software"},
                "rights": [{"id": (cff.get("license") or "mit").lower()}],
            }
        }
        (RELEASE / "zenodo_metadata.json").write_text(json.dumps(zenodo, indent=2), encoding="utf-8")
    except Exception as e:
        print(f"  (warn) could not parse CITATION.cff: {e}")

    print(f"  wrote {len(manifest)} file hashes to release/manifest.json")
    print("  release/ ready: manifest.json, requirements.txt, zenodo_metadata.json")


def upload_zenodo():
    import urllib.request, urllib.error
    token = os.environ.get("ZENODO_TOKEN")
    if not token:
        print("  (skip upload) set ZENODO_TOKEN to auto-upload via the Zenodo API")
        return
    meta = json.loads((RELEASE / "zenodo_metadata.json").read_text(encoding="utf-8"))
    api = "https://zenodo.org/api/records"

    def call(method, url, data=None, ctype="application/json"):
        r = urllib.request.Request(url, method=method, data=data)
        r.add_header("Authorization", f"Bearer {token}")
        if ctype:
            r.add_header("Content-Type", ctype)
        try:
            with urllib.request.urlopen(r) as resp:
                body = resp.read()
                return json.loads(body.decode()) if body else {}
        except urllib.error.HTTPError as e:
            print(f"  !! Zenodo HTTP {e.code}: {e.read().decode()[:400]}")
            sys.exit(1)

    print("  [upload] creating draft record ...")
    dep = call("POST", api, data=json.dumps(meta).encode())
    dep_id = dep["id"]
    print(f"  [upload] record id = {dep_id}")

    files = [f for f in sorted(RELEASE.iterdir()) if f.is_file() and f.name != "zenodo_metadata.json"]
    call("POST", f"{api}/{dep_id}/draft/files", data=json.dumps([{"key": f.name} for f in files]).encode())
    for f in files:
        print(f"  [upload] uploading {f.name} ...")
        call("PUT", f"{api}/{dep_id}/draft/files/{f.name}/content", data=f.read_bytes(), ctype="application/octet-stream")

    print("  [upload] publishing ...")
    pub = call("POST", f"{api}/{dep_id}/draft/actions/publish")
    doi = pub.get("doi") or (pub.get("metadata", {}) or {}).get("doi", "?")
    print(f"\n  DONE. DOI = {doi}")
    print(f"  -> replace '10.5281/zenodo.XXXXXXX' with {doi} in the manuscript + cover letter")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--upload", action="store_true", help="attempt Zenodo upload (needs ZENODO_TOKEN)")
    args = ap.parse_args()
    build_manuscript()
    snapshot()
    if args.upload:
        upload_zenodo()
    print("\nDone. Manuscript rebuilt, audit clean, release/ packaged.")
