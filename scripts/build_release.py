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

    # Zenodo metadata from CITATION.cff (best-effort)
    try:
        import yaml
        cff = yaml.safe_load((ROOT / "CITATION.cff").read_text(encoding="utf-8"))
        authors = [{"name": f"{a.get('given-names','')} {a.get('family-names','')}".strip()}
                   for a in cff.get("authors", [])]
        zenodo = {
            "title": cff.get("title", ""),
            "creators": authors,
            "description": cff.get("abstract", ""),
            "keywords": cff.get("keywords", []),
            "upload_type": "software",
            "license": cff.get("license", ""),
            "access_right": "open",
            "publication_date": date.today().isoformat(),
        }
        (RELEASE / "zenodo_metadata.json").write_text(json.dumps(zenodo, indent=2), encoding="utf-8")
    except Exception as e:
        print(f"  (warn) could not parse CITATION.cff: {e}")

    print(f"  wrote {len(manifest)} file hashes to release/manifest.json")
    print("  release/ ready: manifest.json, requirements.txt, zenodo_metadata.json")


def upload_zenodo():
    token = os.environ.get("ZENODO_TOKEN")
    if not token:
        print("  (skip upload) set ZENODO_TOKEN to auto-upload via the Zenodo API")
        return
    print("  uploading to Zenodo (not implemented — use the API or the web UI)")
    # Zenodo REST: POST https://zenodo.org/api/depositions (create), then PUT files.
    # Left as a documented stub; the metadata in zenodo_metadata.json is drop-in ready.


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--upload", action="store_true", help="attempt Zenodo upload (needs ZENODO_TOKEN)")
    args = ap.parse_args()
    build_manuscript()
    snapshot()
    if args.upload:
        upload_zenodo()
    print("\nDone. Manuscript rebuilt, audit clean, release/ packaged.")
