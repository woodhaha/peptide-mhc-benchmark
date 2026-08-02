#!/usr/bin/env python3
"""Falsification Experiment C (Judge 3): is "non-canonical collapse" real under
a proper null and broader motif?

Judge 3 argued: (a) P9=Y is a legitimate A*03:01/A3-supertype anchor, so the
"non-canonical" positives are partly motif-artifact; (b) uniform-random negatives
let the model win on protein-likeness, not MHC signal.

Fix: (a) broaden the motif dicts to literature-validated anchor sets, re-run the
per-allele vs pan comparison so "canonical" means the broader motif and the
remaining "non-canonical" is genuinely outside known anchors; (b) use
proteome-frequency negatives in the non-canonical evaluator.

This script only REPORTS the impact of motif broadening on the canonical/
non-canonical partition (counts), leaving the numbers for the paper discussion.

Usage: python src/falsify_noncannonical_null.py   (base python; no TF)
"""
import csv, json, urllib.request, urllib.parse, collections

# Literature-validated anchor sets (broader than the 2-residue simplification)
# A3-supertype (A*03/11/31/33/68): P2 aliphatic, P9 basic OR Y/H tolerated.
BROAD_MOTIFS = {
    "HLA-A*02:01": {2: "LMIVATQ", 9: "LVIM"},      # Q2 tolerated in some A2 ligands
    "HLA-A*01:01": {3: "DE", 9: "YF"},
    "HLA-A*03:01": {2: "LVIM", 9: "KRYH"},          # A3-supertype: basic + Y/H
    "HLA-B*07:02": {2: "P", 9: "LIMFV"},            # V/F tolerated at P9
    "HLA-B*44:03": {2: "ED", 9: "FWY"},             # Y tolerated
}
NARROW_MOTIFS = {
    "HLA-A*02:01": {2: "LMIVAT", 9: "LVIM"},
    "HLA-A*01:01": {3: "DE", 9: "YF"},
    "HLA-A*03:01": {2: "LVIM", 9: "KR"},
    "HLA-B*07:02": {2: "P", 9: "LIMF"},
    "HLA-B*44:03": {2: "ED", 9: "FWYL"},
}
BASE = "https://query-api.iedb.org/mhc_export"

def fetch_positive(allele):
    q = urllib.parse.urlencode({
        "select": "epitope__name", "mhc_restriction__name": f"eq.{allele}",
        "epitope__object_type": "eq.Linear peptide",
        "assay__qualitative_measurement": "eq.Positive",
    }, safe="*")
    with urllib.request.urlopen(urllib.request.Request(BASE+"?"+q, headers={"Accept":"application/json"}), timeout=120) as r:
        rows = json.load(r)
    return [x["epitope__name"] for x in rows if x.get("epitope__name") and len(x["epitope__name"])==9]

def canon(pep, m):
    return all(pep[int(p)-1] in r for p, r in m.items())

def main():
    import json
    for allele in NARROW_MOTIFS:
        peps = fetch_positive(allele)
        n_narrow = sum(1 for p in peps if not canon(p, NARROW_MOTIFS[allele]) and max(p.count(a) for a in set(p)) < 4)
        n_broad = sum(1 for p in peps if not canon(p, BROAD_MOTIFS[allele]) and max(p.count(a) for a in set(p)) < 4)
        print(f"{allele}: positives={len(peps)}, 'non-canonical' narrow={n_narrow}, broad={n_broad}, "
              f"reclassified_by_broader_motif={n_narrow-n_broad}")

if __name__ == "__main__":
    main()
