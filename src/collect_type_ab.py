#!/usr/bin/env python3
"""Collect large-sample IEDB positives for 5 alleles; classify canonical / Type A / Type B.

Type classification by motif:
  - canonical: passes BROAD motif (validated supertype anchors)
  - Type A: fails BROAD but passes NARROW-inverse test = tolerated weak anchor
    (in broad but outside narrow); i.e. near-canonical, learnable
  - Type B: fails BOTH narrow and broad at anchor positions = energy-compensated /
    genuinely non-canonical candidate

Writes data/type_ab_dataset.json with all peptides + classification + PMID evidence.

Usage (base python): python src/collect_type_ab.py
"""
import json, urllib.request, urllib.parse, collections

ALLELES = ["HLA-A*02:01", "HLA-A*01:01", "HLA-A*03:01", "HLA-B*07:02", "HLA-B*44:03"]
SAFE = lambda a: a.replace("*", "_").replace(":", "")
BASE = "https://query-api.iedb.org/mhc_export"

# Narrow (2-residue simplification) and broad (validated supertype anchors)
NARROW = {
    "HLA-A*02:01": {2: "LMIVAT", 9: "LVIM"},
    "HLA-A*01:01": {3: "DE", 9: "YF"},
    "HLA-A*03:01": {2: "LVIM", 9: "KR"},
    "HLA-B*07:02": {2: "P", 9: "LIMF"},
    "HLA-B*44:03": {2: "ED", 9: "FWYL"},
}
BROAD = {
    "HLA-A*02:01": {2: "LMIVATQ", 9: "LVIM"},
    "HLA-A*01:01": {3: "DE", 9: "YF"},
    "HLA-A*03:01": {2: "LVIM", 9: "KRYH"},
    "HLA-B*07:02": {2: "P", 9: "LIMFV"},
    "HLA-B*44:03": {2: "ED", 9: "FWYL"},
}

def fetch_all(allele):
    """Fetch positive 9-mer linear epitopes for an allele (server default page)."""
    q = urllib.parse.urlencode({
        "select": "epitope__name,reference__pmid,assay__method",
        "mhc_restriction__name": f"eq.{allele}",
        "epitope__object_type": "eq.Linear peptide",
        "assay__qualitative_measurement": "eq.Positive",
    }, safe="*")
    req = urllib.request.Request(BASE+"?"+q, headers={"Accept": "application/json"})
    with urllib.request.urlopen(req, timeout=180) as r:
        return json.load(r)

def classify(pep, narrow, broad):
    """canonical = passes narrow OR broad; Type A = broad-only tolerated anchor;
    Type B = fails both narrow and broad at some anchor."""
    if len(pep) != 9: return None
    if all(pep[int(p)-1] in r for p, r in narrow.items()): return "canonical"
    if all(pep[int(p)-1] in r for p, r in broad.items()): return "canonical"
    broad_ok = {k: set(broad[k]) for k in broad}
    narrow_ok = {k: set(narrow[k]) for k in narrow}
    if any(pep[int(k)-1] in broad_ok[k] - narrow_ok.get(k, set()) for k in broad_ok):
        return "type_a"
    return "type_b"

def main():
    out = {"alleles": {}}
    for allele in ALLELES:
        rows = fetch_all(allele)
        by = collections.defaultdict(list)
        for x in rows:
            n = x.get("epitope__name")
            if not n or len(n) != 9: continue
            if max(n.count(a) for a in set(n)) >= 4: continue  # skip homopolymer-ish
            by[n].append(x)
        entries = []
        counts = {"canonical": 0, "type_a": 0, "type_b": 0}
        for pep, xs in by.items():
            c = classify(pep, NARROW[allele], BROAD[allele])
            if c is None: continue
            counts[c] += 1
            pmids = sorted({z.get("reference__pmid") for z in xs if z.get("reference__pmid")})
            entries.append({"peptide": pep, "type": c, "n_refs": len(pmids),
                            "pmids": ";".join(str(p) for p in pmids[:5])})
        out["alleles"][allele] = {"counts": counts, "total": len(entries), "entries": entries}
        print(f"{allele}: canonical={counts['canonical']}, TypeA={counts['type_a']}, TypeB={counts['type_b']}, total={len(entries)}")
    with open("data/type_ab_dataset.json", "w", encoding="utf-8") as f:
        json.dump(out, f, indent=1)
    print("saved data/type_ab_dataset.json")

if __name__ == "__main__":
    main()
