#!/usr/bin/env python3
"""Prepare BigMHC input from type_ab_dataset.json.

Extracts 14,488 positives (5 alleles) into per-allele peptide CSVs
for BigMHC predict.py, plus a meta file for sensitivity computation.
"""
import json
import csv
from pathlib import Path

BASE = Path(__file__).resolve().parent.parent
DATA = BASE / "data" / "type_ab_dataset.json"
OUT = BASE / "bigmhc_input"


def main():
    d = json.load(open(DATA, encoding="utf-8"))
    OUT.mkdir(exist_ok=True)
    meta = []
    for allele, a in d["alleles"].items():
        # BigMHC fuzzy-matches allele names; strip to A0201 form for filename safety
        allele_short = allele.replace("HLA-", "").replace("*", "").replace(":", "")
        rows = [[e["peptide"]] for e in a["entries"]]
        with open(OUT / f"{allele_short}.csv", "w", newline="") as f:
            csv.writer(f).writerows(rows)
        for e in a["entries"]:
            meta.append([e["peptide"], allele, e["type"]])
        print(f"{allele}: {len(rows)} peptides -> {allele_short}.csv")
    with open(OUT / "meta.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["peptide", "allele", "type"])
        w.writerows(meta)
    print(f"meta.csv: {len(meta)} rows")


if __name__ == "__main__":
    main()
