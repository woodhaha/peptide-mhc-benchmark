#!/usr/bin/env python3
"""Verify each of the 49 benchmark 'A*02:01' epitopes against the live IEDB API.

Classifies each peptide as:
  A2-POS     -> has A*02:01/A2 Positive (or Positive-*) MHC record
  OTHER      -> has records but none A*02:01 positive
  NO-RECORD  -> zero IEDB records
Saves the full allele breakdown to iedb_epitope_verification.csv.
"""
import csv, json, time, urllib.request
from collections import Counter

BASE = "https://query-api.iedb.org/mhc_search"

def query(pep):
    url = f"{BASE}?linear_sequence=eq.{pep}&select=mhc_allele_name,qualitative_measure"
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
    with urllib.request.urlopen(req, timeout=30) as r:
        return json.loads(r.read())

def is_a2(a):
    return (a or "") in ("HLA-A*02:01", "HLA-A2", "HLA-A*0201")

def is_pos(q):
    return (q or "").startswith("Positive")

rows = list(csv.DictReader(open("02_Data/cleaned/iedb_benchmark_results_v2.csv", encoding="utf-8")))
pos = [r for r in rows if r["true_label"] == "POS"]

results = []
for r in pos:
    pep = r["peptide"]
    try:
        recs = query(pep)
    except Exception as e:
        results.append({"peptide": pep, "status": "ERROR", "n": 0, "a2_pos": 0, "alleles": str(e)})
        print(f"{pep:12} ERROR {e}")
        time.sleep(0.4)
        continue
    a2_pos = sum(1 for x in recs if is_a2(x.get("mhc_allele_name")) and is_pos(x.get("qualitative_measure")))
    alleles = Counter(x.get("mhc_allele_name", "?") for x in recs)
    status = "A2-POS" if a2_pos > 0 else ("OTHER" if recs else "NO-RECORD")
    results.append({"peptide": pep, "status": status, "n": len(recs), "a2_pos": a2_pos,
                    "alleles": "; ".join(f"{k}:{v}" for k, v in alleles.most_common(4))})
    print(f"{pep:12} {status:9} n={len(recs):3} a2pos={a2_pos:2}  {dict(alleles.most_common(3))}")
    time.sleep(0.4)

with open("02_Data/cleaned/iedb_epitope_verification.csv", "w", newline="", encoding="utf-8") as f:
    w = csv.DictWriter(f, fieldnames=["peptide", "status", "n", "a2_pos", "alleles"])
    w.writeheader()
    w.writerows(results)

from collections import Counter as C
st = C(r["status"] for r in results)
print("\n=== SUMMARY ===")
print(dict(st))
print("A2-verified positives:", st.get("A2-POS", 0), "/", len(results))
