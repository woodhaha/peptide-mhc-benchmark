#!/usr/bin/env python3
"""Convert author-year citations to PLOS numbered [n] citations, order of appearance.

Operates on manuscript_bioinformatics_v2.md. Prints the mapping and writes the
converted manuscript to a NEW file (does not overwrite the source).
"""
import re, sys

SRC = "06_Bioinformatics_Submission/manuscript_bioinformatics_v2.md"
DST = "06_Bioinformatics_Submission/manuscript_plos_numbered.md"

# (surname, year, second-author-surname-or-None) in reference-list order
REFS = [
    ("Albert", 2023, None), ("Backert", 2015, "Kohlbacher"), ("Bernett", 2024, None),
    ("Cheng", 2025, None), ("Collins", 2024, None), ("Falk", 1991, None),
    ("Fasoulis", 2024, None), ("Harndahl", 2012, None), ("Henikoff", 1992, "Henikoff"),
    ("Hermann", 2024, None), ("Jappe", 2020, None), ("Kapoor", 2023, "Narayanan"),
    ("Lin", 2023, None), ("Lundegaard", 2011, None), ("Malekzadeh", 2020, None),
    ("Mares", 2025, None), ("Nielsen", 2003, None), ("O'Donnell", 2020, None),
    ("Preibisch", 2026, None), ("Quandt", 2018, None), ("Rammensee", 1995, None),
    ("Reynisson", 2020, None), ("Rives", 2021, None), ("Ruppert", 1993, None),
    ("Sarkizova", 2020, None), ("Schuh", 2026, None), ("Tadros", 2025, None),
    ("Vita", 2019, None), ("Vo", 2025, None), ("Wohlwend", 2025, None),
    ("Wu", 2020, None), ("Xu", 2026, None),
]

def cit_re(surname, year, second):
    pat = re.escape(surname)
    if second:
        pat += r" and " + re.escape(second)
    else:
        pat += r" et al\."
    pat += r", " + str(year)
    return re.compile(pat)

CITS = [(i, cit_re(s, y, s2)) for i, (s, y, s2) in enumerate(REFS)]

text = open(SRC, encoding="utf-8").read()
body, refs_block = text.split("## References", 1)

order = []  # ref indices in order of first appearance

def process_group(inner):
    parts = [p.strip() for p in inner.split(";")]
    nums = []
    for p in parts:
        hit = None
        for ref_idx, rgx in CITS:
            if rgx.fullmatch(p):
                hit = ref_idx
                break
        if hit is None:
            return None  # not a pure-citation group
        if hit not in order:
            order.append(hit)
        nums.append(hit)
    nummap = {r: i + 1 for i, r in enumerate(order)}
    return "[" + ",".join(str(nummap[r]) for r in nums) + "]"

def repl(m):
    res = process_group(m.group(1))
    if res is None:
        return m.group(0)  # keep original parens for non-citations
    return res  # "[n1,n2]"

body2 = re.sub(r"\(([^()]*)\)", repl, body)

# reorder references: split the refs block into individual entries
ref_entries = [l for l in refs_block.strip().split("\n") if l.strip() and l.strip() != "## References"]
# ref_entries[0] is empty (leading newline) — filter empties already done
ref_entries = [l for l in ref_entries if l.strip()]
assert len(ref_entries) == 32, f"expected 32 refs, got {len(ref_entries)}"

reordered = [ref_entries[i] for i in order]
numbered = [f"{n}. {r}" for n, r in enumerate(reordered, 1)]

out = body2 + "## References\n\n" + "\n".join(numbered) + "\n"

open(DST, "w", encoding="utf-8").write(out)

# report
print("=== citation -> number mapping (order of appearance) ===")
for n, i in enumerate(order, 1):
    s, y, s2 = REFS[i]
    print(f"  [{n}] {s}{' and '+s2 if s2 else ' et al.'}, {y}")
print(f"\nTotal unique refs cited: {len(order)} / 32")
uncited = [i for i in range(32) if i not in order]
print(f"Uncited refs: {uncited}")
print(f"Wrote {DST}")
