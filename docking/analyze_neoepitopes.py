#!/usr/bin/env python
"""Neoepitope anchor analysis for HLA-A*02:01 docking"""

peptides = [
    ("p53 R248W neo",  "MNWRPILTI", "M", "I", True,  True,  "+0.407", "CREATED (NB->WB)"),
    ("p53 R248W wt",   "MNRRPILTI", "M", "I", True,  True,  "WT",      "WT non-binder"),
    ("KRAS G12V neo",  "YKLVVVGAV", "K", "V", False, True,  "+0.475", "CREATED (NB->WB)"),
    ("KRAS G12V wt",   "YKLVVVGAG", "K", "G", False, False, "WT",      "WT non-binder"),
    ("KRAS G12V enh",  "LVVVGAVGV", "V", "V", True,  True,  "+0.307", "ENHANCED (WB->SB)"),
    ("KRAS G12V enh_wt","LVVVGAGGV","V", "V", True,  True,  "WT",      "WT weak binder"),
    ("KRAS G12C enh",  "LVVVGACGV", "V", "V", True,  True,  "+0.307", "ENHANCED (WB->SB)"),
]

print("=" * 80)
print("  NEOEPITOPE DOCKING ANALYSIS SUMMARY")
print("=" * 80)
print()
print("Docking template: 1DUZ (HLA-A*02:01 + LLFGYPVYV, 1.80 A)")
print("Anchor motif: P2 = B-pocket (L/M/I/V), P9 = F-pocket (V/L)")
print()

header = f"{'Peptide':<22s} {'P2':<4s} {'P9':<4s} {'P2ok':<6s} {'P9ok':<6s} {'dScore':>7s} {'Effect':<20s}"
print(header)
print("-" * len(header))
for name, seq, p2, p9, p2ok, p9ok, dscore, effect in peptides:
    print(f"{name:<22s} {p2:<4s} {p9:<4s} {str(p2ok):<6s} {str(p9ok):<6s} {dscore:>7s} {effect:<20s}")

print()
print("=" * 80)
print("  STRUCTURAL INSIGHTS")
print("=" * 80)

print()
print("1. p53 R248W neoepitope (MNWRPILTI):")
print("   P2=M (canonical, fits deep B-pocket)")
print("   P9=I (acceptable at F-pocket)")
print("   Mutation R->W at P7: removes positive charge, adds bulky aromatic")
print("   Expected: better burial of hydrophobic W vs exposed charged R at P7")

print()
print("2. KRAS G12V neoepitope (YKLVVVGAV):")
print("   P2=K -- NON-CANONICAL: Lys+ in hydrophobic B-pocket")
print("   P9=V (canonical)")
print("   Mutation G->A at P8: subtle change from Gly to Ala")
print("   Key question: Can K+ at P2 be compensated by E63 in the groove?")
print("   If docking confirms binding, this is a true neoepitope from anchor adaptation")

print()
print("3. KRAS G12V enhanced (LVVVGAVGV):")
print("   P2=V, P9=V: BOTH fully canonical -- BEST anchor profile")
print("   Mutation G->A at P6: flexible Gly to hydrophobic Ala")
print("   Highest chance of strong binding among all candidates")

print()
print("=" * 80)
print("  B-POCKET CONTACTS (from 1DUZ template)")
print("=" * 80)
print("  P2 (LEU) burial: GLU63(2.9A), LYS66(3.0A), TYR99(3.4A),")
print("                   MET45(3.6A), TYR7(3.6A), VAL67(3.6A)")
print()
print("  Interpretation for KRAS K-at-P2:")
print("  - GLU63 is 2.9A from P2 -- could form salt bridge with K!")
print("  - This would be a novel anchor mode: K-E63 electrostatic compensation")
print("  - If true, G12V neoepitope represents a NON-CANONICAL binder")

print()
print("=" * 80)
print("  F-POCKET CONTACTS (from 1DUZ template)")
print("=" * 80)
print("  P9 (VAL) burial: THR143(2.7A), LYS146(2.7A), TYR84(2.7A),")
print("                   ASP77(3.0A), THR80(3.5A)")

print()
print("=" * 80)
print("  DOCKING PRIORITY RANKING")
print("=" * 80)
print()
print("P1: KRAS G12V neo (YKLVVVGAV) vs WT (YKLVVVGAG)")
print("    Tests non-canonical K-at-P2 binding. Highest clinical relevance")
print()
print("P2: p53 R248W neo (MNWRPILTI) vs WT (MNRRPILTI)")
print("    Tests R->W conformational switch at P7. Most frequent p53 mutant")
print()
print("P3: KRAS G12V enh (LVVVGAVGV) vs WT (LVVVGAGGV)")
print("    Best anchor profile. G->A subtle change at P6")
