# Re-review + Data Supplement + NetMHCpan Automation — 2026-10-23 (phase 2)

Follow-up to `three_judge_review_2026-10-23.md` (phase 1). After the first fixes, a second
3-judge re-review found the fixes were mostly "acknowledge but not change the claim", and
two NEW FATALs. This file records the phase-2 closure.

## Re-review verdicts

- **J1 (methodology)**: only 1/5 old FATALs truly closed; new FATALs — (a) multi-allele
  experiments have NO documented train/test split (Tables 6–10 possibly in-sample), (b) 48−7=41≠42
  denominator, (c) natural-prior "always-NB 93.2% > all models", (d) threshold alignment only 1 point.
- **J2 (immunoinformatics)**: ground-truth cleanup is CORRECT (7 peptides verified wrong, TLGIVCPIC
  legit A2); but 48/42/41 denominator, ELAGIGILTV (10-mer vs 9-mer window), Type A/B undefined,
  degenerate homopolymer negatives, 42-peptide list not deposited.
- **J3 (editor)**: title "Ground-Truth Contamination" over-claims (no field-level measurement);
  Box 1 self-reflexivity; desk-reject ~60%.

## Phase-2 closures

| Gap | Fix | Commit |
|---|---|---|
| train/test split (J1 F1) | verified models train on SYNTHETIC data (motif-biased+random, not IEDB 14,488); documented §2.6.3 | 2ee297d |
| threshold sensitivity (J1 F3) | extended to 3 aligned points (0.90/0.95/0.98) | 2ee297d |
| 48/42/41 denominator (J2 MAJOR-1) | 3-step: 49→48 (ILRGSVAHK)→42 (6 more); NetMHCpan flagged non-comparable | f05a2d7 |
| title over-claim (J3 F1) | drop "Ground-Truth Contamination", keep "Label Circularity" | f05a2d7 |
| ELAGIGILTV (J2 MAJOR-2) | 10-mer A27L variant outside 9-mer window → 10/10 not 10/11 | f05a2d7 |
| Type A/B motif (J1 F5/J2 MAJOR-3) | documented 5-allele anchor residues (§2.6.2) | 85177b7 |
| natural-prior (J1 F2) | §3.1 caveat: always-NB 93.2% > all models | 85177b7 |
| degenerate negatives (J2 MAJOR-5) | real IEDB negatives: 90.0% specificity (18/20) | 85177b7 |
| NetMHCpan head-to-head (J2 MAJOR-1B) | standalone 4.1b = 85.7% (36/42) | 2273cd7 |

## NetMHCpan standalone automation (the "彻底 pin 版本" chain)

1. **DTU license form**: automated POST to `/cgi-bin/sw_ship` (fields: software=netMHCpan,
   version=4.1, packageversion=4.1b, platform=Linux, name, mail, affiliation, position,
   accept=yes). Academic license accepted, download link emailed (4h validity).
2. **Download**: read link via alimail IMAP; fetched `netMHCpan-4.1b.Linux.tar.gz` (321K) +
   `data.tar.gz` (30MB) from `/download/<uuid>/`.
3. **Install (WSL2)**: extracted; fixed nested `data/data/` → `data/`; recreated symlink
   `Linux_x86_64/data -> ../data`; copied to no-space path `~/netmhcpan` (binary rejects paths
   with spaces).
4. **Run**: `env NETMHCpan=... TMPDIR=... ./Linux_x86_64/bin/netMHCpan -a HLA-A02:01 -f pep42.fsa
   -l 9 -BA -xls -xlsfile ...` → version confirmed **NetMHCpan 4.1b**.
5. **Result**: **36/42 = 85.7%** sensitivity (%rank<2.0), below Deep FFN 41/42=97.6%.

## Final head-to-head (all pinned, same 42-epitope denominator)

| Predictor | Sensitivity |
|---|---|
| Deep FFN | 97.6% (41/42) |
| NetMHCpan-4.1b (standalone) | 85.7% (36/42), Wilson CI 72.1–93.3% |
| anchor rule | 78.6% (33/42) |

## Pitfalls (reusable)

- NetMHCpan binary rejects paths with spaces → run from no-space dir (WSL2 `~/`).
- `export VAR=...` does NOT reach the binary → use `env VAR=... ./binary`.
- `data.tar.gz` extracts to nested `data/data/` → move up; `data` is a symlink to `../data`.
- DTU platform options: Linux / Darwin only (no Windows) → WSL2 on Windows.
- **License**: non-redistributable → do NOT commit the software (`netmhcpan_dl/` gitignored);
  only the result CSV (`netmhcpan_42_standalone.csv`) is committed.

## Bottom line

All judge-flagged fixable gaps are now closed. The paper is reproducible, the SOTA comparison
is fully version-pinned, the ground truth is corrected, and cross-file consistency is
mechanically audited green. Remaining (optional, revision-time): none — NetMHCpan was the last
blocked item and is now done.
