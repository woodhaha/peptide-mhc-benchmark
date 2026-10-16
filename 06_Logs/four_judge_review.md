# Four-Judge Review — 3 Rounds (2026-10-16)

Panel review of `06_Bioinformatics_Submission/manuscript_bioinformatics_v2.md`,
revise-per-feedback looped 3×. (Precedent: `three_judge_review.md`, 2026-08-02.)

## The four tribunals
- **Judge 1 — Methodology & Statistics**: benchmark design, metrics, confounding, seeds, CIs.
- **Judge 2 — Immunoinformatics**: MHC-I biology, anchors, IEDB, label circularity, predictor comparison.
- **Judge 3 — Structural biophysics**: docking, salt bridge, MD, electrostatics.
- **Judge 4 — Editor & framing**: novelty, claims-vs-evidence, redundancy, title/abstract, consistency.

---

## Round 1 — findings → revisions

**Judge 1 (methodology)**
1. PSSM label source (82.6%) compared to MHCflurry on a *separate 15,000-peptide pool, single seed* — confounded, yet Table 2 presents them side-by-side. → §3.2 rewritten to lead with label-source result and state the confound inline.
2. B\*07:02 non-canonical sens inconsistent: §3.6 "0.60" (5-allele pan) vs Table 6 "0.64" (15-allele pan16). → Clarified both panels in §3.6.
3. Post-hoc domain filter inflating specificity to 100%. → Softened: 75.0% reported as primary, unfiltered metric.

**Judge 2 (immunoinformatics)**
4. "Five false positives were homopolymers … (correct binding predictions …)" — internally contradictory (homopolymers do not bind). → Rewrote as biologically implausible false positives.
5. "NetMHCpan-4.1 achieved only 83.3% … consistent with label circularity" — over-reach (NetMHCpan trains on different MS-elution labels). → Reframed as indicative; pointed to 62.24% label agreement as the direct circularity evidence.

**Judge 3 (structural)**
6. Salt bridge reported via Cα distances; the actual Nζ–OE2 distance missing. → Stated explicitly that Nζ–OE2 was not computed; Cα distances are a coarse proxy.

**Judge 4 (editor)**
7. **Author Contributions section stale** — did not include Yunshi Huang (new first/co-first author) or Xiangyang Xue (co-corresponding), and kept the pre-reorder roles. → Rewrote for all 7 authors.
8. §3.2 titled "Label-source controls" but body re-stated the architecture ranking. → Rewrote to focus on label source.
9. README "IEDB sensitivity 93.9% (46/49)" stale vs manuscript "95.8% (46/48)" (post ILRGSVAHK exclusion). → Fixed README.
10. Uncited epidemiological claim ("ranks second in colorectal, oesophageal…"). → Softened to IARC TP53 database citation.

---

## Round 2 — findings → revisions

1. §3.7 cites "broad-motif retraining" results (A\*03:01 flat 0.634; B\*07:02 Type B 0.680→0.341) but Methods §2.6 never declares the experiment. → Added to §2.6.
2. Repo contains 100 ns MD trajectories + MM-GBSA, but manuscript reports only energy minimisation — undeclared analysis gap. → §2.4 now states MD is an exploratory companion analysis reported separately.
3. §2.3 meta-commentary "satisfying the journal's independent-test-set requirement". → Removed (editorializing).

---

## Round 3 — convergence check

- Grep confirmed zero stale numbers in manuscript: `46/49`, `93.9`, `5,088`, `65.8`, `94.8`, `MHCflurry 2.2.0`, `Yongyu Bai`, "Deep FFN best" all absent.
- Author list (L3), corresponding-author line (L14), and Author Contributions (L230) all consistent for 7 authors.
- Key numbers verified consistent: 5,190 / 46/48 / 95.8% / 34.5% / 0.819 / B\*07:02 0.60 (5-allele) + 0.64 (pan16).

**All three rounds concluded.** Remaining Discussion-vs-Results redundancy (architecture ranking re-stated) is stylistic and accepted; no further substantive edits.
