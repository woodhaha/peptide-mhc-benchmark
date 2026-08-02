# Companion Note — First Draft Complete

**Date**: 2026-08-02 | **Stage**: Authoring (ml-paper-writing) | **Target**: Bioinformatics

## What happened this session

Research was CONCLUDED (10 validated findings). This session moved the project from
**research → authoring**: the ml-paper-writing skill was installed, and the full
first draft of the multi-allele extension companion note was written and verified.

## Deliverables

| Artifact | Path |
|----------|------|
| First draft (IMRaD, ~2,000 words) | `04_Manuscript/companion_multi_allele/companion_note_draft.md` |
| Figure 1 (4-panel findings) | `04_Manuscript/companion_multi_allele/fig1_findings_panel.png` |

## The one-sentence contribution

> **One rebalanced pan-specific Deep FFN replaces a per-allele panel (parity or
> better on 5 alleles, 0.819 pooled on 15), while non-canonical (Type B) binders
> plateau at ~0.65 across every predictor class including structure-based MUNIS —
> the energy-compensated tail is prediction-resistant to all current methods.**

## Three-part story (per ml-paper-writing narrative principle)

1. **Positive, replicable**: rebalanced pan = per-allele parity, 5-seed robust
   (0.92–1.00). Corrects the naive "per-allele wins" benchmark artifact.
2. **Negative, cross-predictor**: Type B plateau 0.607/0.668/0.647 across
   DFFN/MHCflurry/MUNIS — structure does **not** rescue what sequence misses.
3. **Caveat that strengthens the structural story**: P2-Lys/Arg/His all reach
   Glu63 (3.83 Å); the hard tail is prediction-resistant, yet structure predicts
   it binds — arguing for structure as a *complement*, pending experimental tests.

## Number verification (all 15 claims)

| Claim | Source | Status |
|-------|--------|:------:|
| counts 11,225/151/3,112/14,488 | type_ab_dataset.json | ✅ exact |
| DFFN 0.979/0.477/0.607 | summarize_type_ab.py re-run | ✅ exact |
| MHC 0.986/0.669/0.668 | summarize_type_ab.py re-run | ✅ exact |
| MUNIS 0.950/0.647 | munis csv recompute | ✅ exact |
| 5-seed means ± SD | multiseed_results.json | ✅ exact |
| pan16 0.819 (307/375) | three_judge_review | ✅ |
| temporal 0.824/0.936/0.973 | plot_findings_panel.py | ✅ |
| Glu63 3.83 Å; K/R/H reach | p2_virtual_scan.json | ✅ |

Citations reused from the verified manuscript reference list — **no hallucinated
references**.

## Next steps (human input welcome)

1. **Review the draft** — `04_Manuscript/companion_multi_allele/companion_note_draft.md`
2. **Framing check** — emphasized "rebalancing removes the per-allele artifact" as
   the positive; adjust if you prefer different emphasis
3. **Figures** — 4-panel panel exists; may want separate publication-quality figures
   per finding
4. **Reference [3] MUNIS** (Nat Biotechnol 2025) and **[4] Preibisch** (bioRxiv 2026
   preprint) are placeholders inherited from the main manuscript — verify before
   submission
5. Journal = **Bioinformatics** per your saved recommendation

## Prior session's fix (for the record)

Submitted cover letter referenced a **private** repo — corrected to the public
`woodhaha/peptide-mhc-benchmark` and regenerated the DOCX (commit `f640339`).

---

## Update: draft revised after review (2026-08-02, commit ccb9133)

- **[10] Robinson/IMGT-HLA verified** via CrossRef — real reference, not hallucinated.
- Metadata fixed: 2,239 words / 1 table / **2 figures** (was stale ~5,800/5/4).
- `[companion submission]` self-refs → formal citation [15].
- "pre-registered thresholds" → honest "fixed thresholds" (post-hoc 5th-pct canonical).
- **New Figure 2**: rebalanced pan 5-seed mean ± SD — the core positive result.
  Per-allele bars deliberately excluded (apples-to-oranges with the n=25 benchmark).
- Figure 1 (findings panel) now cited in text (1A–1D) with standalone captions.
