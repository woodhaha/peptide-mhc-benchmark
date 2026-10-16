# K-E63 Salt Bridge Structural Analysis — Results Summary (CORRECTED)

> Updated: 2026-10-16 · Supersedes the 2026-06-24 version, which overstated the salt-bridge
> stabilisation. The corrected conclusion is **qualitative only**, consistent with
> manuscript §3.5 and `05_Submission/supplementary/MD_salt_bridge_supplementary.md`.

## Corrected conclusion

A charged P2-Lys in the KRAS G12V neoepitope `YKLVVVGAV` could reach B-pocket Glu63 and
form a salt bridge, but the evidence is **qualitative** and the interaction is **not a
stable, persistent contact** in molecular dynamics. The hypothesis requires a solvated
free-energy calculation and experimental validation.

## What was done

| Method | Result | Interpretation |
|--------|--------|----------------|
| Comparative B-pocket geometry (8 HLA-A\*02:01 structures) | Glu63 at 4.5–4.9 Å from P2 Cα | Conserved acidic B-pocket |
| Energy minimisation (OpenMM, amber14) | P2 Cα–Glu63 OE2 = 3.7 Å | Within Lys side-chain reach (7.6 Å) |
| Coulombic estimate (un-solvated) | ~−30 kcal/mol | **Overstates** — a buried salt bridge contributes only a few kcal/mol after desolvation |
| Molecular dynamics (100 ns, 3 replicas) | K–E63 mean 7.4 Å; < 3 Å in ~14% of frames | Salt bridge is **transient**, not persistent |
| MM-GBSA (final 50 ns) | ΔG = −19.6 ± 46.3 kcal/mol | Too noisy for a precise estimate |
| HDOCKlite docking (10 poses) | K–E63 = 15–47 Å in all poses | **Inconclusive** — peptide placed outside the groove; HDOCKlite is not suited to peptide-MHC docking |

## What was wrong in the prior version

- ~~"−29.6 kcal/mol, ~10× stronger than canonical P2 burial"~~ — an un-solvated Coulombic
  estimate that ignored desolvation penalties; removed.
- ~~"First documented case of a charged P2 anchor"~~ — overstatement; the configuration is
  hypothetical and unvalidated.

## Status

The structural argument is restricted to minimisation and qualitative electrostatics
(manuscript §3.5). The 100-ns MD is reported in Supplementary Note S. Neither the Coulombic
estimate nor the HDOCKlite docking supports a **quantitative** salt-bridge claim.
