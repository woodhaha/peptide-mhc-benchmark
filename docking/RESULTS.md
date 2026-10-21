# K-E63 Salt Bridge Structural Analysis — Results Summary (REVISED)

> Updated: 2026-10-21 · The 2026-10-16 version incorrectly labelled an MD analysis of the
> template peptide (LLFGYPVYV, P2=Leu) as a K-E63 salt-bridge analysis. The MD-derived
> numbers (mean 7.4 Å, <3 Å ~14%, MM-GBSA −19.6±46.3) measured Glu63-to-P2-Leu contacts
> on the wrong peptide and are **withdrawn**. See `docking/sb_contacts.py`
> (docstring: "E63-P2Leu interface").

## Corrected conclusion

A charged P2-Lys in the KRAS G12V neoepitope `YKLVVVGAV` could reach B-pocket Glu63 and
form a salt bridge, but this is a **qualitative, untested structural conjecture**. No
molecular-dynamics quantification on the correct peptide exists: the 100-ns trajectories
simulated the template peptide (P2=Leu), and the K-at-P2 trajectory is ~1.2 ns and
unequilibrated. The hypothesis requires a properly-built K-at-P2 model, a solvated
free-energy calculation, and experimental validation.

## What was done (reliable)

| Method | Result | Interpretation |
|--------|--------|----------------|
| Comparative B-pocket geometry (8 HLA-A\*02:01 structures) | Glu63 at 4.5–4.9 Å from P2 Cα | Conserved acidic B-pocket |
| Energy minimisation (OpenMM, amber14) | P2 Cα–Glu63 OE2 = 3.7 Å | Within Lys side-chain reach (7.6 Å) |
| Coulombic estimate (un-solvated) | ~−30 kcal/mol | **Overstates** — desolvation not included |

## Withdrawn (incorrect peptide / unavailable)

- ~~Molecular dynamics (100 ns, 3 replicas): "K–E63 mean 7.4 Å; <3 Å ~14%"~~ — actually
  Glu63–P2-Leu on the template peptide `LLFGYPVYV` (see `docking/sb_contacts.py`,
  docstring "E63-P2Leu interface").
- ~~MM-GBSA (final 50 ns): ΔG = −19.6 ± 46.3 kcal/mol~~ — same wrong-peptide trajectories.
- The K-at-P2 trajectory (`md_kras_k_at_p2_traj.dcd`) is ~1.2 ns (510–1710 ps) with
  uncontrolled temperature (344–509 K) and the `1DUZ_K2_mutant.pdb` lysine lacks a side
  chain (only N/CA/CB/O), so no Nζ–OE2 distance exists.

## Status

The structural argument in manuscript §3.5 is restricted to comparative B-pocket geometry
and qualitative electrostatics; no MD quantification is claimed.
