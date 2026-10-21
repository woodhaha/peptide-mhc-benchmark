# Supplementary Note S — WITHDRAWN (2026-10-21)

This note previously reported a molecular-dynamics analysis of the "K-E63 salt bridge".
It has been **withdrawn**: the 100-ns trajectories in `md results/salt_bridge/` simulated
the 1DUZ template peptide `LLFGYPVYV` (P2 = Leu), not the K-at-P2 mutant. The
Glu63-to-P2-Leu distances were mislabelled as "K-E63" salt-bridge contacts.

- The K-at-P2 trajectory (`docking/md_kras_k_at_p2_traj.dcd`) is ~1.2 ns (510–1710 ps)
  with uncontrolled temperature (344–509 K) — unequilibrated and unusable.
- `docking/1DUZ_K2_mutant.pdb` has a lysine without a side chain (only N/CA/CB/O), so no
  Nζ–OE2 distance exists.

The manuscript §3.5 therefore presents the salt bridge only as a **qualitative, untested
structural conjecture** from comparative B-pocket geometry and qualitative electrostatics.
See `docking/RESULTS.md` for the full audit trail.
