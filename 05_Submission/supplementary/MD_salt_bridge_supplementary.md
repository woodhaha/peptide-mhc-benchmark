# Supplementary Note S — Molecular Dynamics of the K-E63 Salt Bridge

This note reports the molecular-dynamics (MD) analysis referenced in the main text
(§2.4), which is deliberately excluded from the primary structural argument and
presented here for completeness.

## S.1 Methods

The KRAS G12V neoepitope `YKLVVVGAV` was modelled into the HLA-A\*02:01 peptide-binding
groove using the 1DUZ template (template peptide `LLFGYPVYV` replaced in silico). The
complex was solvated and equilibrated, then simulated for 100 ns with OpenMM 8.5.2 and
the amber14 force field. Three independent replicas (seeds 42 and 99) were run to assess
convergence. MM-GBSA binding energies were computed over the final 50 ns (200 frames).

## S.2 Results

**System stability.** The complex remained stable over 100 ns: the MHC α-chain,
β2-microglobulin, and peptide backbone RMSDs converged to 0.15, 0.09, and 0.07 Å,
respectively (Fig. S-MD2). Whole-complex backbone RMSD reached 6.1 Å, reflecting
flexible regions outside the peptide-binding groove.

**The K-E63 salt bridge is transient, not persistent.** The P2-Lys–Glu63 minimum
distance fluctuated between 3.6 and 5.9 Å across the three replicas (mean ≈ 4.1 Å),
with transient close contacts (< 3 Å) in ~14% of frames in the primary trajectory
(Fig. S-MD3; Fig. S-MD1). The salt bridge therefore forms intermittently but does not
persist as a stable interaction over the simulated timescale.

**Binding energy.** MM-GBSA over the final 50 ns gave ΔG = −19.6 ± 46.3 kcal/mol
(mean ± SD, n = 200), dominated by gas-phase electrostatics (ΔE = −19.5 kcal/mol) with a
negligible solvation term. The large variance precludes a precise binding estimate
(Fig. S-MD4).

**Hydrogen bonding.** 3,574 unique hydrogen bonds were observed across the trajectory,
with 14 persistent pairs (occupancy > 30%).

## S.3 Interpretation

The MD is consistent with the conservative framing in the main text (§3.5): a charged
P2-Lys can transiently reach Glu63, but the K-E63 salt bridge is not a stable,
persistent interaction over the simulated timescale. These data support treating the
salt-bridge hypothesis as qualitative, pending a more rigorous free-energy calculation
(e.g., alchemical transformation or umbrella sampling) and experimental validation
(e.g., X-ray crystallography or cryo-EM of the YKLVVVGAV–HLA-A\*02:01 complex).

## S.4 Formal docking (HDOCKlite)

Formal peptide–MHC docking was attempted with HDOCKlite (rigid-body protein–protein
docking) for the KRAS G12V neoepitope and controls. All ten returned poses placed the
peptide outside the binding groove (P2-Lys–Glu63 ≥ 15 Å), reflecting the unsuitability of
a rigid-body docking tool for peptide–MHC docking; these poses were therefore not used for
the structural argument. The template-based geometry and minimisation (§S.2, main text
§3.5) remain the basis of the structural hypothesis.

## Figures

- **Fig. S-MD1** — 100-ns summary panel (system overview).
- **Fig. S-MD2** — Backbone RMSD over 100 ns (system stability).
- **Fig. S-MD3** — P2-Lys–Glu63 distance time series (transient salt-bridge contacts).
- **Fig. S-MD4** — MM-GBSA binding-energy distribution.
