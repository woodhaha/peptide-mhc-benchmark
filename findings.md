# Findings — Multi-Allele Extension

## Current Understanding (2026-08-02)

The single-allele manuscript (A*02:01, submitted to Bioinformatics) established:
**label quality > architecture** for peptide-MHC binding prediction. The multi-allele
extension asks whether the Deep FFN + BLOSUM62 recipe generalizes across HLA alleles.

**Key result so far:** The Deep FFN generalizes cleanly to other alleles when (a) the
evaluation benchmark is properly curated and (b) training uses fixed deterministic
splits with mixed biased+random data. 3/3 previously-"failed" alleles now pass.

## Patterns and Insights

1. **"Model failures" were benchmark artifacts.** A*03:01's benchmark was
   contaminated with A*02:01 epitopes; B*44:03's collapsed to 5 near-duplicate
   peptides (F1 math cap ~0.63). Neither reflected model capability. *Lesson: verify
   benchmark curation before trusting a negative result.*
2. **Training variance is a reproducibility trap.** Unseeded random data generation
   produced 22-89% sensitivity swings across runs for A*01:01. Seed-42 fixed splits
   eliminated it. *Lesson: always save training data to disk with a fixed seed.*
3. **Mixed training data matters.** Pure motif-biased data causes catastrophic
   overfitting (from prior run-001). 50/50 biased+random gives the model both the
   motif signal and the ability to discriminate.

## Lessons and Constraints

- BLOSUM62 canonical-motif epitopes are the reliable eval signal; non-canonical
  binders are a separate (harder) problem already documented for A*02:01.
- mhcflurry's Class1AffinityPredictor loads in base anaconda (torch OK) but not in
  the r-tensorflow venv (torch c10.dll broken). Split workflow: label in base, train
  in venv.
- Old `.h5`/`.keras` models saved under keras 2.x have a `quantization_config` kwarg
  that keras 3.12 rejects — retrain rather than patch-load.

## Update (2026-08-02, mission complete)

**5/5 alleles PASS** with the Deep FFN + BLOSUM62 recipe + curated IEDB benchmarks:
sens 100%, spec 97-99%, F1 0.94-0.98 for all. The prior run-001 "failures" were
benchmark-curation and training-variance artifacts.

**Non-canonical robustness confirmed (per-allele models):** models retain 0.88-1.0
sensitivity on IEDB positives that violate canonical anchor motifs — they do not
simply memorize anchors. Precision 0.81-0.96 (lowest for A*02:01, its established
hard case). This directly supports the manuscript's claim that non-canonical binders
(K-E63 salt bridge) are the genuine hard case, not a failure of cross-allele
generalization.

**Step B — pan-specific model loses non-canonical robustness:** A single allele-aware
model (BLOSUM62 + 5-dim one-hot) passes all 5 canonical benchmarks but non-canonical
sensitivity collapses — worst A*03:01 (1.0 → 0.48, 14/25 missed vs 0/25 for the
per-allele model). Missed epitopes are P9 = H/Y/R variants. The allele one-hot is
weak conditioning; the shared pathway compromises to an "average MHC-I" motif.
**Design guidance: per-allele models for neoantigen/non-canonical discovery; pan
model only for canonical screening at scale.**

**Negative result — PSSM conditioning makes it worse:** Replacing the one-hot with a
20-dim anchor-preference PSSM drops non-canonical sens in 4/5 alleles (A*03:01 0.48 →
0.36). The PSSM encodes the canonical-motif bias the shared model over-relies on.
**Mechanism settled: the limitation is shared parameters, not conditioning strength.**
Any pan architecture learns the average MHC-I motif; per-allele models are required
for non-canonical capacity.

## CORRECTION (2026-08-02, three-judge review)

The 3-part story above was largely FALSIFIED by independent adversarial review +
falsification experiments:

1. **Canonical "100% / 97-99% / F1 0.94-0.98"** — n=25 Wilson CI floor is 86.7%; a
   zero-parameter anchor-checker (bind iff P2/P9 match motif) achieves sens 1.0/F1
   0.85-1.0 on all 5 alleles. The canonical benchmark is a motif tautology, not
   evidence of capability. Only A*02:01 with protein-background negatives shows the
   model beating the rule (F1 0.943 vs 0.847).
2. **"Per-allele > pan for non-canonical"** — the A*03:01 collapse (1.0→0.48) was a
   10.6× data-size artifact. Pan trained on full data recovers 0.48→1.0. The
   comparison is UNTESTED, not resolved.
3. **"Shared params limit pan" / "PSSM worse"** — falsified; it was data imbalance.
4. **"Non-canonical" partition** — partly a motif-dict artifact: A*03:01 P9=Y/H is a
   legitimate A3-supertype anchor (640 positives reclassified); B*07:02 P9=V/F (299).
   True non-canonical (Type B, e.g. A*02:01 P2-K salt bridge) is different from
   dict-coverage bug (Type A).

**Surviving kernel (defensible):** with full data, a pan model reaches parity with
per-allele models on non-canonical binders; per-allele superiority is not supported.
Prior run-001 "failures" were curation/variance artifacts (still valid).

**Before publishable:** protein-background negatives, motifs from independent sources,
multi-seed, external labels to break MHCflurry circularity, report CIs. The manuscript's
K-E63 salt-bridge claim is NOT tested by this work.

## Resolution (2026-08-02) — the honest, defensible positive result

After correction + rebalancing, the multi-allele thread converges on a strong claim:

**A single rebalanced pan-specific Deep FFN (B*07:02 upsampled 4x) reaches 0.84-0.96
non-canonical sensitivity across all 5 alleles — parity or better than per-allele
models.** The B*07:02 exception (0.60) was data scarcity, not an architectural limit.

The full arc:
1. Naive story ("5/5 pass, per-allele>pan, PSSM-dilution") — FALSIFIED by 3-judge
   review + falsification experiments (canonical = motif tautology; pan collapse =
   data artifact; non-canonical partition partly motif-dict artifact).
2. Corrected story: per-allele and pan-full parity in 3/5, B*07:02 over-anchoring.
3. Resolution: rebalancing fixes B*07:02; rebalanced pan beats or matches per-allele
   on all 5 alleles.

**Key lesson:** the first pass over-claimed. Independent adversarial review +
falsification + honest re-measurement produced a narrower but far stronger claim:
one model, properly fed, generalizes across 5 alleles on canonical AND non-canonical
binders.

## Seed-robustness confirmed (2026-08-02)

5-seed rebalanced pan (0,1,42,123,2024), non-canonical sens mean±SD:
A*02:01 0.98±0.02, A*01:01 0.92±0.04, A*03:01 0.98±0.02, B*07:02 0.92±0.05,
B*44:03 1.00±0.00. Test acc 0.952-0.956 (SD ~0.002). The single-seed criticism is
closed — the rebalanced-pan-generalizes claim is robust to seed.

## Novel resolution: Type A vs Type B non-canonical binders (2026-08-02)

The multi-allele work's real value is NOT the pan-vs-per-allele comparison (that's
field-consensus replication). It's the **Type A / Type B non-canonical distinction**:

- **Type A (tolerated weak anchors):** P9-Y/H/V/F variants. A well-trained model
  (rebalanced pan, 5 alleles) catches these at 0.92-0.98 — learnable, near-canonical.
- **Type B (energy-compensated):** e.g. the manuscript's A*02:01 P2-K (KRAS G12V
  salt bridge), P9-K/F/C. NOT rescued by 2.6x more data (13k model scores them
  0.000-0.254). Genuinely prediction-resistant.

**Manuscript implication:** the submitted manuscript's claim that non-canonical
binders are "systematically missed" is CONFIRMED, but sharpened to Type B only.
Ordinary non-canonical binders are actually well-handled by well-trained sequence
models. This is novel external validation FOR the manuscript, not a contradiction.

## Large-sample Type A/B test — REVISED (2026-08-02, n=14,488)

Collected ALL IEDB positive 9-mers for 5 alleles (14,488: canonical 11,225,
Type A 151, Type B 3,112), scored with both predictors (Deep FFN + MHCflurry).

| Type | n | Deep FFN sens | MHCflurry sens |
|------|--:|:-------------:|:--------------:|
| canonical | 11,225 | 0.979 [0.976,0.982] | 0.986 [0.984,0.988] |
| type_a | 151 | 0.477 [0.399,0.556] | 0.669 [0.590,0.739] |
| type_b | 3,112 | 0.607 [0.590,0.624] | 0.668 [0.652,0.685] |

**Revisions to the small-n story (CRITICAL):**
1. Non-canonical penalty is real but MODERATE (~0.98 -> ~0.6), not catastrophic.
2. Type A is NOT easy — DFF catches only 0.477, worse than Type B. The n=25 "Type A
   at 0.92-0.98" was cherry-picked top-evidence epitopes + pan model; at scale with
   narrow-biased training, broad-only tolerated anchors (P9-Y/H) are out-of-distribution
   and score low.
3. Type B is heterogeneous. The manuscript's 3 FNs + KRAS P2-K are the HARD TAIL that
   even MHCflurry misses (3.8-13.8 uM) — genuinely prediction-resistant, cross-predictor.

**Robust, manuscript-relevant facts:** canonical binders caught ~0.98 by both predictors;
motif-failing peptides caught ~0.6-0.67 (moderate penalty); the extreme non-canonical
tail (energy-compensated, e.g. P2-K salt bridge) is prediction-resistant for both.
Manuscript should frame its claim as "extreme non-canonical tail prediction-resistant,"
not "non-canonical binders systematically missed."

## MUNIS full-population result — structure does NOT rescue non-canonical binders (2026-08-02)

Scored all 3,112 Type B + 3,112 canonical with MUNIS (AlphaFold-Multimer pMHC structure
predictor) on AutoDL RTX 6000D:

| Predictor | Type B sens (pooled) |
|-----------|:--------------------:|
| Deep FFN (sequence) | 0.607 |
| MHCflurry (sequence) | 0.668 |
| **MUNIS (structure)** | **0.647** |

**The three-tier pipeline claim is NOT supported.** Structure-based prediction catches
Type B non-canonical binders at the same rate as sequence models (~0.65). The n=1
QYDPVAALF rescue was not representative.

**Honest scientific value:** this directly challenges the field's assumption that
"structure will solve what sequence misses." Non-canonical binding (Type B, incl. the
P2-K salt-bridge class) is genuinely hard for ALL predictor classes tested. The
manuscript should re-frame its three-tier pipeline as a hypothesis needing experimental
validation.

## Open Questions

- Does a broader-motif-trained model (include Type A residues in training data) recover
  the type_a 0.48?
- Would a *fine-tuned-on-pMHC* structure model (not generic AF2-multimer) do better?
  (MUNIS is trained on pMHC, so this is a probe of current limits.)
- External-label validation (unseen by MHCflurry) still the general publishability gate.
