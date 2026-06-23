# Decision Log — Run 001

**Mission:** Replace BLOSUM62 with PLM embeddings for peptide-MHC binding prediction
**Max runtime:** 4 hours
**Started:** 2026-06-19

---

## Iteration 1: ESM-2 t6 (8M) Mean Pooling ❌
- **Approach:** Extract 320-dim ESM-2 t6 embeddings via mean pooling, train Deep FFN
- **Result:** ❌ FAILED
- **Metrics:** Test accuracy 65.9% vs BLOSUM62 baseline 91.9%
- **Diagnosis:** ESM-2 t6 (8M params) is too small to capture peptide-MHC binding features. Mean pooling over 9 positions loses critical anchor residue information. BLOSUM62's position-specific substitution scores are more directly relevant for this task.
- **Decision:** Small PLM embeddings are not a drop-in replacement for BLOSUM62. Larger models (t12/t30) or position-aware encoding may help.

---

## Mission Status: FAILED (Inconclusive)

ESM-2 t6 (smallest variant) tested. Larger variants (t12 35M/480-dim, t30 150M/640-dim) not tested due to computational constraints. The hypothesis that PLM embeddings could replace BLOSUM62 remains open — a negative result with the smallest model doesn't rule out larger ones.

### Recommendation
- Try ESM-2 t12 (35M) or t30 (150M) for more expressive embeddings
- Consider per-position embeddings instead of mean pooling
- Consider fine-tuning ESM-2 on binding data rather than using frozen embeddings
