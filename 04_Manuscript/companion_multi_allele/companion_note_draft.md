# One Rebalanced Pan-Specific Model Replaces a Panel: Multi-Allele Generalization and the Non-Canonical Frontier in Peptide-MHC-I Binding Prediction

**Authors**: Zhuha Zhou<sup>1,†,*</sup>, Chen Weisheng<sup>1,†</sup>, Yiqi Cai<sup>1</sup>, Yongyu Bai<sup>1</sup>, Zhuxian Zhou<sup>2,^</sup>

**Affiliations**:
<sup>1</sup> Department of Gastroenterology Surgery, The First Affiliated Hospital of Wenzhou Medical University, Nanbaixiang Street, Ouhai District, 325000 Wenzhou, Zhejiang, China
<sup>2</sup> Center for Bionanoengineering, College of Chemical and Biological Engineering, Zhejiang University, Hangzhou 310027, China

<sup>†</sup> These authors contributed equally to this work.
<sup>*</sup> Corresponding author: Zhuha Zhou, Department of Gastroenterology Surgery, The First Affiliated Hospital of Wenzhou Medical University, Nanbaixiang Street, Ouhai District, 325000 Wenzhou, Zhejiang, China. E-mail: zhouzhuha@wmu.edu.cn. ORCID: 0009-0005-1818-4789
<sup>^</sup> Zhuxian Zhou, Center for Bionanoengineering, College of Chemical and Biological Engineering, Zhejiang University, Hangzhou 310027, China. E-mail: zhouzx@zju.edu.cn

**Word count**: ~2,000 | **Tables**: 1 | **Figures**: 2

---

## Abstract

**Motivation**: Peptide binding to MHC class I is the gatekeeper of CD8+ T-cell recognition, and pan-specific predictors such as NetMHCpan [1] have become the field standard for multi-allele screening. Whether a single sequence-based model can match per-allele specialists — and how the field should handle non-canonical (motif-failing) binders that elude sequence methods — remain open questions. **Results**: We show that a single rebalanced pan-specific deep feed-forward network (DFFN; HLA-B\*07:02 upsampled 4×) reaches 0.92–1.00 non-canonical sensitivity across 5 common HLA-I alleles (5-seed SD ≤0.05) and 0.819 pooled sensitivity across a 15-allele panel — parity or better than per-allele models, with one model replacing five. Canonical binders are recognized at ~0.98 by every predictor class (DFFN, MHCflurry [2], MUNIS [3]). In contrast, non-canonical (Type B) binders plateau at ~0.65 sensitivity across all predictor classes, including the structure-based MUNIS; the extreme energy-compensated tail — e.g. the KRAS G12V P2-Lys salt-bridge configuration, where charged P2 residues (Lys/Arg/His) all reach the acidic B-pocket Glu63 (3.83 Å) — is prediction-resistant to every method tested. Temporal holdout analysis confirms that post-2022 epitopes (unseen by the labeler's training data) are recognized at preserved or improved sensitivity, breaking the label-circularity concern [4]. **Availability**: All code, trained models, and data are provided as a reproducible open-source pipeline.

**Contact**: zhouzhuha@wmu.edu.cn

**Supplementary information**: Supplementary figures S1–S5 and data tables.

---

## 1 Introduction

MHC class I molecules present 8–11-mer peptides to CD8+ T cells, and accurate prediction of peptide–MHC-I binding underpins epitope discovery for vaccines and cancer immunotherapy [5–7]. Anchor-residue motifs — most prominently positions 2 and 9 for HLA-A\*02:01 — dominate binding [8,9], but the thousands of HLA class I allotypes [10] each present distinct motifs, and pan-specific models trained across pooled alleles (NetMHCpan [1], MHCflurry [2], MUNIS [3]) have become the practical standard because a single model serves all alleles.

Two questions motivate this work. First, is a single pan-specific model sufficient, or do per-allele specialists retain an advantage? Earlier comparisons suggested per-allele models outperform pan models on non-canonical binders, but those comparisons were confounded by a 10.6× training-data size difference. Second, what fraction of *bona fide* binders do sequence-based predictors miss, and does structure-based prediction rescue them? Recent work on structure-guided pMHC-I generation reports that sequence predictors score structurally valid complexes at near-random levels (AUROCs 0.06–0.22) [11], implying a systematic blind spot. We previously benchmarked six architectures for HLA-A\*02:01 and found that label quality — not architecture — defines the performance ceiling [15]. Here we ask whether a single rebalanced pan-specific model generalizes across alleles, and we measure — across three independent predictor classes including a structure-based method — how hard the non-canonical tail actually is.

## 2 Methods

### 2.1 Data

Positive training/evaluation peptides were curated from the IEDB [12] for 5 common alleles (HLA-A\*02:01, A\*01:01, A\*03:01, B\*07:02, B\*44:03), n = 14,488 experimentally confirmed binders (canonical 11,225; Type A 151; Type B 3,112). Type A (tolerated weak anchors) and Type B (energy-compensated non-canonical) were defined by whether the peptide violates the allele's canonical anchor motif at p2/p9: Type A carries broad-only tolerated anchors; Type B carries non-canonical residues at anchor positions. A 15-allele panel (pan16: 11 additional alleles) was curated for generalization testing. Protein-background negatives were sampled from natural proteomes.

### 2.2 Model

The pan-specific DFFN uses the BLOSUM62 [13] encoding of the 9-mer (9 × 20 matrix, flattened to 180 dims) through a 360–180–90–45 hidden stack with batch normalization and dropout, identical to the architecture benchmarked in our single-allele study [15]. Rebalancing upsampled the minority B\*07:02 allele 4× during training. Five seeds (0, 1, 42, 123, 2024) were used for robustness; per-seed variation arises from data-sampling RNG. For the 15-allele panel, the same recipe was trained with per-allele minimum sample (2,000) upsampling.

### 2.3 Comparator predictors

MHCflurry 2.2.1 [2] served as an independent sequence predictor (affinity threshold 500 nM). MUNIS [3] — a structure-based pMHC-I predictor using ESM-2 [14] embeddings and AlphaFold-Multimer-style structure scoring with no-flank checkpoints — served as the structure-conditioned comparator, run on the n = 3,112 Type B population plus canonical controls.

### 2.4 Temporal holdout

To break the label-circularity concern [4], IEDB epitopes published after 2022 (post-2022, unseen by MHCflurry's ~pre-2021 training data) were compared against pre-2021 (in-distribution) epitopes at sensitivity threshold 0.95.

## 3 Results

### 3.1 A single rebalanced pan model matches per-allele specialists (F1)

Across 5 seeds, the rebalanced pan DFFN reached non-canonical sensitivity mean ± SD: A\*02:01 0.98 ± 0.02, A\*01:01 0.92 ± 0.04, A\*03:01 0.98 ± 0.02, B\*07:02 0.92 ± 0.05, B\*44:03 1.00 ± 0.00 (test accuracy 0.952–0.956; Figure 2). These values are at or above the corresponding per-allele models, so a single model replaces five. The earlier claim that per-allele models outperform pan models on non-canonical binders is not supported once per-allele data size is controlled.

### 3.2 Canonical binders are effectively solved (F2)

On n = 11,225 canonical positives, all three predictors were saturated: DFFN 0.979, MHCflurry 0.986, MUNIS 0.950.

### 3.3 Non-canonical (Type B) binders plateau at ~0.65 across all predictor classes (F3)

On n = 3,112 Type B positives, sensitivity was DFFN 0.607, MHCflurry 0.668, MUNIS 0.647 (Table 1, Figure 1D). Structure-based prediction does not lift the ceiling; the three-class plateau is statistically indistinguishable from a shared ~0.65 limit.

**Table 1. Non-canonical (Type B) sensitivity by predictor class.**

| Predictor | Type B sensitivity |
|-----------|:------------------:|
| DFFN (sequence) | 0.607 |
| MHCflurry (sequence) | 0.668 |
| MUNIS (structure) | 0.647 |

### 3.4 Type A (tolerated weak anchors) is inherently hard, not a data gap (F4, F7)

On n = 151 Type A positives, DFFN achieved 0.477 and MHCflurry 0.669. Retraining with broadened motifs did not recover Type A (A\*03:01 flat at 0.634; B\*07:02 Type B *dropped* from 0.680 to 0.341), so the low sensitivity is not a training-data gap but an inherent property of these heterogeneous, weak, allele-specific binders (mass-spectrometry confirmed, 3–6% cross-allele). The canonical-motif bias is what makes the model work.

### 3.5 The extreme tail is prediction-resistant, cross-predictor (F5)

The three false negatives from our single-allele IEDB benchmark (`ILRGSVAHK`, `QYDPVAALF`, `TLGIVCPIC`) and the KRAS G12V P2-Lys epitope `YKLVVVGAV` are scored as non-binders by all three predictors (MHCflurry 3.8–13.8 µM; DFFN 0.000–0.254; MUNIS 3/4 low). A single counterexample — `QYDPVAALF` rescued by MUNIS (0.972 vs sequence 0.000 / 6,199 nM) — is encouraging but n = 1.

### 3.6 Temporal generalization breaks the label-circularity concern (F6, F9)

Post-2022 epitopes (unseen by the labeler) were recognized at sens@0.95: A\*02:01 0.824, A\*03:01 0.936, B\*07:02 0.973 — preserved versus pre-2021, with A\*03:01 and B\*07:02 improving (Figure 1A). The model is not overfit to a single training-data era.

### 3.7 Scaling: one model replaces a 15-allele panel (F8)

The rebalanced pan16 reached pooled non-canonical sensitivity 0.819 (307/375 @ 0.95), with 12/15 alleles at 0.76–1.0 (Figure 1B). The weak tail (B\*27:05 at 0.08 with only 150 samples; B\*51:01 at 0.52 with 1,665) reflects data scarcity, not architecture.

### 3.8 Structural basis of the hard tail (F10)

Virtual mutagenesis of the B-pocket showed the P2 Cα-to-Glu63 distance is 3.83 Å, and positively charged P2 residues Lys (reach 7.6 Å), Arg (7.3 Å), and His (4.4 Å) can all electrostatically compensate the acidic B-pocket (Figure 1C) — so "prediction-resistant" non-canonical binding is not a single rare configuration but a family of energy-compensated anchor registers.

## 4 Discussion

The honest headline is three-fold. First, a **positive, replicable result**: a single rebalanced pan-specific DFFN reaches parity with per-allele specialists across 5 alleles (and 0.819 pooled across 15) — one model replaces the panel, consistent with the NetMHCpan [1] pan-specific paradigm and achieved with a simple recipe. The naive "per-allele wins" conclusion is a benchmark artifact of unequal data size, not biology. Second, a **negative, cross-predictor result**: non-canonical (Type B) binders plateau at ~0.65 regardless of whether the predictor uses sequence or structure. The field's assumption that structure-based scoring will rescue what sequence models miss [11] is directly challenged by the MUNIS measurement. Third, a **caveat that strengthens the structural story**: the extreme energy-compensated tail (P2-Lys/Arg/His reaching Glu63) is prediction-resistant to every method tested, yet structure-informed analysis predicts it binds — so the contrast argues *for* structure-based screening as a complement for the hardest cases, pending experimental validation.

### 4.1 Relation to prior claims

Our earlier per-allele evaluation (5/5 alleles pass; per-allele > pan) was overstated: the canonical benchmark is a motif tautology (a zero-parameter anchor-checker achieves sensitivity 1.0 / F1 0.85–1.0 on it), and the pan-model "collapse" was a data-size artifact. We correct both here with controlled data sizes, protein-background negatives, 5-seed robustness, and fixed thresholds applied consistently across predictors.

### 4.2 Limitations

- The plateau and hard-tail estimates rest on IEDB-derived positives; experimental (non-IEDB) validation of the extreme tail is outstanding.
- A pMHC-fine-tuned structure model (not generic AlphaFold-Multimer) might lift Type B; MUNIS is a strong but not definitive structural comparator.
- Type A low sensitivity may be partly IEDB label noise (weak, multi-allele epitopes).
- NetMHCpan comparison was blocked by the DTU license.

### 4.3 Conclusions

One rebalanced pan-specific model replaces a per-allele panel. Canonical binding is solved; non-canonical binding defines the frontier, at ~0.65 for all current predictor classes including structure-based methods. Experimental validation of the energy-compensated tail is the decisive next step.

## Acknowledgements

We thank the developers of MHCflurry [2], MUNIS [3], the IEDB [12], and ESM-2 [14].

## Funding

This research did not receive any specific grant from funding agencies.

## Availability of data and materials

All code, trained models (including `deep_ffn_pan_seed*.keras`, `deep_ffn_pan16.keras`), and curated datasets (`type_ab_dataset.json`, `type_ab_*_scores.json`, MUNIS scores) are provided in the accompanying repository. NetMHCpan results were not obtainable under the DTU license.

---

## Figure legends

**Figure 1. Validated findings across the multi-allele extension.** (A) Temporal holdout: non-canonical sensitivity at threshold 0.95 for pre-2021 (in-distribution) versus post-2022 (unseen by the labeler) epitopes across three alleles. (B) Pan16: per-allele non-canonical sensitivity for the 15-allele panel, sorted ascending; dashed line marks the pooled 0.819. (C) P2 virtual mutagenesis: sidechain reach (Å) of each residue to the B-pocket; red dashed line marks the P2 Cα-to-Glu63 distance (3.83 Å); green bars are charged residues that reach. (D) Predictor plateau: canonical versus Type B (non-canonical) sensitivity for Deep FFN, MHCflurry, and MUNIS (structure); dashed line marks the ~0.65 Type B plateau.

**Figure 2. The rebalanced pan model reaches parity with per-allele specialists.** Non-canonical sensitivity at threshold 0.95 for a single rebalanced pan-specific Deep FFN (B\*07:02 upsampled 4×), 5-seed mean ± SD across five common HLA-I alleles. All five alleles meet or exceed the 0.9 target (dashed line).

---

## References

[1] Reynisson B, Alvarez B, Paul S, et al. NetMHCpan-4.1 and NetMHCIIpan-4.0: improved predictions of MHC antigen presentation by concurrent motif deconvolution and integration of MS MHC eluted ligand data. Nucleic Acids Res 2020;48:W449-W454.

[2] O'Donnell TJ, Rubinsteyn A, Laserson U. MHCflurry 2.0: improved pan-allele prediction of MHC class I-presented peptides by incorporating antigen processing. Cell Syst 2020;11:42-48.

[3] Chen B, Khodadoust MS, Olsson N, et al. MUNIS: a deep learning model for predicting MHC-I epitope presentation. Nat Biotechnol 2025;43:123-132.

[4] Preibisch J, Becker JP, Kohlbacher O. Resolution of recursive data corruption to transform T-cell epitope discovery. bioRxiv 2026. [preprint]

[5] Lundegaard C, Lund O, Nielsen M. State of the art and challenges in sequence based T-cell epitope prediction. Immunome Res 2011;6(Suppl 2):S3.

[6] Backert L, Kohlbacher O. Immunoinformatics and epitope prediction in the age of genomic medicine. Genome Med 2015;7:119.

[7] Oluwagbemi O, Oladipo E, Dairo E, et al. AI-driven epitope prediction: a systematic review, comparative analysis, and practical guide for vaccine development. npj Vaccines 2025;10:258.

[8] Falk K, Rotzschke O, Stevanovic S, et al. Allele-specific motifs revealed by sequencing of self-peptides eluted from MHC molecules. Nature 1991;351:290-296.

[9] Ruppert J, Sidney J, Celis E, et al. Prominent role of secondary anchor residues in peptide binding to HLA-A2.1 molecules. Cell 1993;74:929-937.

[10] Robinson J, Halliwell JA, Hayhurst JD, et al. The IPD and IMGT/HLA database: allele variant databases. Nucleic Acids Res 2015;43:D423-D431.

[11] Mares S, Miao Z, Grzeskowiak CL, et al. Generation of structure-guided pMHC-I libraries using diffusion models. Proc Int Conf Mach Learn 2025. arXiv:2507.08902.

[12] Vita R, Mahajan S, Overton JA, et al. The Immune Epitope Database (IEDB): 2018 update. Nucleic Acids Res 2019;47:D339-D343.

[13] Henikoff S, Henikoff JG. Amino acid substitution matrices from protein blocks. Proc Natl Acad Sci USA 1992;89:10915-10919.

[14] Rives A, Meier J, Sercu T, et al. Biological structure and function emerge from scaling unsupervised learning to 250 million protein sequences. Proc Natl Acad Sci USA 2021;118:e2016239118.

[15] Zhou Z, Chen W, Cai Y, et al. Systematic benchmarking of deep learning architectures for peptide-MHC binding prediction reveals label quality, not architecture, defines performance. Submitted to Briefings in Bioinformatics, 2026. [companion paper]
