# Benchmark Hygiene Reveals the Non-Canonical Frontier of Peptide-MHC-I Binding Prediction

**Authors**: Yunshi Huang<sup>1,†</sup>, Weisheng Chen<sup>2,†</sup>, Yiqi Cai<sup>2</sup>, Xiaopeng Zhuang<sup>2</sup>, Zhuxian Zhou<sup>3</sup>, Xiangyang Xue<sup>4,*</sup>, Zhuha Zhou<sup>2,*</sup>

**Word count**: ~3,000 (main text) | **Tables**: 5 | **Figures**: 9

**Affiliations**:  
<sup>1</sup> Department of General Surgery, The First Affiliated Hospital of Wenzhou Medical University, Wenzhou, China  
<sup>2</sup> Department of Gastroenterology Surgery, The First Affiliated Hospital of Wenzhou Medical University, Nanbaixiang Street, Ouhai District, 325000 Wenzhou, Zhejiang, China  
<sup>3</sup> Center for Bionanoengineering, College of Chemical and Biological Engineering, Zhejiang University, Hangzhou 310027, China  
<sup>4</sup> Zhejiang Key Laboratory of Intelligent Cancer Biomarker Discovery and Translation, The First Affiliated Hospital of Wenzhou Medical University, Wenzhou, China  

<sup>†</sup> These authors contributed equally to this work.  
<sup>*</sup> Corresponding authors: Xiangyang Xue, wzxxy@wmu.edu.cn, ORCID: 0000-0002-8670-2016; Zhuha Zhou, zhouzhuha@wmu.edu.cn, ORCID: 0009-0005-1818-4789

**ORCIDs**: Yunshi Huang 0009-0009-8852-1479; Weisheng Chen 0000-0002-5367-9446; Yiqi Cai 0000-0002-4014-6544; Xiaopeng Zhuang 0009-0003-4289-0734; Zhuxian Zhou 0000-0002-7104-9915; Xiangyang Xue 0000-0002-8670-2016; Zhuha Zhou 0009-0005-1818-4789.  

---

## Abstract

Peptide-MHC-I binding predictors are routinely trained on computationally labelled data, yet benchmark scores that reproduce a labeler's decision boundary are rarely separated from true binding accuracy. We present a controlled benchmark of six deep-learning architectures for HLA-A\*02:01 under multi-seed uncertainty, and demonstrate how benchmark composition manufactures false conclusions: two apparent findings dissolved under re-analysis, each traced to a benchmark-composition artifact. The genuinely hard problem is non-canonical binding: the trained model recovers binders that a zero-parameter anchor rule misses. We validate against 42 re-verified IEDB epitopes (97.6% sensitivity), and show that two high-fidelity predictors that rank peptides almost identically (Spearman ρ = 0.96) nevertheless give the same hard label to only 62% of peptides, a discrepancy that reflects a ~4-fold offset in their percentile cutpoints and largely vanishes when the cutpoints are aligned. A benchmark-hygiene checklist distils the countermeasures. Code, models, and data are available at github.com/woodhaha/peptide-mhc-benchmark.

## Author summary

Predicting which peptides bind to MHC class I molecules — the first step in T-cell immune recognition — underpins cancer-vaccine and immunotherapy design. Many machine-learning tools now make these predictions, but it is hard to tell whether a tool's reported accuracy reflects genuine biological insight or merely the quirks of how it was trained and evaluated. We benchmarked six neural-network architectures for the most common human MHC allele (HLA-A\*02:01) and audited the benchmark itself. Two findings emerged. First, benchmark scores are highly sensitive to where the "binder" threshold is drawn: two leading predictors, despite near-identical rankings, give the same hard label to only 62% of peptides; the disagreement traces to a roughly four-fold difference in where each draws its binder threshold and mostly disappears once the thresholds are aligned. Second, the genuinely difficult frontier is non-canonical binding — peptides that violate the classic anchor-position rules, which simple motif rules miss but the trained model recovers. We distil these lessons into a reusable benchmark-hygiene checklist so that others can avoid reporting evaluation artifacts as biological findings.

---

## 1 Introduction

MHC class I molecules present 8–11-mer peptides to CD8+ T cells, and accurate prediction of peptide–MHC-I binding underpins epitope discovery for vaccine design and cancer immunotherapy (Falk et al., 1991; Lundegaard et al., 2011). For HLA-A\*02:01 — the most prevalent class I allele globally — binding is dominated by two primary anchor positions: position 2 (p2, hydrophobic: Leu/Met/Ile/Val) and position 9 (p9, C-terminal, Val/Leu preferred) (Ruppert et al., 1993; Rammensee et al., 1995). Computational predictors have evolved from position-specific scoring matrices to neural networks and pan-specific tools such as NetMHCpan (Reynisson et al., 2020), MHCflurry (O'Donnell et al., 2020), BigMHC (Albert et al., 2023), and MUNIS (Wohlwend et al., 2025), increasingly using protein language model embeddings (Rives et al., 2021).

Several gaps remain. A particularly under-examined one is **cutpoint sensitivity**: because >55% of IEDB entries are computationally rather than experimentally labelled (Preibisch et al., 2026), a predictor's apparent accuracy can track the decision boundary of whichever labeler trained it rather than true binding — and two labelers that agree on peptide ranking may still disagree on hard labels wherever their percentile thresholds differ. The impact of training-label quality is seldom quantified, and general guidance on avoiding data leakage (Bernett et al., 2024) has not been instantiated as a peptide-MHC-specific set of tripwires. Head-to-head architecture comparisons under fully controlled conditions (identical training data, encoding, evaluation) are uncommon. And the translation of binding prediction to real-world epitope discovery is rarely demonstrated in a unified, reproducible pipeline.

To isolate the effects of architecture choice, label quality, and encoding, we present a controlled benchmark. We compare six architectures under multi-seed uncertainty and find that architecture effects are small and largely indistinguishable after accounting for seed and sampling uncertainty, confirm that random labels fall to chance, and evaluate ESM-2 embeddings against BLOSUM62 with a mean-pooled ablation. We further validate against 42 A\*02:01-verified IEDB epitopes, and scan proteins and cancer mutations as a translational application. We report this as a benchmark following the criteria of the PLOS Computational Biology benchmark category (Peters et al., 2018).

## 2 Methods

### 2.1 Training data and labels

100,000 random 9-mer peptides were generated by uniform sampling of the 20 standard amino acids. Each was labelled with MHCflurry 2.2.1 against HLA-A\*02:01 and assigned to Strong Binder (SB, percentile rank <0.5%), Weak Binder (WB, 0.5%≤rank<2.0%), or Non-Binder (NB, rank ≥2.0%) following the netMHCpan convention (Reynisson et al., 2020). The resulting distribution (1.7% SB, 5.1% WB, 93.2% NB) was balanced by downsampling to 1,730 per class (n=5,190), partitioned 90/10 into training (n=4,671) and held-out test (n=519) sets with stratified sampling (seed 42). Two comparator labelling strategies were implemented: a PSSM encoding crystallographic A\*02:01 anchor preferences, and random synthetic labels as a baseline control. A held-out benchmark of 42 experimentally validated IEDB epitopes (HLA-A\*02:01-restricted, with positive binding/T-cell assay evidence, re-verified against the IEDB on 2026-10-23) plus 20 homopolymer negatives (≤2 unique amino acids) was curated from the IEDB (query-api.iedb.org); label-corpus overlap was determined by sequence-identity match against the IEDB entries in MHCflurry's training corpus.

### 2.2 Peptide encoding

Each 9-mer was encoded with the BLOSUM62 matrix (Henikoff and Henikoff, 1992), min–max normalized, producing a 9×20 matrix flattened to 180 dimensions for fully connected networks. ESM-2 embeddings (Lin et al., 2023) were extracted at t6 (2,880-dim per-peptide; 320 per position × 9) and t12 (4,320-dim; 480 per position × 9) depths — the two smallest ESM-2 checkpoints (8M and 35M parameters); a mean-pooled t6 variant (320-dim) served as an ablation control destroying positional information.

### 2.3 Model architectures

Six classifiers were implemented: (i) FFN (Nielsen et al., 2003): 180→90→3, 49K params; (ii) Deep FFN: 180→360→180→90→45→3 with batch normalisation and dropout, 152K params; (iii) CNN: Conv2D(32,3×3) + FFN body, 1.05M params; (iv) LSTM(64)→dense(32)→3, 24K params; (v) ResNet-style stem + 3 residual blocks, 325K params; (vi) Random Forest (100 trees). For ESM-2 embeddings an FFN classifier with a 256→128→64 stack was trained. All networks used cross-entropy loss with Adam or RMSprop (lr 0.001), 150 epochs, batch size 50, 20% internal validation split, early stopping (patience 10), and learning-rate reduction on plateau. Five-fold stratified cross-validation was performed for the FFN. The held-out test set (519 peptides) was disjoint from training and used for all reported metrics. All models were trained and evaluated on a single CPU (TensorFlow 2.x on native Windows, without GPU acceleration); confidence intervals are reported wherever sample sizes permit.

### 2.4 Protein scanning, mutation analysis, and structure

Ten therapeutically relevant proteins (MART-1, gp100, tyrosinase, NY-ESO-1, WT1, p53, KRAS, influenza M1, CMV pp65, SARS-CoV-2 spike RBD) were scanned with a sliding 9-mer window (stride 1; 3,357 windows). Sixteen recurrent hotspot mutations (7 p53, 9 KRAS) were analyzed by comparing wild-type and mutant binding scores; mutations were classified as Created/Enhanced/Destroyed/Unchanged. The KRAS G12V neoepitope YKLVVVGAV carries a non-canonical lysine at P2.

### 2.5 Feature extraction

Learned features were extracted from the penultimate dense layer via truncated models. Hierarchical clustering (Ward's method, Euclidean distance) was applied to representations of top epitope candidates; feature–binding-score correlation was assessed by Pearson's r with a 1,000-iteration permutation test to control for chance correlations. t-SNE visualisation of the 45-dimensional penultimate layer was used to inspect class separation, and integrated gradients identified position-specific anchor contributions.

## 3 Results

The results are organised to isolate the failure mode of interest. Sections 3.1–3.3 vary the three variables that could confound a benchmark — architecture, label source, and encoding — and show that each has only a small effect at this data scale. This controlled setup is what allows §3.4 to attribute the large label-agreement gap between two high-fidelity predictors to percentile-cutpoint placement rather than to architecture choice, label quality, or feature encoding.

### 3.1 Training data and architecture comparison

MHCflurry labelling of 100,000 random 9-mers yielded 1.7% SB, 5.1% WB, 93.2% NB; after balancing, 5,190 peptides (1,730/class). Across five training seeds, the CNN achieved the highest accuracy (91.5% ± 0.8%; Table 1, Figure 1), modestly ahead of the Deep FFN (90.6% ± 1.2%), the FFN (89.8% ± 1.1%), and the LSTM (89.7% ± 0.9%) — the top four architectures within ~2 percentage points, a range within the combined seed-plus-binomial uncertainty (n=519; 95% CI CNN [0.886, 0.943] vs Deep FFN [0.872, 0.940]) — while ResNet (87.6% ± 1.8%) and Random Forest (79.7% ± 0.4%) were clearly worse. Two external predictors frame this comparison (Table 1): MHCflurry 2.2.1 itself scores 100% on the labels it produced — an expected ceiling, not evidence of superiority — while the independent state-of-the-art predictor NetMHCpan-4.1, which was not trained on these labels, agrees with them on only 62.24% (323/519). The 8.5-point gap between the labeler (100%) and the best finite-sample approximation (91.5% for the CNN) is information loss from finite samples; the far larger ~29-point gap to NetMHCpan (62.24%) reflects primarily cutpoint placement rather than label fidelity (see §3.4 and Discussion). The WB class was consistently the most difficult (F1 0.752–0.880). These class-balanced accuracies overstate real-world utility: under the uniform-random-peptide prior (1.7% SB, 5.1% WB, 93.2% NB), a trivial always-non-binder classifier reaches 93.2% accuracy, above every model's balanced accuracy, so natural-prior precision/PPV should be reported (Box 1, item 5).

**Table 1.** Model performance on the MHCflurry-labelled test set (mean ± SD across five training seeds, n=519).

| Model | Accuracy (%) | Parameters |
|-------|:-----------:|:----------:|
| CNN | **91.5 ± 0.8** | 1,053,863 |
| Deep FFN | 90.6 ± 1.2 | 152,823 |
| FFN (Nielsen) | 89.8 ± 1.1 | 49,143 |
| LSTM | 89.7 ± 0.9 | 23,939 |
| ResNet | 87.6 ± 1.8 | 324,931 |
| Random Forest | 79.7 ± 0.4 | -- |
| MHCflurry 2.2.1 (labeler) | 100 | — |
| NetMHCpan-4.1 (independent SOTA) | 62.24 | — |

† Accuracy is agreement with the MHCflurry test labels (n=519). MHCflurry scores 100% because it produced the labels — an expected ceiling, not evidence of superiority. NetMHCpan-4.1, which was not trained on these labels, agrees with them on only 62.24%, quantifying cutpoint sensitivity (see Discussion).

The modest spread among the top four architectures reflects the positional nature of peptide–MHC binding: anchor positions p2/p9 dominate, so architectures that preserve positional information perform comparably. ResNet and the Random Forest were clearly worse, consistent with the absence of translation invariance in a 9×20 encoding and, for the Random Forest, the high dimensionality of the BLOSUM62 encoding.

### 3.2 Label-source controls

Label source, not architecture, sets the practical ceiling of what a model can learn (Table 2, Figure 2). Training on MHCflurry labels reached 90.6% accuracy, whereas a deterministic PSSM label source reached 82.6% (on a separate 15,000-peptide pool, so not a clean head-to-head comparison) and an i.i.d. random-label null fell to chance (34.5 ± 1.7%, close to the 33.3% expected for three balanced classes). The random-label result confirms the models learn genuine signal rather than noise; the PSSM–MHCflurry gap, though confounded by pool size, suggests MHCflurry's labels are at least as learnable as a hand-crafted scoring matrix.

**Table 2.** Deep FFN performance across label sources (accuracy, mean ± SD across five seeds where noted).

| Label source | Accuracy (%) | Notes |
|-------------|:-----------:|-------|
| DFFN on MHCflurry labels | 90.6 ± 1.2 | 5,190 peptides, 5 seeds |
| PSSM (biophysics) | 82.6 | separate 15,000-peptide pool, single seed |
| Random (i.i.d.) | 34.5 ± 1.7 | chance-level null, 5 seeds |

### 3.3 Encoding

On the Deep FFN (the architecture used for the encoding comparisons), BLOSUM62 generalized stably under cross-validation (92.1 ± 0.4%; per-fold 92.0/92.6/91.7/92.4/91.8%). ESM-2 t6 per-position embeddings scored 93.3% on a single split but collapsed to 82.2 ± 1.0% under CV — the 2,880-dimensional features overfit the 4,671-sample training set (Table 3, Figure 7). Mean-pooled t6 embeddings collapsed to 65.9%, confirming that positional information is essential. The t12 variant (4,320-dim) overfit similarly (90.9–91.5%). The single-split 93.3% result reflects a favourable partition rather than a genuine encoding advantage at this training scale. ESM-2's representational advantage materialises only at large scale: the 650M checkpoint (11,520-dim per-peptide) reached 92.4 ± 0.5% under CV, matching BLOSUM62 (92.1 ± 0.4%) but at 64× the input dimensionality. At the data scale typical of single-allele studies, compact biologically motivated encodings such as BLOSUM62 remain the more cost-effective choice, while the small ESM-2 checkpoints overfit.

**Table 3.** Classification accuracy by encoding strategy.

| Encoding | Dim | Accuracy (%) | Notes |
|----------|:---:|:-----------:|-------|
| BLOSUM62 | 180 | 92.1 ± 0.4 (CV) | Stable CV generalisation |
| ESM-2 t6 mean pool | 320 | 65.9 | Positional info destroyed |
| ESM-2 t12 per-pos | 4,320 | 90.9 | Overfitting |
| ESM-2 t12 + L2 | 4,320 | 91.5 | Partially mitigated |
| ESM-2 t6 per-pos | 2,880 | 93.3 (82.2 ± 1.0 CV) | Single-split optimistic |
| ESM-2 t33 (650M) per-pos | 11,520 | 92.4 ± 0.5 (CV) | Matches BLOSUM62 at 64× dim |

### 3.4 Held-out validation and translation

The Deep FFN achieved 97.6% sensitivity (41/42; Wilson 95% CI 87.7–99.6%), 75.0% specificity (15/20), ROC AUC 0.970 against IEDB epitopes (Table 5, Figure 3). Of the original 49 A\*02:01-restricted epitopes, one (`ILRGSVAHK`) was removed on re-curation for allele mis-attribution (HLA-A\*03:01, previously mis-attributed to B\*08:01), leaving 48. IEDB re-verification (query-api.iedb.org, accessed 2026-10-23) then excluded six further peptides lacking A\*02:01 restriction — `LLFDRFENL` and `LLHHAFDSL` (no IEDB records), `YVKQNTLKL` (HLA-B\*07:02), `QYDPVAALF` and `CMTWNQMNL` (HLA-A\*24:02), and `VIYQYMDDL` (no positive A\*02:01 record) — giving the final 42-epitope benchmark. The single remaining false negative, `TLGIVCPIC` (HPV16 E7), carries a non-canonical p9 anchor (C). An independent predictor, NetMHCpan-4.1b (DTU standalone, run locally), achieved 85.7% sensitivity (36/42; Wilson 95% CI 72.1–93.3%) on the corrected 42-epitope set, statistically indistinguishable from the Deep FFN (41/42 = 97.6%; exact McNemar p = 0.0625, five discordant pairs) — though with only five discordant pairs this comparison is underpowered and establishes non-inferiority, not equivalence. NetMHCpan's pan-specific training on MS-elution-based labels makes this an indicative rather than clean comparison; the more direct evidence of cutpoint sensitivity is the 62.24% label agreement (323/519 test peptides) reported in the Discussion. To contextualise the benchmark, we include a zero-parameter anchor-checker baseline (a rule classifying any 9-mer with canonical p2/p9 residues as a binder): it achieves sensitivity 0.786 (33/42) and specificity 0.900 (18/20), below the trained model (0.976), because nine genuine binders carry non-canonical p2/p9 anchors. It recovers eight of these nine (Wilson 95% CI 52–98%) — the single residual miss, `TLGIVCPIC` (p9 = C), defines the hard non-canonical tail — suggesting — on a sample of only nine peptides, too small for strong inference — that part of the informative signal may lie beyond the anchor motif. The five false positives were homopolymers (≤2 unique amino acids) with canonical anchors — biologically implausible binders whose low sequence complexity the classifier fails to penalise. A domain-knowledge filter (≤2 unique amino acids → non-binder) removed all five, but we report 75.0% as the primary, unfiltered specificity to avoid post-hoc inflation. As a secondary negative set, 20 experimentally measured A\*02:01 negatives (IEDB qualitative_measure = Negative) yielded 90.0% specificity (18/20), confirming that the homopolymer set is a conservative lower bound on realistic specificity. To characterise the labels independently of the labeler's own decision boundary, we additionally scored the Deep FFN against a large label-independent set of experimentally measured affinities (IEDB binding-affinity corpus, accessed 2026-08-02) — 3,438 unique 9-mers with measured IC50 <500 nM and 3,086 with IC50 >5,000 nM (the 500–5,000 nM grey zone was excluded to yield unambiguous binders and non-binders), restricted to direct binding assays (competitive fluorescence/radioactivity) and excluding MS-elution-derived, predictor-dependent entries. The MHCflurry-trained Deep FFN recognised 94.4% of the experimental binders and rejected 74.5% of the non-binders (threshold 0.5; ROC AUC 0.942); the rebalanced per-allele model reached 89.8% sensitivity and 89.1% specificity (AUC 0.945). To characterise the labelers directly, we scored MHCflurry 2.2.1 itself on the 6,524 affinities (AUC 0.950) alongside NetMHCpan-4.1b (0.937 EL-rank, 0.960 BA-rank) — both high, though 91.6% of these peptides appear in MHCflurry's training corpus, so these values are largely in-corpus. The two predictors also rank peptides almost identically (Spearman ρ = 0.96), yet their hard labels agree on only 62.24% of the 519 test peptides (Cohen's κ = 0.44). The disagreement is asymmetric — NetMHCpan downgrades MHCflurry's weak binders to non-binders (73% of WB) and strong to weak (30%) — and largely vanishes under the natural class prior (96% agreement), localising the 62% figure to cutpoint placement rather than label definition or fidelity. To quantify this cutpoint sensitivity, we varied the percentile-rank threshold and recomputed the two predictors' binary binder agreement (Figure 8). Agreement is non-monotonic, dipping to 73.4% at the standard 2% threshold, where MHCflurry calls 69.6% of the balanced test set binders versus 43.0% for NetMHCpan. The two rank distributions are offset roughly four-fold — MHCflurry's top-2% threshold corresponds to NetMHCpan's top-8% — and aligning the cutpoints recovers 96.1% agreement. The predictors thus share predictive skill but differ in calibration (where the binder threshold is drawn), a decoupling recently formalised as a 'free knob' in threshold-based evaluation (Munim et al., 2026). These two effects are independent: the 62% three-class agreement reflects the balanced prior, whereas the 73% binary agreement at 2% reflects the cutpoint offset; correcting either effect separately restores ~96% agreement. This pattern generalises across objectives (Figure 9): two presentation predictors — MUNIS and BigMHC, which add antigen processing to binding — correlate with each other at ρ = 0.90 but with the two binding predictors at only |ρ| ≈ 0.72–0.78, below the binding predictors' mutual 0.96. Cutpoint sensitivity is therefore specific to predictors sharing a binding-affinity objective, while cross-objective disagreement reflects a genuinely different task. Protein scanning of 10 proteins (3,357 9-mer windows) recovered all ten in-scope known epitopes (100%; Figure 4), including immunodominant `GILGFVFTL` (M1), `NLVPMVATV` (pp65), and `SLLMWITQC` (NY-ESO-1); the MART-1 clinical variant `ELAGIGILTV` is a 10-mer A27L anchor-modified peptide absent from the wild-type MART-1 sequence, and therefore lies outside the 9-mer scanning space. Mutation scanning of 16 hotspot mutations identified 7 epitope-altering events (Table 4, Figure 5). p53 R248W was predicted to create the neoepitope `MNWRPILTI` (the R→W substitution at P3 alters a solvent-exposed, TCR-facing position, so any immunogenicity would arise from creating a non-self sequence rather than improved MHC binding, a mechanism our binding-only model cannot directly assess); KRAS G12V was predicted to create `YKLVVVGAV` via G→V at p9 creating a canonical C-terminal anchor. p53 R248W is among the most frequent TP53 missense mutations across cancers (~3.5% of somatic TP53 mutations, IARC), and mutation-specific T-cell responses have been documented in patients and in HLA-A2/DRB1-transgenic mice (Malekzadeh et al., 2020; Quandt et al., 2018); KRAS G12C enhanced an existing epitope (WB→SB), suggesting synergy with pharmacological G12C inhibition. Three mutations (KRAS G13D, KRAS A146T, p53 R249S) destroyed existing epitopes, representing potential immune-evasion mechanisms.

**Table 4.** Epitope-altering cancer hotspot mutations.

| Mutation | WT peptide | Mutant peptide | Effect | Δ score |
|----------|-----------|----------------|--------|:-----:|
| p53 R248W | MNRRPILTI | MNWRPILTI | CREATED | +0.41 |
| KRAS G12V | YKLVVVGAG | YKLVVVGAV | CREATED | +0.48 |
| KRAS G12V | LVVVGAGGV | LVVVGAVGV | ENHANCED | +0.31 |
| KRAS G12C | LVVVGAGGV | LVVVGACGV | ENHANCED | +0.31 |
| KRAS G13D | GVGKSALTI | DVDKSALTI | DESTROYED | −0.47 |
| KRAS A146T | GIPFIETSA | GIPFIETST | DESTROYED | −0.42 |
| p53 R249S | GMNRRPILT | GMNRSPILT | DESTROYED | −0.16 |

**Table 5.** IEDB benchmark performance.

| Metric | Value |
|--------|:-----:|
| Sensitivity | 97.6% (41/42) |
| NetMHCpan-4.1b sensitivity | 85.7% (36/42) |
| Specificity | 75.0% (15/20) |
| Precision | 89.1% |
| F1 | 0.932 |
| ROC AUC | 0.970 |

*NetMHCpan-4.1b standalone (DTU academic license, run locally under WSL2).

### 3.5 Learned feature representations

For interpretability, we extracted the Deep FFN penultimate dense layer (Figure 6). t-SNE showed partial class separation, with strong binders forming a distinct cluster and weak binders a transition zone; integrated gradients attributed position-specific anchor contributions (p9 dominated strong binders, p2/p3 weak binders). Because these features are the penultimate layer of a model trained to predict binding, their correlation with the binding score is expected and is reported for completeness rather than as an independent finding.

## 4 Discussion

### Methodological warnings

Two apparent findings in this project dissolved under controlled re-analysis, each illustrating how benchmark composition can manufacture false conclusions. First, a claimed 'label quality dominates architecture' result collapsed once the degenerate random-label null (which performed at chance, ~33%) was removed — the realistic label effect is smaller than the architecture range. Second, a reported 'simple architecture beats complex ones' ranking was an artifact of stale model checkpoints: re-training with fixed seeds showed a CNN instead ranked above the Deep FFN (91.5% vs 90.6%, higher in four of five seeds). We recommend that peptide-MHC benchmarks report multi-seed estimates with confidence intervals, protein-background negatives, frozen data snapshots, and a state-of-the-art comparator such as NetMHCpan, to avoid these failure modes. The same failure modes recur beyond peptide-MHC: data leakage from protein language-model pretraining distorts downstream task performance whenever train/test splits are not pretraining-aware (Hermann et al., 2024), and leakage is now recognised as a field-wide reproducibility problem affecting hundreds of papers across seventeen fields (Kapoor and Narayanan, 2023). Dedicated guidance now exists for avoiding leakage in biological machine learning (Bernett et al., 2024) and for reporting clinical prediction models (Collins et al., 2024); a recent audit of biomolecular benchmarks found cross-split contamination and label conflicts pervasive across the field (Schuh et al., 2026); and label errors in test sets have been shown to destabilise machine-learning benchmarks more broadly (Northcutt et al., 2021). Box 1 instantiates these general principles for peptide-MHC benchmarking specifically.

**Box 1. Benchmark hygiene recommendations.** The two artifacts translate into reusable practices for immunoinformatics and machine-learning benchmarks: (1) use a true i.i.d. random-label null (chance ≈ 33% for three balanced classes) rather than a degenerate control; (2) freeze data snapshots, model checkpoints, and label versions, and report multi-seed mean ± SD with confidence intervals; (3) report label provenance and, alongside hard labels, the continuous-score rank correlation with an independent state-of-the-art comparator (NetMHCpan), to separate cutpoint sensitivity from fidelity; (4) report performance under a realistic class prior (e.g. proteome-scan prevalence), not only a class-balanced prior, since a class-balanced benchmark overstates real-world utility. Each failure mode reproduced a false conclusion that was caught only by controlled re-analysis.

Architecture choice has a small and statistically inconclusive effect on peptide-MHC binding prediction. Across five seeds, a CNN (91.5%) ranked above a Deep FFN (90.6%) by ~0.9 percentage points (higher in four of five seeds; a small margin within seed-to-seed variability), while a Random Forest was clearly worse (79.7%), for a maximum inter-architecture range of ~12 percentage points. The top four architectures fell within ~2 percentage points, indicating that data quality and the decision boundary matter more than architectural complexity for this task. Random-label controls performed at chance (~33%), confirming that the observed performance reflects genuine signal rather than noise overfitting.

A PSSM-derived label source yielded lower accuracy than MHCflurry labels (82.6% vs 90.6%, though on a separate peptide pool and so not a clean comparison). This does not support the hypothesis that deterministic biophysical labels are inherently superior to ML-derived labels; rather, it suggests that MHCflurry's labels are at least as learnable as a simple PSSM. We emphasise that training on any computational label — including MHCflurry — re-deposits predictor bias into downstream models, and recommend that benchmarks report label provenance and the WB-class F1 as a standard diagnostic.

An independent predictor, NetMHCpan-4.1, agreed with the MHCflurry labels on only 62.24% (323/519) of the held-out test peptides, versus 90.6% for a model trained on those labels. This ~28-point gap could reflect low label fidelity, low rank agreement, or different cutpoints. Scoring NetMHCpan on the same 6,524 measured affinities gives 0.937 AUC (§3.4), and its ranks correlate with MHCflurry's at Spearman ρ = 0.96, so neither fidelity nor ranking explains the gap; it localises to the placement of percentile cutpoints (NetMHCpan systematically downgrading weak binders).

Two predictors (MHCflurry and NetMHCpan) rank peptides almost identically (Spearman ρ = 0.96) and share a high AUC on experimental affinities (~0.95), yet agree on only 62% of hard labels. This separation of predictive skill (ranking, AUC) from calibration (where the binder threshold is drawn) is a general property of threshold-based evaluation, recently formalised as a 'free knob' in model benchmarking (Munim et al., 2026); two models with equal predictive skill can disagree substantially at any fixed threshold (Sanz-Guerrero and von der Wense, 2026). The two labelers' near-identical ranking is consistent with prior systematic benchmarking (Paul et al., 2020; Zhao and Sher, 2018); the contribution here is quantifying that hard-label agreement is nonetheless only 62%, and localising that gap to percentile-cutpoint placement rather than label fidelity or definition. Benchmarks that report accuracy against a single computational labeler therefore measure that labeler's specific cutpoints, and overstate generalisation whenever labelers draw different thresholds. The countermeasures in Box 1 are allele-agnostic practices — a random-label null, multi-seed reporting, a rank-correlation comparator, and natural-prior reporting apply to any allele — whereas the ~4-fold offset we quantify is an A\*02:01-specific instantiation; whether other alleles show a smaller or larger offset is precisely what a per-allele balanced benchmark is designed to answer.

BLOSUM62 proved more robust than ESM-2 at single-allele scales: the 2,880-dimensional per-position embeddings overfit the 4,671-sample training set, and the single-split 93.3% reflects a favourable partition rather than a genuine advantage. The shallower 6-layer ESM-2 variant outperformed the deeper 12-layer model, and a PSSM label source (82.6%) underperformed MHCflurry labels. The language-model route may still be viable with continued pre-training on HLA-associated peptides (Mares et al., 2025).

The non-canonical tail — where a zero-parameter anchor rule fails but the trained model still recovers eight of nine such binders — extends to non-coding, fusion-derived, and post-translationally modified peptides, for which current predictors remain under-powered (Vo et al., 2025). A new wave of transformer- and PLM-based predictors continues to improve canonical prediction (Xu et al., 2026; Cheng et al., 2025), but whether any of them lifts the non-canonical plateau remains to be seen.

### Translational implications

p53 R248W, which created `MNWRPILTI`, is among the most frequent TP53 missense mutations (~3.5% of somatic mutations, IARC); its vaccine potential requires experimental confirmation (IFN-γ ELISpot). KRAS G12C, which enhanced an existing epitope (WB→SB), is targeted by sotorasib and adagrasib, suggesting combination immunotherapy. We note that binding prediction alone does not establish immunogenicity: peptide–MHC stability predicts immunogenicity better than affinity (Harndahl et al., 2012; Backert and Kohlbacher, 2015), and high-throughput thermostability profiling now quantifies this at immunopeptidome scale (Jappe et al., 2020; Fasoulis et al., 2024).

### Limitations

Training labels are MHCflurry-derived rather than experimental affinities, so our single-allele benchmark is itself subject to the cutpoint sensitivity we describe; we quantify this by cross-referencing the 42 A\*02:01-verified IEDB epitopes against MHCflurry's training corpus; 40 of the 42 overlapped with MHCflurry's training corpus, and the two non-overlapping epitopes (`YILEETSVM`, `YLQPRTFLL`) were nonetheless correctly classified, providing limited evidence of generalisation beyond the training distribution. The single-allele design was deliberate — holding the allele constant isolates label-source variation from allele-driven confounds. The cutpoint-sensitivity quantification (the ~4-fold offset) is accordingly specific to the A\*02:01 benchmark. Testing whether it transfers is not a matter of re-scoring the same peptides for other alleles: those peptides are A\*02:01-selected and mostly non-binders elsewhere (rank agreement drops to Spearman ρ = 0.73 and 0.82 for A\*01:01 and A\*03:01 versus 0.96 for A\*02:01). A proper cross-allele test therefore requires a balanced benchmark per allele — which is itself the methodological point: per-allele balanced benchmarks, not pooled multi-allele ones, are the correct unit of analysis; a systematic per-allele survey remains future work. All translational predictions and the cancer neoepitope candidates require experimental confirmation (e.g., IFN-γ ELISpot or tetramer staining). Antigen-processing steps (proteasomal cleavage, TAP transport) are not modelled.

## 5 Conclusion

Architecture mattered little: the CNN ranked above the Deep FFN in four of five seeds, the top four architectures fell within ~2 percentage points, and only the Random Forest lagged (79.7%). Random-label controls confirmed chance-level performance, ruling out noise overfitting. Canonical binding is largely captured by the anchor motif, but the trained model (0.98 sensitivity) recovers non-canonical binders that a zero-parameter anchor rule (0.79) misses. Beyond peptide-MHC, the two artifacts we corrected — a degenerate null and stale checkpoints — illustrate a general benchmark-composition failure mode that reproduces false conclusions across machine-learning benchmarks. Benchmark hygiene reveals that two high-fidelity labelers, near-identical in ranking (ρ = 0.96), still agree on only 62% of hard labels because their percentile cutpoints differ ~4-fold; aligning the cutpoints restores agreement. The Box 1 checklist distils these lessons into a reusable safeguard for other benchmarks. Experimental validation of the non-canonical tail is the decisive next step.

## Figure legends

**Figure 1.** Model performance comparison across six architectures on the MHCflurry-labelled test set (accuracy, 5-seed mean ± SD).

**Figure 2.** Architecture comparison (Table 1) and label-source controls: MHCflurry (90.6%), PSSM (82.6%), random labels at chance (~33%).

**Figure 3.** IEDB benchmark: ROC curve (AUC 0.970), confusion matrix, per-epitope scores.

**Figure 4.** Protein epitope scanning: SB/WB/NB distribution across 10 proteins and top epitope candidates.

**Figure 5.** Cancer hotspot mutation analysis: delta scores and wild-type/mutant pairs for 7 epitope-altering mutations.

**Figure 6.** t-SNE of the Deep FFN 45-dimensional penultimate-layer representation, coloured by binding class (SB/WB/NB).

**Figure 7.** ESM-2 embedding comparison across five encoding strategies.

**Figure 8.** Cutpoint sensitivity (the calibration "free knob"). (a) Binary binder agreement between MHCflurry and NetMHCpan as a function of a shared percentile-rank threshold; agreement dips at the standard 2% cutpoint. (b) Agreement with MHCflurry's 2% binder call as NetMHCpan's threshold is varied, showing an ~4-fold calibration offset (MHCflurry 2% ≈ NetMHCpan 8%).

**Figure 9.** Cross-predictor rank-correlation matrix (|Spearman ρ|). Predictors sharing an objective correlate highly (binding–binding 0.96; presentation–presentation 0.90), whereas cross-objective correlation is moderate (0.72–0.78), reflecting the additional antigen-processing component of presentation predictors.

## Data and code availability

All code, trained models, feature matrices, prediction results, and the experimental-label benchmark (`experimental_label_scores.json`: 6,524 peptides — 3,438 binders with measured IC50 <500 nM and 3,086 non-binders with IC50 >5,000 nM) are available at github.com/woodhaha/peptide-mhc-benchmark and archived at Zenodo (DOI: 10.5281/zenodo.XXXXXXX). The Deep FFN is distributed as the `pmhcbenchmark` Python package, with a `pmhc-predict` command-line tool and a standalone Windows GUI (`pmhc_predict.exe`); exported weights (`pmhc_models.npz`) enable tensorflow-free inference. The single-allele benchmark-hygiene analysis constitutes the complete study. A manuscript-consistency checker (`peptide_audit.py` with `audit_config.json`) automates the reference-order, factual-tripwire, and figure-inventory checks used in this audit.

## Author contributions

Conceptualization: Yunshi Huang, Weisheng Chen, Zhuha Zhou. Methodology: Yunshi Huang, Weisheng Chen. Formal analysis: Yunshi Huang, Weisheng Chen, Yiqi Cai. Investigation: Yunshi Huang, Weisheng Chen, Yiqi Cai. Validation: Xiaopeng Zhuang, Zhuxian Zhou. Writing – original draft: Yunshi Huang, Weisheng Chen, Yiqi Cai. Writing – review & editing: Xiaopeng Zhuang, Zhuxian Zhou, Xiangyang Xue, Zhuha Zhou. Supervision: Xiangyang Xue, Zhuha Zhou.

## Competing interests

The authors declare no competing interests.

## Funding

This research did not receive specific grant funding.

## Ethics statement

This study used only publicly available, de-identified data from the Immune Epitope Database (IEDB) and public protein databases; no human subjects, animals, or patient-derived samples were involved, and no ethics approval or informed consent was required.

## AI use disclosure

AI-assisted tools (Claude Code, Anthropic) were used for code development and manuscript editing. All scientific content, interpretations, and conclusions were authored and approved by the human authors, who take full responsibility for the manuscript's content. No AI tools were listed as authors.

## References

Albert BA, Yang Y, Shao XM, et al. Deep neural networks predict class I major histocompatibility complex epitope presentation and transfer learn neoepitope immunogenicity. Nat Mach Intell 2023;5:861-872.
Backert L, Kohlbacher O. Immunoinformatics and epitope prediction in the age of genomic medicine. Genome Med 2015;7:119.
Bernett J, Blumenthal DB, Grimm DG, et al. Guiding questions to avoid data leakage in biological machine learning applications. Nat Methods 2024. doi:10.1038/s41592-024-02362-y.
Cheng X, Mei H, Chen P, Wu H, Liu R, Lei Y, Wang P. ESMpHLA: evolutionary scale model-based deep learning prediction of HLA class I binding peptides. HLA 2025;105:e70263.
Collins GS, Moons KGM, Dhiman P, et al. TRIPOD+AI statement: updated guidance for reporting clinical prediction models that use regression or machine learning methods. BMJ 2024;385:e078378.
Falk K, Rotzschke O, Stevanovic S, et al. Allele-specific motifs revealed by sequencing of self-peptides eluted from MHC molecules. Nature 1991;351:290-296.
Fasoulis R, Rigo MM, Antunes DA, et al. Transfer learning improves pMHC kinetic stability and immunogenicity predictions. ImmunoInformatics 2024. doi:10.1016/j.immuno.2023.100030.
Harndahl M, Rasmussen M, Roder G, et al. Peptide-MHC class I stability is a better predictor than peptide affinity of CTL immunogenicity. Eur J Immunol 2012;42:1405-1416.
Henikoff S, Henikoff JG. Amino acid substitution matrices from protein blocks. Proc Natl Acad Sci USA 1992;89:10915-10919.
Hermann L, Fiedler T, Nguyen HA, et al. Beware of data leakage from protein LLM pretraining. Proc Mach Learn Res 2024;261:106-116.
Jappe EC, Garde C, Ramarathinam SH, et al. Thermostability profiling of MHC-bound peptides: a new dimension in immunopeptidomics and aid for immunotherapy design. Nat Commun 2020. doi:10.1038/s41467-020-20166-4.
Kapoor S, Narayanan A. Leakage and the reproducibility crisis in machine-learning-based science. Patterns 2023;4:100804.
Lin Z, Akin H, Rao R, et al. Evolutionary-scale prediction of atomic-level protein structure with a language model. Science 2023;379:1123-1130.
Lundegaard C, Lund O, Nielsen M. State of the art and challenges in sequence based T-cell epitope prediction. Immunome Res 2011;6(Suppl 2):S3.
Malekzadeh P, Yossef R, Cafri G, et al. Antigen experienced T cells from peripheral blood recognize p53 neoantigens. Clin Cancer Res 2020. doi:10.1158/1078-0432.CCR-19-1874.
Mares SE, Espinoza Weinberger A, Ioannidis NM. Continued domain-specific pre-training of protein language models for pMHC-I binding prediction. Proc Mach Learn Res 2025;311:304-325.
Munim MH, Saiem BA, Sany A-A, Hashem T. A free knob: decoupling calibration and predictive skill in threshold-based evaluation. arXiv 2026;2609.33457.
Nielsen M, Lundegaard C, Worning P, et al. Reliable prediction of T-cell epitopes using neural networks with novel sequence representations. Protein Sci 2003;12:1007-1017.
Northcutt CG, Athalye A, Mueller J. Pervasive label errors in test sets destabilize machine learning benchmarks. Proc Conf Neural Inf Process Syst Datasets Benchmarks Track 2021.
O'Donnell TJ, Rubinsteyn A, Laserson U. MHCflurry 2.0: improved pan-allele prediction of MHC class I-presented peptides by incorporating antigen processing. Cell Syst 2020;11:42-48.
Paul S, Croft NP, Purcell AW, Tscharke DC, Sette A, Nielsen M, Peters B. Benchmarking predictions of MHC class I restricted T cell epitopes in a comprehensively studied model system. PLoS Comput Biol 2020;16:e1007757.
Peters B, Brenner SE, Wang E, Slonim D, Kann MG. Putting benchmarks in their rightful place: The heart of computational biology. PLoS Comput Biol 2018;14:e1006494.
Preibisch G, Tyrolski M, Kucharski P, et al. Resolution of recursive data corruption to transform T-cell epitope discovery. bioRxiv 2026.
Quandt J, Schlude C, Bartoschek M, et al. Long-peptide vaccination with driver gene mutations in p53 and Kras induces cancer mutation-specific effector as well as regulatory T cell responses. Oncoimmunology 2018;7:e1500671.
Rammensee HG, Friede T, Stevanovic S. MHC ligands and peptide motifs: first listing. Immunogenetics 1995;41:178-228.
Reynisson B, Alvarez B, Paul S, et al. NetMHCpan-4.1 and NetMHCIIpan-4.0: improved predictions of MHC antigen presentation by concurrent motif deconvolution and integration of MS MHC eluted ligand data. Nucleic Acids Res 2020;48:W449-W454.
Rives A, Meier J, Sercu T, et al. Biological structure and function emerge from scaling unsupervised learning to 250 million protein sequences. Proc Natl Acad Sci USA 2021;118:e2016239118.
Ruppert J, Sidney J, Celis E, et al. Prominent role of secondary anchor residues in peptide binding to HLA-A2.1 molecules. Cell 1993;74:929-937.
Sanz-Guerrero M, von der Wense K. Calibration as a first-class criterion in LLM evaluation. arXiv 2026;2609.26489.
Schuh MG, Daniluk A, Sieber SA. Auditing widely used biomolecular benchmarks reveals systematic data inconsistencies. Chem Sci 2026;17:17455-17467.
Vo TH, McNeela E, O'Donovan O, et al. Artificial intelligence and the evolving landscape of immunopeptidomics. Proteomics Clin Appl 2025. doi:10.1002/prca.70018.
Wohlwend J, Nathan A, Barzilay R, et al. Deep learning enhances the prediction of HLA class I-presented CD8+ T cell epitopes in foreign pathogens. Nat Mach Intell 2025. doi:10.1038/s42256-024-00971-y.
Xu H, Ni Y, Chai Z, Cui X, Lei X, Liu L, Shan J, Qian C. TransBindpMHCI: a transformer-based model for pan-specific MHC-I peptide binding prediction. BMC Bioinformatics 2026;27:117.
Zhao W, Sher X. Systematically benchmarking peptide-MHC binding predictors: From synthetic to naturally processed epitopes. PLoS Comput Biol 2018;14:e1006457.
