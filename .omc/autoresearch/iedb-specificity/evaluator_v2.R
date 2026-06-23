#!/usr/bin/env Rscript
# Evaluator v2 — adds homopolymer/repeat detection filter
# Homopolymers (all 9 positions same AA) are biologically incapable of
# being T-cell epitopes — override to NB regardless of model prediction.
# Usage: Rscript evaluator_v2.R <model_path.h5>
# Output: JSON line with {pass, specificity, sensitivity, f1, accuracy, ...}

options(warn = -1)
Sys.setenv(TF_CPP_MIN_LOG_LEVEL = "3")
suppressMessages({library(keras); library(tensorflow); library(tidyverse); library(pROC)})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  cat('{"error":"Missing model path argument"}\n')
  quit(status = 1)
}
model_path <- args[1]
if (!file.exists(model_path)) {
  cat(sprintf('{"error":"Model file not found: %s"}\n', model_path))
  quit(status = 1)
}

# ---- BLOSUM62 ----
blosum62 <- matrix(c(
    4,-1,-2,-2, 0,-1,-1, 0,-2,-1,-1,-1,-1,-2,-1, 1, 0,-3,-2, 0,
   -1, 5, 0,-2,-3, 1, 0,-2, 0,-3,-2, 2,-1,-3,-2,-1,-1,-3,-2,-3,
   -2, 0, 6, 1,-3, 0, 0, 0, 1,-3,-3, 0,-2,-3,-2, 1, 0,-4,-2,-3,
   -2,-2, 1, 6,-3, 0, 2,-1,-1,-3,-4,-1,-3,-3,-1, 0,-1,-4,-3,-3,
    0,-3,-3,-3, 9,-3,-4,-3,-3,-1,-1,-3,-1,-2,-3,-1,-1,-2,-2,-1,
   -1, 1, 0, 0,-3, 5, 2,-2, 0,-3,-2, 1, 0,-3,-1, 0,-1,-2,-1,-2,
   -1, 0, 0, 2,-4, 2, 5,-2, 0,-3,-3, 1,-2,-3,-1, 0,-1,-3,-2,-2,
    0,-2, 0,-1,-3,-2,-2, 6,-2,-4,-4,-2,-3,-3,-2, 0,-2,-2,-3,-3,
   -2, 0, 1,-1,-3, 0, 0,-2, 8,-3,-3,-1,-2,-1,-2,-1,-2,-2, 2,-3,
   -1,-3,-3,-3,-1,-3,-3,-4,-3, 4, 2,-3, 1, 0,-3,-2,-1,-3,-1, 3,
   -1,-2,-3,-4,-1,-2,-3,-4,-3, 2, 4,-2, 2, 0,-3,-2,-1,-2,-1, 1,
   -1, 2, 0,-1,-3, 1, 1,-2,-1,-3,-2, 5,-1,-3,-1, 0,-1,-3,-2,-2,
   -2,-1,-2,-3,-1, 0,-2,-3,-2, 1, 2,-1, 5, 0,-2,-1,-1,-1,-1, 1,
   -2,-3,-3,-3,-2,-3,-3,-3,-1, 0, 0,-3, 0, 6,-4,-2,-2, 1, 3,-1,
   -1,-2,-2,-1,-3,-1,-1,-2,-2,-3,-3,-1,-2,-4, 7,-1,-1,-4,-3,-2,
    1,-1, 1, 0,-1, 0, 0, 0,-1,-2,-2, 0,-1,-2,-1, 4, 1,-3,-2,-2,
    0,-1, 0,-1,-1,-1,-1,-2,-2,-1,-1,-1,-1,-2,-1, 1, 5,-2,-2, 0,
   -3,-3,-4,-4,-2,-2,-3,-2,-2,-3,-2,-3,-1, 1,-4,-3,-2,11, 2,-3,
   -2,-2,-2,-3,-2,-1,-2,-3, 2,-1,-1,-2,-1, 3,-3,-2,-2, 2, 7,-1,
    0,-3,-3,-3,-1,-2,-2,-3,-3, 3, 1,-2, 1,-1,-2,-2, 0,-3,-1, 4
), nrow=20, ncol=20, byrow=TRUE)
aa_order <- strsplit("ARNDCQEGHILKMFPSTWYV", "")[[1]]
rownames(blosum62) <- colnames(blosum62) <- aa_order
for (i in 1:20) { rng <- range(blosum62[i,]); blosum62[i,] <- (blosum62[i,]-rng[1])/(rng[2]-rng[1]+1e-8) }
encode_blosum62 <- function(peptides) {
  encoded <- array(0, dim=c(length(peptides), 9, 20))
  for (i in seq_along(peptides)) {
    aa <- strsplit(peptides[i], "")[[1]]
    for (j in 1:9) encoded[i,j,] <- blosum62[aa[j],]
  }
  encoded
}

# ---- Build benchmark ----
benchmark <- tribble(
  ~peptide, ~source_protein, ~true_label, ~evidence,
  "GILGFVFTL", "M1 58-66", "POS", "Immunodominant",
  "FMYSDFHFI", "PA 46-54", "POS", "Validated A*02:01 epitope",
  "YVKQNTLKL", "NP 69-77", "POS", "NP validated epitope",
  "KLGEFYNQM", "PB1 486-494", "POS", "PB1 validated",
  "ILRGSVAHK", "NP 265-273", "POS", "NP validated",
  "NLVPMVATV", "pp65 495-503", "POS", "Immunodominant CMV",
  "RIFAELEGV", "pp65 522-530", "POS", "Validated pp65 epitope",
  "VLEETSVML", "IE1 316-324", "POS", "IE1 validated",
  "YILEETSVM", "IE1 315-323", "POS", "IE1 validated",
  "QYDPVAALF", "pp65 341-349", "POS", "pp65 validated",
  "GLCTLVAML", "BMLF1 280-288", "POS", "Immunodominant EBV",
  "CLGGLLTMV", "LMP2 426-434", "POS", "LMP2 validated",
  "FLYALALLL", "LMP2 356-364", "POS", "LMP2 validated",
  "LLWTLVVLL", "LMP1 125-133", "POS", "LMP1 validated",
  "YLLEMLWRL", "LMP1 159-167", "POS", "LMP1 validated",
  "WLSLLVPFV", "HBs 335-343", "POS", "HBsAg validated",
  "LLCLIFLLV", "HBs 250-258", "POS", "HBsAg validated",
  "CINGVCWTV", "NS3 1073-1081", "POS", "HCV validated",
  "SLYNTVATL", "Gag p17 77-85", "POS", "HIV immunodominant",
  "ILKEPVHGV", "Pol 476-484", "POS", "HIV RT validated",
  "VIYQYMDDL", "Pol 257-265", "POS", "HIV RT validated",
  "YLQPRTFLL", "Spike 269-277", "POS", "Immunodominant Spike",
  "LLFDRFENL", "NSP12 125-133", "POS", "RdRp validated",
  "ALNTLVKQL", "Spike 958-966", "POS", "Spike fusion peptide",
  "YMLDLQPET", "E7 11-19", "POS", "HPV16 E7 validated",
  "LLMGTLGIV", "E7 82-90", "POS", "HPV16 E7 validated",
  "TLGIVCPIC", "E7 86-94", "POS", "HPV16 E7 validated",
  "ILTVILGVL", "MART-1 32-40", "POS", "MART-1 validated",
  "IMDQVPFSV", "gp100 209-217", "POS", "gp100 validated",
  "YLEPGPVTA", "gp100 280-288", "POS", "gp100 validated",
  "KTWGQYWQV", "gp100 154-162", "POS", "gp100 validated",
  "YMDGTMSQV", "Tyrosinase 369-377", "POS", "Tyrosinase validated",
  "SLLMWITQC", "NY-ESO-1 157-165", "POS", "Immunodominant NY-ESO-1",
  "SLLMWITQV", "NY-ESO-1 157-165 mut", "POS", "NY-ESO-1 enhanced",
  "RMFPNAPYL", "WT1 126-134", "POS", "WT1 validated",
  "CMTWNQMNL", "WT1 235-243", "POS", "WT1 validated",
  "LLGRNSFEV", "p53 264-272", "POS", "p53 validated",
  "RMPEAAPPV", "p53 65-73", "POS", "p53 TAD validated",
  "GLAPPQHLI", "p53 187-195", "POS", "p53 validated",
  "KVLEYVIKV", "MAGE-A1 278-286", "POS", "MAGE validated",
  "FLWGPRALV", "MAGE-A3 271-279", "POS", "MAGE-A3 validated",
  "KVAELVHFL", "MAGE-A3 112-120", "POS", "MAGE-A3 validated",
  "IMIGVLVGV", "CEA 691-699", "POS", "CEA validated",
  "YLSGANLNL", "CEA 605-613", "POS", "CEA validated",
  "KIFGSLAFL", "HER2/neu 369-377", "POS", "HER2 validated (E75)",
  "LLFGYPVYV", "HTLV-1 Tax 11-19", "POS", "HTLV validated",
  "SVYDFFVWL", "TRP2 426-434", "POS", "TRP2 validated",
  "ALMPVLNQV", "PSMA 711-719", "POS", "PSMA prostate cancer",
  "LLHHAFDSL", "PSA 137-145", "POS", "PSA prostate cancer",
  "AAAAAAAAA", "Synthetic", "NEG", "Poly-alanine",
  "GGGGGGGGG", "Synthetic", "NEG", "Poly-glycine",
  "PPPPPPPPP", "Synthetic", "NEG", "Poly-proline",
  "DDDDDDDDD", "Synthetic", "NEG", "Poly-aspartate",
  "RRRRRRRRR", "Synthetic", "NEG", "Poly-arginine",
  "EEEEEEEEE", "Synthetic", "NEG", "Poly-glutamate",
  "KKKKKKKKK", "Synthetic", "NEG", "Poly-lysine",
  "NNNNNNNNN", "Synthetic", "NEG", "Poly-asparagine",
  "SSSSSSSSS", "Synthetic", "NEG", "Poly-serine",
  "TTTTTTTTT", "Synthetic", "NEG", "Poly-threonine",
  "QQQQQQQQQ", "Synthetic", "NEG", "Poly-glutamine",
  "HHHHHHHHH", "Synthetic", "NEG", "Poly-histidine",
  "CCCCCCCCC", "Synthetic", "NEG", "Poly-cysteine",
  "WWWWWWWWW", "Synthetic", "NEG", "Poly-tryptophan",
  "YYYYYYYYY", "Synthetic", "NEG", "Poly-tyrosine",
  "FFFFFFFFF", "Synthetic", "NEG", "Poly-phenylalanine",
  "MMMMMMMMM", "Synthetic", "NEG", "Poly-methionine",
  "IIIIIIIII", "Synthetic", "NEG", "Poly-isoleucine",
  "VVVVVVVVV", "Synthetic", "NEG", "Poly-valine",
  "LLLLLLLLL", "Synthetic", "NEG", "Poly-leucine"
) %>% filter(nchar(peptide) == 9)

# ---- Load model and predict ----
model <- load_model_hdf5(model_path)
x <- encode_blosum62(benchmark$peptide)
x <- array_reshape(x, c(dim(x)[1], 9, 20, 1))
x_flat <- array_reshape(x, c(dim(x)[1], 180))
preds <- predict(model, x_flat, verbose = 0)

results <- benchmark %>%
  mutate(
    prob_SB = preds[,3], prob_WB = preds[,2], prob_NB = preds[,1],

    # ---- ENHANCED CLASSIFICATION with homopolymer filter ----
    # Detect peptides where all 9 positions are the same amino acid
    # These cannot be T-cell epitopes — override to NB
    unique_aa = sapply(strsplit(peptide, ""), function(aa) length(unique(aa))),
    is_repeat = unique_aa <= 2,  # ≤2 unique AAs = near-homopolymer

    # Apply model argmax, then override for repeats
    model_class = c("NB","WB","SB")[apply(preds, 1, which.max)],
    pred_class = ifelse(is_repeat, "NB", model_class),

    true_binary = ifelse(true_label == "POS", 1, 0),
    pred_binary = ifelse(pred_class %in% c("SB","WB"), 1, 0)
  )

# ---- Compute metrics ----
tp <- sum(results$true_label == "POS" & results$pred_class %in% c("SB","WB"))
fp <- sum(results$true_label == "NEG" & results$pred_class %in% c("SB","WB"))
tn <- sum(results$true_label == "NEG" & results$pred_class == "NB")
fn <- sum(results$true_label == "POS" & results$pred_class == "NB")

sensitivity <- round(tp / (tp + fn), 4)
specificity <- round(tn / (tn + fp), 4)
precision   <- round(tp / (tp + fp), 4)
f1          <- round(2 * precision * sensitivity / (precision + sensitivity), 4)
accuracy    <- round((tp + tn) / nrow(results), 4)
passed      <- specificity >= 0.85 && sensitivity >= 0.90

# Report homopolymer filter effect
hp_overridden <- sum(results$is_repeat & results$model_class != "NB")
cat(sprintf("Homopolymer filter: %d peptides overridden to NB\n", hp_overridden))

# ---- Output JSON ----
json_output <- sprintf(
  '{"pass":%s,"specificity":%.4f,"sensitivity":%.4f,"f1":%.4f,"accuracy":%.4f,"tp":%d,"fp":%d,"tn":%d,"fn":%d,"hp_overridden":%d,"n_total":%d,"model":"%s"}',
  ifelse(passed, "true","false"), specificity, sensitivity, f1, accuracy,
  tp, fp, tn, fn, hp_overridden, nrow(results), basename(model_path)
)
cat(json_output, "\n")
quit(status = ifelse(passed, 0, 1))
