# Benchmark against known IEDB-validated HLA-A*02:01 epitopes
options(warn = -1)
Sys.setenv(TF_CPP_MIN_LOG_LEVEL = "3")
suppressMessages({library(keras); library(tensorflow); library(tidyverse); library(pROC)})
setwd('D:/Researching/Peptide epitope')

cat("Loading model...\n")
model <- load_model_hdf5("models/FFN_Deep.h5")

# ---- BLOSUM62 encoder ----
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

# ---- IEDB-validated HLA-A*02:01 epitope benchmark set ----
# Sources: IEDB, literature, known immunodominant epitopes
# Positive: experimentally validated T-cell epitopes
# Negative: known non-binders or very unlikely binders

benchmark <- tribble(
  ~peptide,        ~source_protein,         ~true_label, ~evidence,
  # ====== VIRAL: INFLUENZA ======
  "GILGFVFTL",    "M1 58-66",             "POS", "Immunodominant, 100+ studies",
  "FMYSDFHFI",    "PA 46-54",             "POS", "Validated A*02:01 epitope",
  "ILGFVFTLTV",   "M1 59-68",             "POS", "Extended M1 epitope",  # 10mer - will exclude
  "YVKQNTLKL",    "NP 69-77",             "POS", "NP validated epitope",
  "KLGEFYNQM",    "PB1 486-494",          "POS", "PB1 validated",
  "ILRGSVAHK",    "NP 265-273",           "POS", "NP validated (also B*08:01)",

  # ====== VIRAL: CMV ======
  "NLVPMVATV",    "pp65 495-503",         "POS", "Immunodominant CMV",
  "RIFAELEGV",    "pp65 522-530",         "POS", "Validated pp65 epitope",
  "VLEETSVML",    "IE1 316-324",          "POS", "IE1 validated",
  "YILEETSVM",    "IE1 315-323",          "POS", "IE1 validated",
  "QYDPVAALF",    "pp65 341-349",         "POS", "pp65 validated",

  # ====== VIRAL: EBV ======
  "GLCTLVAML",    "BMLF1 280-288",        "POS", "Immunodominant EBV",
  "CLGGLLTMV",    "LMP2 426-434",         "POS", "LMP2 validated",
  "FLYALALLL",    "LMP2 356-364",         "POS", "LMP2 validated",
  "LLWTLVVLL",    "LMP1 125-133",         "POS", "LMP1 validated",
  "YLLEMLWRL",    "LMP1 159-167",         "POS", "LMP1 validated",

  # ====== VIRAL: HBV ======
  "FLPSDFFPSV",   "HBc 18-27",            "POS", "HBV core validated",
  "WLSLLVPFV",    "HBs 335-343",          "POS", "HBsAg validated",
  "LLVPFVQWFV",   "HBs 338-347",          "POS", "HBsAg validated",  # 10mer
  "LLCLIFLLV",    "HBs 250-258",          "POS", "HBsAg validated",

  # ====== VIRAL: HCV ======
  "CINGVCWTV",    "NS3 1073-1081",        "POS", "HCV validated",
  "KLVALGINAV",   "NS3 1406-1415",        "POS", "HCV validated",
  "LLFNILGGWV",   "NS3 1428-1437",        "POS", "HCV validated",  # 10mer

  # ====== VIRAL: HIV ======
  "SLYNTVATL",    "Gag p17 77-85",        "POS", "HIV immunodominant",
  "ILKEPVHGV",    "Pol 476-484",          "POS", "HIV RT validated",
  "VIYQYMDDL",    "Pol 257-265",          "POS", "HIV RT validated",
  "FLGKIWPSHK",   "Gag p24 255-264",      "POS", "HIV Gag",  # 10mer

  # ====== VIRAL: SARS-CoV-2 ======
  "YLQPRTFLL",    "Spike 269-277",        "POS", "Immunodominant Spike",
  "LLFDRFENL",    "NSP12 125-133",        "POS", "RdRp validated",
  "TLACFVLAAV",   "ORF1ab 2238-2247",     "POS", "Validated CoV epitope",
  "ALNTLVKQL",    "Spike 958-966",        "POS", "Spike fusion peptide",

  # ====== VIRAL: HPV ======
  "YMLDLQPET",    "E7 11-19",             "POS", "HPV16 E7 validated",
  "LLMGTLGIV",    "E7 82-90",             "POS", "HPV16 E7 validated",
  "TLGIVCPIC",    "E7 86-94",             "POS", "HPV16 E7 validated",

  # ====== TUMOR: Melanoma ======
  "ELAGIGILTV",   "MART-1 26-35",         "POS", "Clinical anchor-modified variant",
  "ILTVILGVL",    "MART-1 32-40",         "POS", "MART-1 validated",
  "IMDQVPFSV",    "gp100 209-217",        "POS", "gp100 validated",
  "YLEPGPVTA",    "gp100 280-288",        "POS", "gp100 validated",
  "KTWGQYWQV",    "gp100 154-162",        "POS", "gp100 validated",
  "YMDGTMSQV",    "Tyrosinase 369-377",   "POS", "Tyrosinase validated",

  # ====== TUMOR: Cancer/Testis ======
  "SLLMWITQC",    "NY-ESO-1 157-165",     "POS", "Immunodominant NY-ESO-1",
  "SLLMWITQV",    "NY-ESO-1 157-165 mut", "POS", "NY-ESO-1 enhanced variant",

  # ====== TUMOR: WT1/p53 ======
  "RMFPNAPYL",    "WT1 126-134",          "POS", "WT1 validated",
  "CMTWNQMNL",    "WT1 235-243",          "POS", "WT1 validated",
  "LLGRNSFEV",    "p53 264-272",          "POS", "p53 validated",
  "RMPEAAPPV",    "p53 65-73",            "POS", "p53 TAD validated",
  "GLAPPQHLI",    "p53 187-195",          "POS", "p53 validated",

  # ====== TUMOR: MAGE/CEA/HER2 ======
  "KVLEYVIKV",    "MAGE-A1 278-286",      "POS", "MAGE validated",
  "FLWGPRALV",    "MAGE-A3 271-279",      "POS", "MAGE-A3 validated",
  "KVAELVHFL",    "MAGE-A3 112-120",      "POS", "MAGE-A3 validated",
  "IMIGVLVGV",    "CEA 691-699",          "POS", "CEA validated",
  "YLSGANLNL",    "CEA 605-613",          "POS", "CEA validated",
  "KIFGSLAFL",    "HER2/neu 369-377",     "POS", "HER2 validated (E75)",

  # ====== TUMOR: Other ======
  "LLFGYPVYV",    "HTLV-1 Tax 11-19",     "POS", "HTLV validated",
  "VMAGVGSPYV",   "TRP2 180-188",         "POS", "TRP2 melanoma",
  "SVYDFFVWL",    "TRP2 426-434",         "POS", "TRP2 validated",
  "LLFGLALIEV",   "Survivin 96-104",      "POS", "Survivin validated",
  "ELTLGEFLKL",   "Survivin 95-104",      "POS", "Survivin validated",  # 10mer
  "YLQLVFGIEV",   "Telomerase 540-548",   "POS", "hTERT validated",
  "ALMPVLNQV",    "PSMA 711-719",         "POS", "PSMA prostate cancer",
  "LLHHAFDSL",    "PSA 137-145",          "POS", "PSA prostate cancer",
  "FLTPKKLQCV",   "PSA 154-163",          "POS", "PSA validated",  # 10mer
  "VISNDVCAQV",   "PSA 248-257",          "POS", "PSA validated",  # 10mer

  # ====== NEGATIVE CONTROLS (known non-binders) ======
  "AAAAAAAAA",    "Synthetic",            "NEG", "Poly-alanine",
  "GGGGGGGGG",    "Synthetic",            "NEG", "Poly-glycine",
  "PPPPPPPPP",    "Synthetic",            "NEG", "Poly-proline",
  "DDDDDDDDD",    "Synthetic",            "NEG", "Poly-aspartate",
  "RRRRRRRRR",    "Synthetic",            "NEG", "Poly-arginine",
  "EEEEEEEEE",    "Synthetic",            "NEG", "Poly-glutamate",
  "KKKKKKKKK",    "Synthetic",            "NEG", "Poly-lysine",
  "NNNNNNNNN",    "Synthetic",            "NEG", "Poly-asparagine",
  "SSSSSSSSS",    "Synthetic",            "NEG", "Poly-serine",
  "TTTTTTTTT",    "Synthetic",            "NEG", "Poly-threonine",
  "QQQQQQQQQ",    "Synthetic",            "NEG", "Poly-glutamine",
  "HHHHHHHHH",    "Synthetic",            "NEG", "Poly-histidine",
  "CCCCCCCCC",    "Synthetic",            "NEG", "Poly-cysteine",
  "WWWWWWWWW",    "Synthetic",            "NEG", "Poly-tryptophan",
  "YYYYYYYYY",    "Synthetic",            "NEG", "Poly-tyrosine",
  "FFFFFFFFF",    "Synthetic",            "NEG", "Poly-phenylalanine",
  "MMMMMMMMM",    "Synthetic",            "NEG", "Poly-methionine (anchor-only, no peptide context)",
  "IIIIIIIII",    "Synthetic",            "NEG", "Poly-isoleucine",
  "VVVVVVVVV",    "Synthetic",            "NEG", "Poly-valine (anchor-only)",
  "LLLLLLLLL",    "Synthetic",            "NEG", "Controversial: predicted SB but unnatural"
)

# Filter to 9-mers only
benchmark <- benchmark %>% filter(nchar(peptide) == 9)
cat(sprintf("Benchmark set: %d peptides (%d POS, %d NEG)\n",
    nrow(benchmark), sum(benchmark$true_label=="POS"), sum(benchmark$true_label=="NEG")))

# ---- Predict ----
cat("Predicting...\n")
x <- encode_blosum62(benchmark$peptide)
x <- array_reshape(x, c(dim(x)[1], 9, 20, 1))
x_flat <- array_reshape(x, c(dim(x)[1], 180))
preds <- predict(model, x_flat, verbose=0)

results <- benchmark %>%
  mutate(
    prob_SB  = preds[,3],
    prob_WB  = preds[,2],
    prob_NB  = preds[,1],
    pred_class = c("NB","WB","SB")[apply(preds,1,which.max)],
    binding_score = prob_SB + prob_WB * 0.5,
    p2 = substr(peptide,2,2),
    p9 = substr(peptide,9,9),
    # Binary classification for ROC: POS=1, NEG=0
    true_binary = ifelse(true_label == "POS", 1, 0),
    pred_binary = prob_SB + prob_WB  # SB or WB = positive call
  )

# ---- Performance metrics ----
cat("\n")
cat("==============================================================\n")
cat("  IEDB BENCHMARK RESULTS\n")
cat("==============================================================\n\n")

# Confusion matrix at default threshold (pred_class)
tp <- sum(results$true_label == "POS" & results$pred_class %in% c("SB","WB"))
fp <- sum(results$true_label == "NEG" & results$pred_class %in% c("SB","WB"))
tn <- sum(results$true_label == "NEG" & results$pred_class == "NB")
fn <- sum(results$true_label == "POS" & results$pred_class == "NB")

sensitivity <- tp / (tp + fn)
specificity <- tn / (tn + fp)
precision   <- tp / (tp + fp)
f1 <- 2 * precision * sensitivity / (precision + sensitivity)
accuracy <- (tp + tn) / nrow(results)

cat(sprintf("Sensitivity (Recall): %.1f%% (%d/%d)\n", sensitivity*100, tp, tp+fn))
cat(sprintf("Specificity:          %.1f%% (%d/%d)\n", specificity*100, tn, tn+fp))
cat(sprintf("Precision (PPV):      %.1f%%\n", precision*100))
cat(sprintf("F1 Score:             %.3f\n", f1))
cat(sprintf("Accuracy:             %.1f%%\n", accuracy*100))

# ROC AUC
roc_obj <- roc(results$true_binary, results$pred_binary)
cat(sprintf("ROC AUC:              %.3f\n", auc(roc_obj)))

# ---- Per-category breakdown ----
cat("\n--- Sensitivity by Category ---\n")
results %>%
  mutate(category = case_when(
    grepl("Influenza|CMV|EBV|HBV|HCV|HIV|SARS|HPV|HTLV", source_protein) ~ "Viral",
    grepl("MART|gp100|Tyros|TRP|Melanoma", source_protein) ~ "Melanoma",
    grepl("NY-ESO|MAGE|WT1|p53|CEA|HER2|Surviv|hTERT|Telom|PSA|PSMA", source_protein) ~ "Tumor (other)",
    TRUE ~ "Negative Control"
  )) %>%
  filter(true_label == "POS") %>%
  group_by(category) %>%
  summarise(
    n = n(),
    found = sum(pred_class %in% c("SB","WB")),
    sens = sprintf("%.1f%%", 100*found/n),
    .groups = "drop"
  ) %>%
  print(n=10)

# ---- False negatives ----
cat("\n--- False Negatives (missed epitopes) ---\n")
results %>%
  filter(true_label == "POS", pred_class == "NB") %>%
  select(peptide, source_protein, prob_SB, prob_WB, prob_NB, p2, p9) %>%
  mutate(prob_SB=round(prob_SB*100,1), prob_WB=round(prob_WB*100,1)) %>%
  print(n=20, width=120)

# ---- False positives ----
cat("\n--- False Positives (wrongly predicted as binder) ---\n")
results %>%
  filter(true_label == "NEG", pred_class %in% c("SB","WB")) %>%
  select(peptide, source_protein, pred_class, prob_SB, prob_WB, p2, p9) %>%
  mutate(prob_SB=round(prob_SB*100,1), prob_WB=round(prob_WB*100,1)) %>%
  print(n=20, width=120)

# ---- Save ----
write_csv(results, "data/iedb_benchmark_results.csv")
cat(sprintf("\nSaved: data/iedb_benchmark_results.csv (%d peptides)\n", nrow(results)))
