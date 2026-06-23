# Scan real protein sequences for HLA-A*02:01 epitopes
# Uses best trained model (Deep FFN, MHCflurry, 91.9% accuracy)

options(warn = -1)
Sys.setenv(TF_CPP_MIN_LOG_LEVEL = "3")
suppressMessages({
  library(keras)
  library(tensorflow)
  library(tidyverse)
})

setwd('D:/Researching/Peptide epitope')

# ---- 1. Load model ----
cat("Loading Deep FFN model (MHCflurry-trained, 91.9%)...\n")
model <- load_model_hdf5("models/FFN_Deep.h5")
cat(sprintf("Model ready: %s parameters\n", format(count_params(model), big.mark = ",")))

# ---- 2. BLOSUM62 encoder ----
encode_blosum62 <- function(peptides) {
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
  ), nrow = 20, ncol = 20, byrow = TRUE)
  aa_order <- strsplit("ARNDCQEGHILKMFPSTWYV", "")[[1]]
  rownames(blosum62) <- colnames(blosum62) <- aa_order
  for (i in 1:20) {
    rng <- range(blosum62[i, ])
    blosum62[i, ] <- (blosum62[i, ] - rng[1]) / (rng[2] - rng[1] + 1e-8)
  }
  encoded <- array(0, dim = c(length(peptides), 9, 20))
  for (i in seq_along(peptides)) {
    aa <- strsplit(peptides[i], "")[[1]]
    for (j in 1:9) encoded[i, j, ] <- blosum62[aa[j], ]
  }
  encoded
}

# ---- 3. Scan function ----
scan_protein <- function(protein_seq, protein_name, model) {
  # Generate all 9-mer windows
  seq_len <- nchar(protein_seq)
  n_peptides <- seq_len - 8
  peptides <- sapply(1:n_peptides, function(i) substr(protein_seq, i, i + 8))

  cat(sprintf("Scanning %s (%d aa) → %d 9-mer peptides\n",
              protein_name, seq_len, n_peptides))

  # Encode and predict
  x <- encode_blosum62(peptides)
  x <- array_reshape(x, c(dim(x)[1], 9, 20, 1))
  x_flat <- array_reshape(x, c(dim(x)[1], 180))

  preds <- predict(model, x_flat, verbose = 0)

  # Build results
  results <- tibble(
    protein   = protein_name,
    start_pos = 1:n_peptides,
    end_pos   = (1:n_peptides) + 8,
    peptide   = peptides,
    p2        = substr(peptides, 2, 2),
    p9        = substr(peptides, 9, 9),
    prob_SB   = preds[, 3],
    prob_WB   = preds[, 2],
    prob_NB   = preds[, 1],
    pred_class = c("NB", "WB", "SB")[apply(preds, 1, which.max)]
  ) %>%
    # Binding score weighted toward strong binders
    mutate(
      binding_score = prob_SB * 1.0 + prob_WB * 0.5,
      has_canonical_anchors = p2 %in% c("L","M","I","V","A","T","Q") &
                               p9 %in% c("V","L","I","A","M","T")
    ) %>%
    arrange(desc(binding_score))

  return(results)
}

# ---- 4. Target proteins ----
# Protein 1: MART-1 / Melan-A (melanoma antigen)
# Known A*02:01 epitopes:
#   EAAGIGILTV (pos 26-35, natural) — verified binder
#   ELAGIGILTV (pos 26-35, anchor-modified) — enhanced binder
#   ILTVILGVL  (pos 32-40) — verified binder
MART1 <- "MPREDAHFIYGYPKKGHGHSYTTAEEAAGIGILTVILGVLLLIGCWYCRRRNGYRALMDKSLHVGTQCALTRRCPQEGFDHRDSKVSLQEKNCEPVVPNAPPAYEKLSAEQSPPPYSP"

# Protein 2: SARS-CoV-2 Spike RBD (receptor binding domain)
# Known A*02:01 epitopes:
#   YLQPRTFLL (pos 269-277) — verified immunodominant
#   RLFRKSNLK (pos 454-462, polybasic) — reported binder
#   FIAGLIAIV (pos 1220-1228) — fusion peptide region
SPIKE_RBD <- "NITNLCPFGEVFNATRFASVYAWNRKRISNCVADYSVLYNSASFSTFKCYGVSPTKLNDLCFTNVYADSFVIRGDEVRQIAPGQTGKIADYNYKLPDDFTGCVIAWNSNNLDSKVGGNYNYLYRLFRKSNLKPFERDISTEIYQAGSTPCNGVEGFNCYFPLQSYGFQPTNGVGYQPYRVVVLSFELLHAPATVCGPKKSTNLVKNKCVNFNFNGLTGTGVLTESNKKFLPFQQFGRDIDDTTDAVRDPQTLEILDITPCSFGGVSVITPGTNTSNQVAVLYQDVNCTEVPVAIHADQLTPTWRVYSTGSNVFQTRAGCLIGAEHVNNSYECDIPIGAGICASYQTQTNSPRRAR"

# Protein 3: p53 (tumor suppressor)
# Known A*02:01 epitopes:
#   LLGRNSFEV (pos 265-273) — validated
#   RMPEAAPPV (pos 65-73)
#   GLAPPQHLIRV (pos 187-197, but 10mer)
P53 <- "MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGPDEAPRMPEAAPPVAPAPAAPTPAAPAPAPSWPLSSSVPSQKTYQGSYGFRLGFLHSGTAKSVTCTYSPALNKMFCQLAKTCPVQLWVDSTPPPGTRVRAMAIYKQSQHMTEVVRRCPHHERCSDSDGLAPPQHLIRVEGNLRVEYLDDRNTFRHSVVVPYEPPEVGSDCTTIHYNYMCNSSCMGGMNRRPILTIITLEDSSGNLLGRNSFEVRVCACPGRDRRTEEENLRKKGEPHHELPPGSTKRALPNNTSSSPQPKKKPL"

# ---- 5. Run scans ----
cat("\n")
cat("==============================================================\n")
cat("  PROTEIN EPITOPE SCAN — HLA-A*02:01\n")
cat("  Model: Deep FFN (MHCflurry, 91.9% accuracy)\n")
cat("==============================================================\n\n")

results_list <- list()
for (seq_data in list(
  list(seq = MART1,     name = "MART-1/Melan-A"),
  list(seq = SPIKE_RBD, name = "SARS-CoV-2 Spike RBD"),
  list(seq = P53,       name = "p53 Tumor Suppressor")
)) {
  res <- scan_protein(seq_data$seq, seq_data$name, model)
  results_list[[seq_data$name]] <- res
}

all_results <- bind_rows(results_list)

# ---- 6. Report top epitope candidates ----
cat("\n")
cat("==============================================================\n")
cat("  TOP 10 EPITOPE CANDIDATES PER PROTEIN\n")
cat("==============================================================\n\n")

for (prot in names(results_list)) {
  cat(sprintf("--- %s ---\n", prot))
  res <- results_list[[prot]]
  top <- head(res, 10)

  # Format output
  top %>%
    mutate(
      score = sprintf("%.3f", binding_score),
      `P(SB)` = sprintf("%.1f%%", prob_SB * 100),
      `P(WB)` = sprintf("%.1f%%", prob_WB * 100),
      Anchors = ifelse(has_canonical_anchors, "✓", "—"),
      pos = sprintf("%d-%d", start_pos, end_pos)
    ) %>%
    select(pos, peptide, pred_class, score, `P(SB)`, `P(WB)`, Anchors) %>%
    print(n = 10, width = 100)
  cat("\n")
}

# ---- 7. Known epitope validation ----
cat("==============================================================\n")
cat("  KNOWN EPITOPE VALIDATION\n")
cat("==============================================================\n\n")

known_epitopes <- tribble(
  ~protein,        ~peptide,      ~description,
  "MART-1",        "EAAGIGILTV", "Natural MART-1 26-35 (verified A*02:01 binder)",
  "MART-1",        "ELAGIGILTV", "Anchor-modified MART-1 (enhanced binder, clinical)",
  "MART-1",        "ILTVILGVL",  "MART-1 32-40 (verified binder)",
  "Spike RBD",     "YLQPRTFLL",  "Spike 269-277 (immunodominant A*02:01 epitope)",
  "Spike RBD",     "RLFRKSNLK",  "Spike 454-462 (polybasic, reported binder)",
  "p53",           "LLGRNSFEV",  "p53 265-273 (validated A*02:01 epitope)"
)

cat("Encoding and predicting known epitopes...\n")
x_known <- encode_blosum62(known_epitopes$peptide)
x_known <- array_reshape(x_known, c(dim(x_known)[1], 9, 20, 1))
x_known_flat <- array_reshape(x_known, c(dim(x_known)[1], 180))
preds_known <- predict(model, x_known_flat, verbose = 0)

known_results <- known_epitopes %>%
  mutate(
    prob_SB    = round(preds_known[, 3] * 100, 1),
    prob_WB    = round(preds_known[, 2] * 100, 1),
    prob_NB    = round(preds_known[, 1] * 100, 1),
    prediction = c("NB", "WB", "SB")[apply(preds_known, 1, which.max)],
    pct_rank   = "",
    validated  = ""
  )

# Map to protein scan results to find percentile rank
for (i in 1:nrow(known_results)) {
  prot_res <- results_list[[
    ifelse(grepl("MART", known_results$protein[i]), "MART-1/Melan-A",
    ifelse(grepl("Spike", known_results$protein[i]), "SARS-CoV-2 Spike RBD",
    "p53 Tumor Suppressor"))
  ]]
  # Find position in ranked results
  idx <- which(prot_res$peptide == known_results$peptide[i])
  if (length(idx) > 0) {
    known_results$pct_rank[i] <- sprintf("Top %.1f%%", idx / nrow(prot_res) * 100)
    known_results$validated[i] <- ifelse(
      prot_res$pred_class[idx] %in% c("SB", "WB"), "✓ CONFIRMED", "✗ MISSED")
  }
}

known_results %>%
  select(protein, peptide, prediction, prob_SB, prob_WB, pct_rank, validated, description) %>%
  print(n = 10, width = 130)

# ---- 8. Save results ----
write_csv(all_results, "data/protein_epitope_scan.csv")
cat(sprintf("\nFull scan saved: data/protein_epitope_scan.csv (%d peptides)\n",
            nrow(all_results)))

# Summary stats
sb_count <- sum(all_results$pred_class == "SB")
wb_count <- sum(all_results$pred_class == "WB")
cat(sprintf("\nScan summary (%d total peptides):\n", nrow(all_results)))
cat(sprintf("  Strong Binders:  %d (%.1f%%)\n", sb_count, 100*sb_count/nrow(all_results)))
cat(sprintf("  Weak Binders:    %d (%.1f%%)\n", wb_count, 100*wb_count/nrow(all_results)))
cat(sprintf("  Non-Binders:     %d (%.1f%%)\n",
            nrow(all_results) - sb_count - wb_count,
            100*(nrow(all_results)-sb_count-wb_count)/nrow(all_results)))
