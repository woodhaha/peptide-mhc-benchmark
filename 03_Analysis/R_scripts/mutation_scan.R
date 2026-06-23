# Cancer hotspot mutation scanning -- neoantigen prediction
# Scan p53 and KRAS mutations for epitope-creating/destroying events
options(warn = -1)
Sys.setenv(TF_CPP_MIN_LOG_LEVEL = "3")
suppressMessages({library(keras); library(tensorflow); library(tidyverse)})
setwd('D:/Researching/Peptide epitope')

cat("Loading model...\n")
model <- load_model_hdf5("models/FFN_Deep.h5")

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

predict_binding <- function(peptides, model) {
  x <- encode_blosum62(peptides)
  x <- array_reshape(x, c(dim(x)[1], 9, 20, 1))
  x_flat <- array_reshape(x, c(dim(x)[1], 180))
  preds <- predict(model, x_flat, verbose = 0)
  tibble(
    peptide   = peptides,
    prob_SB   = preds[,3],
    prob_WB   = preds[,2],
    prob_NB   = preds[,1],
    score     = preds[,3] + preds[,2] * 0.5,
    pred      = c("NB","WB","SB")[apply(preds,1,which.max)]
  )
}

# ---- Protein sequences ----
P53_SEQ <- "MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGPDEAPRMPEAAPPVAPAPAAPTPAAPAPAPSWPLSSSVPSQKTYQGSYGFRLGFLHSGTAKSVTCTYSPALNKMFCQLAKTCPVQLWVDSTPPPGTRVRAMAIYKQSQHMTEVVRRCPHHERCSDSDGLAPPQHLIRVEGNLRVEYLDDRNTFRHSVVVPYEPPEVGSDCTTIHYNYMCNSSCMGGMNRRPILTIITLEDSSGNLLGRNSFEVRVCACPGRDRRTEEENLRKKGEPHHELPPGSTKRALPNNTSSSPQPKKKPL"

KRAS_SEQ <- "MTEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFLCVFAINNTKSFEDIHHYREQIKRVKDSEDVPMVLVGNKCDLPSRTVDTKQAQDLARSYGIPFIETSAKTRQGVDDAFYTLVREIRKHKEKMSKDGKKKKKKSKTKCVIM"

get_window_peptides <- function(seq, pos, window = 8) {
  # Returns 9 9-mers spanning pos (pos-8 to pos, each windowed)
  starts <- (pos - window):pos
  starts <- starts[starts >= 1 & (starts + 8) <= nchar(seq)]
  sapply(starts, function(s) substr(seq, s, s + 8))
}

# ---- Hotspot mutation definitions ----
# Format: position, wt_aa, mut_aa, description
p53_mutations <- tribble(
  ~pos, ~wt, ~mut, ~description,
  175,  "R", "H", "p53 R175H -- structural mutant, LFS hotspot",
  245,  "G", "S", "p53 G245S -- DNA contact, Li-Fraumeni hotspot",
  248,  "R", "W", "p53 R248W -- DNA contact, most frequent",
  249,  "R", "S", "p53 R249S -- aflatoxin-associated, HCC",
  273,  "R", "H", "p53 R273H -- DNA contact, second most frequent",
  282,  "R", "W", "p53 R282W -- DNA contact, LFS hotspot",
  220,  "Y", "C", "p53 Y220C -- structural cavity mutant"
)

kras_mutations <- tribble(
  ~pos, ~wt, ~mut, ~description,
  12,   "G", "D", "KRAS G12D -- pancreatic/colon, most common",
  12,   "G", "V", "KRAS G12V -- pancreatic/lung",
  12,   "G", "C", "KRAS G12C -- lung adenocarcinoma (druggable)",
  12,   "G", "R", "KRAS G12R -- pancreatic",
  13,   "G", "D", "KRAS G13D -- colon",
  61,   "Q", "H", "KRAS Q61H -- melanoma",
  61,   "Q", "L", "KRAS Q61L -- NSCLC",
  61,   "Q", "R", "KRAS Q61R -- pancreatic",
  146,  "A", "T", "KRAS A146T -- colon"
)

# ---- Run mutation scan ----
cat("\n")
cat("==============================================================\n")
cat("  CANCER HOTSPOT MUTATION SCAN\n")
cat("  Neoantigen Epitope Prediction -- HLA-A*02:01\n")
cat("==============================================================\n\n")

all_results <- list()

for (protein_name in c("p53", "KRAS")) {
  seq <- if (protein_name == "p53") P53_SEQ else KRAS_SEQ
  mutations <- if (protein_name == "p53") p53_mutations else kras_mutations

  cat(sprintf("=== %s (%d mutations) ===\n\n", protein_name, nrow(mutations)))

  for (i in 1:nrow(mutations)) {
    m <- mutations[i,]
    # Verify WT residue
    actual_wt <- substr(seq, m$pos, m$pos)
    if (actual_wt != m$wt) {
      cat(sprintf("  SKIP %s: expected %s at pos %d, found %s\n",
                  m$description, m$wt, m$pos, actual_wt))
      next
    }

    # Get window peptides for WT
    wt_peps <- get_window_peptides(seq, m$pos)
    # Create mutant sequence and get window peptides
    mut_seq <- seq
    substr(mut_seq, m$pos, m$pos) <- m$mut
    mut_peps <- get_window_peptides(mut_seq, m$pos)

    # Predict for both
    wt_pred  <- predict_binding(wt_peps, model)
    mut_pred <- predict_binding(mut_peps, model)

    # Compare
    comparison <- wt_pred %>%
      rename_with(~paste0("wt_", .x), -peptide) %>%
      bind_cols(mut_pred %>% select(-peptide) %>% rename_with(~paste0("mut_", .x))) %>%
      mutate(
        peptide     = wt_peps,
        wt_peptide  = wt_peps,
        mut_peptide = mut_peps,
        delta_score = mut_score - wt_score,
        effect = case_when(
          wt_pred == "NB" & mut_pred %in% c("SB","WB") ~ "CREATED (neoepitope)",
          wt_pred %in% c("SB","WB") & mut_pred == "NB" ~ "DESTROYED",
          wt_pred %in% c("SB","WB") & mut_score > wt_score + 0.1 ~ "ENHANCED",
          mut_pred %in% c("SB","WB") & mut_score < wt_score - 0.1 ~ "WEAKENED",
          TRUE ~ "unchanged"
        ),
        pos_rel = paste0("p", 1:9),
        mutation = m$description
      )

    # Keep only informative peptides (where mutation affects the window)
    comparison <- comparison %>%
      mutate(
        mut_pos_in_peptide = m$pos - (which(comparison$pos_rel == "p1"):(which(comparison$pos_rel == "p1")+8))[1] + 1
      )

    all_results[[paste0(protein_name, "_", m$pos, m$wt, ">", m$mut)]] <- comparison

    # Print significant changes
    hits <- comparison %>% filter(effect != "unchanged")
    if (nrow(hits) > 0) {
      cat(sprintf("  %s (pos %d: %s>%s)\n", m$description, m$pos, m$wt, m$mut))
      for (j in 1:nrow(hits)) {
        h <- hits[j,]
        cat(sprintf("    %s: %s > %s  |  score %.3f > %.3f  |  %s\n",
                    h$wt_peptide, h$wt_pred, h$mut_pred,
                    h$wt_score, h$mut_score, h$effect))
      }
      cat("\n")
    } else {
      cat(sprintf("  %s: no epitope change\n", m$description))
    }
  }
  cat("\n")
}

# ---- Combine and save ----
full_results <- bind_rows(all_results)
write_csv(full_results, "data/mutation_scan_results.csv")

# ---- Summary ----
cat("==============================================================\n")
cat("  MUTATION SCAN SUMMARY\n")
cat("==============================================================\n\n")

summary <- full_results %>%
  filter(effect != "unchanged") %>%
  select(mutation, wt_peptide, mut_peptide, wt_pred, mut_pred,
         wt_score, mut_score, delta_score, effect, pos_rel)

cat("Significant epitope-altering mutations:\n\n")

summary %>%
  mutate(
    wt_score  = round(wt_score, 3),
    mut_score = round(mut_score, 3),
    delta     = round(delta_score, 3)
  ) %>%
  arrange(desc(abs(delta_score))) %>%
  select(mutation, wt_peptide, mut_peptide, wt_pred, mut_pred, delta, effect) %>%
  print(n = 50, width = 130)

# Neoantigen candidates (created epitopes)
neo <- summary %>% filter(effect == "CREATED (neoepitope)")
cat(sprintf("\nNeoepitopes created: %d\n", nrow(neo)))
if (nrow(neo) > 0) {
  neo %>% select(mutation, mut_peptide, mut_pred, mut_score) %>% print(n=20)
}

# Destroyed epitopes
destroyed <- summary %>% filter(effect == "DESTROYED")
cat(sprintf("\nEpitopes destroyed: %d\n", nrow(destroyed)))
if (nrow(destroyed) > 0) {
  destroyed %>% select(mutation, wt_peptide, wt_pred, wt_score) %>% print(n=20)
}

cat(sprintf("\nFull results saved: data/mutation_scan_results.csv (%d rows)\n",
            nrow(full_results)))
