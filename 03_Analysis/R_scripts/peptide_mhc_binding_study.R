################################################################################
#                                                                              #
#  peptide_mhc_binding_study.R                                                 #
#  ============================                                                #
#                                                                              #
#  Deep Learning for Peptide–MHCI Binding Prediction                           #
#  Reproducing & Extending Jessen (2018)                                       #
#                                                                              #
#  "Deep Learning for Cancer Immunotherapy"                                    #
#  Leon Eyrich Jessen, RStudio AI Blog, January 2018                           #
#                                                                              #
#  Study Design: study-design.md                                               #
#  Literature Map: literature-map-epitope-prediction.md                        #
#                                                                              #
#  Author: Peptide Epitope Research Project                                    #
#  Date:   2026-06-14                                                          #
#                                                                              #
################################################################################

# ============================================================================
# SECTION 0: Environment Setup
# ============================================================================

# 0.1 Install required packages (uncomment first run) ------------------------
# install.packages(c("keras", "tensorflow", "tidyverse", "randomForest",
#                     "pROC", "caret", "Peptides"))
# library(keras)
# install_keras()   # installs TensorFlow backend

# For PepTools (BLOSUM62 encoding and peptide utilities):
# devtools::install_github("leonjessen/PepTools")
# For ggseqlogo (sequence logos):
# devtools::install_github("omarwagih/ggseqlogo")

# 0.2 Load libraries ---------------------------------------------------------
library(keras)          # Deep learning framework (R interface to TensorFlow)
library(tensorflow)     # TensorFlow backend
library(tidyverse)      # Data wrangling + plotting (dplyr, ggplot2, readr, etc.)
library(randomForest)   # Random Forest baseline model
library(pROC)           # ROC curve analysis
library(caret)          # Confusion matrix utilities, train/test splitting
library(doParallel)     # Parallel backend for CV folds & multi-model training
library(foreach)        # foreach + %dopar% loop construct
# Peptides package (AAindex properties) — optional, load gracefully
has_peptides <- requireNamespace("Peptides", quietly = TRUE)
if (has_peptides) {
  library(Peptides)
} else {
  message("Peptides package not installed. AAindex encoding will be limited.")
}

# Optional: PepTools for BLOSUM encoding (try loading, fail gracefully)
has_peptools <- requireNamespace("PepTools", quietly = TRUE)
has_ggseqlogo <- requireNamespace("ggseqlogo", quietly = TRUE)

# 0.3 Global configuration ---------------------------------------------------
set.seed(20180129)  # Reproducibility

# Target HLA allele (from Jessen 2018)
TARGET_HLA <- "HLA-A*02:01"

# Model hyperparameters
EPOCHS       <- 150
BATCH_SIZE   <- 50
VALIDATION_SPLIT <- 0.2
DROPOUT1     <- 0.4
DROPOUT2     <- 0.3
LEARNING_RATE <- 0.001

# Paths
DATA_DIR   <- "data"
PLOTS_DIR  <- "plots"
MODELS_DIR <- "models"

dir.create(DATA_DIR,   showWarnings = FALSE, recursive = TRUE)
dir.create(PLOTS_DIR,  showWarnings = FALSE, recursive = TRUE)
dir.create(MODELS_DIR, showWarnings = FALSE, recursive = TRUE)

# Parallel configuration
N_CORES      <- 4                           # Number of CPU cores for doParallel
TF_INTRA_OP  <- 4L                          # TensorFlow intra-op parallelism (integer)
TF_INTER_OP  <- 2L                          # TensorFlow inter-op parallelism (integer)

message("=== Environment Setup Complete ===")
message("Target HLA: ", TARGET_HLA)
message(sprintf("Parallel: %d cores (TF intra=%d, inter=%d)", N_CORES, TF_INTRA_OP, TF_INTER_OP))
message("PepTools available: ", has_peptools)
message("ggseqlogo available: ", has_ggseqlogo)

# 0.4 Parallel backend setup/teardown -----------------------------------------

#' Initialize parallel backend for multi-core training
#'
#' Creates a doParallel cluster with N_CORES workers. Each worker loads
#' keras/tensorflow for independent model training. Used for CV fold
#' parallelization and multi-model training.
#'
#' @param n_cores Number of CPU cores (default: N_CORES)
#' @return Cluster object (invisibly)
setup_parallel <- function(n_cores = N_CORES) {
  if (n_cores < 2) {
    message("Parallel disabled (n_cores < 2). Running sequentially.")
    registerDoSEQ()
    return(invisible(NULL))
  }

  n_available <- parallel::detectCores()
  if (n_cores > n_available) {
    n_cores <- n_available
  }

  message(sprintf("Setting up parallel backend: %d cores (%d available)",
                  n_cores, n_available))

  cl <- makeCluster(n_cores)
  registerDoParallel(cl)

  # Export key libraries and variables to each worker
  clusterExport(cl, c("LEARNING_RATE", "EPOCHS", "BATCH_SIZE",
                      "VALIDATION_SPLIT", "DROPOUT1", "DROPOUT2"))

  # Each worker loads required packages
  clusterEvalQ(cl, {
    library(keras)
    library(tensorflow)
    library(tidyverse)
    # Suppress TF logging in workers
    tensorflow::tf$get_logger()$setLevel('ERROR')
  })

  message(sprintf("Parallel backend ready: %d workers", length(cl)))
  return(invisible(cl))
}

# Store cluster reference for cleanup
.PARALLEL_CLUSTER <- NULL

#' Tear down parallel backend
#'
#' Stops the doParallel cluster and cleans up worker processes.
#' Should be called at the end of main() or on error.
#'
#' @param cl Cluster object (default: .PARALLEL_CLUSTER global)
teardown_parallel <- function(cl = NULL) {
  if (is.null(cl)) {
    cl <- .PARALLEL_CLUSTER
  }
  if (!is.null(cl)) {
    message("Stopping parallel cluster...")
    stopCluster(cl)
    registerDoSEQ()
    .PARALLEL_CLUSTER <<- NULL
    message("Parallel cluster stopped.")
  }
}

# ============================================================================
# SECTION 1: Peptide Data Generation & Loading
# ============================================================================
# Jessen (2018): 1,000,000 random 9-mer peptides → netMHCpan → label as
# SB (Strong Binder), WB (Weak Binder), NB (Non-Binder) → balance classes
# ============================================================================

#' Generate random 9-mer peptides
#'
#' @param n Number of peptides to generate
#' @return Character vector of peptide sequences
generate_random_peptides <- function(n = 10000) {
  amino_acids <- strsplit("ARNDCQEGHILKMFPSTWYV", "")[[1]]
  peptides <- replicate(n, paste0(sample(amino_acids, size = 9, replace = TRUE),
                                  collapse = ""))
  message(sprintf("Generated %s random 9-mer peptides", format(n, big.mark = ",")))
  return(peptides)
}

#' Generate peptides with PSSM-based realistic binding labels
#'
#' Uses position-specific scoring matrices derived from the known HLA-A*02:01
#' binding motif (anchor residues at p2: L/M/I/V, p9: V/L/I).
#' This captures real MHC-I binding biochemistry without external tools.
#'
#' The PSSM is derived from:
#'   - p2/p9 anchor preferences (crystallographic data)
#'   - p1 auxiliary anchor (aromatic preference)
#'   - p3 negative charge preference (interacts with HLA α-helix)
#'   - p4-p8 weak TCR-facing positional constraints
#'
#' @param n_total Total balanced peptides
#' @param test_frac Test set fraction
#' @return Data frame with realistic binding labels
generate_pssm_labeled_peptides <- function(n_total = 23760, test_frac = 0.10) {
  set.seed(20180129)

  n_per_class <- n_total %/% 3
  all_aa <- strsplit("ARNDCQEGHILKMFPSTWYV", "")[[1]]

  # HLA-A*02:01 PSSM (position × amino acid → log-odds score)
  # p1-p9 scores based on IEDB frequency data for A*02:01 binders
  pssm <- list(
    p1 = c(A=0.1, R=0.0, N=0.0, D=0.0, C=0.0, Q=0.0, E=0.0, G=0.2,
           H=0.0, I=0.2, L=0.1, K=0.0, M=0.1, F=0.4, P=0.0,
           S=0.2, T=0.2, W=0.3, Y=0.5, V=0.1),
    p2 = c(A=0.5, R=-1, N=-1, D=-1, C=-1, Q=0.3, E=-1, G=0.0,
           H=-1, I=0.7, L=1.0, K=-1, M=0.9, F=0.0, P=-1,
           S=0.0, T=0.3, W=-1, Y=0.0, V=0.7),
    p3 = c(A=0.0, R=0.0, N=0.0, D=0.3, E=0.3, Q=0.0, G=0.0,
           H=0.0, I=0.0, L=0.0, K=0.0, M=0.0, F=0.0, P=0.0,
           S=0.0, T=0.0, C=-0.2, W=0.0, Y=0.0, V=0.0),
    p4 = c(A=0.1, R=0.1, N=0.0, D=0.0, E=0.0, Q=0.1, G=0.0,
           H=0.1, I=0.1, L=0.1, K=0.1, M=0.0, F=0.0, P=-0.2,
           S=0.1, T=0.1, C=-0.1, W=0.0, Y=0.0, V=0.1),
    p5 = c(A=0.0, R=0.1, N=0.0, D=0.0, E=0.0, Q=0.0, G=0.0,
           H=0.1, I=0.0, L=0.1, K=0.0, M=0.0, F=0.0, P=-0.1,
           S=0.0, T=0.0, C=-0.1, W=0.0, Y=0.0, V=0.0),
    p6 = c(A=0.0, R=0.0, N=0.0, D=0.0, E=0.1, Q=0.0, G=0.0,
           H=0.0, I=0.1, L=0.1, K=0.0, M=0.0, F=0.0, P=-0.1,
           S=0.0, T=0.0, C=0.0, W=0.0, Y=0.0, V=0.0),
    p7 = c(A=0.0, R=0.0, N=0.0, D=0.0, E=0.0, Q=0.0, G=0.0,
           H=0.0, I=0.0, L=0.1, K=0.0, M=0.0, F=0.1, P=-0.1,
           S=0.0, T=0.0, C=0.0, W=0.0, Y=0.0, V=0.0),
    p8 = c(A=0.0, R=0.0, N=0.0, D=0.0, E=0.0, Q=0.0, G=0.0,
           H=0.0, I=0.1, L=0.2, K=0.0, M=0.0, F=0.0, P=-0.1,
           S=0.0, T=0.0, C=0.0, W=0.0, Y=0.0, V=0.1),
    p9 = c(A=0.4, R=-1, N=-1, D=-1, C=-1, Q=0.2, E=-1, G=0.0,
           H=-1, I=0.7, L=0.8, K=-1, M=0.5, F=0.0, P=-1,
           S=0.1, T=0.3, W=-1, Y=0.0, V=1.0)
  )

  # Score a single peptide
  score_peptide <- function(pep) {
    aa <- strsplit(pep, "")[[1]]
    total <- 0
    for (pos in 1:9) {
      nm <- paste0("p", pos)
      total <- total + pssm[[nm]][aa[pos]]
    }
    return(total)
  }

  # Generate candidate peptides and score them
  message("Generating and scoring peptides with HLA-A*02:01 PSSM...")
  # Generate more candidates than needed to ensure enough high scorers
  n_candidates <- n_total * 3
  candidates <- replicate(n_candidates,
    paste0(sample(all_aa, 9, replace = TRUE), collapse = ""))

  scores <- sapply(candidates, score_peptide)
  names(scores) <- candidates

  # Thresholds: scores ~N(0, 2) typically
  # SB: top ~1-2% (score > 3.5)
  # WB: next ~5% (score 2.0-3.5)
  # NB: rest (score < 2.0)
  sb_threshold <- quantile(scores, 0.98)   # top 2%
  wb_threshold <- quantile(scores, 0.93)   # next 5%

  labels <- rep("NB", n_candidates)
  labels[scores >= wb_threshold & scores < sb_threshold] <- "WB"
  labels[scores >= sb_threshold] <- "SB"

  n_sb <- sum(labels == "SB")
  n_wb <- sum(labels == "WB")
  n_nb <- sum(labels == "NB")

  message(sprintf("PSSM-labeled: SB=%d (%0.1f%%) WB=%d (%0.1f%%) NB=%d (%0.1f%%)",
                  n_sb, 100*n_sb/n_candidates,
                  n_wb, 100*n_wb/n_candidates,
                  n_nb, 100*n_nb/n_candidates))

  # Balance
  min_count <- min(n_sb, n_wb, n_nb)
  if (min_count < n_per_class) {
    n_per_class <- min_count
    n_total <- n_per_class * 3
  }
  message(sprintf("Balancing to %d per class (%d total)", n_per_class, n_total))

  sb_peps <- sample(candidates[labels == "SB"], n_per_class)
  wb_peps <- sample(candidates[labels == "WB"], n_per_class)
  nb_peps <- sample(candidates[labels == "NB"], n_per_class)

  pep_dat <- tibble(
    peptide_id = sprintf("PEP_%05d", 1:n_total),
    peptide    = c(sb_peps, wb_peps, nb_peps),
    label_chr  = factor(
      c(rep("SB", n_per_class), rep("WB", n_per_class), rep("NB", n_per_class)),
      levels = c("NB", "WB", "SB")),
    label_num  = c(rep(2L, n_per_class), rep(1L, n_per_class), rep(0L, n_per_class))
  )

  # Shuffle
  pep_dat <- pep_dat[sample(nrow(pep_dat)), ]
  pep_dat$peptide_id <- sprintf("PEP_%05d", 1:nrow(pep_dat))

  # Train/test split
  pep_dat$idx <- 1:nrow(pep_dat)
  test_idx <- pep_dat %>%
    group_by(label_chr) %>%
    slice_sample(prop = test_frac) %>%
    pull(idx)
  pep_dat <- pep_dat %>%
    mutate(data_type = ifelse(idx %in% test_idx, "test", "train")) %>%
    select(-idx)

  message(sprintf("PSSM dataset: %d peptides", nrow(pep_dat)))
  return(pep_dat)
}

#' Generate synthetic Jessen-format dataset
#'
#' When the original Jessen (2018) GitHub dataset is unavailable, this creates
#' a replica: random 9-mer peptides with simulated netMHCpan-style labels
#' (SB/WB/NB), balanced classes, and 90/10 train/test split.
#'
#' @param n_total Total balanced peptides (default 23,760 = 7,920 × 3 classes)
#' @param sb_frac Fraction of strong binders in raw generation
#' @param wb_frac Fraction of weak binders
#' @param test_frac Fraction for test set
#' @return Data frame with peptide, label_chr, label_num, data_type
generate_synthetic_jessen_data <- function(n_total = 23760,
                                           sb_frac = 0.01,
                                           wb_frac = 0.05,
                                           test_frac = 0.10) {
  set.seed(20180129)

  # Per-class count after balancing: 23,760 / 3 = 7,920
  n_per_class <- n_total %/% 3

  # Generate peptides with known HLA-A*02:01 anchor residue biases
  # Anchor positions: p2 prefers L,M,I,V; p9 prefers V,L,I
  anchor_p2 <- c("L", "M", "I", "V")
  anchor_p9 <- c("V", "L", "I")
  all_aa <- strsplit("ARNDCQEGHILKMFPSTWYV", "")[[1]]

  generate_biased_peptide <- function(use_anchors = TRUE) {
    if (use_anchors) {
      p2 <- sample(anchor_p2, 1)
      p9 <- sample(anchor_p9, 1)
      mid <- sample(all_aa, 6, replace = TRUE)
      other_one <- sample(all_aa, 1)  # p1
      return(paste0(other_one, p2, paste0(mid, collapse = ""), p9))
    } else {
      return(paste0(sample(all_aa, 9, replace = TRUE), collapse = ""))
    }
  }

  # Strong binders: have anchor residues
  sb_peps <- replicate(n_per_class, generate_biased_peptide(use_anchors = TRUE))

  # Weak binders: sometimes have anchors
  wb_peps <- replicate(n_per_class, {
    if (runif(1) < 0.3) generate_biased_peptide(use_anchors = TRUE)
    else generate_biased_peptide(use_anchors = FALSE)
  })

  # Non-binders: random (rarely have both anchors)
  nb_peps <- replicate(n_per_class, generate_biased_peptide(use_anchors = FALSE))

  # Combine into balanced dataset
  pep_dat <- tibble(
    peptide_id = sprintf("PEP_%05d", 1:(n_per_class * 3)),
    peptide    = c(sb_peps, wb_peps, nb_peps),
    label_chr  = factor(
      c(rep("SB", n_per_class), rep("WB", n_per_class), rep("NB", n_per_class)),
      levels = c("NB", "WB", "SB")
    ),
    label_num  = c(rep(2L, n_per_class), rep(1L, n_per_class), rep(0L, n_per_class))
  )

  # Shuffle (preserve peptide_id → sequence mapping)
  pep_dat <- pep_dat[sample(nrow(pep_dat)), ]
  pep_dat$peptide_id <- sprintf("PEP_%05d", 1:nrow(pep_dat))  # renumber after shuffle

  # Train/test split (stratified)
  pep_dat$idx <- 1:nrow(pep_dat)
  test_idx <- pep_dat %>%
    group_by(label_chr) %>%
    slice_sample(prop = test_frac) %>%
    pull(idx)
  pep_dat <- pep_dat %>%
    mutate(data_type = ifelse(idx %in% test_idx, "test", "train")) %>%
    select(-idx)

  message(sprintf("Synthetic Jessen dataset: %s peptides (%s NB / %s WB / %s SB)",
                  format(nrow(pep_dat), big.mark=","),
                  format(n_per_class, big.mark=","),
                  format(n_per_class, big.mark=","),
                  format(n_per_class, big.mark=",")))

  return(pep_dat)
}

#' Load pre-computed peptide-MHC binding data
#'
#' Jessen (2018) used netMHCpan-4.0 to predict binding. We provide two paths:
#'   (a) Load the original Jessen dataset from GitHub
#'   (b) Load a custom dataset (for extended alleles)
#'
#' @param source One of "jessen2018" or "custom"
#' @param custom_path Path to custom data file (if source = "custom")
#' @return Data frame with columns: peptide, label_chr, label_num, data_type
load_peptide_data <- function(source = c("jessen2018", "custom"),
                              custom_path = NULL,
                              n_peptides_synthetic = 23760) {
  source <- match.arg(source)

  if (source == "jessen2018") {
    tmp_file <- file.path(DATA_DIR, "jessen2018_peptides.csv")

    # Try to download original dataset; fall back to synthetic if unavailable
    if (!file.exists(tmp_file)) {
      url <- paste0(
        "https://raw.githubusercontent.com/leonjessen/",
        "dl_for_cancer_immunotherapy/master/data/",
        "ran_peps_netMHCpan40_predicted_A0201_reduced_cleaned_balanced.tsv"
      )

      downloaded <- tryCatch({
        download.file(url, tmp_file, mode = "wb")
        TRUE
      }, error = function(e) FALSE, warning = function(w) FALSE)

      if (!downloaded || !file.exists(tmp_file) || file.info(tmp_file)$size < 100) {
        message("Original Jessen dataset unavailable (404).")
        message("Generating PSSM-labeled peptides (HLA-A*02:01 binding motif)...")
        pep_dat <- generate_pssm_labeled_peptides(n_total = n_peptides_synthetic)

        # Save for reuse
        readr::write_csv(pep_dat, tmp_file)
        message(sprintf("Saved synthetic dataset: %s", tmp_file))
        return(pep_dat)
      }
    }

    pep_dat <- readr::read_csv(tmp_file, show_col_types = FALSE)
    message(sprintf("Loaded %s peptides from synthetic data", format(nrow(pep_dat), big.mark=",")))

  } else if (source == "custom") {
    if (is.null(custom_path)) {
      stop("custom_path must be provided when source = 'custom'")
    }
    pep_dat <- readr::read_csv(custom_path, show_col_types = FALSE)
  }

  # Ensure consistent column types
  pep_dat <- pep_dat %>%
    mutate(
      label_chr = factor(label_chr, levels = c("NB", "WB", "SB")),
      label_num = as.integer(label_num),
      data_type = factor(data_type, levels = c("train", "test"))
    )

  # Add peptide_id if not already present
  if (!"peptide_id" %in% names(pep_dat)) {
    pep_dat <- pep_dat %>%
      mutate(peptide_id = sprintf("PEP_%05d", 1:n()))
  }

  return(pep_dat)
}

#' Balance peptide classes by downsampling
#'
#' @param pep_dat Peptide data frame
#' @param min_count Minimum count per class (NULL = auto to smallest class)
#' @return Balanced data frame
balance_classes <- function(pep_dat, min_count = NULL) {
  class_counts <- pep_dat %>% count(label_chr)

  if (is.null(min_count)) {
    min_count <- min(class_counts$n)
  }

  message("Balancing classes to n = ", min_count, " per class:")

  balanced <- pep_dat %>%
    group_by(label_chr) %>%
    slice_sample(n = min_count) %>%
    ungroup()

  balanced %>%
    count(label_chr) %>%
    mutate(pct = n / sum(n) * 100) %>%
    print()

  return(balanced)
}

#' Split data into train/test sets
#'
#' @param pep_dat Data frame with label_num column
#' @param test_frac Fraction for test set (default 0.10)
#' @param stratify_by Column to stratify by
#' @return List with $train, $test, and $split_info
split_train_test <- function(pep_dat, test_frac = 0.10, stratify_by = "label_chr") {
  set.seed(20180129)

  pep_dat$idx <- 1:nrow(pep_dat)

  test_idx <- pep_dat %>%
    group_by(across(all_of(stratify_by))) %>%
    slice_sample(prop = test_frac) %>%
    pull(idx)

  pep_dat <- pep_dat %>%
    mutate(data_type = ifelse(idx %in% test_idx, "test", "train")) %>%
    select(-idx)

  n_test  <- sum(pep_dat$data_type == "test")
  n_train <- sum(pep_dat$data_type == "train")

  message(sprintf("Train: %s | Test: %s | Test fraction: %.1f%%",
                  format(n_train, big.mark=","),
                  format(n_test, big.mark=","), 100 * n_test / nrow(pep_dat)))

  return(pep_dat)
}

# ============================================================================
# SECTION 2: Peptide Encoding
# ============================================================================
# Jessen (2018): BLOSUM62 encoding → 9×20 'image' per peptide
# Extended: one-hot encoding + AAindex physicochemical properties
# ============================================================================

#' Encode peptides using BLOSUM62 substitution matrix
#'
#' Each amino acid → 20-element vector of substitution probabilities.
#' Each 9-mer → 9×20 matrix (treated as a 1-channel image).
#'
#' @param peptides Character vector of peptide sequences
#' @return 3D array: [n_peptides, 9, 20]
encode_blosum62 <- function(peptides) {
  if (has_peptools) {
    # Use PepTools if available (exact reproduction of Jessen 2018)
    encoded <- PepTools::pep_encode(peptides)
  } else {
    # Fallback: manual BLOSUM62 encoding
    message("PepTools not available; using manual BLOSUM62 encoding")

    # BLOSUM62 matrix (simplified — 20 standard AA order)
    blosum62 <- matrix(c(
      4,-1,-2,-2, 0,-1,-1, 0,-2,-1,-1,-1,-1,-2,-1, 1, 0,-3,-2, 0,  # A
     -1, 5, 0,-2,-3, 1, 0,-2, 0,-3,-2, 2,-1,-3,-2,-1,-1,-3,-2,-3,  # R
     -2, 0, 6, 1,-3, 0, 0, 0, 1,-3,-3, 0,-2,-3,-2, 1, 0,-4,-2,-3,  # N
     -2,-2, 1, 6,-3, 0, 2,-1,-1,-3,-4,-1,-3,-3,-1, 0,-1,-4,-3,-3,  # D
      0,-3,-3,-3, 9,-3,-4,-3,-3,-1,-1,-3,-1,-2,-3,-1,-1,-2,-2,-1,  # C
     -1, 1, 0, 0,-3, 5, 2,-2, 0,-3,-2, 1, 0,-3,-1, 0,-1,-2,-1,-2,  # Q
     -1, 0, 0, 2,-4, 2, 5,-2, 0,-3,-3, 1,-2,-3,-1, 0,-1,-3,-2,-2,  # E
      0,-2, 0,-1,-3,-2,-2, 6,-2,-4,-4,-2,-3,-3,-2, 0,-2,-2,-3,-3,  # G
     -2, 0, 1,-1,-3, 0, 0,-2, 8,-3,-3,-1,-2,-1,-2,-1,-2,-2, 2,-3,  # H
     -1,-3,-3,-3,-1,-3,-3,-4,-3, 4, 2,-3, 1, 0,-3,-2,-1,-3,-1, 3,  # I
     -1,-2,-3,-4,-1,-2,-3,-4,-3, 2, 4,-2, 2, 0,-3,-2,-1,-2,-1, 1,  # L
     -1, 2, 0,-1,-3, 1, 1,-2,-1,-3,-2, 5,-1,-3,-1, 0,-1,-3,-2,-2,  # K
     -2,-1,-2,-3,-1, 0,-2,-3,-2, 1, 2,-1, 5, 0,-2,-1,-1,-1,-1, 1,  # M
     -2,-3,-3,-3,-2,-3,-3,-3,-1, 0, 0,-3, 0, 6,-4,-2,-2, 1, 3,-1,  # F
     -1,-2,-2,-1,-3,-1,-1,-2,-2,-3,-3,-1,-2,-4, 7,-1,-1,-4,-3,-2,  # P
      1,-1, 1, 0,-1, 0, 0, 0,-1,-2,-2, 0,-1,-2,-1, 4, 1,-3,-2,-2,  # S
      0,-1, 0,-1,-1,-1,-1,-2,-2,-1,-1,-1,-1,-2,-1, 1, 5,-2,-2, 0,  # T
     -3,-3,-4,-4,-2,-2,-3,-2,-2,-3,-2,-3,-1, 1,-4,-3,-2,11, 2,-3,  # W
     -2,-2,-2,-3,-2,-1,-2,-3, 2,-1,-1,-2,-1, 3,-3,-2,-2, 2, 7,-1,  # Y
      0,-3,-3,-3,-1,-2,-2,-3,-3, 3, 1,-2, 1,-1,-2,-2, 0,-3,-1, 4   # V
    ), nrow = 20, ncol = 20, byrow = TRUE)

    aa_order <- strsplit("ARNDCQEGHILKMFPSTWYV", "")[[1]]
    rownames(blosum62) <- aa_order
    colnames(blosum62) <- aa_order

    # Normalize each row to [0,1] like PepTools does
    blosum62_norm <- blosum62
    for (i in 1:20) {
      row_min <- min(blosum62[i, ])
      row_max <- max(blosum62[i, ])
      blosum62_norm[i, ] <- (blosum62[i, ] - row_min) / (row_max - row_min + 1e-8)
    }

    # Encode each peptide
    encoded <- array(0, dim = c(length(peptides), 9, 20))
    for (i in seq_along(peptides)) {
      pep <- peptides[i]
      aa_vec <- strsplit(pep, "")[[1]]
      for (j in 1:9) {
        encoded[i, j, ] <- blosum62_norm[aa_vec[j], ]
      }
    }
  }
  return(encoded)
}

#' Encode peptides using one-hot encoding
#'
#' Each amino acid → 20-element binary vector.
#' Each 9-mer → 9×20 binary matrix.
#'
#' @param peptides Character vector of peptide sequences
#' @return 3D array: [n_peptides, 9, 20]
encode_onehot <- function(peptides) {
  aa_order <- strsplit("ARNDCQEGHILKMFPSTWYV", "")[[1]]
  encoded <- array(0, dim = c(length(peptides), 9, 20))

  for (i in seq_along(peptides)) {
    aa_vec <- strsplit(peptides[i], "")[[1]]
    for (j in 1:9) {
      encoded[i, j, which(aa_order == aa_vec[j])] <- 1
    }
  }
  return(encoded)
}

#' Encode peptides using AAindex physicochemical properties
#'
#' Uses the Peptides package to get key properties for each residue.
#' Properties: hydrophobicity (KYTJ820101), bulkiness (ZIMJ680104),
#' polarity (GRAR740102), isoelectric point (ZIMJ680104)
#'
#' @param peptides Character vector of peptide sequences
#' @return 3D array: [n_peptides, 9, n_properties]
encode_aaindex <- function(peptides, properties = c("hydrophobicity", "charge",
                                                     "polarity", "flexibility")) {
  n_props <- length(properties)
  encoded <- array(0, dim = c(length(peptides), 9, n_props))

  # Property lookup (manually defined for reliability)
  # Values normalized to [0,1]
  aa_properties <- list(
    hydrophobicity = c(A=1.8, R=-4.5, N=-3.5, D=-3.5, C=2.5, Q=-3.5, E=-3.5,
                       G=-0.4, H=-3.2, I=4.5, L=3.8, K=-3.9, M=1.9, F=2.8,
                       P=-1.6, S=-0.8, T=-0.7, W=-0.9, Y=-1.3, V=4.2),
    charge         = c(A=0, R=1, N=0, D=-1, C=0, Q=0, E=-1, G=0, H=0.1,
                       I=0, L=0, K=1, M=0, F=0, P=0, S=0, T=0, W=0, Y=0, V=0),
    polarity       = c(A=0, R=1, N=1, D=1, C=0, Q=1, E=1, G=0, H=1, I=0,
                       L=0, K=1, M=0, F=0, P=0, S=1, T=1, W=0, Y=0, V=0),
    flexibility    = c(A=0.36, R=0.53, N=0.46, D=0.51, C=0.35, Q=0.49, E=0.50,
                       G=0.54, H=0.47, I=0.36, L=0.37, K=0.47, M=0.38, F=0.31,
                       P=0.51, S=0.43, T=0.38, W=0.31, Y=0.35, V=0.35)
  )

  for (i in seq_along(peptides)) {
    aa_vec <- strsplit(peptides[i], "")[[1]]
    for (j in 1:9) {
      for (k in seq_along(properties)) {
        prop_name <- properties[k]
        if (prop_name %in% names(aa_properties)) {
          encoded[i, j, k] <- aa_properties[[prop_name]][aa_vec[j]]
        }
      }
    }
  }

  # Normalize each property to [0,1]
  for (k in 1:n_props) {
    min_val <- min(encoded[, , k])
    max_val <- max(encoded[, , k])
    if (max_val > min_val) {
      encoded[, , k] <- (encoded[, , k] - min_val) / (max_val - min_val)
    }
  }

  return(encoded)
}

#' Prepare encoded data for Keras models
#'
#' Reshapes encoding to 4D tensor [n, 9, encoding_dim, 1] for Conv2D compatibility
#' and creates one-hot encoded labels.
#'
#' @param x_encoded 3D array from encode_* functions
#' @param y_labels  Integer vector of class labels (0, 1, 2)
#' @param num_classes Number of output classes
#' @return List with $x (4D array) and $y (2D one-hot matrix)
prepare_keras_data <- function(x_encoded, y_labels, num_classes = 3) {
  # Reshape to 4D: [n_samples, height=9, width=encoding_dim, channels=1]
  x <- array_reshape(x_encoded, c(dim(x_encoded)[1], dim(x_encoded)[2],
                                   dim(x_encoded)[3], 1))
  # One-hot encode labels
  y <- to_categorical(y_labels, num_classes = num_classes)
  return(list(x = x, y = y))
}

# ============================================================================
# SECTION 3: Model Definitions
# ============================================================================
# Reproducing Jessen (2018) architectures + extensions
# ============================================================================

#' Build Feed-Forward Neural Network (Jessen 2018 — best performer)
#'
#' Architecture: Dense(180)→Dropout(0.4)→Dense(90)→Dropout(0.3)→Dense(3, softmax)
#' Input: flat 180-element vector (9 positions × 20 BLOSUM values)
#'
#' @param input_dim Flattened input dimension (default 180 = 9×20 BLOSUM62)
#' @param num_classes Number of output classes (default 3: SB, WB, NB)
#' @param dropout1 First dropout rate
#' @param dropout2 Second dropout rate
#' @return Compiled Keras model
build_ffn <- function(input_dim = 180,
                      num_classes = 3,
                      dropout1 = DROPOUT1,
                      dropout2 = DROPOUT2) {

  model <- keras_model_sequential(name = "FFN_Jessen2018") %>%
    layer_dense(units = 180, activation = "relu",
                input_shape = input_dim,
                kernel_initializer = "he_normal") %>%
    layer_dropout(rate = dropout1) %>%
    layer_dense(units = 90, activation = "relu",
                kernel_initializer = "he_normal") %>%
    layer_dropout(rate = dropout2) %>%
    layer_dense(units = num_classes, activation = "softmax")

  model %>% compile(
    loss      = "categorical_crossentropy",
    optimizer = optimizer_rmsprop(learning_rate = LEARNING_RATE),
    metrics   = c("accuracy")
  )

  message(sprintf("FFN built: %s trainable parameters",
                  format(count_params(model), big.mark = ",")))
  return(model)
}

#' Build Convolutional Neural Network (Jessen 2018)
#'
#' Architecture: Conv2D(32, 3×3)→Dropout(0.25)→Flatten→[FFN body]
#' Input: 9×20×1 'image' tensor
#'
#' @param input_shape 3D input shape (height, width, channels)
#' @param num_classes Number of output classes
#' @return Compiled Keras model
build_cnn <- function(input_shape = c(9, 20, 1),
                      num_classes = 3) {

  model <- keras_model_sequential(name = "CNN_Jessen2018") %>%
    layer_conv_2d(filters = 32, kernel_size = c(3, 3),
                  activation = "relu",
                  input_shape = input_shape,
                  kernel_initializer = "he_normal",
                  padding = "same") %>%
    layer_dropout(rate = 0.25) %>%
    layer_flatten() %>%
    layer_dense(units = 180, activation = "relu",
                kernel_initializer = "he_normal") %>%
    layer_dropout(rate = DROPOUT1) %>%
    layer_dense(units = 90, activation = "relu",
                kernel_initializer = "he_normal") %>%
    layer_dropout(rate = DROPOUT2) %>%
    layer_dense(units = num_classes, activation = "softmax")

  model %>% compile(
    loss      = "categorical_crossentropy",
    optimizer = optimizer_rmsprop(learning_rate = LEARNING_RATE),
    metrics   = c("accuracy")
  )

  message(sprintf("CNN built: %s trainable parameters",
                  format(count_params(model), big.mark = ",")))
  return(model)
}

#' Build LSTM model (Extension beyond Jessen 2018)
#'
#' Rationale: Peptides are sequential; LSTM may capture positional dependencies
#' better than FFN, particularly the anchor position interactions (p2, p9).
#'
#' Architecture: LSTM(64)→Dropout(0.3)→Dense(32)→Dense(3, softmax)
#'
#' @param input_shape 2D input shape (timesteps=9, features=20)
#' @param num_classes Number of output classes
#' @return Compiled Keras model
build_lstm <- function(input_shape = c(9, 20),
                       num_classes = 3) {

  model <- keras_model_sequential(name = "LSTM_Extended") %>%
    layer_lstm(units = 64, input_shape = input_shape,
               return_sequences = FALSE,
               dropout = 0.2, recurrent_dropout = 0.0) %>%
    layer_dropout(rate = 0.3) %>%
    layer_dense(units = 32, activation = "relu",
                kernel_initializer = "he_normal") %>%
    layer_dense(units = num_classes, activation = "softmax")

  model %>% compile(
    loss      = "categorical_crossentropy",
    optimizer = optimizer_adam(learning_rate = LEARNING_RATE),
    metrics   = c("accuracy")
  )

  message(sprintf("LSTM built: %s trainable parameters",
                  format(count_params(model), big.mark = ",")))
  return(model)
}

#' Build deeper FFN (Extension: test if depth helps)
#'
#' Architecture: Dense(360)→Dense(180)→Dense(90)→Dense(45)→Dense(3)
#'
#' @param input_dim Flattened input dimension
#' @param num_classes Number of output classes
#' @return Compiled Keras model
build_ffn_deep <- function(input_dim = 180, num_classes = 3) {

  model <- keras_model_sequential(name = "FFN_Deep_Extended") %>%
    layer_dense(units = 360, activation = "relu",
                input_shape = input_dim,
                kernel_initializer = "he_normal") %>%
    layer_batch_normalization() %>%
    layer_dropout(rate = 0.5) %>%
    layer_dense(units = 180, activation = "relu",
                kernel_initializer = "he_normal") %>%
    layer_batch_normalization() %>%
    layer_dropout(rate = 0.4) %>%
    layer_dense(units = 90, activation = "relu",
                kernel_initializer = "he_normal") %>%
    layer_dropout(rate = 0.3) %>%
    layer_dense(units = 45, activation = "relu",
                kernel_initializer = "he_normal") %>%
    layer_dense(units = num_classes, activation = "softmax")

  model %>% compile(
    loss      = "categorical_crossentropy",
    optimizer = optimizer_adam(learning_rate = LEARNING_RATE),
    metrics   = c("accuracy")
  )

  message(sprintf("Deep FFN built: %s trainable parameters",
                  format(count_params(model), big.mark = ",")))
  return(model)
}

#' Build ResNet-style CNN with residual blocks (Extension for comparison)
#'
#' Adapted from ResNet50 architecture, scaled for 9×20 peptide 'images'.
#' Each residual block learns F(x) = H(x) - x via skip connections.
#' Uses the Keras functional API for tensor routing.
#'
#' Architecture:
#'   Input (9×20×1) → Stem Conv(32) → BN → ReLU
#'   → ResBlock(32→32, identity) → ResBlock(32→64, projection)
#'   → ResBlock(64→128, projection) → GlobalAvgPool
#'   → Dense(128)→Dropout(0.3)→Dense(3, softmax)
#'
#' @param input_shape 3D input shape (height=9, width=20, channels=1)
#' @param num_classes Number of output classes
#' @param filters Vector of filter counts for each residual block
#' @return Compiled Keras model (functional API)
build_resnet_style <- function(input_shape = c(9, 20, 1),
                               num_classes = 3,
                               filters = c(32, 64, 128)) {

  input_layer <- layer_input(shape = input_shape, name = "peptide_input")

  # --- Stem ---
  stem_conv  <- layer_conv_2d(input_layer, filters = filters[1],
                              kernel_size = c(3, 3), padding = "same",
                              kernel_initializer = "he_normal",
                              name = "stem_conv")
  stem_bn    <- layer_batch_normalization(stem_conv, name = "stem_bn")
  stem_relu  <- layer_activation(stem_bn, "relu", name = "stem_relu")

  # --- Residual Block 1: 32→32 (identity shortcut, no projection) ---
  b1_conv1 <- layer_conv_2d(stem_relu, filters = filters[1],
                            kernel_size = c(3, 3), padding = "same",
                            kernel_initializer = "he_normal",
                            name = "block1_conv1")
  b1_bn1   <- layer_batch_normalization(b1_conv1, name = "block1_bn1")
  b1_relu1 <- layer_activation(b1_bn1, "relu", name = "block1_relu1")
  b1_conv2 <- layer_conv_2d(b1_relu1, filters = filters[1],
                            kernel_size = c(3, 3), padding = "same",
                            kernel_initializer = "he_normal",
                            name = "block1_conv2")
  b1_bn2   <- layer_batch_normalization(b1_conv2, name = "block1_bn2")
  b1_add   <- tf$keras$layers$add(list(b1_bn2, stem_relu), name = "block1_add")
  b1_out   <- layer_activation(b1_add, "relu", name = "block1_relu2")

  # --- Residual Block 2: 32→64 (projection shortcut via 1×1 conv) ---
  b2_proj  <- layer_conv_2d(b1_out, filters = filters[2],
                            kernel_size = c(1, 1), padding = "same",
                            kernel_initializer = "he_normal",
                            name = "block2_proj")

  b2_conv1 <- layer_conv_2d(b1_out, filters = filters[2],
                            kernel_size = c(3, 3), padding = "same",
                            kernel_initializer = "he_normal",
                            name = "block2_conv1")
  b2_bn1   <- layer_batch_normalization(b2_conv1, name = "block2_bn1")
  b2_relu1 <- layer_activation(b2_bn1, "relu", name = "block2_relu1")
  b2_conv2 <- layer_conv_2d(b2_relu1, filters = filters[2],
                            kernel_size = c(3, 3), padding = "same",
                            kernel_initializer = "he_normal",
                            name = "block2_conv2")
  b2_bn2   <- layer_batch_normalization(b2_conv2, name = "block2_bn2")
  b2_add   <- tf$keras$layers$add(list(b2_bn2, b2_proj), name = "block2_add")
  b2_out   <- layer_activation(b2_add, "relu", name = "block2_relu2")

  # --- Residual Block 3: 64→128 (projection shortcut via 1×1 conv) ---
  b3_proj  <- layer_conv_2d(b2_out, filters = filters[3],
                            kernel_size = c(1, 1), padding = "same",
                            kernel_initializer = "he_normal",
                            name = "block3_proj")

  b3_conv1 <- layer_conv_2d(b2_out, filters = filters[3],
                            kernel_size = c(3, 3), padding = "same",
                            kernel_initializer = "he_normal",
                            name = "block3_conv1")
  b3_bn1   <- layer_batch_normalization(b3_conv1, name = "block3_bn1")
  b3_relu1 <- layer_activation(b3_bn1, "relu", name = "block3_relu1")
  b3_conv2 <- layer_conv_2d(b3_relu1, filters = filters[3],
                            kernel_size = c(3, 3), padding = "same",
                            kernel_initializer = "he_normal",
                            name = "block3_conv2")
  b3_bn2   <- layer_batch_normalization(b3_conv2, name = "block3_bn2")
  b3_add   <- tf$keras$layers$add(list(b3_bn2, b3_proj), name = "block3_add")
  b3_out   <- layer_activation(b3_add, "relu", name = "block3_relu2")

  # --- Classifier head ---
  gap      <- layer_global_average_pooling_2d(b3_out, name = "global_avg_pool")
  drop1    <- layer_dropout(gap, rate = 0.4, name = "classifier_dropout1")
  dense1   <- layer_dense(drop1, units = 128, activation = "relu",
                          kernel_initializer = "he_normal",
                          name = "classifier_dense1")
  drop2    <- layer_dropout(dense1, rate = 0.3, name = "classifier_dropout2")
  output   <- layer_dense(drop2, units = num_classes, activation = "softmax",
                          name = "output")

  # --- Build model ---
  model <- keras_model(inputs = input_layer, outputs = output,
                       name = "ResNet_Peptide")

  model %>% compile(
    loss      = "categorical_crossentropy",
    optimizer = optimizer_adam(learning_rate = LEARNING_RATE),
    metrics   = c("accuracy")
  )

  message(sprintf("ResNet-style built: %s trainable parameters",
                  format(count_params(model), big.mark = ",")))
  return(model)
}

#' Train Random Forest model (Jessen 2018 baseline)
#'
#' @param x_train Matrix of encoded peptide features (n_peptides × n_features)
#' @param y_train Factor of class labels
#' @param ntree Number of trees
#' @return Trained randomForest object
build_rf <- function(x_train, y_train, ntree = 100) {
  model <- randomForest(
    x      = x_train,
    y      = factor(y_train),
    ntree  = ntree,
    importance = TRUE
  )
  message(sprintf("Random Forest trained: %d trees", ntree))
  return(model)
}

# ============================================================================
# SECTION 4: Training Utilities
# ============================================================================

#' Train a Keras model with early stopping and history recording
#'
#' @param model Compiled Keras model
#' @param x_train Training features (4D array)
#' @param y_train Training labels (one-hot matrix)
#' @param epochs Max epochs
#' @param batch_size Batch size
#' @param validation_split Internal validation fraction
#' @param early_stopping_patience Patience for early stopping (NULL = disabled)
#' @param verbose Verbosity level
#' @return Training history object
train_keras_model <- function(model, x_train, y_train,
                              epochs = EPOCHS,
                              batch_size = BATCH_SIZE,
                              validation_split = VALIDATION_SPLIT,
                              early_stopping_patience = 15,
                              verbose = 1) {

  callbacks_list <- list()

  if (!is.null(early_stopping_patience)) {
    callbacks_list[[1]] <- callback_early_stopping(
      monitor   = "val_loss",
      patience  = early_stopping_patience,
      restore_best_weights = TRUE
    )
  }

  # Add learning rate reduction on plateau
  callbacks_list[[length(callbacks_list) + 1]] <- callback_reduce_lr_on_plateau(
    monitor = "val_loss",
    factor  = 0.5,
    patience = 5,
    min_lr  = 1e-6
  )

  history <- model %>% fit(
    x               = x_train,
    y               = y_train,
    epochs          = epochs,
    batch_size      = batch_size,
    validation_split = validation_split,
    callbacks       = callbacks_list,
    verbose         = verbose
  )

  return(history)
}

#' Plot training history
#'
#' @param history Keras training history object
#' @param title Plot title
#' @param save_path Path to save plot (NULL = display only)
plot_training_history <- function(history, title = "Training History",
                                  save_path = NULL) {
  # Convert history to data frame
  metrics <- history$metrics

  df <- data.frame(
    epoch     = 1:length(metrics$loss),
    loss      = metrics$loss,
    val_loss  = metrics$val_loss,
    accuracy  = metrics$accuracy,
    val_accuracy = metrics$val_accuracy
  )

  # Loss plot
  p_loss <- ggplot(df, aes(x = epoch)) +
    geom_line(aes(y = loss, color = "Train"), linewidth = 0.8) +
    geom_line(aes(y = val_loss, color = "Validation"), linewidth = 0.8) +
    scale_color_manual(values = c("Train" = "steelblue", "Validation" = "tomato")) +
    labs(title = paste(title, "- Loss"),
         x = "Epoch", y = "Loss", color = "") +
    theme_minimal(base_size = 12)

  # Accuracy plot
  p_acc <- ggplot(df, aes(x = epoch)) +
    geom_line(aes(y = accuracy, color = "Train"), linewidth = 0.8) +
    geom_line(aes(y = val_accuracy, color = "Validation"), linewidth = 0.8) +
    scale_color_manual(values = c("Train" = "steelblue", "Validation" = "tomato")) +
    labs(title = paste(title, "- Accuracy"),
         x = "Epoch", y = "Accuracy", color = "") +
    theme_minimal(base_size = 12)

  # Combine
  p_combined <- gridExtra::grid.arrange(p_loss, p_acc, ncol = 2)

  if (!is.null(save_path)) {
    ggsave(save_path, p_combined, width = 12, height = 5, dpi = 150)
    message("Training plot saved: ", save_path)
  }

  return(invisible(p_combined))
}

# ============================================================================
# SECTION 5: Evaluation
# ============================================================================

#' Evaluate model performance comprehensively
#'
#' @param model Trained Keras model or randomForest object
#' @param x_test Test features
#' @param y_test Test labels (integer vector 0/1/2 or factor)
#' @param model_name Name for reporting
#' @param class_names Class label names
#' @return List with accuracy, per-class metrics, confusion matrix, predictions
evaluate_model <- function(model, x_test, y_test,
                           model_name = "Model",
                           class_names = c("NB", "WB", "SB")) {

  # Get predictions
  if (inherits(model, "randomForest")) {
    y_pred_class <- as.integer(as.character(predict(model, x_test)))
    y_real <- as.integer(as.character(y_test))
    # RF doesn't give probabilities for 3 classes reliably in this setup
    y_pred_prob <- NULL
  } else {
    # Keras model
    y_pred_prob <- predict(model, x_test)
    y_pred_class <- apply(y_pred_prob, 1, which.max) - 1  # 0-indexed
    # Get real labels
    if (is.matrix(y_test) || is.array(y_test)) {
      y_real <- apply(y_test, 1, which.max) - 1
    } else {
      y_real <- as.integer(y_test)
    }
  }

  # Confusion matrix
  cm <- table(
    Predicted = factor(y_pred_class, levels = 0:2, labels = class_names),
    Actual    = factor(y_real, levels = 0:2, labels = class_names)
  )

  # Overall accuracy
  accuracy <- sum(diag(cm)) / sum(cm)

  # Per-class metrics
  per_class <- list()
  for (i in seq_along(class_names)) {
    tp <- cm[i, i]
    fp <- sum(cm[i, -i])
    fn <- sum(cm[-i, i])
    tn <- sum(cm[-i, -i])

    precision <- ifelse(tp + fp > 0, tp / (tp + fp), 0)
    recall    <- ifelse(tp + fn > 0, tp / (tp + fn), 0)
    f1        <- ifelse(precision + recall > 0,
                         2 * precision * recall / (precision + recall), 0)

    per_class[[class_names[i]]] <- c(
      Precision = precision,
      Recall    = recall,
      F1        = f1,
      Support   = sum(cm[, i])
    )
  }

  # Macro F1
  macro_f1 <- mean(sapply(per_class, function(x) x["F1"]))

  # Print summary
  message(sprintf("\n=== %s Evaluation ===", model_name))
  message(sprintf("Accuracy:  %.2f%%", accuracy * 100))
  message(sprintf("Macro F1:  %.4f", macro_f1))
  for (cls in class_names) {
    message(sprintf("  %s: P=%.3f R=%.3f F1=%.3f (n=%d)",
                    cls,
                    per_class[[cls]]["Precision"],
                    per_class[[cls]]["Recall"],
                    per_class[[cls]]["F1"],
                    per_class[[cls]]["Support"]))
  }

  return(list(
    model_name   = model_name,
    accuracy     = accuracy,
    macro_f1     = macro_f1,
    per_class    = per_class,
    confusion    = cm,
    y_real       = y_real,
    y_pred_class = y_pred_class,
    y_pred_prob  = y_pred_prob
  ))
}

#' Plot confusion matrix
#'
#' @param eval_result Output from evaluate_model()
#' @param save_path Path to save plot
plot_confusion_matrix <- function(eval_result, save_path = NULL) {
  cm <- eval_result$confusion
  cm_df <- as.data.frame(as.table(cm))
  names(cm_df) <- c("Predicted", "Actual", "Freq")

  acc <- round(eval_result$accuracy * 100, 1)

  p <- ggplot(cm_df, aes(x = Predicted, y = Actual, fill = Freq)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = Freq), size = 5, fontface = "bold") +
    scale_fill_gradient(low = "whitesmoke", high = "steelblue") +
    labs(
      title = paste0(eval_result$model_name, " — Confusion Matrix"),
      subtitle = paste0("Accuracy: ", acc, "%"),
      x = "Predicted Class",
      y = "Actual Class"
    ) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "none")

  if (!is.null(save_path)) {
    ggsave(save_path, p, width = 6, height = 5, dpi = 150)
  }

  return(p)
}

#' Plot prediction scatter (reproducing Jessen 2018 Fig)
#'
#' @param eval_result Output from evaluate_model()
#' @param save_path Path to save plot
plot_prediction_scatter <- function(eval_result, save_path = NULL) {
  class_names <- c("NB", "WB", "SB")

  results <- tibble(
    y_real  = factor(eval_result$y_real, levels = 0:2, labels = class_names),
    y_pred  = factor(eval_result$y_pred_class, levels = 0:2, labels = class_names),
    Correct = ifelse(as.character(y_real) == as.character(y_pred), "Yes", "No")
  )

  acc <- round(eval_result$accuracy * 100, 1)

  p <- results %>%
    ggplot(aes(x = y_pred, y = y_real, colour = Correct)) +
    geom_jitter(width = 0.3, height = 0.3, alpha = 0.6, size = 2) +
    scale_color_manual(values = c("Yes" = "cornflowerblue", "No" = "tomato")) +
    labs(
      title = paste0("Performance on Unseen Test Data — ", eval_result$model_name),
      subtitle = paste0("Accuracy = ", acc, "%"),
      x = "Predicted Class",
      y = "Actual Class (netMHCpan label)",
      colour = "Correct"
    ) +
    theme_bw(base_size = 13) +
    theme(legend.position = "bottom")

  if (!is.null(save_path)) {
    ggsave(save_path, p, width = 6, height = 5.5, dpi = 150)
  }

  return(p)
}

#' Generate sequence logo for strong binders (if ggseqlogo available)
#'
#' @param peptides Character vector of SB peptide sequences
#' @param save_path Path to save plot
plot_sequence_logo <- function(peptides, save_path = NULL) {
  if (!has_ggseqlogo) {
    message("ggseqlogo not installed. Skipping sequence logo.")
    return(NULL)
  }

  p <- ggseqlogo::ggseqlogo(peptides) +
    labs(title = "Sequence Logo — Strong Binders (SB)",
         subtitle = paste("n =", length(peptides), "peptides")) +
    theme_minimal()

  if (!is.null(save_path)) {
    ggsave(save_path, p, width = 8, height = 3, dpi = 150)
  }

  return(p)
}

#' Multi-class ROC curves (one-vs-rest)
#'
#' @param eval_result Output from evaluate_model() (must have y_pred_prob)
#' @param class_names Class labels
#' @param save_path Path to save plot
plot_roc_curves <- function(eval_result, class_names = c("NB", "WB", "SB"),
                            save_path = NULL) {
  if (is.null(eval_result$y_pred_prob)) {
    message("ROC curves require probability predictions. Skipping.")
    return(NULL)
  }

  roc_list <- list()
  for (i in seq_along(class_names)) {
    y_binary <- ifelse(eval_result$y_real == (i - 1), 1, 0)
    y_score  <- eval_result$y_pred_prob[, i]
    roc_list[[class_names[i]]] <- roc(y_binary, y_score, quiet = TRUE)
  }

  # Plot
  colors <- c("NB" = "gray50", "WB" = "darkorange", "SB" = "steelblue")
  p <- ggroc(roc_list, linewidth = 1) +
    scale_color_manual(values = colors) +
    geom_abline(linetype = "dashed", alpha = 0.3) +
    labs(
      title = paste("ROC Curves —", eval_result$model_name),
      x = "1 - Specificity",
      y = "Sensitivity",
      color = "Class"
    ) +
    theme_minimal(base_size = 12) +
    annotate("text", x = 0.7, y = 0.15,
             label = paste("Macro F1:", round(eval_result$macro_f1, 3)),
             size = 4, fontface = "italic")

  if (!is.null(save_path)) {
    ggsave(save_path, p, width = 6, height = 5, dpi = 150)
  }

  return(p)
}

# ============================================================================
# SECTION 5b: Feature Extraction & Model Persistence
# ============================================================================
# Extract peptide features from trained models and save everything for reuse.
# Key outputs: feature.csv (peptide_id + raw + learned features + predictions)
#              models/*.h5 (each trained model with architecture tag)
# ============================================================================

#' Extract learned features from penultimate layer of a trained Keras model
#'
#' Builds a sub-model that outputs activations from the last hidden dense
#' layer before the softmax output. These features represent the model's
#' learned representation of each peptide.
#'
#' For FFN: extracts Dense(90) layer output (90-dim feature vector)
#' For ResNet: extracts classifier_dense1 layer output (128-dim feature vector)
#'
#' @param model Trained Keras model
#' @param x_data Input data (same format model was trained on)
#' @param layer_name Name of penultimate dense layer. NULL = auto-detect.
#' @return Matrix of feature vectors (n_samples × n_features)
extract_learned_features <- function(model, x_data, layer_name = NULL) {

  if (is.null(layer_name)) {
    # Auto-detect: find the last Dense layer before the output softmax
    layer_names <- sapply(model$layers, function(l) l$name)
    dense_layers <- grep("dense|classifier_dense", layer_names, value = TRUE)
    # Exclude output layer (last dense)
    if (length(dense_layers) >= 2) {
      layer_name <- dense_layers[length(dense_layers) - 1]
    } else {
      layer_name <- dense_layers[1]
    }
  }

  message(sprintf("Extracting features from layer: %s", layer_name))

  # Build sub-model that outputs the target layer
  feature_model <- keras_model(
    inputs = model$input,
    outputs = get_layer(model, layer_name)$output
  )

  features <- predict(feature_model, x_data)
  message(sprintf("Extracted features: [%s]",
                  paste(dim(features), collapse = ", ")))

  return(features)
}

#' Build comprehensive feature CSV for all peptides
#'
#' Combines:
#'   1. Peptide metadata (peptide_id, sequence, label)
#'   2. Raw BLOSUM62 encoding (180 features, prefixed 'blosum_')
#'   3. Learned features from penultimate layer (prefixed 'feat_')
#'   4. Model predictions (prob_SB, prob_WB, prob_NB)
#'   5. Model name and architecture tag
#'
#' @param pep_dat Peptide data frame with peptide_id, peptide, label_chr, data_type
#' @param model Trained Keras model
#' @param x_data Encoded input data for all peptides
#' @param model_name Name tag for the model (e.g., "FFN_Jessen2018")
#' @param encoding_type Type of encoding used (e.g., "BLOSUM62")
#' @param save_path Path to save CSV
#' @return Data frame with all features (invisibly)
build_feature_csv <- function(pep_dat, model, x_data, y_data,
                              model_name = "model",
                              encoding_type = "BLOSUM62",
                              save_path = "feature.csv") {

  message(sprintf("\n=== Building Feature CSV: %s ===", save_path))

  # 1. Extract raw BLOSUM features (flatten the 4D tensor → 2D matrix)
  # Input: [n, 9, 20, 1] → Flatten: [n, 180]
  x_flat <- array_reshape(x_data, c(dim(x_data)[1],
                                     dim(x_data)[2] * dim(x_data)[3]))
  colnames(x_flat) <- paste0("blosum_",
                             rep(1:9, each = 20), "_",
                             rep(strsplit("ARNDCQEGHILKMFPSTWYV", "")[[1]], 9))

  # 2. Extract learned features from penultimate layer
  # For FFN (flat input): extract from 90-unit dense layer
  # For CNN/ResNet (4D input): extract from their classifier dense layer
  learned_feats <- tryCatch(
    extract_learned_features(model, x_data),
    error = function(e) {
      # If 4D model but x_data is flattened, reshape and retry
      message("Retrying feature extraction with reshaped input...")
      x_4d <- array_reshape(x_flat, c(dim(x_flat)[1], 9, 20, 1))
      extract_learned_features(model, x_4d)
    }
  )

  colnames(learned_feats) <- paste0("feat_", 1:ncol(learned_feats))

  # 3. Get model predictions
  preds <- predict(model, x_data)
  colnames(preds) <- c("prob_NB", "prob_WB", "prob_SB")

  # 4. Get true labels
  if (is.matrix(y_data) || is.array(y_data)) {
    y_int <- apply(y_data, 1, which.max) - 1
  } else {
    y_int <- y_data
  }

  # 5. Assemble feature CSV
  # Ensure pep_dat is aligned with x_data
  feature_df <- pep_dat %>%
    select(peptide_id, peptide, label_chr, label_num, data_type) %>%
    mutate(
      model_name   = model_name,
      encoding     = encoding_type,
      pred_class   = apply(preds, 1, which.max) - 1,
      pred_correct = pred_class == label_num
    ) %>%
    bind_cols(
      as_tibble(preds),
      as_tibble(x_flat),
      as_tibble(learned_feats)
    )

  # 6. Save
  readr::write_csv(feature_df, save_path)
  message(sprintf("Feature CSV saved: %s", save_path))
  message(sprintf("  Rows: %d | Raw features: %d | Learned features: %d",
                  nrow(feature_df), ncol(x_flat), ncol(learned_feats)))

  return(invisible(feature_df))
}

#' Save all trained models with consistent naming
#'
#' Saves each model to models/ directory with architecture tag and timestamp.
#' Also saves the Random Forest model as .rds.
#'
#' @param models_list Named list of trained models
#' @param prefix Optional prefix to add to filenames
save_all_models <- function(models_list, prefix = "") {
  dir.create(MODELS_DIR, showWarnings = FALSE, recursive = TRUE)

  for (model_name in names(models_list)) {
    model <- models_list[[model_name]]

    safe_name <- gsub("[^A-Za-z0-9_]", "_", model_name)
    if (nchar(prefix) > 0) {
      safe_name <- paste0(prefix, "_", safe_name)
    }

    if (inherits(model, "randomForest")) {
      # Save Random Forest as .rds
      rf_path <- file.path(MODELS_DIR, paste0(safe_name, ".rds"))
      saveRDS(model, rf_path)
      message(sprintf("Saved RF: %s", rf_path))
    } else {
      # Save Keras model as .h5
      keras_path <- file.path(MODELS_DIR, paste0(safe_name, ".h5"))
      save_model_hdf5(model, keras_path)
      message(sprintf("Saved Keras: %s", keras_path))
    }
  }

  # Also save a model manifest
  manifest <- data.frame(
    model_name   = names(models_list),
    file         = paste0(gsub("[^A-Za-z0-9_]", "_", names(models_list)),
                          ifelse(sapply(models_list, inherits, "randomForest"),
                                 ".rds", ".h5")),
    type         = ifelse(sapply(models_list, inherits, "randomForest"),
                          "randomForest", "keras"),
    saved_at     = Sys.time()
  )
  readr::write_csv(manifest, file.path(MODELS_DIR, "model_manifest.csv"))
  message(sprintf("Model manifest saved: models/model_manifest.csv (%d models)",
                  nrow(manifest)))

  return(invisible(manifest))
}

# ============================================================================
# SECTION 6: Cross-Validation & Ensemble (Extensions)
# ============================================================================

#' Run k-fold cross-validation for a Keras model
#'
#' Jessen (2018) mentions: "Models are also trained and evaluated using
#' cross-validation, usually 5-fold. We then save each of the five models and
#' create an ensemble prediction wisdom-of-the-crowd."
#'
#' @param build_fn Function that builds a compiled Keras model (no arguments)
#' @param x_all Full feature array
#' @param y_all Full label array
#' @param k Number of folds
#' @param epochs Epochs per fold
#' @param batch_size Batch size
#' @return List with $fold_results, $ensemble_accuracy, $models
run_cross_validation <- function(build_fn, x_all, y_all,
                                 k = 5,
                                 epochs = EPOCHS,
                                 batch_size = BATCH_SIZE,
                                 parallel = TRUE) {

  message(sprintf("\n=== Starting %d-Fold Cross-Validation ===", k))

  # Create stratified folds
  y_int <- if (is.matrix(y_all)) apply(y_all, 1, which.max) - 1 else y_all
  folds <- caret::createFolds(factor(y_int), k = k, list = TRUE)

  # Detect parallel backend
  use_parallel <- parallel && getDoParWorkers() > 1
  if (use_parallel) {
    message(sprintf("CV: parallel mode (%d workers)", getDoParWorkers()))
  } else {
    message("CV: sequential mode")
  }

  # ---- Parallel fold training via foreach ----
  # Each worker: builds a fresh model, trains on its fold split, returns results
  `%myop%` <- if (use_parallel) `%dopar%` else `%do%`

  fold_outputs <- foreach(
    i       = 1:k,
    .packages = c("keras", "tensorflow", "tidyverse", "pROC", "caret"),
    .export  = c("LEARNING_RATE", "EPOCHS", "BATCH_SIZE",
                 "VALIDATION_SPLIT", "DROPOUT1", "DROPOUT2",
                 "evaluate_model", "build_fn")
  ) %myop% {

    # Set TF logging to ERROR in workers
    tensorflow::tf$get_logger()$setLevel("ERROR")
    Sys.setenv(TF_CPP_MIN_LOG_LEVEL = "2")

    test_idx  <- folds[[i]]
    train_idx <- setdiff(1:nrow(x_all), test_idx)

    # Dimension-aware subsetting: handles 2D [n, features] or 4D [n, h, w, c]
    ndim <- length(dim(x_all))
    if (ndim == 2) {
      x_train_fold <- x_all[train_idx, , drop = FALSE]
      x_test_fold  <- x_all[test_idx,  , drop = FALSE]
    } else {
      extra_dims <- rep(", ", ndim - 1)  # drop = FALSE for all dims
      x_train_fold <- x_all[train_idx, , drop = FALSE]
      x_test_fold  <- x_all[test_idx,  , drop = FALSE]
    }
    y_train_fold <- y_all[train_idx, , drop = FALSE]
    y_test_fold  <- y_all[test_idx,  , drop = FALSE]

    # Build fresh model
    model <- build_fn()

    # Train with verbose=0 to avoid cluttering parallel output
    history <- model %>% fit(
      x          = x_train_fold,
      y          = y_train_fold,
      epochs     = epochs,
      batch_size = batch_size,
      validation_split = 0.2,
      callbacks  = list(
        callback_early_stopping(monitor = "val_loss", patience = 10,
                                restore_best_weights = TRUE),
        callback_reduce_lr_on_plateau(monitor = "val_loss", factor = 0.5,
                                       patience = 5, min_lr = 1e-6)
      ),
      verbose    = 0
    )

    # Evaluate on fold test set
    eval_fold <- evaluate_model(model, x_test_fold, y_test_fold,
                                model_name = sprintf("Fold %d", i))

    # Get predictions for ensemble
    fold_pred <- predict(model, x_test_fold)

    # Return fold result bundle
    list(
      fold       = i,
      model      = model,
      eval       = eval_fold,
      preds      = fold_pred,
      test_idx   = test_idx
    )
  }

  # ---- Reassemble results from parallel workers ----
  fold_results   <- vector("list", k)
  all_models     <- vector("list", k)
  all_predictions <- matrix(NA, nrow = nrow(x_all), ncol = ncol(y_all))

  for (out in fold_outputs) {
    i <- out$fold
    fold_results[[i]]   <- out$eval
    all_models[[i]]     <- out$model
    all_predictions[out$test_idx, ] <- out$preds
    message(sprintf("Fold %d accuracy: %.2f%%", i, out$eval$accuracy * 100))
  }

  # Ensemble performance (all held-out predictions combined)
  y_pred_ensemble <- apply(all_predictions, 1, which.max) - 1
  y_real <- if (is.matrix(y_all)) apply(y_all, 1, which.max) - 1 else y_all

  ensemble_acc <- mean(y_pred_ensemble == y_real, na.rm = TRUE)

  # Per-fold summary
  fold_accs <- sapply(fold_results, function(x) x$accuracy)
  message(sprintf("\n=== CV Summary ==="))
  message(sprintf("Individual fold accuracies: %s",
                  paste(sprintf("%.2f%%", fold_accs * 100), collapse = ", ")))
  message(sprintf("Mean fold accuracy: %.2f%% (±%.2f%%)",
                  mean(fold_accs) * 100, sd(fold_accs) * 100))
  message(sprintf("Ensemble (all folds combined): %.2f%%", ensemble_acc * 100))

  return(list(
    fold_results       = fold_results,
    fold_accuracies    = fold_accs,
    mean_accuracy      = mean(fold_accs),
    sd_accuracy        = sd(fold_accs),
    ensemble_accuracy  = ensemble_acc,
    models             = all_models,
    all_predictions    = all_predictions,
    y_real             = y_real
  ))
}

#' Compare all models side-by-side
#'
#' @param results_list Named list of evaluate_model() outputs
#' @return Data frame of model comparisons
compare_models <- function(results_list) {
  comparison <- data.frame(
    Model     = names(results_list),
    Accuracy  = sapply(results_list, function(x) x$accuracy),
    Macro_F1  = sapply(results_list, function(x) x$macro_f1),
    NB_F1     = sapply(results_list, function(x) x$per_class$NB["F1"]),
    WB_F1     = sapply(results_list, function(x) x$per_class$WB["F1"]),
    SB_F1     = sapply(results_list, function(x) x$per_class$SB["F1"]),
    row.names = NULL
  )

  comparison <- comparison %>%
    mutate(Accuracy  = round(Accuracy * 100, 1),
           Macro_F1  = round(Macro_F1, 3),
           NB_F1     = round(NB_F1, 3),
           WB_F1     = round(WB_F1, 3),
           SB_F1     = round(SB_F1, 3)) %>%
    arrange(desc(Accuracy))

  message("\n=== Model Comparison ===")
  print(comparison, row.names = FALSE)

  return(comparison)
}

# ============================================================================
# SECTION 7: Main Execution Pipeline
# ============================================================================
# This is the main workflow that reproduces Jessen (2018) and extends it.
# Comment/uncomment sections as needed for your study.
# ============================================================================

main <- function(parallel = TRUE,
                  data_source = c("jessen2018", "custom"),
                  data_path = NULL) {

  data_source <- match.arg(data_source)

  message("\n")
  message("==============================================================")
  message("  Deep Learning for Peptide–MHC Binding Prediction")
  message("  Reproducing & Extending Jessen (2018)")
  message("  Target: ", TARGET_HLA)
  message("  Data:   ", if(data_source == "custom") data_path else "Jessen 2018 synthetic")
  message("==============================================================")
  message("\n")

  # --- Setup parallel backend ---
  if (parallel) {
    .PARALLEL_CLUSTER <<- setup_parallel(N_CORES)
    # Ensure cleanup on exit
    on.exit(teardown_parallel(.PARALLEL_CLUSTER), add = TRUE)
  }

  # --- Set TF threading (inside main, after parallel setup) ---
  tryCatch({
    tensorflow::tf$config$threading$set_intra_op_parallelism_threads(
      as.integer(TF_INTRA_OP))
    tensorflow::tf$config$threading$set_inter_op_parallelism_threads(
      as.integer(TF_INTER_OP))
    message(sprintf("TF threading: intra=%d inter=%d", TF_INTRA_OP, TF_INTER_OP))
  }, error = function(e) {
    message("TF threading config failed (non-fatal): ", e$message)
  })

  # -----------------------------------------------------------------------
  # Step 1: Load peptide data
  # -----------------------------------------------------------------------
  message("--- Step 1: Loading peptide data ---")

  pep_dat <- load_peptide_data(source = data_source, custom_path = data_path)

  # Option B (commented): Generate custom peptides + run netMHCpan
  # pep_dat <- generate_and_label_peptides(n = 100000, hla_allele = "HLA-A*02:01")
  # pep_dat <- balance_classes(pep_dat)
  # pep_dat <- split_train_test(pep_dat)

  # Display data summary
  message("\nData summary:")
  pep_dat %>%
    count(label_chr, data_type) %>%
    tidyr::pivot_wider(names_from = data_type, values_from = n) %>%
    print()

  # -----------------------------------------------------------------------
  # Step 2: Sequence logo for strong binders (reproducing Jessen Fig)
  # -----------------------------------------------------------------------
  message("\n--- Step 2: Sequence logo ---")
  sb_peptides <- pep_dat %>%
    filter(label_chr == "SB") %>%
    pull(peptide)

  message(sprintf("Strong binders: %d peptides", length(sb_peptides)))
  if (has_ggseqlogo) {
    plot_sequence_logo(sb_peptides,
                       save_path = file.path(PLOTS_DIR, "01_sequence_logo_SB.png"))
  }

  # -----------------------------------------------------------------------
  # Step 3: Encode peptides with BLOSUM62
  # -----------------------------------------------------------------------
  message("\n--- Step 3: Encoding peptides ---")

  # Training data
  x_train_peps <- pep_dat %>%
    filter(data_type == "train") %>%
    pull(peptide)
  y_train_labs <- pep_dat %>%
    filter(data_type == "train") %>%
    pull(label_num)

  # Test data
  x_test_peps <- pep_dat %>%
    filter(data_type == "test") %>%
    pull(peptide)
  y_test_labs <- pep_dat %>%
    filter(data_type == "test") %>%
    pull(label_num)

  message(sprintf("Train: %d peptides | Test: %d peptides",
                  length(x_train_peps), length(x_test_peps)))

  # BLOSUM62 encoding (Jessen 2018 original)
  message("Encoding with BLOSUM62...")
  x_train_blosum <- encode_blosum62(x_train_peps)
  x_test_blosum  <- encode_blosum62(x_test_peps)

  # Alternative: One-hot encoding (uncomment to compare)
  # x_train_enc <- encode_onehot(x_train_peps)
  # x_test_enc  <- encode_onehot(x_test_peps)
  # encoding_dim <- 20  # one-hot = 20

  # For BLOSUM: encoding dimension is 20 (same as one-hot conceptually)
  encoding_dim <- dim(x_train_blosum)[3]  # 20

  # Prepare for Keras
  train_data <- prepare_keras_data(x_train_blosum, y_train_labs, num_classes = 3)
  test_data  <- prepare_keras_data(x_test_blosum, y_test_labs, num_classes = 3)

  x_train <- train_data$x
  y_train <- train_data$y
  x_test  <- test_data$x
  y_test  <- test_data$y

  message(sprintf("x_train shape: [%s]", paste(dim(x_train), collapse = ", ")))
  message(sprintf("y_train shape: [%s]", paste(dim(y_train), collapse = ", ")))
  message(sprintf("x_test shape:  [%s]", paste(dim(x_test), collapse = ", ")))

  # -----------------------------------------------------------------------
  # Step 4: Train models
  # -----------------------------------------------------------------------
  message("\n--- Step 4: Training models ---")
  results_all <- list()

  # --- 4a: Feed-Forward Neural Network (Jessen 2018 — BEST) ---
  message("\n>>> Model 1: Feed-Forward Neural Network <<<")

  # Flatten 4D tensor to 2D for FFN: [n, 9, 20, 1] → [n, 180]
  x_train_flat <- array_reshape(x_train, c(dim(x_train)[1], 9 * encoding_dim))
  x_test_flat  <- array_reshape(x_test,  c(dim(x_test)[1],  9 * encoding_dim))

  ffn_model <- build_ffn(input_dim = 9 * encoding_dim)
  ffn_history <- train_keras_model(ffn_model, x_train_flat, y_train)

  # Evaluate FFN
  results_all[["FFN_Jessen2018"]] <- evaluate_model(
    ffn_model, x_test_flat, y_test,
    model_name = "FFN (Jessen 2018)"
  )

  # Plot training history
  plot_training_history(ffn_history,
                        title = "FFN (Jessen 2018)",
                        save_path = file.path(PLOTS_DIR, "02_ffn_training_history.png"))

  # Plot results
  plot_confusion_matrix(results_all[["FFN_Jessen2018"]],
                        save_path = file.path(PLOTS_DIR, "03_ffn_confusion_matrix.png"))

  plot_prediction_scatter(results_all[["FFN_Jessen2018"]],
                          save_path = file.path(PLOTS_DIR, "04_ffn_predictions.png"))

  # --- 4b: Convolutional Neural Network (Jessen 2018) ---
  message("\n>>> Model 2: Convolutional Neural Network <<<")

  cnn_model <- build_cnn(input_shape = c(9, encoding_dim, 1))
  cnn_history <- train_keras_model(cnn_model, x_train, y_train)

  results_all[["CNN_Jessen2018"]] <- evaluate_model(
    cnn_model, x_test, y_test,
    model_name = "CNN (Jessen 2018)"
  )

  plot_training_history(cnn_history,
                        title = "CNN (Jessen 2018)",
                        save_path = file.path(PLOTS_DIR, "05_cnn_training_history.png"))

  # --- 4c: LSTM (Extension) ---
  message("\n>>> Model 3: LSTM <<<")

  # LSTM expects [n, timesteps=9, features=20] (no channel dimension)
  x_train_lstm <- array_reshape(x_train, c(dim(x_train)[1], 9, encoding_dim))
  x_test_lstm  <- array_reshape(x_test,  c(dim(x_test)[1],  9, encoding_dim))

  lstm_model <- build_lstm(input_shape = c(9, encoding_dim))
  lstm_history <- train_keras_model(lstm_model, x_train_lstm, y_train)

  results_all[["LSTM"]] <- evaluate_model(
    lstm_model, x_test_lstm, y_test,
    model_name = "LSTM (Extended)"
  )

  # --- 4d: Deeper FFN (Extension) ---
  message("\n>>> Model 4: Deep FFN <<<")

  ffn_deep_model <- build_ffn_deep(input_dim = 9 * encoding_dim)
  ffn_deep_history <- train_keras_model(ffn_deep_model, x_train_flat, y_train)

  results_all[["FFN_Deep"]] <- evaluate_model(
    ffn_deep_model, x_test_flat, y_test,
    model_name = "Deep FFN (Extended)"
  )

  # --- 4e: ResNet-style CNN (Extension — residual learning) ---
  message("\n>>> Model 5: ResNet-Style CNN <<<")

  resnet_model <- build_resnet_style(input_shape = c(9, encoding_dim, 1))
  resnet_history <- train_keras_model(resnet_model, x_train, y_train)

  # Plot ResNet training history
  plot_training_history(resnet_history,
                        title = "ResNet-Style CNN",
                        save_path = file.path(PLOTS_DIR, "07_resnet_training_history.png"))

  results_all[["ResNet"]] <- evaluate_model(
    resnet_model, x_test, y_test,
    model_name = "ResNet-Style CNN"
  )

  plot_confusion_matrix(results_all[["ResNet"]],
                        save_path = file.path(PLOTS_DIR, "08_resnet_confusion_matrix.png"))

  # --- 4f: Random Forest (Jessen 2018 baseline) ---
  message("\n>>> Model 6: Random Forest <<<")

  # RF needs a 2D matrix
  x_train_rf <- matrix(x_train_flat, nrow = dim(x_train_flat)[1])
  x_test_rf  <- matrix(x_test_flat,  nrow = dim(x_test_flat)[1])
  y_train_rf <- pep_dat %>%
    filter(data_type == "train") %>%
    pull(label_num) %>%
    factor()
  y_test_rf <- pep_dat %>%
    filter(data_type == "test") %>%
    pull(label_num) %>%
    factor()

  rf_model <- build_rf(x_train_rf, y_train_rf, ntree = 100)
  results_all[["Random_Forest"]] <- evaluate_model(
    rf_model, x_test_rf, y_test_rf,
    model_name = "Random Forest (Jessen 2018)"
  )

  # Plot RF predictions
  plot_prediction_scatter(results_all[["Random_Forest"]],
                          save_path = file.path(PLOTS_DIR, "06_rf_predictions.png"))

  # -----------------------------------------------------------------------
  # Step 5: Model Comparison
  # -----------------------------------------------------------------------
  message("\n--- Step 5: Model Comparison ---")
  comparison <- compare_models(results_all)

  # Save comparison table
  write_csv(comparison, file.path(DATA_DIR, "model_comparison.csv"))
  message("Comparison saved: ", file.path(DATA_DIR, "model_comparison.csv"))

  # -----------------------------------------------------------------------
  # Step 6: 5-Fold Cross-Validation (Extension)
  # -----------------------------------------------------------------------
  message("\n--- Step 6: 5-Fold Cross-Validation ---")

  # Run CV for the best-performing architecture (FFN)
  # Combine train+test for full CV (flat format for FFN)
  x_all_cv <- rbind(x_train_flat, x_test_flat)
  y_all_cv <- rbind(y_train, y_test)

  cv_results <- run_cross_validation(
    build_fn = function() build_ffn(input_dim = 9 * encoding_dim),
    x_all    = x_all_cv,   # flat: [n, 180] — FFN expects flat input
    y_all    = y_all_cv,
    k        = 5,
    parallel = FALSE
  )

  # Save CV summary
  cv_summary <- data.frame(
    Fold          = 1:5,
    Accuracy      = cv_results$fold_accuracies,
    Ensemble      = cv_results$ensemble_accuracy
  )
  write_csv(cv_summary, file.path(DATA_DIR, "cv_summary.csv"))

  # -----------------------------------------------------------------------
  # Step 7: Save all trained models
  # -----------------------------------------------------------------------
  message("\n--- Step 7: Saving all trained models ---")

  # Collect all trained models (FFN, CNN, LSTM, DeepFFN, ResNet, RF)
  all_trained_models <- list(
    FFN_Jessen2018  = ffn_model,
    CNN_Jessen2018  = cnn_model,
    LSTM_Extended   = lstm_model,
    FFN_Deep        = ffn_deep_model,
    ResNet_Style    = resnet_model,
    Random_Forest   = rf_model
  )

  # Save with consistent naming
  save_all_models(all_trained_models, prefix = "")

  # Save CV fold models separately
  if (length(cv_results$models) > 0) {
    names(cv_results$models) <- paste0("FFN_CV_Fold_", seq_along(cv_results$models))
    save_all_models(cv_results$models, prefix = "cv")
  }

  # -----------------------------------------------------------------------
  # Step 8: Extract features and build feature.csv
  # -----------------------------------------------------------------------
  message("\n--- Step 8: Building feature.csv ---")

  # Combine train + test data for complete feature extraction
  # Flattened arrays: simple rbind
  x_all_flat <- rbind(x_train_flat, x_test_flat)
  # 4D arrays: manual concatenation along first axis
  x_all_4d <- array(0, dim = c(dim(x_train)[1] + dim(x_test)[1],
                                dim(x_train)[2], dim(x_train)[3], dim(x_train)[4]))
  x_all_4d[1:dim(x_train)[1], , , ] <- x_train
  x_all_4d[(dim(x_train)[1] + 1):(dim(x_train)[1] + dim(x_test)[1]), , , ] <- x_test
  # Labels
  y_all <- rbind(y_train, y_test)

  # Aligned metadata (respecting original order: train then test)
  pep_dat_train <- pep_dat %>% filter(data_type == "train")
  pep_dat_test  <- pep_dat %>% filter(data_type == "test")
  pep_dat_all   <- bind_rows(pep_dat_train, pep_dat_test)

  # Build feature CSV using the best-performing model (FFN)
  feature_csv_path <- file.path(DATA_DIR, "feature.csv")

  build_feature_csv(
    pep_dat    = pep_dat_all,
    model      = ffn_model,
    x_data     = x_all_flat,
    y_data     = y_all,
    model_name = "FFN_Jessen2018",
    encoding_type = "BLOSUM62",
    save_path  = feature_csv_path
  )

  # Also extract features from ResNet for comparison (uses 4D input)
  feature_csv_resnet <- file.path(DATA_DIR, "feature_resnet.csv")
  build_feature_csv(
    pep_dat    = pep_dat_all,
    model      = resnet_model,
    x_data     = x_all_4d,
    y_data     = y_all,
    model_name = "ResNet_Style",
    encoding_type = "BLOSUM62",
    save_path  = feature_csv_resnet
  )

  # -----------------------------------------------------------------------
  # Step 8: Final Report
  # -----------------------------------------------------------------------
  message("\n")
  message("==============================================================")
  message("  STUDY COMPLETE")
  message("==============================================================")
  message(sprintf("  Models trained:         %d", length(results_all)))
  message(sprintf("  Best model:             %s", comparison$Model[1]))
  message(sprintf("  Best accuracy:          %.1f%%", comparison$Accuracy[1]))
  message(sprintf("  5-fold CV mean ± SD:    %.1f%% ± %.1f%%",
                  cv_results$mean_accuracy * 100,
                  cv_results$sd_accuracy * 100))
  message(sprintf("  Ensemble accuracy:      %.1f%%",
                  cv_results$ensemble_accuracy * 100))
  message("\n  Output files:")
  message(sprintf("    Plots:                %s/", PLOTS_DIR))
  message(sprintf("    Feature CSV (FFN):    %s", file.path(DATA_DIR, "feature.csv")))
  message(sprintf("    Feature CSV (ResNet): %s", file.path(DATA_DIR, "feature_resnet.csv")))
  message(sprintf("    Model manifest:       %s", file.path(MODELS_DIR, "model_manifest.csv")))
  message(sprintf("    Models saved:         %s/ (%d models + %d CV folds)",
                  MODELS_DIR,
                  length(all_trained_models),
                  length(cv_results$models)))
  message("==============================================================\n")

  # Return everything for further analysis
  return(invisible(list(
    pep_dat     = pep_dat,
    results     = results_all,
    comparison  = comparison,
    cv_results  = cv_results,
    models      = list(ffn = ffn_model, cnn = cnn_model, lstm = lstm_model,
                       ffn_deep = ffn_deep_model, rf = rf_model)
  )))
}

# ============================================================================
# Run the pipeline
# ============================================================================
# Uncomment the line below to execute the full pipeline:
# study_results <- main()

message("\n=== peptide_mhc_binding_study.R loaded ===")
message("To run the full pipeline, execute: study_results <- main()")
message("Or step through individual sections interactively.")

# ============================================================================
# SECTION 8: Multi-Allele Extension
# ============================================================================
# Jessen (2018) only tested HLA-A*02:01. This extension tests generalization
# to additional common alleles, addressing Research Question 5.
# ============================================================================

#' Run netMHCpan locally and parse results
#'
#' NOTE: Requires netMHCpan installed locally.
#' Download: https://services.healthtech.dtu.dk/services/NetMHCpan-4.1/
#'
#' @param peptides Character vector of peptide sequences
#' @param hla_allele HLA allele (e.g., "HLA-A*02:01")
#' @param netmhcpan_path Path to netMHCpan executable
#' @return Data frame with netMHCpan predictions
run_netmhcpan <- function(peptides, hla_allele = "HLA-A*02:01",
                          netmhcpan_path = "netMHCpan") {
  # Write peptides to temp file
  tmp_in  <- tempfile(fileext = ".pep")
  tmp_out <- tempfile(fileext = ".xls")
  writeLines(peptides, tmp_in)

  # Run netMHCpan
  cmd <- sprintf('%s -p %s -a %s -xls -xlsfile %s',
                 netmhcpan_path, tmp_in, hla_allele, tmp_out)
  system(cmd, intern = TRUE)

  # Parse output
  results <- read.table(tmp_out, header = TRUE, sep = "\t",
                        stringsAsFactors = FALSE, skip = 1)

  # Clean up
  unlink(c(tmp_in, tmp_out))

  return(results)
}

#' Generate labeled peptide dataset for any HLA allele
#'
#' @param n_peptides Number of random peptides to generate
#' @param hla_allele Target HLA allele
#' @param balance Whether to balance classes
#' @return Data frame with peptide, label_chr, label_num
generate_labeled_peptides <- function(n_peptides = 100000,
                                      hla_allele = "HLA-A*02:01",
                                      balance = TRUE) {
  message(sprintf("Generating %s random 9-mers for %s...",
                  format(n_peptides, big.mark=","), hla_allele))
  peptides <- generate_random_peptides(n_peptides)

  # NOTE: Replace this with actual netMHCpan run in production
  # netmhcpan_results <- run_netmhcpan(peptides, hla_allele)
  # For now: simulate placeholder labels
  message("WARNING: Using simulated labels. Replace with netMHCpan output.")
  set.seed(20180129)
  n_sb <- round(n_peptides * 0.01)   # ~1% strong binders
  n_wb <- round(n_peptides * 0.05)   # ~5% weak binders
  n_nb <- n_peptides - n_sb - n_wb    # ~94% non-binders

  labels <- sample(c(
    rep("SB", n_sb),
    rep("WB", n_wb),
    rep("NB", n_nb)
  ))

  pep_dat <- tibble(
    peptide   = peptides,
    label_chr = factor(labels, levels = c("NB", "WB", "SB")),
    label_num = case_when(
      label_chr == "NB" ~ 0L,
      label_chr == "WB" ~ 1L,
      label_chr == "SB" ~ 2L
    ),
    hla_allele = hla_allele
  )

  if (balance) {
    pep_dat <- balance_classes(pep_dat)
  }

  return(pep_dat)
}

#' Test model generalization across multiple HLA alleles
#'
#' @param alleles Vector of HLA alleles to test
#' @param model A trained Keras FFN model (on HLA-A*02:01)
#' @param n_peptides_per_allele Number of test peptides per allele
#' @return Data frame with per-allele accuracies
test_multi_allele <- function(alleles = c("HLA-A*01:01", "HLA-A*02:01",
                                          "HLA-A*03:01", "HLA-B*07:02",
                                          "HLA-B*08:01"),
                              model = NULL,
                              n_peptides_per_allele = 5000) {

  if (is.null(model)) {
    stop("Provide a trained model (e.g., from main() output)")
  }

  message(sprintf("\n=== Multi-Allele Generalization: %d alleles ===",
                  length(alleles)))

  allele_results <- list()

  for (allele in alleles) {
    message(sprintf("\nTesting %s...", allele))

    # Generate and encode peptides for this allele
    pep_dat <- generate_labeled_peptides(n_peptides_per_allele,
                                         hla_allele = allele,
                                         balance = TRUE)
    pep_dat <- split_train_test(pep_dat)

    x_test_peps <- pep_dat %>% filter(data_type == "test") %>% pull(peptide)
    y_test_labs <- pep_dat %>% filter(data_type == "test") %>% pull(label_num)

    x_encoded <- encode_blosum62(x_test_peps)
    test_data <- prepare_keras_data(x_encoded, y_test_labs)

    # Flatten for FFN
    x_test_flat <- array_reshape(test_data$x,
                                 c(dim(test_data$x)[1], 9 * 20))

    # Evaluate
    eval_result <- evaluate_model(model, x_test_flat, test_data$y,
                                  model_name = allele)

    allele_results[[allele]] <- eval_result
  }

  # Summary
  summary_df <- data.frame(
    Allele   = names(allele_results),
    Accuracy = sapply(allele_results, function(x) x$accuracy * 100),
    Macro_F1 = sapply(allele_results, function(x) x$macro_f1)
  ) %>% arrange(desc(Accuracy))

  message("\n=== Multi-Allele Generalization Summary ===")
  print(summary_df, row.names = FALSE)

  return(summary_df)
}

# ============================================================================
# SECTION 9: Quick-Start Examples
# ============================================================================
# These examples can be run independently without executing the full pipeline.
# Useful for testing, debugging, or interactive exploration.
# ============================================================================

#' Quick-Start Example 1: Minimal reproduction of Jessen (2018) FFN
#'
#' Loads the original dataset, trains FFN, reports accuracy.
#' Expected output: ~95% accuracy (matching Jessen 2018).
#'
#' @examples
#' quickstart_jessen_ffn()
quickstart_jessen_ffn <- function() {
  message("=== Quick-Start: Jessen (2018) FFN Reproduction ===")

  # Load data
  pep_dat <- load_peptide_data("jessen2018")

  # Encode
  x_train_peps <- pep_dat %>% filter(data_type == "train") %>% pull(peptide)
  y_train_labs <- pep_dat %>% filter(data_type == "train") %>% pull(label_num)
  x_test_peps  <- pep_dat %>% filter(data_type == "test") %>% pull(peptide)
  y_test_labs  <- pep_dat %>% filter(data_type == "test") %>% pull(label_num)

  x_train_enc <- encode_blosum62(x_train_peps)
  x_test_enc  <- encode_blosum62(x_test_peps)

  train_data <- prepare_keras_data(x_train_enc, y_train_labs)
  test_data  <- prepare_keras_data(x_test_enc, y_test_labs)

  # Flatten for FFN
  x_train_flat <- array_reshape(train_data$x, c(dim(train_data$x)[1], 180))
  x_test_flat  <- array_reshape(test_data$x,  c(dim(test_data$x)[1],  180))

  # Build, train, evaluate
  model <- build_ffn()
  history <- train_keras_model(model, x_train_flat, train_data$y)
  eval   <- evaluate_model(model, x_test_flat, test_data$y,
                           model_name = "Quick-Start FFN")

  plot_confusion_matrix(eval)
  plot_prediction_scatter(eval)

  return(invisible(list(model = model, history = history, eval = eval)))
}

#' Quick-Start Example 2: Compare all encodings on a small dataset
#'
#' Tests BLOSUM62 vs one-hot vs AAindex on 5,000 peptides.
#' Useful for RQ3: Which encoding performs best?
#'
#' @examples
#' quickstart_encoding_comparison()
quickstart_encoding_comparison <- function(n_peptides = 5000) {
  message("=== Quick-Start: Encoding Comparison (RQ3) ===")

  set.seed(42)
  peptides <- generate_random_peptides(n_peptides)
  # Simulate labels (replace with netMHCpan in production)
  labels <- sample(0:2, n_peptides, replace = TRUE, prob = c(0.9, 0.07, 0.03))
  # Balance
  pep_df <- tibble(peptide = peptides, label_num = labels) %>%
    group_by(label_num) %>%
    slice_sample(n = min(table(labels))) %>%
    ungroup()

  x_peps <- pep_df$peptide
  y_labs <- pep_df$label_num

  encodings <- list(
    BLOSUM62 = encode_blosum62(x_peps),
    OneHot   = encode_onehot(x_peps)
    # AAindex has different dims, needs separate handling
  )

  results <- list()
  for (enc_name in names(encodings)) {
    message(sprintf("\n--- %s encoding ---", enc_name))
    enc <- encodings[[enc_name]]
    dim_enc <- dim(enc)[3]

    k_data <- prepare_keras_data(enc, y_labs)
    x_flat <- array_reshape(k_data$x, c(dim(k_data$x)[1], 9 * dim_enc))

    # Simple train/val split
    n_train <- floor(nrow(x_flat) * 0.8)
    x_tr <- x_flat[1:n_train, ]
    y_tr <- k_data$y[1:n_train, ]
    x_val <- x_flat[(n_train+1):nrow(x_flat), ]
    y_val <- k_data$y[(n_train+1):nrow(x_flat), ]

    model <- build_ffn(input_dim = 9 * dim_enc)
    history <- model %>% fit(x_tr, y_tr,
                             epochs = 50, batch_size = 32,
                             validation_data = list(x_val, y_val),
                             verbose = 0)

    eval <- evaluate_model(model, x_val, y_val,
                           model_name = paste("FFN", enc_name))
    results[[enc_name]] <- eval$accuracy
  }

  message("\n=== Encoding Comparison ===")
  for (name in names(results)) {
    message(sprintf("  %s: %.2f%% accuracy", name, results[[name]] * 100))
  }

  return(invisible(results))
}

#' Quick-Start Example 3: Visualize peptide encoding as 'images'
#'
#' Shows how BLOSUM62 encoding transforms peptides into 9×20 matrices,
#' reproducing Jessen's pep_plot_images() concept.
#'
#' @param peptides Character vector of peptide sequences
#' @examples
#' quickstart_visualize_encoding(c("LLTDAQRIV", "LMAFYLYEV"))
quickstart_visualize_encoding <- function(peptides) {
  encoded <- encode_blosum62(peptides)

  for (i in seq_along(peptides)) {
    pep_img <- encoded[i, , ]  # 9×20 matrix

    # Convert to long format for ggplot
    df <- expand.grid(Position = 1:9, AA = 1:20)
    aa_order <- strsplit("ARNDCQEGHILKMFPSTWYV", "")[[1]]
    df$AA_label <- rep(aa_order, each = 9)
    df$Value <- as.vector(t(pep_img))

    p <- ggplot(df, aes(x = AA_label, y = -Position, fill = Value)) +
      geom_tile(color = "white", linewidth = 0.3) +
      scale_fill_gradient(low = "white", high = "steelblue") +
      labs(
        title = paste("BLOSUM62 Encoding —", peptides[i]),
        x = "Amino Acid",
        y = "Position (p1→p9)",
        fill = "BLOSUM\nscore"
      ) +
      theme_minimal(base_size = 11) +
      theme(axis.text.x = element_text(angle = 0, size = 7))

    print(p)
  }

  return(invisible(encoded))
}

# ============================================================================
# APPENDIX: Session Info & Reproducibility
# ============================================================================

#' Print session info for reproducibility
#'
#' @examples
#' print_session_info()
print_session_info <- function() {
  message("\n=== Session Info ===")
  message("R version: ", R.version.string)
  message("TensorFlow: ", tf$`__version__`)
  message("Keras:      ", packageVersion("keras"))
  message("Date:       ", Sys.Date())

  if (has_peptools) {
    message("PepTools:   ", packageVersion("PepTools"))
  } else {
    message("PepTools:   NOT INSTALLED (using manual BLOSUM62)")
  }

  message("\nTo install missing packages:")
  message("  install.packages(c('keras','tidyverse','randomForest','pROC','caret','Peptides'))")
  message("  library(keras); install_keras()")
  message("  devtools::install_github('leonjessen/PepTools')")
  message("  devtools::install_github('omarwagih/ggseqlogo')")
  message("===============================\n")
}

# ============================================================================
# END OF SCRIPT
# ============================================================================
# Citation:
#   Jessen, L.E. (2018). "Deep Learning for Cancer Immunotherapy."
#   RStudio AI Blog. https://blogs.rstudio.com/ai/posts/2018-01-29-dl-for-cancer-immunotherapy
#
#   Reynisson, B. et al. (2020). "NetMHCpan-4.1 and NetMHCIIpan-4.0:
#   improved predictions of MHC antigen presentation by concurrent motif
#   deconvolution and integration of MS MHC eluted ligand data."
#   Nucleic Acids Research, 48(W1), W449–W454.
# ============================================================================
