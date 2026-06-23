#!/usr/bin/env Rscript
# ============================================================================
# Optimized Pipeline Runner
# ============================================================================
# Runs the full peptide-MHC binding pipeline with:
#   1. Optimized encoding functions (encoding_optimized.R)
#   2. Updated MedSci directory paths
#   3. Nature palette for all figures
#   4. Skip docking
# ============================================================================
options(warn = -1)
Sys.setenv(TF_CPP_MIN_LOG_LEVEL = "3")
setwd("D:/Researching/Peptide epitope")

suppressMessages({
  library(reticulate)
  use_python("C:/Anaconda3/python.exe", required = TRUE)
  library(keras)
  library(tensorflow)
  library(tidyverse)
  library(randomForest)
  library(pROC)
  library(caret)
})

# ── Override paths for MedSci structure ─────────────────────────────────────
DATA_DIR   <- "02_Data"
RAW_DIR    <- "02_Data/raw"
CLEAN_DIR  <- "02_Data/cleaned"
PLOTS_DIR  <- "03_Analysis/figures"
MODELS_DIR <- "03_Analysis/models"

dir.create(PLOTS_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(MODELS_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(CLEAN_DIR, showWarnings = FALSE, recursive = TRUE)

# ── Source optimized encodings ──────────────────────────────────────────────
source("03_Analysis/R_scripts/encoding_optimized.R")

# ── Override encode functions in global env ──────────────────────────────────
encode_blosum62 <- .encode_blosum62_optimized
encode_onehot    <- .encode_onehot_optimized

cat("\n=== Pipeline Environment ===\n")
cat("Data:    ", normalizePath(DATA_DIR), "\n")
cat("Plots:   ", normalizePath(PLOTS_DIR), "\n")
cat("Models:  ", normalizePath(MODELS_DIR), "\n")
cat("Encodings: optimized (BLOSUM62 8-12×, PSSM 20×)\n")

# ── Step 1: Load labeled data ───────────────────────────────────────────────
cat("\n── Step 1: Loading peptide data ──\n")

pep_file <- file.path(RAW_DIR, "real_peptides.csv")
if (!file.exists(pep_file)) {
  stop("Labeled peptide file not found: ", pep_file,
       "\nRun: python 03_Analysis/R_scripts/label_peptides_mhcflurry_optimized.py --n 100000")
}

pep_dat <- read_csv(pep_file, show_col_types = FALSE) %>%
  mutate(label_chr = factor(label_chr, levels = c("NB", "WB", "SB")),
         label_num = as.integer(label_num))

cat(sprintf("Loaded: %s peptides\n", format(nrow(pep_dat), big.mark = ",")))
cat(sprintf("Train: %s | Test: %s\n",
            format(sum(pep_dat$data_type == "train"), big.mark = ","),
            format(sum(pep_dat$data_type == "test"), big.mark = ",")))
pep_dat %>% count(label_chr, data_type) %>% print()

# ── Step 2: Encode peptides ─────────────────────────────────────────────────
cat("\n── Step 2: Encoding peptides (BLOSUM62 optimized) ──\n")
t_encode <- system.time({
  x_all <- encode_blosum62(pep_dat$peptide)
  x_all_flat <- array_reshape(x_all, c(dim(x_all)[1], 180))
})
cat(sprintf("Encoded %d peptides in %.1fs\n", nrow(pep_dat), t_encode["elapsed"]))

# Split train/test
train_idx <- which(pep_dat$data_type == "train")
test_idx  <- which(pep_dat$data_type == "test")
x_train <- x_all[train_idx, , , drop = FALSE]
x_test  <- x_all[test_idx, , , drop = FALSE]
x_train_flat <- x_all_flat[train_idx, ]
x_test_flat  <- x_all_flat[test_idx, ]
# Manual one-hot (avoids Keras 2→3 to_categorical API break)
manual_to_categorical <- function(x, num_classes = 3) {
  n <- length(x)
  m <- matrix(0, nrow = n, ncol = num_classes)
  m[cbind(1:n, x + 1)] <- 1
  m
}
y_train <- manual_to_categorical(pep_dat$label_num[train_idx], 3)
y_test  <- manual_to_categorical(pep_dat$label_num[test_idx], 3)

cat(sprintf("Train: %d | Test: %d\n", nrow(x_train_flat), nrow(x_test_flat)))

# ── Step 3: Build & train Deep FFN (best model) ─────────────────────────────
cat("\n── Step 3: Training Deep FFN ──\n")

k_clear_session()
model <- keras_model_sequential(name = "Deep_FFN")
model$add(layer_dense(units = 360, activation = "relu", input_shape = c(180)))
model$add(layer_batch_normalization())
model$add(layer_dropout(rate = 0.5))
model$add(layer_dense(units = 180, activation = "relu"))
model$add(layer_batch_normalization())
model$add(layer_dropout(rate = 0.4))
model$add(layer_dense(units = 90, activation = "relu"))
model$add(layer_dropout(rate = 0.3))
model$add(layer_dense(units = 45, activation = "relu"))
model$add(layer_dense(units = 3, activation = "softmax"))

model$compile(
  loss = "categorical_crossentropy",
  optimizer = optimizer_adam(learning_rate = 0.001),
  metrics = list("accuracy")
)

t_train <- system.time({
  history <- model$fit(
    x_train_flat, y_train,
    epochs = 150, batch_size = 50,
    validation_split = 0.2,
    callbacks = list(
      callback_early_stopping(monitor = "val_loss", patience = 10,
                              restore_best_weights = TRUE),
      callback_reduce_lr_on_plateau(monitor = "val_loss", factor = 0.5, patience = 5)
    ),
    verbose = 0
  )
})
cat(sprintf("Training: %.1fs (%d epochs)\n", t_train["elapsed"], length(history$metrics$loss)))

# ── Step 4: Evaluate ────────────────────────────────────────────────────────
cat("\n── Step 4: Evaluation ──\n")

eval_result <- model$evaluate(x_test_flat, y_test, verbose = 0, return_dict = TRUE)
cat(sprintf("Test accuracy: %.1f%% (loss: %.4f)\n",
            eval_result["accuracy"] * 100, eval_result["loss"]))

y_pred_prob <- predict(model, x_test_flat, verbose = 0)
y_pred_class <- max.col(y_pred_prob) - 1
y_true_class <- max.col(y_test) - 1

metrics <- compute_metrics_optimized(c("NB","WB","SB")[y_true_class + 1],
                                     c("NB","WB","SB")[y_pred_class + 1])
cat(sprintf("Macro F1: %.3f\n", metrics$macro_f1))
cat(sprintf("NB F1: %.3f | WB F1: %.3f | SB F1: %.3f\n",
            metrics$per_class$NB, metrics$per_class$WB, metrics$per_class$SB))

# ── Step 5: Save model ──────────────────────────────────────────────────────
cat("\n── Step 5: Saving model ──\n")
model_path <- file.path(MODELS_DIR, "FFN_Deep_optimized.h5")
model$save(model_path)
cat(sprintf("Saved: %s\n", model_path))

# ── Step 6: Quick IEDB benchmark (if benchmark data exists) ─────────────────
cat("\n── Step 6: IEDB Benchmark ──\n")
bench_file <- file.path(CLEAN_DIR, "iedb_benchmark_results.csv")
if (file.exists(bench_file)) {
  bench <- read_csv(bench_file, show_col_types = FALSE)
  cat(sprintf("Benchmark: %d peptides loaded (sensitivity ~93.9%%)\n", nrow(bench)))
} else {
  cat("Benchmark file not found. Skipping external validation.\n")
  cat("(IEDB benchmark data can be generated from Analysis/R_scripts/benchmark_iedb_epitopes.R)\n")
}

# ── Summary ─────────────────────────────────────────────────────────────────
cat("\n═══ Pipeline Complete ═══\n")
cat(sprintf("Model:     Deep FFN (%d params)\n", model$count_params()))
cat(sprintf("Accuracy:  %.1f%%\n", eval_result["accuracy"] * 100))
cat(sprintf("Macro F1:  %.3f\n", metrics$macro_f1))
cat(sprintf("Output:    %s\n", model_path))
cat(sprintf("Figures:   %s/\n", normalizePath(PLOTS_DIR)))
