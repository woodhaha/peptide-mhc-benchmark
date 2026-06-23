#!/usr/bin/env Rscript
# Retrain Deep FFN with class weights to penalize NB→SB/WB false positives
# Usage: Rscript retrain_class_weighted.R <nb_weight> <output_path.h5>
# Example: Rscript retrain_class_weighted.R 3.0 models/ffn_cw_nb3.h5

options(warn = -1)
Sys.setenv(TF_CPP_MIN_LOG_LEVEL = "3")
suppressMessages({library(keras); library(tensorflow); library(tidyverse)})

args <- commandArgs(trailingOnly = TRUE)
nb_weight <- ifelse(length(args) >= 1, as.numeric(args[1]), 3.0)
output_path <- ifelse(length(args) >= 2, args[2], "models/ffn_class_weighted.h5")

cat(sprintf("Class-Weighted Retraining: NB weight = %.1f\n", nb_weight))
cat(sprintf("Output: %s\n", output_path))

# ---- Constants ----
EPOCHS <- 100
BATCH_SIZE <- 64
LEARNING_RATE <- 0.001
PATIENCE <- 15

# ---- Load training data ----
cat("Loading training data...\n")
setwd("D:/Researching/Peptide epitope")
pep_dat <- read_csv("Data/raw/real_peptides.csv", show_col_types = FALSE)

# Encode with BLOSUM62
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

# Encode all peptides
x <- encode_blosum62(pep_dat$peptide)
x_flat <- array_reshape(x, c(dim(x)[1], 180))

# Labels
y <- to_categorical(pep_dat$label_num, 3)
colnames(y) <- c("NB", "WB", "SB")

# Train/test split (same as original: 90/10 stratified)
set.seed(20260615)
idx <- 1:nrow(pep_dat)
test_idx <- c()
for (cl in 0:2) {
  class_idx <- which(pep_dat$label_num == cl)
  test_idx <- c(test_idx, sample(class_idx, round(length(class_idx) * 0.1)))
}
train_idx <- setdiff(idx, test_idx)

x_train <- x_flat[train_idx, ]
y_train <- y[train_idx, ]
x_test  <- x_flat[test_idx, ]
y_test  <- y[test_idx, ]

cat(sprintf("Train: %d  Test: %d\n", nrow(x_train), nrow(x_test)))

# Class frequencies
class_counts <- table(pep_dat$label_num[train_idx])
cat(sprintf("Train classes: NB=%d WB=%d SB=%d\n",
    class_counts[1], class_counts[2], class_counts[3]))

# ---- Build Deep FFN (same architecture as original) ----
input_layer <- layer_input(shape = c(180))
x <- input_layer %>%
  layer_dense(units = 360, activation = "relu", kernel_initializer = "he_normal") %>%
  layer_batch_normalization() %>%
  layer_dropout(rate = 0.4) %>%
  layer_dense(units = 180, activation = "relu", kernel_initializer = "he_normal") %>%
  layer_batch_normalization() %>%
  layer_dropout(rate = 0.4) %>%
  layer_dense(units = 90, activation = "relu", kernel_initializer = "he_normal") %>%
  layer_dropout(rate = 0.3) %>%
  layer_dense(units = 45, activation = "relu", kernel_initializer = "he_normal") %>%
  layer_dense(units = 3, activation = "softmax")

model <- keras_model(inputs = input_layer, outputs = x)

model %>% compile(
  loss      = "categorical_crossentropy",
  optimizer = optimizer_adam(learning_rate = LEARNING_RATE),
  metrics   = c("accuracy")
)

cat(sprintf("Model built: %s parameters\n", format(count_params(model), big.mark = ",")))

# ---- Train with class weights ----
class_weight <- list("0" = nb_weight, "1" = 1.0, "2" = 1.0)
cat(sprintf("Class weights: NB=%.1f, WB=1.0, SB=1.0\n", nb_weight))

early_stop <- callback_early_stopping(monitor = "val_loss", patience = PATIENCE,
                                       restore_best_weights = TRUE)
reduce_lr <- callback_reduce_lr_on_plateau(monitor = "val_loss", factor = 0.5,
                                            patience = 5, min_lr = 1e-6)

cat("Training...\n")
history <- model %>% fit(
  x_train, y_train,
  epochs = EPOCHS,
  batch_size = BATCH_SIZE,
  validation_split = 0.1,
  class_weight = class_weight,
  callbacks = list(early_stop, reduce_lr),
  verbose = 1
)

# Quick eval on test set
test_preds <- predict(model, x_test, verbose = 0)
test_class <- apply(test_preds, 1, which.max) - 1
true_class <- apply(y_test, 1, which.max) - 1
test_acc <- mean(test_class == true_class)
cat(sprintf("Test accuracy: %.4f\n", test_acc))

# ---- Save model ----
dir.create(dirname(output_path), showWarnings = FALSE, recursive = TRUE)
save_model_hdf5(model, output_path)
cat(sprintf("Model saved: %s\n", output_path))
