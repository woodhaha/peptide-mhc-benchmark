#!/usr/bin/env Rscript
# Retrain Deep FFN with homopolymer NB examples added to training data
# Usage: Rscript retrain_with_homopolymers.R <output_path.h5>
options(warn = -1)
Sys.setenv(TF_CPP_MIN_LOG_LEVEL = "3")
suppressMessages({library(keras); library(tensorflow); library(tidyverse)})

args <- commandArgs(trailingOnly = TRUE)
output_path <- ifelse(length(args) >= 1, args[1], "models/ffn_homopoly.h5")
cat(sprintf("Data-Augmented Retraining: adding homopolymers as NB\nOutput: %s\n", output_path))

EPOCHS <- 100; BATCH_SIZE <- 64; LEARNING_RATE <- 0.001; PATIENCE <- 15
setwd("D:/Researching/Peptide epitope")

# ---- Load training data ----
cat("Loading training data...\n")
pep_dat <- read_csv("Data/raw/real_peptides.csv", show_col_types = FALSE)

# ---- Add homopolymers as explicit NB examples ----
homopolymers <- c("AAAAAAAAA","GGGGGGGGG","PPPPPPPPP","DDDDDDDDD","RRRRRRRRR",
                   "EEEEEEEEE","KKKKKKKKK","NNNNNNNNN","SSSSSSSSS","TTTTTTTTT",
                   "QQQQQQQQQ","HHHHHHHHH","CCCCCCCCC","WWWWWWWWW","YYYYYYYYY",
                   "FFFFFFFFF","MMMMMMMMM","IIIIIIIII","VVVVVVVVV","LLLLLLLLL")

cat(sprintf("Adding %d homopolymers as NB training examples...\n", length(homopolymers)))
homo_df <- tibble(
  peptide_id = paste0("HOMO_", 1:length(homopolymers)),
  peptide    = homopolymers,
  label_chr  = factor(rep("NB", length(homopolymers)), levels = c("NB","WB","SB")),
  label_num  = rep(0L, length(homopolymers)),
  data_type  = rep("train", length(homopolymers))
)

# Also add to train set — ensure homopolymers are in training
pep_dat$data_type <- "train"
pep_aug <- bind_rows(pep_dat, homo_df)
cat(sprintf("Training data: %d → %d peptides (with homopolymers)\n", nrow(pep_dat), nrow(pep_aug)))

# ---- BLOSUM62 encode ----
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

# Encode all peptides (including homopolymers)
x <- encode_blosum62(pep_aug$peptide)
x_flat <- array_reshape(x, c(dim(x)[1], 180))
y <- to_categorical(pep_aug$label_num, 3)

# Split: 90/10 for non-homopolymer data, keep all homopolymers in training
set.seed(20260615)
homo_idx <- which(grepl("^HOMO_", pep_aug$peptide_id))
non_homo_idx <- setdiff(1:nrow(pep_aug), homo_idx)

test_nh <- sample(non_homo_idx, round(length(non_homo_idx) * 0.1))
train_nh <- setdiff(non_homo_idx, test_nh)
train_idx <- sort(c(train_nh, homo_idx))  # all homopolymers in train
test_idx <- sort(test_nh)

x_train <- x_flat[train_idx,]; y_train <- y[train_idx,]
x_test  <- x_flat[test_idx,];  y_test  <- y[test_idx,]

cat(sprintf("Train: %d (incl %d homopolymers)  Test: %d\n",
    nrow(x_train), length(intersect(train_idx, homo_idx)), nrow(x_test)))

# ---- Build Deep FFN ----
input_layer <- layer_input(shape = c(180))
xx <- input_layer %>%
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
model <- keras_model(inputs = input_layer, outputs = xx)

model %>% compile(
  loss = "categorical_crossentropy",
  optimizer = optimizer_adam(learning_rate = LEARNING_RATE),
  metrics = c("accuracy")
)
cat(sprintf("Model: %s params\n", format(count_params(model), big.mark = ",")))

# Moderate class weight
class_weight <- list("0" = 3.0, "1" = 1.0, "2" = 1.0)

early_stop <- callback_early_stopping(monitor = "val_loss", patience = PATIENCE,
                                       restore_best_weights = TRUE)
reduce_lr <- callback_reduce_lr_on_plateau(monitor = "val_loss", factor = 0.5,
                                            patience = 5, min_lr = 1e-6)

cat("Training with homopolymer augmentation...\n")
history <- model %>% fit(
  x_train, y_train,
  epochs = EPOCHS, batch_size = BATCH_SIZE,
  validation_split = 0.1,
  class_weight = class_weight,
  callbacks = list(early_stop, reduce_lr),
  verbose = 1
)

test_preds <- predict(model, x_test, verbose = 0)
test_class <- apply(test_preds, 1, which.max) - 1
true_class <- apply(y_test, 1, which.max) - 1
cat(sprintf("Test accuracy: %.4f\n", mean(test_class == true_class)))

dir.create(dirname(output_path), showWarnings = FALSE, recursive = TRUE)
save_model_hdf5(model, output_path)
cat(sprintf("Model saved: %s\n", output_path))
