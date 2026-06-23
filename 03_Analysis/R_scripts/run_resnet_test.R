# ResNet comparison test — minimal output
options(warn = -1)

# Suppress TF logging
Sys.setenv(TF_CPP_MIN_LOG_LEVEL = "2")

setwd('D:/Researching/Peptide epitope')
suppressMessages(source('peptide_mhc_binding_study.R'))

cat('\n=== ResNet vs FFN Comparison ===\n')

# Data prep
pep_dat <- load_peptide_data("jessen2018")
x_train_peps <- pep_dat %>% filter(data_type == "train") %>% pull(peptide)
y_train_labs <- pep_dat %>% filter(data_type == "train") %>% pull(label_num)
x_test_peps  <- pep_dat %>% filter(data_type == "test") %>% pull(peptide)
y_test_labs  <- pep_dat %>% filter(data_type == "test") %>% pull(label_num)

x_train_enc <- encode_blosum62(x_train_peps)
x_test_enc  <- encode_blosum62(x_test_peps)
train_data <- prepare_keras_data(x_train_enc, y_train_labs)
test_data  <- prepare_keras_data(x_test_enc, y_test_labs)
x_train_flat <- array_reshape(train_data$x, c(dim(train_data$x)[1], 180))
x_test_flat  <- array_reshape(test_data$x,  c(dim(test_data$x)[1], 180))

results <- list()

# FFN baseline
cat('\nTraining FFN...\n')
ffn <- build_ffn(input_dim = 180)
capture.output(ffn_hist <- train_keras_model(ffn, x_train_flat, train_data$y, epochs = 50, verbose = 0))
results[["FFN_Jessen"]] <- evaluate_model(ffn, x_test_flat, test_data$y, "FFN (Jessen)")

# ResNet
cat('Training ResNet...\n')
resnet <- build_resnet_style(input_shape = c(9, 20, 1))
capture.output(resnet_hist <- train_keras_model(resnet, train_data$x, train_data$y, epochs = 50, verbose = 0))
results[["ResNet"]] <- evaluate_model(resnet, test_data$x, test_data$y, "ResNet-Style CNN")

# Head-to-head
cat('\n================================================\n')
cat('  MODEL COMPARISON: FFN vs ResNet-Style CNN\n')
cat('  (Synthetic data, 50 epochs, 23,760 peptides)\n')
cat('================================================\n')
for (name in names(results)) {
  r <- results[[name]]
  cat(sprintf('\n  [%s]\n', name))
  cat(sprintf('  Accuracy : %5.1f%%\n', r$accuracy * 100))
  cat(sprintf('  Macro F1 : %5.3f\n', r$macro_f1))
  for (cls in names(r$per_class)) {
    m <- r$per_class[[cls]]
    cat(sprintf('    %-3s: P=%5.3f  R=%5.3f  F1=%5.3f\n', cls, m["Precision"], m["Recall"], m["F1"]))
  }
}
cat('\nDone.\n')
