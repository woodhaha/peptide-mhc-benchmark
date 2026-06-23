options(warn=-1); Sys.setenv('TF_CPP_MIN_LOG_LEVEL'='3')
suppressMessages({library(keras); library(tensorflow); library(tidyverse)})
setwd('D:/Researching/Peptide epitope')

# Load model
model <- load_model_hdf5('Analysis/FFN_Deep.h5')

# Read existing benchmark CSV
results <- read_csv('Data/cleaned/iedb_benchmark_results.csv', show_col_types=FALSE)

# Add homopolymer filter column
results <- results %>% mutate(
  unique_aa = sapply(strsplit(peptide,''), function(aa) length(unique(aa))),
  is_repeat = unique_aa <= 2,
  model_class = pred_class,
  pred_class_v2 = ifelse(is_repeat, 'NB', model_class)
)

# Write updated CSVs
write_csv(results, 'Data/cleaned/iedb_benchmark_results_v2.csv')
write_csv(results, 'Manuscript/data/iedb_benchmark_results_v2.csv')

# Summary
tp_v2 <- sum(results$true_label == 'POS' & results$pred_class_v2 %in% c('SB','WB'))
fp_v2 <- sum(results$true_label == 'NEG' & results$pred_class_v2 %in% c('SB','WB'))
tn_v2 <- sum(results$true_label == 'NEG' & results$pred_class_v2 == 'NB')
fn_v2 <- sum(results$true_label == 'POS' & results$pred_class_v2 == 'NB')
cat(sprintf('Updated benchmark: Sens=%.1f%% Spec=%.1f%% FP=%d TN=%d\n',
    100*tp_v2/(tp_v2+fn_v2), 100*tn_v2/(tn_v2+fp_v2), fp_v2, tn_v2))
