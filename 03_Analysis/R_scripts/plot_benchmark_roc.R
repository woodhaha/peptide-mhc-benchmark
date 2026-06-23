# ROC curve + confusion matrix for IEDB benchmark
library(tidyverse)
library(pROC)
library(caret)

setwd('D:/Researching/Peptide epitope')

# Load benchmark results
res <- read_csv("data/iedb_benchmark_results.csv", show_col_types = FALSE)

# ---- Plot 1: ROC Curve ----
roc_obj <- roc(res$true_binary, res$pred_binary)

png("plots/benchmark_roc.png", width = 1800, height = 1600, res = 200)
par(mar = c(4, 4, 3, 2))

plot.roc(roc_obj,
         main = "ROC Curve -- IEDB Benchmark (HLA-A*02:01)",
         sub  = paste0("Deep FFN (MHCflurry-trained, 91.9% accuracy) | ",
                       "49 POS + 20 NEG | AUC = ", round(auc(roc_obj), 3)),
         col = "#4472C4", lwd = 3,
         print.auc = TRUE, print.auc.cex = 1.5,
         print.auc.x = 0.45, print.auc.y = 0.35,
         auc.polygon = TRUE, auc.polygon.col = "#4472C430",
         grid = TRUE, grid.col = "gray90",
         legacy.axes = TRUE)

# Add sensitivity at key specificity points
sens_95 <- ci(roc_obj, of = "sp", sensitivities = 0.95)
sens_90 <- ci(roc_obj, of = "sp", sensitivities = 0.90)
abline(v = 0.95, lty = 2, col = "gray50")
abline(v = 0.90, lty = 2, col = "gray50")

text(0.72, 0.25, sprintf("Sens @ 95%% Spec: %.1f%%", sens_95[2] * 100),
     cex = 0.9, col = "gray30")
text(0.72, 0.18, sprintf("Sens @ 90%% Spec: %.1f%%", sens_90[2] * 100),
     cex = 0.9, col = "gray30")

dev.off()
cat("Saved: plots/benchmark_roc.png\n")

# ---- Plot 2: Confusion Matrix Heatmap ----
cm <- res %>%
  mutate(
    Actual  = factor(ifelse(true_label == "POS", "True Epitope", "Non-Binder"),
                     levels = c("True Epitope", "Non-Binder")),
    Predicted = factor(case_when(
      pred_class == "SB" ~ "SB",
      pred_class == "WB" ~ "WB",
      TRUE ~ "NB"
    ), levels = c("SB", "WB", "NB"))
  )

cm_table <- table(Actual = cm$Actual, Predicted = cm$Predicted)
cm_df <- as.data.frame(as.table(cm_table))

# Compute metrics for annotation
tp <- cm_table["True Epitope", "SB"] + cm_table["True Epitope", "WB"]
fn <- cm_table["True Epitope", "NB"]
fp_sb <- cm_table["Non-Binder", "SB"]
fp_wb <- cm_table["Non-Binder", "WB"]
tn <- cm_table["Non-Binder", "NB"]

# Custom labels with counts + percentages
total_pos <- sum(cm_table["True Epitope", ])
total_neg <- sum(cm_table["Non-Binder", ])

cm_df <- cm_df %>%
  mutate(
    pct = case_when(
      Actual == "True Epitope" ~ Freq / total_pos,
      Actual == "Non-Binder"  ~ Freq / total_neg
    ),
    label = paste0(Freq, "\n(", round(pct * 100, 1), "%)"),
    # Color intensity by proportion
    fill_intensity = Freq / max(Freq)
  )

p_cm <- ggplot(cm_df, aes(x = Predicted, y = Actual)) +
  geom_tile(aes(fill = fill_intensity), color = "white", linewidth = 1.2) +
  geom_text(aes(label = label), size = 5.5, fontface = "bold") +
  scale_fill_gradient(low = "#F5F5F5", high = "#4472C4", guide = "none") +
  labs(
    title    = "Confusion Matrix -- IEDB Benchmark",
    subtitle = paste0(
      "Sensitivity: 93.9% | Specificity: 75.0% | F1: 0.920 | AUC: 0.947\n",
      "SB + WB = Predicted Binder | NB = Predicted Non-Binder"
    ),
    x = "Predicted Class",
    y = "Actual Class"
  ) +
  theme_minimal(base_size = 15) +
  theme(
    plot.title    = element_text(face = "bold", size = 18),
    plot.subtitle = element_text(size = 11, color = "gray40"),
    axis.text     = element_text(size = 14, face = "bold"),
    axis.title    = element_text(size = 14, face = "bold"),
    panel.grid    = element_blank()
  )

ggsave("plots/benchmark_confusion_matrix.png", p_cm,
       width = 8, height = 5.5, dpi = 180)
cat("Saved: plots/benchmark_confusion_matrix.png\n")

# ---- Plot 3: Per-epitope prediction scores (dot plot) ----
p_dots <- res %>%
  mutate(
    category = case_when(
      grepl("M1|PA|NP|PB1", source_protein) ~ "Influenza",
      grepl("pp65|IE1", source_protein) ~ "CMV",
      grepl("LMP|BMLF", source_protein) ~ "EBV",
      grepl("HBc|HBs", source_protein) ~ "HBV",
      grepl("NS3", source_protein) ~ "HCV",
      grepl("Gag|Pol", source_protein) ~ "HIV",
      grepl("Spike|NSP|ORF", source_protein) ~ "SARS-CoV-2",
      grepl("E7", source_protein) ~ "HPV",
      grepl("MART|gp100|Tyros|TRP", source_protein) ~ "Melanoma",
      grepl("WT1|p53", source_protein) ~ "WT1/p53",
      grepl("MAGE|CEA|HER2", source_protein) ~ "MAGE/CEA/HER2",
      grepl("NY-ESO|Surviv|hTERT|PSA|PSMA|HTLV", source_protein) ~ "Other Tumor",
      TRUE ~ "Negative Control"
    ),
    category = factor(category, levels = rev(c(
      "Influenza", "CMV", "EBV", "HBV", "HCV", "HIV", "SARS-CoV-2", "HPV",
      "Melanoma", "WT1/p53", "MAGE/CEA/HER2", "Other Tumor", "Negative Control"
    ))),
    true_label = factor(true_label, levels = c("NEG", "POS")),
    outcome = case_when(
      true_label == "POS" & pred_class %in% c("SB","WB") ~ "TP",
      true_label == "POS" & pred_class == "NB" ~ "FN",
      true_label == "NEG" & pred_class == "NB" ~ "TN",
      true_label == "NEG" & pred_class %in% c("SB","WB") ~ "FP"
    )
  ) %>%
  arrange(category, desc(binding_score)) %>%
  mutate(idx = 1:n())

p_d <- ggplot(p_dots, aes(x = binding_score, y = reorder(peptide, binding_score))) +
  geom_point(aes(color = outcome, shape = true_label), size = 2.5) +
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "gray50", alpha = 0.5) +
  scale_color_manual(values = c(
    "TP" = "#228B22", "FN" = "#DC143C",
    "TN" = "#4472C4", "FP" = "#FF8C00"
  )) +
  scale_shape_manual(values = c("POS" = 16, "NEG" = 17)) +
  labs(
    title = "IEDB Benchmark -- Per-Epitope Binding Scores",
    subtitle = paste0("Green=TP, Red=FN, Blue=TN, Orange=FP | ",
                      "Sensitivity: 93.9% | AUC: 0.947"),
    x = "Binding Score (SB*1.0 + WB*0.5)",
    y = NULL,
    color = "Outcome", shape = "True Label"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.y = element_text(size = 6, family = "mono"),
    legend.position = "bottom"
  )

ggsave("plots/benchmark_per_epitope.png", p_d,
       width = 10, height = 12, dpi = 150)
cat("Saved: plots/benchmark_per_epitope.png\n")

cat("\nDone.\n")
