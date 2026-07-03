# ============================================================
# Peptide Epitope Prediction — Nature Palette Figure Regeneration
# All manuscript figures + supplementary, PNG + PDF output
# Palette: NPG (Nature Publishing Group) via ggsci
# ============================================================
options(warn = -1)
Sys.setenv(TF_CPP_MIN_LOG_LEVEL = "3")
suppressMessages({
  library(tidyverse)
  library(ggsci)
  library(pROC)
  library(gridExtra)
  library(grid)
  library(ggrepel)
})

setwd('D:/Researching/Peptide epitope')
dir.create("03_Analysis/figures", showWarnings = FALSE, recursive = TRUE)
outdir <- "03_Analysis/figures"

# ── Nature NPG Palette ──────────────────────────────────────
npg_red   <- "#E64B35"
npg_blue  <- "#4DBBD5"
npg_green <- "#00A087"
npg_navy  <- "#3C5488"
npg_pink  <- "#F39B7F"
npg_grey  <- "#8491B4"
npg_orange <- "#DC0000"
npg_brown <- "#7E6148"
npg_purple <- "#B09C85"

nature_2 <- c(npg_navy, npg_red)        # 2-category: blue+red
nature_3 <- c(npg_navy, npg_red, npg_green)  # 3-category
nature_4 <- c(npg_navy, npg_red, npg_green, npg_orange)
nature_6 <- c(npg_navy, npg_red, npg_green, npg_orange, npg_blue, npg_purple)

# Helper: save both PNG + PDF
save_both <- function(p, name, w = 10, h = 6, dpi = 300) {
  ggsave(file.path(outdir, paste0(name, ".png")), p, width = w, height = h, dpi = dpi)
  ggsave(file.path(outdir, paste0(name, ".pdf")), p, width = w, height = h, device = "pdf")
  cat(sprintf("  ✓ %s.{png,pdf}\n", name))
}

# ── Clean old figures ───────────────────────────────────────
old_files <- list.files(outdir, pattern = "\\.(png|pdf)$", full.names = TRUE)
if (length(old_files) > 0) {
  file.remove(old_files)
  cat(sprintf("Deleted %d old figure files\n", length(old_files)))
}

# ============================================================
# FIGURE 1: Model Performance Comparison
# ============================================================
cat("\n── Figure 1: Model Comparison ──\n")

mhcflurry <- tribble(
  ~Model, ~Accuracy, ~Macro_F1, ~NB_F1, ~WB_F1, ~SB_F1,
  "MHCflurry 2.2.0†", 100.0, 1.000, 1.000, 1.000, 1.000,
  "Deep FFN",        91.9, 0.921, 0.969, 0.880, 0.913,
  "FFN (Jessen)",    90.9, 0.911, 0.969, 0.875, 0.891,
  "CNN",             90.0, 0.901, 0.963, 0.860, 0.882,
  "ResNet",          84.3, 0.847, 0.907, 0.792, 0.843,
  "LSTM",            83.3, 0.836, 0.935, 0.752, 0.822,
  "Random Forest",   81.1, 0.814, 0.927, 0.723, 0.793
) %>% mutate(Model = factor(Model, levels = rev(c(
  "MHCflurry 2.2.0†", "Deep FFN", "FFN (Jessen)", "CNN", "ResNet", "LSTM", "Random Forest"
))))

# Panel A: Accuracy bar chart
# MHCflurry bar in grey with dashed border; neural models in nature_7 colors
nature_7 <- c("grey60", nature_6)  # grey for MHCflurry, then 6 nature colors
mhcflurry_flag <- mhcflurry$Model == "MHCflurry 2.2.0†"

p1a <- ggplot(mhcflurry, aes(x = Accuracy, y = Model, fill = Model)) +
  geom_col(width = 0.6, alpha = 0.9) +
  geom_text(aes(label = ifelse(Model == "MHCflurry 2.2.0†",
                               sprintf("%.1f%%†", Accuracy),
                               sprintf("%.1f%%", Accuracy))),
            hjust = -0.15, size = 3.8, fontface = "bold") +
  scale_fill_manual(values = setNames(nature_7, levels(mhcflurry$Model)), guide = "none") +
  labs(title = "Model Accuracy — HLA-A*02:01 Peptide-MHC Binding",
       subtitle = "MHCflurry 2.2.0 labels, BLOSUM62 encoding, 5,088 balanced 9-mers\n†MHCflurry evaluated on its own training distribution (label circularity — see Section 4.6)",
       x = "Accuracy (%)", y = NULL) +
  xlim(0, 105) +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold"),
        panel.grid.major.y = element_blank())
save_both(p1a, "Figure1_model_comparison", w = 10, h = 5)

# Panel B: Per-class F1 for Deep FFN
deep_f1 <- mhcflurry %>% filter(Model == "Deep FFN") %>%
  select(NB_F1, WB_F1, SB_F1) %>%
  pivot_longer(everything(), names_to = "Class", values_to = "F1") %>%
  mutate(Class = factor(Class, levels = c("SB_F1", "WB_F1", "NB_F1"),
                        labels = c("Strong\nBinder", "Weak\nBinder", "Non\nBinder")))

p1b <- ggplot(deep_f1, aes(x = Class, y = F1, fill = Class)) +
  geom_col(width = 0.5, alpha = 0.9) +
  geom_text(aes(label = sprintf("%.3f", F1)), vjust = -0.5, size = 5, fontface = "bold") +
  scale_fill_manual(values = c(npg_red, npg_orange, npg_navy), guide = "none") +
  labs(title = "Per-Class F1 — Deep FFN (Best Model)",
       subtitle = "Weak Binder class is the primary locus of performance variation",
       y = "F1 Score", x = NULL) +
  ylim(0, 1.05) +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold"))
save_both(p1b, "Figure1_per_class_f1", w = 7, h = 5)

# ============================================================
# FIGURE 2: Label Quality Effect
# ============================================================
cat("\n── Figure 2: Label Quality ──\n")

pssm <- tribble(
  ~Model, ~Accuracy,
  "Deep FFN",        94.8,
  "FFN (Jessen)",    93.2,
  "CNN",             93.2,
  "ResNet",          88.4,
  "LSTM",            87.8,
  "Random Forest",   87.3
) %>% mutate(Model = factor(Model, levels = rev(c(
  "Deep FFN", "FFN (Jessen)", "CNN", "ResNet", "LSTM", "Random Forest"
))))

gap_df <- pssm %>% rename(PSSM = Accuracy) %>%
  left_join(mhcflurry %>% select(Model, Accuracy) %>% rename(MHCflurry = Accuracy), by = "Model") %>%
  mutate(Gap = PSSM - MHCflurry)

p2 <- ggplot(gap_df, aes(x = Gap, y = Model)) +
  geom_col(fill = npg_red, alpha = 0.8, width = 0.55) +
  geom_text(aes(label = sprintf("+%.1f%%", Gap)), hjust = -0.15,
            size = 4.5, fontface = "bold", color = npg_red) +
  labs(title = "Label Quality Dominates Architecture",
       subtitle = expression("PSSM − MHCflurry accuracy gap | Δ" ~ "3.0 pp > any architectural improvement"),
       x = "Accuracy Gap (pp)", y = NULL) +
  xlim(0, 8) +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold"),
        panel.grid.major.y = element_blank())
save_both(p2, "Figure2_label_comparison", w = 8, h = 4.5)

# ============================================================
# FIGURE 3: IEDB Benchmark
# ============================================================
cat("\n── Figure 3: IEDB Benchmark ──\n")

res <- read_csv("02_Data/cleaned/iedb_benchmark_results.csv", show_col_types = FALSE)

# ROC curve (base R plot → wrap in grid)
roc_obj <- roc(res$true_binary, res$pred_binary)

png(file.path(outdir, "Figure3_iedb_benchmark.png"), width = 2000, height = 1800, res = 250)
par(mar = c(4, 4, 3, 2), bg = "white")
plot.roc(roc_obj,
         main = "IEDB Benchmark — HLA-A*02:01",
         sub  = sprintf("Deep FFN | 49 POS + 20 NEG | AUC = %.3f | Sensitivity = 93.9%%", auc(roc_obj)),
         col = npg_navy, lwd = 3,
         print.auc = TRUE, print.auc.cex = 1.5,
         print.auc.x = 0.45, print.auc.y = 0.35,
         auc.polygon = TRUE, auc.polygon.col = paste0(npg_navy, "30"),
         grid = TRUE, grid.col = "gray90",
         legacy.axes = TRUE)
dev.off()
# PDF version
pdf(file.path(outdir, "Figure3_iedb_benchmark.pdf"), width = 8, height = 7)
par(mar = c(4, 4, 3, 2))
plot.roc(roc_obj, main = "IEDB Benchmark — HLA-A*02:01",
         sub = sprintf("Deep FFN | 49 POS + 20 NEG | AUC = %.3f | Sensitivity = 93.9%%", auc(roc_obj)),
         col = npg_navy, lwd = 3, print.auc = TRUE, print.auc.cex = 1.5,
         print.auc.x = 0.45, print.auc.y = 0.35,
         auc.polygon = TRUE, auc.polygon.col = paste0(npg_navy, "30"),
         grid = TRUE, grid.col = "gray90", legacy.axes = TRUE)
dev.off()
cat("  ✓ Figure3_iedb_benchmark.{png,pdf}\n")

# Confusion matrix
cm <- res %>%
  mutate(
    Actual = factor(ifelse(true_label == "POS", "True Epitope", "Non-Binder"),
                    levels = c("True Epitope", "Non-Binder")),
    Predicted = factor(case_when(
      pred_class == "SB" ~ "SB", pred_class == "WB" ~ "WB", TRUE ~ "NB"
    ), levels = c("SB", "WB", "NB"))
  )
cm_table <- table(Actual = cm$Actual, Predicted = cm$Predicted)
cm_df <- as.data.frame(as.table(cm_table)) %>%
  mutate(
    pct = case_when(
      Actual == "True Epitope" ~ Freq / sum(cm_table["True Epitope", ]),
      Actual == "Non-Binder"  ~ Freq / sum(cm_table["Non-Binder", ])
    ),
    label = paste0(Freq, "\n(", round(pct * 100, 1), "%)")
  )

p3_cm <- ggplot(cm_df, aes(x = Predicted, y = Actual)) +
  geom_tile(aes(fill = Freq), color = "white", linewidth = 1.2) +
  geom_text(aes(label = label), size = 5.5, fontface = "bold") +
  scale_fill_gradient(low = "#F0F4FA", high = npg_navy, guide = "none") +
  labs(title = "Confusion Matrix — IEDB Benchmark",
       subtitle = "Sensitivity: 93.9% | Specificity: 75.0% | F1: 0.920 | SB+WB = Predicted Binder") +
  theme_minimal(base_size = 15) +
  theme(plot.title = element_text(face = "bold", size = 18),
        plot.subtitle = element_text(size = 11, color = "gray40"),
        axis.text = element_text(size = 14, face = "bold"),
        panel.grid = element_blank())
save_both(p3_cm, "Figure3_confusion_matrix", w = 8, h = 5.5)

# Per-epitope dot plot
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
    true_label = factor(true_label, levels = c("NEG", "POS")),
    outcome = case_when(
      true_label == "POS" & pred_class %in% c("SB","WB") ~ "TP",
      true_label == "POS" & pred_class == "NB" ~ "FN",
      true_label == "NEG" & pred_class == "NB" ~ "TN",
      TRUE ~ "FP"
    ),
    p9 = str_sub(peptide, 9, 9)
  ) %>% arrange(desc(binding_score)) %>% mutate(idx = 1:n())

  fn_annotate <- p_dots %>% filter(outcome == "FN") %>% slice_head(n = 3)
  p_dots_sub <- "Sensitivity: 93.9% | AUC: 0.947 | 49 POS + 20 NEG | FN labelled with non-canonical p9 anchor"

p3_dots <- ggplot(p_dots, aes(x = binding_score, y = reorder(peptide, binding_score))) +
  geom_point(aes(color = outcome), size = 2.5) +
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "gray50", alpha = 0.5) +
  geom_text(data = fn_annotate, aes(label = paste0("p9=", p9)),
            hjust = -0.3, size = 3.2, color = "#E64B35", fontface = "italic") +
  scale_color_manual(values = c("TP" = npg_green, "FN" = npg_red,
                                 "TN" = npg_navy, "FP" = npg_orange)) +
  labs(title = "IEDB Benchmark — Per-Epitope Binding Scores",
       subtitle = p_dots_sub,
       x = "Binding Score (SB x 1.0 + WB x 0.5)", y = NULL, color = "Outcome") +
  xlim(0, 1.15) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"),
        axis.text.y = element_text(size = 6, family = "mono"),
        legend.position = "bottom")
save_both(p3_dots, "Figure3_per_epitope", w = 10, h = 12)

# ============================================================
# FIGURE 4: Protein Epitope Scanning
# ============================================================
cat("\n── Figure 4: Protein Epitope Scanning ──\n")

scan <- read_csv("02_Data/cleaned/protein_epitope_scan_extended.csv", show_col_types = FALSE)

scan <- scan %>%
  rename(bind_class = pred_class) %>%
  mutate(bind_class = factor(bind_class, levels = c("SB", "WB", "NB")))

p4 <- scan %>%
  count(protein, bind_class) %>%
  group_by(protein) %>% mutate(pct = n / sum(n)) %>% ungroup() %>%
  mutate(protein = fct_reorder(protein, pct * (bind_class == "SB"), .fun = sum)) %>%
  ggplot(aes(x = protein, y = pct, fill = bind_class)) +
  geom_col(alpha = 0.9) +
  scale_fill_manual(values = c("SB" = npg_red, "WB" = npg_orange, "NB" = npg_grey),
                    labels = c("SB" = "Strong Binder", "WB" = "Weak Binder", "NB" = "Non-Binder")) +
  scale_y_continuous(labels = scales::percent) +
  labs(title = "Protein Epitope Scanning — 10 Proteins, 3,536 9-mers",
       subtitle = sprintf("HLA-A*02:01 | Deep FFN | Predicted across 10 therapeutically relevant proteins",
                          sum(scan$bind_class == "SB")),
       x = NULL, y = "Proportion", fill = "Predicted Class") +
  coord_flip() +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold"), legend.position = "bottom")
save_both(p4, "Figure4_protein_scan", w = 10, h = 6)

# ============================================================
# FIGURE 5: Mutation Analysis
# ============================================================
cat("\n── Figure 5: Mutation Analysis ──\n")

mut <- read_csv("02_Data/cleaned/mutation_scan_results.csv", show_col_types = FALSE)
delta <- mut %>%
  filter(effect != "unchanged") %>%
  mutate(
    mut_clean = str_replace(mutation, " -- .*", ""),
    protein   = ifelse(str_detect(mut_clean, "p53"), "p53", "KRAS"),
    label     = paste0(mut_clean, "\n", wt_peptide, " → ", mut_peptide),
    effect_f  = factor(effect, levels = c("CREATED (neoepitope)", "ENHANCED", "DESTROYED", "WEAKENED"))
  )

p5 <- ggplot(delta, aes(x = reorder(label, delta_score), y = delta_score, fill = effect_f)) +
  geom_col(width = 0.55, alpha = 0.9) +
  geom_hline(yintercept = 0, linewidth = 0.8, color = "gray40") +
  geom_text(aes(label = sprintf("%+.3f", delta_score),
                y = ifelse(delta_score > 0, delta_score + 0.03, delta_score - 0.03)),
            size = 4, fontface = "bold", color = "gray20") +
  scale_fill_manual(values = c("CREATED (neoepitope)" = npg_green,
                                "ENHANCED" = npg_navy,
                                "DESTROYED" = npg_red,
                                "WEAKENED" = npg_orange)) +
  coord_flip() +
  labs(title = "Mutation Epitope Score Changes — p53 & KRAS Hotspots",
       subtitle = "HLA-A*02:01 binding score | Green=neoepitope, Red=destroyed, Blue=enhanced",
       x = NULL, y = expression(Delta ~ "Binding Score"), fill = "Effect") +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold", size = 15),
        panel.grid.major.y = element_blank(),
        legend.position = "bottom")
save_both(p5, "Figure5_mutation_delta", w = 12, h = 5.5)

# Paired WT vs mutant
paired <- delta %>%
  select(mut_clean, wt_peptide, mut_peptide, wt_score, mut_score, effect_f) %>%
  pivot_longer(c(wt_score, mut_score), names_to = "type", values_to = "score") %>%
  mutate(type = ifelse(type == "wt_score", "Wild-Type", "Mutant"),
         type = factor(type, levels = c("Wild-Type", "Mutant")),
         label = paste0(mut_clean, "\n", wt_peptide))

p5_paired <- ggplot(paired, aes(x = type, y = score, group = label)) +
  geom_line(aes(color = effect_f), linewidth = 1.2, alpha = 0.7) +
  geom_point(aes(color = effect_f, shape = type), size = 4) +
  scale_color_manual(values = c("CREATED (neoepitope)" = npg_green,
                                 "ENHANCED" = npg_navy,
                                 "DESTROYED" = npg_red,
                                 "WEAKENED" = npg_orange)) +
  scale_y_continuous(limits = c(0, 0.9)) +
  facet_wrap(~label, nrow = 2) +
  labs(title = "WT vs Mutant Binding Scores — p53 & KRAS",
       subtitle = "HLA-A*02:01 | Green=Neoepitope, Red=Destroyed",
       x = NULL, y = "Binding Score", color = "Effect", shape = NULL) +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold", size = 15),
        strip.text = element_text(size = 8.5, face = "bold"),
        legend.position = "bottom")
save_both(p5_paired, "Figure5_mutation_paired", w = 14, h = 7)

# ============================================================
# FIGURE 6: Feature Representations
# ============================================================
cat("\n── Figure 6: Feature Representations ──\n")

# BLOSUM62 matrix
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
for (i in 1:20) { rng <- range(blosum62[i,]); blosum62[i,] <- (blosum62[i,] - rng[1]) / (rng[2] - rng[1] + 1e-8) }

encode_blosum62 <- function(peptides) {
  encoded <- array(0, dim = c(length(peptides), 9, 20))
  for (i in seq_along(peptides)) {
    aa <- strsplit(peptides[i], "")[[1]]
    for (j in 1:9) encoded[i, j, ] <- blosum62[aa[j], ]
  }
  encoded
}

# Top epitopes
epitopes <- tribble(
  ~peptide,       ~label,               ~type,
  "GILGFVFTL",   "M1 58-66",           "Confirmed",
  "NLVPMVATV",   "pp65 495-503",       "Confirmed",
  "RMFPNAPYL",   "WT1 126-134",        "Confirmed",
  "LLGRNSFEV",   "p53 264-272",        "Confirmed",
  "IMDQVPFSV",   "gp100 209-217",      "Confirmed",
  "YMNGTMSQV",   "Tyrosinase 369-377", "Confirmed",
  "YLQPRTFLL",   "Spike 269-277",      "Confirmed",
  "SLLMWITQC",   "NY-ESO-1 157-165",   "Confirmed",
  "ALMDKSLHV",   "MART-1 56-64",       "Novel",
  "KIADYNYKL",   "Spike RBD 87-95",    "Novel",
  "RMPEAAPPV",   "p53 65-73",          "Novel",
  "LLTEVETYV",   "M1 3-11",            "Novel",
  "RLLQTGIHV",   "pp65 40-48",         "Novel",
  "SLGEQQYSV",   "WT1 187-195",        "Novel",
  "ILRGSVAHK",   "NP 265-273 (FN)",    "False Negative",
  "QYDPVAALF",   "pp65 341-349 (FN)",  "False Negative",
  "TLGIVCPIC",   "E7 86-94 (FN)",      "False Negative",
  "AAAAAAAAA",   "Poly-A (Neg)",       "Negative Control",
  "LLLLLLLLL",   "Poly-L (FP)",        "False Positive",
  "DDDDDDDDD",   "Poly-D (Neg)",       "Negative Control"
)

x <- encode_blosum62(epitopes$peptide)

# BLOSUM heatmap grid
plot_list <- list()
for (i in 1:nrow(epitopes)) {
  mat <- x[i, , ]
  df <- expand.grid(Pos = factor(9:1, levels = 9:1), AA = factor(aa_order, levels = aa_order))
  df$Value <- as.vector(t(mat[nrow(mat):1, ]))
  type_col <- case_when(
    epitopes$type[i] == "Confirmed" ~ npg_navy,
    epitopes$type[i] == "Novel" ~ npg_green,
    epitopes$type[i] == "False Negative" ~ npg_red,
    epitopes$type[i] == "False Positive" ~ npg_orange,
    TRUE ~ npg_grey
  )
  p <- ggplot(df, aes(x = AA, y = Pos, fill = Value)) +
    geom_tile(color = "white", linewidth = 0.2) +
    scale_fill_gradient(low = "white", high = npg_navy, limits = c(0, 1), name = "BLOSUM") +
    labs(title = sprintf("%s | %s", epitopes$peptide[i], epitopes$label[i]),
         subtitle = epitopes$type[i], x = NULL, y = "Position") +
    theme_minimal(base_size = 8) +
    theme(plot.title = element_text(size = 9, face = "bold", family = "mono",
                                    color = type_col),
          plot.subtitle = element_text(size = 7, color = type_col),
          axis.text.x = element_text(size = 5.5, angle = 90, hjust = 1),
          axis.text.y = element_text(size = 6),
          legend.position = "right",
          legend.key.height = unit(0.6, "cm"),
          panel.grid = element_blank())
  plot_list[[i]] <- p
}

png(file.path(outdir, "Figure6_blosum62_heatmap.png"), width = 4800, height = 4800, res = 250)
grid.arrange(grobs = plot_list, ncol = 4,
             top = textGrob("BLOSUM62 Encoding — Top Epitopes (HLA-A*02:01)",
                            gp = gpar(fontface = "bold", fontsize = 20, col = npg_navy)))
dev.off()
pdf(file.path(outdir, "Figure6_blosum62_heatmap.pdf"), width = 19, height = 19)
grid.arrange(grobs = plot_list, ncol = 4,
             top = textGrob("BLOSUM62 Encoding — Top Epitopes (HLA-A*02:01)",
                            gp = gpar(fontface = "bold", fontsize = 20, col = npg_navy)))
dev.off()
cat("  ✓ Figure6_blosum62_heatmap.{png,pdf}\n")

# ============================================================
# FIGURE 7: ESM-2 Embedding Comparison
# ============================================================
cat("\n── Figure 7: ESM-2 Comparison ──\n")

esm_data <- tribble(
  ~Encoding,                      ~Accuracy, ~Dimensionality, ~Type,
  "BLOSUM62 (baseline)",          91.9,      180,             "Baseline",
  "ESM-2 t6 mean-pooled",        65.9,      320,             "Pooled",
  "ESM-2 t12 per-position + L2", 91.5,      4320,            "Per-position",
  "ESM-2 t12 per-position",      90.9,      4320,            "Per-position",
  "ESM-2 t6 per-position",       93.3,      2880,            "Per-position"
) %>% mutate(
  Encoding = factor(Encoding, levels = rev(c(
    "ESM-2 t6 per-position", "ESM-2 t12 per-position",
    "ESM-2 t12 per-position + L2", "ESM-2 t6 mean-pooled",
    "BLOSUM62 (baseline)"
  ))),
  Type = factor(Type, levels = c("Baseline", "Per-position", "Pooled"))
)

p7 <- ggplot(esm_data, aes(x = Accuracy, y = Encoding, fill = Type)) +
  geom_col(width = 0.5, alpha = 0.9) +
  geom_text(aes(label = sprintf("%.1f%%", Accuracy)),
            hjust = -0.15, size = 4.2, fontface = "bold") +
  geom_vline(xintercept = 91.9, linetype = "dashed", color = npg_grey, linewidth = 0.8) +
  annotate("text", x = 91.9, y = 0.8, label = "BLOSUM62 baseline",
           hjust = -0.05, size = 3.5, color = npg_grey, fontface = "italic") +
  scale_fill_manual(values = c("Baseline" = npg_grey, "Per-position" = npg_navy,
                                "Pooled" = npg_red)) +
  labs(title = "ESM-2 Embedding Strategies — Accuracy Comparison",
       subtitle = "t6 per-position (2,880-dim): +1.4 pp over BLOSUM62 | Mean-pooled collapses to 65.9%",
       x = "Accuracy (%)", y = NULL, fill = "Encoding Type") +
  xlim(0, 100) +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold"),
        panel.grid.major.y = element_blank(),
        legend.position = "bottom")
save_both(p7, "Figure7_esm2_comparison", w = 11, h = 5)

# ============================================================
# SUPPLEMENTARY FIGURES
# ============================================================
cat("\n── Supplementary Figures ──\n")

# S1: 5-fold CV fold-by-fold
cv_data <- tribble(
  ~Source,     ~Fold, ~Accuracy,
  "PSSM",      1,     90.42, "PSSM", 2, 91.71, "PSSM", 3, 91.72, "PSSM", 4, 89.90, "PSSM", 5, 91.82,
  "MHCflurry", 1,     90.86, "MHCflurry", 2, 89.38, "MHCflurry", 3, 89.88, "MHCflurry", 4, 89.29, "MHCflurry", 5, 88.79
)
p_s1 <- ggplot(cv_data, aes(x = factor(Fold), y = Accuracy, color = Source, group = Source)) +
  geom_line(linewidth = 1.2) + geom_point(size = 3.5) +
  scale_color_manual(values = c("PSSM" = npg_navy, "MHCflurry" = npg_red)) +
  labs(title = "5-Fold Cross-Validation — PSSM vs MHCflurry",
       subtitle = "PSSM: μ=91.1%±0.9% | MHCflurry: μ=89.6%±0.8%",
       x = "CV Fold", y = "Accuracy (%)", color = "Label Source") +
  ylim(87, 93) + theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold"), legend.position = "bottom")
save_both(p_s1, "FigureS1_per_class_f1", w = 8, h = 5)

# S2: Per-class F1 across all models
all_f1 <- mhcflurry %>%
  select(Model, NB_F1, WB_F1, SB_F1) %>%
  pivot_longer(-Model, names_to = "Class", values_to = "F1") %>%
  mutate(Class = factor(Class, levels = c("SB_F1", "WB_F1", "NB_F1"),
                        labels = c("Strong Binder", "Weak Binder", "Non-Binder")))
p_s2 <- ggplot(all_f1, aes(x = F1, y = Model, fill = Class)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6, alpha = 0.9) +
  scale_fill_manual(values = c("Strong Binder" = npg_navy, "Weak Binder" = npg_red,
                                "Non-Binder" = npg_green)) +
  labs(title = "Per-Class F1 Across All Architectures",
       subtitle = "Weak Binder class shows the most architecture-dependent variation",
       x = "F1 Score", y = NULL, fill = "Class") +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold"), legend.position = "bottom")
save_both(p_s2, "FigureS2_cv_folds", w = 9, h = 5)

# S3: Protein scan per-class for top 5 proteins
top5_proteins <- scan %>%
  filter(protein %in% c("CMV pp65", "gp100/PMEL", "M1", "Tyrosinase", "NY-ESO-1")) %>%
  count(protein, bind_class)
p_s3 <- ggplot(top5_proteins, aes(x = protein, y = n, fill = bind_class)) +
  geom_col(alpha = 0.9) +
  scale_fill_manual(values = c("SB" = npg_red, "WB" = npg_orange, "NB" = npg_grey)) +
  labs(title = "Top 5 Proteins — Epitope Class Distribution",
       subtitle = "SB=Strong Binder, WB=Weak Binder, NB=Non-Binder",
       x = NULL, y = "Number of 9-mers", fill = "Class") +
  coord_flip() + theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold"), legend.position = "bottom")
save_both(p_s3, "FigureS3_confusion_matrix", w = 8, h = 4.5)

# S4: Mutation protein map
p53_pos <- c(175, 220, 245, 248, 249, 273, 282)
kras_pos <- c(12, 13, 61, 146)
map_df <- bind_rows(
  tibble(pos = p53_pos, protein = "p53",
         label = c("R175H", "Y220C", "G245S", "R248W*", "R249S", "R273H", "R282W"),
         has_effect = c(FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE)),
  tibble(pos = kras_pos, protein = "KRAS",
         label = c("G12V/C/D/R*", "G13D*", "Q61H/L/R", "A146T*"),
         has_effect = c(TRUE, TRUE, FALSE, TRUE))
)
p_s4 <- ggplot(map_df, aes(x = pos, y = protein, color = has_effect)) +
  geom_segment(data = tibble(protein = c("p53","KRAS"), start = c(1,1), end = c(323,188)),
               aes(x = start, xend = end, y = protein, yend = protein),
               color = "gray70", linewidth = 3, inherit.aes = FALSE) +
  geom_point(size = 5) +
  geom_text_repel(aes(label = label), size = 3.2, fontface = "bold", min.segment.length = 0,
                  box.padding = 0.8, point.padding = 0.3, max.overlaps = 20) +
  scale_color_manual(values = c("TRUE" = npg_red, "FALSE" = "gray60"),
                     labels = c("TRUE" = "Epitope-altering", "FALSE" = "Silent")) +
  scale_x_continuous(limits = c(-10, 370)) +
  labs(title = "Mutation Hotspot Map — p53 & KRAS",
       subtitle = "* = epitope-altering | Red alters MHC-I binding prediction",
       x = "Amino Acid Position", y = NULL) +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold"), legend.position = "bottom")
save_both(p_s4, "FigureS4_mutation_paired", w = 12, h = 3.5)

# S5: TCGA/IARC mutation frequency — p53 R248W and KRAS G12V
s5_data <- tribble(
  ~Gene, ~Mutation, ~Cancer_Type, ~Frequency, ~Frequency_Label,
  "p53", "R248W", "Skin SCC",          3.53,  "3.5% of all TP53 mutations",
  "p53", "R248W", "Colorectal",        3.53,  "3.5% of all TP53 mutations",
  "p53", "R248W", "Oesophageal ADC",   3.53,  "3.5% of all TP53 mutations",
  "p53", "R248W", "Breast",            3.53,  "3.5% of all TP53 mutations",
  "p53", "R248W", "Ovary",             3.53,  "3.5% of all TP53 mutations",
  "p53", "R248W", "Pancreas",          3.53,  "3.5% of all TP53 mutations",
  "KRAS", "G12V", "Pancreatic ADC",   31.0,  "28–36% of KRAS-mutant PAAD",
  "KRAS", "G12V", "Lung ADC",         19.0,  "~19% of KRAS-mutant LUAD",
  "KRAS", "G12V", "Colorectal",       17.0,  "~17% of KRAS-mutant CRC",
  "KRAS", "G12V", "Small intestine",  14.0,  "~14% of KRAS-mutant SIC",
  "KRAS", "G12V", "Ampullary",        13.0,  "~13% of KRAS-mutant"
) %>% mutate(
  highlight = case_when(
    Cancer_Type == "Pancreatic ADC" ~ "Neoepitope: YKLVVVGAV",
    Cancer_Type == "Lung ADC" ~ "Neoepitope: YKLVVVGAV",
    TRUE ~ "Background"
  ),
  Gene = factor(Gene, levels = c("p53", "KRAS")),
  y_order = paste(Gene, Cancer_Type, sep = " — ")
)

p_s5 <- ggplot(s5_data, aes(x = Frequency, y = reorder(Cancer_Type, Frequency))) +
  geom_col(aes(fill = highlight, alpha = highlight), width = 0.6) +
  geom_text(aes(label = Frequency_Label), hjust = -0.05, size = 3.5, fontface = "italic") +
  facet_wrap(~Gene, scales = "free_y", ncol = 1) +
  scale_fill_manual(values = c("Neoepitope: YKLVVVGAV" = "#E64B35",
                                "Background" = "#8491B4"), guide = "none") +
  scale_alpha_manual(values = c("Neoepitope: YKLVVVGAV" = 1.0, "Background" = 0.5), guide = "none") +
  labs(title = "Cancer Mutation Frequency — p53 R248W & KRAS G12V",
       subtitle = "R248W is among the top 8 most common TP53 hotspots (IARC R18); G12V is the 2nd most common KRAS mutation in PAAD (TCGA/FMI/GENIE consensus)\nRed bars mark cancers where the neoepitope candidates (MNWRPILTI / YKLVVVGAV) are therapeutically relevant",
       x = "Frequency (% of gene-mutant cases)", y = NULL) +
  xlim(0, 42) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold", size = 14, color = "white"),
        strip.background = element_rect(fill = "#3C5488", color = NA),
        panel.spacing = unit(1.5, "lines"))
save_both(p_s5, "FigureS5_tcga_frequency", w = 10, h = 5.5)

# ============================================================
# DONE
# ============================================================
cat("\n═══ All figures regenerated with Nature NPG palette ═══\n")
cat(sprintf("Output: %s/\n", outdir))
cat("Palette: NPG (Nature Publishing Group) via ggsci\n")
cat("Formats: PNG (300 dpi) + PDF (vector)\n")
