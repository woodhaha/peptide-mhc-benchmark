#!/usr/bin/env Rscript
# Redraw ALL project figures with Nature journal color palette
# Nature colors: Blue=#377EB8, Red=#E41A1C, Green=#4DAF4A, Purple=#984EA3,
#                Orange=#FF7F00, Brown=#A65628, Gray=#999999

library(tidyverse)
library(pROC)
library(caret)
library(ggrepel)

setwd('D:/Researching/Peptide epitope')

# ---- Nature Theme ----
nature_theme <- theme_minimal(base_size = 12) +
  theme(
    text = element_text(family = 'sans', color = '#333333'),
    plot.title = element_text(face = 'bold', size = 14, color = '#222222'),
    plot.subtitle = element_text(size = 10, color = '#666666'),
    axis.title = element_text(size = 11, face = 'bold', color = '#444444'),
    axis.text = element_text(size = 9, color = '#555555'),
    legend.position = 'bottom',
    legend.title = element_text(size = 10, face = 'bold'),
    legend.text = element_text(size = 9),
    panel.grid.major = element_line(color = '#E5E5E5', linewidth = 0.3),
    panel.grid.minor = element_blank(),
    plot.background = element_rect(fill = 'white', color = NA),
    panel.background = element_rect(fill = 'white', color = NA),
    strip.text = element_text(face = 'bold', size = 9),
    strip.background = element_rect(fill = '#F5F5F5', color = NA)
  )

# Nature color palette
NATURE_COLORS <- c(
  'blue'   = '#377EB8',
  'red'    = '#E41A1C',
  'green'  = '#4DAF4A',
  'purple' = '#984EA3',
  'orange' = '#FF7F00',
  'brown'  = '#A65628',
  'gray'   = '#999999',
  'dark'   = '#333333'
)

# For multi-class: PSSM=blue, MHCflurry=orange, etc.
PSSM_COL    <- '#377EB8'  # Nature blue
MHCFL_COL   <- '#FF7F00'  # Nature orange
CREATED_COL <- '#4DAF4A'  # Nature green
DESTROY_COL <- '#E41A1C'  # Nature red
ENHANCE_COL <- '#377EB8'  # Nature blue
WEAKEN_COL  <- '#999999'  # Nature gray

dir.create('figures', showWarnings = FALSE)
dir.create('plots', showWarnings = FALSE)

cat('========================================\n')
cat('  REDRAWING ALL FIGURES - NATURE STYLE\n')
cat('========================================\n\n')

############################################################
# FIGURE 1: Model Accuracy Comparison
############################################################
cat('[1/8] Model comparison...\n')

pssm <- tribble(
  ~Model, ~Accuracy, ~Macro_F1, ~NB_F1, ~WB_F1, ~SB_F1,
  "Deep FFN", 94.8, 0.948, 0.976, 0.925, 0.944,
  "FFN (Jessen)", 93.2, 0.932, 0.970, 0.900, 0.926,
  "CNN", 93.2, 0.932, 0.976, 0.901, 0.920,
  "ResNet", 88.4, 0.881, 0.917, 0.813, 0.914,
  "LSTM", 87.8, 0.879, 0.954, 0.825, 0.858,
  "Random Forest", 87.3, 0.875, 0.954, 0.823, 0.848
) %>% mutate(Source = "PSSM")

mhcflurry <- tribble(
  ~Model, ~Accuracy, ~Macro_F1, ~NB_F1, ~WB_F1, ~SB_F1,
  "Deep FFN", 91.9, 0.921, 0.969, 0.880, 0.913,
  "FFN (Jessen)", 90.9, 0.911, 0.969, 0.875, 0.891,
  "CNN", 90.0, 0.901, 0.963, 0.860, 0.882,
  "ResNet", 84.3, 0.847, 0.907, 0.792, 0.843,
  "LSTM", 83.3, 0.836, 0.935, 0.752, 0.822,
  "Random Forest", 81.1, 0.814, 0.927, 0.723, 0.793
) %>% mutate(Source = "MHCflurry")

cv_data <- tibble(
  Source = rep(c("PSSM", "MHCflurry"), each = 5),
  Fold = rep(1:5, 2),
  Accuracy = c(90.42, 91.71, 91.72, 89.90, 91.82, 90.86, 89.38, 89.88, 89.29, 88.79)
)

both <- bind_rows(pssm, mhcflurry) %>%
  mutate(Model = factor(Model, levels = rev(c("Deep FFN", "FFN (Jessen)", "CNN", "ResNet", "LSTM", "Random Forest"))))

# Plot 1a: Accuracy bars
p1a <- ggplot(both, aes(x = Accuracy, y = Model, fill = Source)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  geom_text(aes(label = sprintf("%.1f%%", Accuracy), group = Source),
            position = position_dodge(width = 0.7), hjust = -0.15, size = 3.2, fontface = 'bold') +
  scale_fill_manual(values = c("PSSM" = PSSM_COL, "MHCflurry" = MHCFL_COL)) +
  labs(title = 'Model Accuracy Comparison', subtitle = 'PSSM vs MHCflurry-derived labels | HLA-A*02:01',
       x = 'Accuracy (%)', y = NULL, fill = 'Label Source') +
  xlim(0, 100) + nature_theme

ggsave('figures/Figure1_model_comparison.pdf', p1a, width = 10, height = 5)

# Plot 1b: Per-class F1
deep_ffn <- both %>% filter(Model == "Deep FFN") %>%
  select(Source, NB_F1, WB_F1, SB_F1) %>%
  pivot_longer(-Source, names_to = "Class", values_to = "F1") %>%
  mutate(Class = factor(Class, levels = c("SB_F1", "WB_F1", "NB_F1"),
                        labels = c("SB", "WB", "NB")))

p1b <- ggplot(deep_ffn, aes(x = Class, y = F1, fill = Source)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  geom_text(aes(label = sprintf("%.3f", F1), group = Source),
            position = position_dodge(width = 0.7), vjust = -0.5, size = 4, fontface = 'bold') +
  scale_fill_manual(values = c("PSSM" = PSSM_COL, "MHCflurry" = MHCFL_COL)) +
  labs(title = 'Per-Class F1: Deep FFN', subtitle = 'Weak Binder class shows largest PSSM-MHCflurry gap',
       y = 'F1 Score', x = NULL, fill = 'Label Source') +
  ylim(0, 1.05) + nature_theme

ggsave('figures/FigureS1_per_class_f1.pdf', p1b, width = 7, height = 4.5)

# Plot 1c: CV folds
p1c <- ggplot(cv_data, aes(x = factor(Fold), y = Accuracy, color = Source, group = Source)) +
  geom_line(linewidth = 1) + geom_point(size = 3) +
  geom_hline(data = cv_data %>% group_by(Source) %>% summarise(m = mean(Accuracy)),
             aes(yintercept = m, color = Source), linetype = 'dashed', alpha = 0.4) +
  scale_color_manual(values = c("PSSM" = PSSM_COL, "MHCflurry" = MHCFL_COL)) +
  labs(title = '5-Fold Cross-Validation', subtitle = 'Deep FFN | PSSM CV=91.1±0.9%, MHCflurry CV=89.6±0.8%',
       x = 'Fold', y = 'Accuracy (%)', color = 'Label Source') +
  nature_theme

ggsave('figures/FigureS2_cv_folds.pdf', p1c, width = 7, height = 4.5)

cat('  -> Figure1, S1, S2 saved\n')

############################################################
# FIGURE 2: IEDB Benchmark (ROC + Confusion Matrix)
############################################################
cat('[2/8] IEDB benchmark...\n')

res <- read_csv("Data/cleaned/iedb_benchmark_results_v2.csv", show_col_types = FALSE)

# ROC
roc_obj <- roc(res$true_binary, res$pred_binary)

pdf("figures/Figure2_benchmark_roc.pdf", width = 8, height = 7)
par(bg = 'white', family = 'sans')
plot.roc(roc_obj,
         main = '', col = NATURE_COLORS['blue'], lwd = 3,
         print.auc = TRUE, print.auc.cex = 1.6, print.auc.col = '#222222',
         print.auc.x = 0.48, print.auc.y = 0.35,
         auc.polygon = TRUE, auc.polygon.col = '#377EB820',
         grid = TRUE, grid.col = '#E5E5E5', grid.lty = 1,
         legacy.axes = TRUE)
title(main = expression(bold('IEDB Benchmark ROC Curve')), line = 2.5, cex.main = 1.8, col.main = '#222222')
title(main = 'Deep FFN | HLA-A*02:01 | 49 epitopes + 20 controls', line = 1, cex.main = 1, col.main = '#666666')
dev.off()

# Confusion matrix
cm <- res %>%
  mutate(Actual = factor(ifelse(true_label == "POS", "True Epitope", "Non-Binder"),
                         levels = c("True Epitope", "Non-Binder")),
         Predicted = factor(pred_class_v2, levels = c("SB", "WB", "NB")))
cm_table <- table(Actual = cm$Actual, Predicted = cm$Predicted)
cm_df <- as.data.frame(as.table(cm_table))
total_pos <- sum(cm_table["True Epitope", ])
total_neg <- sum(cm_table["Non-Binder", ])
cm_df <- cm_df %>% mutate(
  pct = case_when(Actual == "True Epitope" ~ Freq / total_pos, Actual == "Non-Binder" ~ Freq / total_neg),
  label = paste0(Freq, "\n(", round(pct * 100, 1), "%)"),
  fill_intensity = Freq / max(Freq))

p_cm <- ggplot(cm_df, aes(x = Predicted, y = Actual)) +
  geom_tile(aes(fill = fill_intensity), color = 'white', linewidth = 1.2) +
  geom_text(aes(label = label), size = 5.5, fontface = 'bold', color = '#333333') +
  scale_fill_gradient(low = '#F0F4F8', high = NATURE_COLORS['blue'], guide = 'none') +
  labs(title = 'Confusion Matrix — IEDB Benchmark',
       subtitle = 'Sensitivity: 93.9% | Specificity: 75.0% | F1: 0.920\nSB + WB = Predicted Binder | NB = Predicted Non-Binder',
       x = 'Predicted Class', y = 'Actual Class') +
  nature_theme + theme(panel.grid = element_blank())

ggsave('figures/FigureS3_confusion_matrix.pdf', p_cm, width = 7.5, height = 5)
cat('  -> Figure2, S3 saved\n')

############################################################
# FIGURE 3: Feature Maps (BLOSUM62 encodings)
############################################################
cat('[3/8] BLOSUM62 feature maps...\n')

# Simplified BLOSUM62 heatmap for top epitopes
# Using the Nature blue-red diverging palette
blosum62 <- matrix(c(
   4,-1,-2,-2, 0,-1,-1, 0,-2,-1,-1,-1,-1,-2,-1, 1, 0,-3,-2, 0,
  -1, 5, 0,-2,-3, 1, 0,-2, 0,-3,-2, 2,-1,-3,-2,-1,-1,-3,-2,-3,
  -2, 0, 6, 1,-3, 0, 0, 0, 1,-3,-3, 0,-2,-3,-2, 1, 0,-4,-2,-3,
  -2,-2, 1, 6,-3, 0, 2,-1,-1,-3,-4,-1,-3,-3,-1, 0,-1,-4,-3,-3,
   0,-3,-3,-3, 9,-3,-4,-3,-3,-1,-1,-3,-1,-2,-3,-1,-1,-2,-2,-1,
  -1, 1, 0, 0,-3, 5, 2,-2, 0,-3,-2, 1, 0,-3,-1, 0,-1,-2,-1,-2,
  -1, 0, 0, 2,-4, 2, 5,-2, 0,-3,-3, 1,-2,-3,-1, 0,-1,-3,-2,-2,
   0,-2, 0,-1,-3,-2,-2, 6,-2,-4,-4,-2,-3,-3,-2, 0,-2,-2,-3,-3,
  -2, 0, 1,-1,-3, 0, 0,-2, 8,-3,-3,-1,-2,-1,-2,-1,-2, 2, 2,-3,
  -1,-3,-3,-3,-1,-3,-3,-4,-3, 4, 2,-3, 1, 0,-3,-2,-1,-3,-1, 3,
  -1,-2,-3,-4,-1,-2,-3,-4,-3, 2, 4,-2, 2, 0,-3,-2,-1,-2,-1, 1,
  -1, 2, 0,-1,-3, 1, 1,-2,-1,-3,-2, 5,-1,-3,-1, 0,-1,-3,-2,-2,
  -1,-1,-2,-3,-1, 0,-2,-3,-2, 1, 2,-1, 5, 0,-2,-1,-1,-1,-1, 1,
  -2,-3,-3,-3,-2,-3,-3,-3,-1, 0, 0,-3, 0, 6,-4,-2,-2, 1, 3,-1,
  -1,-2,-2,-1,-3,-1,-1,-2,-2,-3,-3,-1,-2,-4, 7,-1,-1,-4,-3,-2,
   1,-1, 1, 0,-1, 0, 0, 0,-1,-2,-2, 0,-1,-2,-1, 4, 1,-3,-2,-2,
   0,-1, 0,-1,-1,-1,-1,-2,-2,-1,-1,-1,-1,-2,-1, 1, 5,-2,-2, 0,
  -3,-3,-4,-4,-2,-2,-3,-2, 2,-3,-2,-3,-1, 1,-4,-3,-2,11, 2,-3,
  -2,-2,-2,-3,-2,-1,-2,-3, 2,-1,-1,-2,-1, 3,-3,-2,-2, 2, 7,-1,
   0,-3,-3,-3,-1,-2,-2,-3,-3, 3, 1,-2, 1,-1,-2,-2, 0,-3,-1, 4
), nrow = 20, ncol = 20, byrow = TRUE)

aas <- c('A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V')
blosum_norm <- t(apply(blosum62, 1, function(x) (x - min(x)) / (max(x) - min(x) + 1e-8)))
blosum_df <- expand.grid(AA = aas, Position = 1:20)
blosum_df$Score <- as.vector(blosum_norm)
blosum_df$Position <- factor(blosum_df$Position)
blosum_df$AA <- factor(blosum_df$AA, levels = rev(aas))

p_blosum <- ggplot(blosum_df, aes(x = Position, y = AA, fill = Score)) +
  geom_tile(color = 'white', linewidth = 0.3) +
  scale_fill_gradient2(low = '#377EB8', mid = 'white', high = '#E41A1C', midpoint = 0.5,
                       name = 'Normalized\nScore') +
  labs(title = 'BLOSUM62 Substitution Matrix', subtitle = 'Row-normalized | Blue=low, Red=high',
       x = 'Amino Acid Position', y = NULL) +
  nature_theme + theme(panel.grid = element_blank(), legend.position = 'right')

ggsave('figures/Figure3_blosum62_heatmap.pdf', p_blosum, width = 8, height = 6)
cat('  -> Figure3 saved\n')

############################################################
# FIGURE 4: Mutation Delta Scores
############################################################
cat('[4/8] Mutation analysis...\n')

mut <- read_csv("Data/cleaned/mutation_scan_results.csv", show_col_types = FALSE)

delta <- mut %>%
  filter(effect != "unchanged") %>%
  mutate(
    mut_clean = str_replace(mutation, " -- .*", ""),
    protein = ifelse(str_detect(mut_clean, "p53"), "p53", "KRAS"),
    label = paste0(mut_clean, "  ", wt_peptide, " > ", mut_peptide),
    effect_f = factor(effect, levels = c("CREATED (neoepitope)", "ENHANCED", "DESTROYED", "WEAKENED")))

# Delta score bar chart
p_mut <- ggplot(delta, aes(x = reorder(label, delta_score), y = delta_score)) +
  geom_col(aes(fill = effect_f), width = 0.55, alpha = 0.9) +
  geom_hline(yintercept = 0, linewidth = 0.6, color = '#999999') +
  geom_text(aes(label = sprintf("%+.3f", delta_score),
                y = ifelse(delta_score > 0, delta_score + 0.03, delta_score - 0.03)),
            size = 3.8, fontface = 'bold', color = '#333333') +
  scale_fill_manual(values = c("CREATED (neoepitope)" = CREATED_COL, "ENHANCED" = ENHANCE_COL,
                                "DESTROYED" = DESTROY_COL, "WEAKENED" = WEAKEN_COL)) +
  coord_flip() +
  labs(title = 'Mutation Epitope Score Changes — p53 & KRAS Hotspots',
       subtitle = 'HLA-A*02:01 | Green=Created, Blue=Enhanced, Red=Destroyed',
       x = NULL, y = expression(Delta*' Binding Score'), fill = 'Effect') +
  nature_theme + theme(panel.grid.major.y = element_blank())

ggsave('figures/Figure4_mutation_delta.pdf', p_mut, width = 12, height = 5)

# Paired WT vs Mutant
paired <- delta %>%
  select(mut_clean, wt_peptide, mut_peptide, wt_score, mut_score, effect_f) %>%
  pivot_longer(c(wt_score, mut_score), names_to = "type", values_to = "score") %>%
  mutate(type = ifelse(type == "wt_score", "Wild-Type", "Mutant"),
         type = factor(type, levels = c("Wild-Type", "Mutant")),
         label = paste0(mut_clean, "\n", wt_peptide))

p_paired <- ggplot(paired, aes(x = type, y = score, group = label)) +
  geom_line(aes(color = effect_f), linewidth = 1, alpha = 0.6) +
  geom_point(aes(color = effect_f, shape = type), size = 3.5) +
  scale_color_manual(values = c("CREATED (neoepitope)" = CREATED_COL, "ENHANCED" = ENHANCE_COL,
                                 "DESTROYED" = DESTROY_COL)) +
  scale_shape_manual(values = c("Wild-Type" = 1, "Mutant" = 16)) +
  scale_y_continuous(limits = c(0, 0.9)) +
  facet_wrap(~label, nrow = 2) +
  labs(title = 'WT vs Mutant Binding Scores', subtitle = 'HLA-A*02:01 | p53 & KRAS hotspots',
       x = NULL, y = 'Binding Score', color = 'Effect', shape = NULL) +
  nature_theme + theme(panel.spacing = unit(1, 'lines'))

ggsave('figures/FigureS4_mutation_paired.pdf', p_paired, width = 14, height = 6.5)

cat('  -> Figure4, S4 saved\n')

############################################################
# FIGURE 5: Protein Mutation Map
############################################################
cat('[5/8] Protein mutation map...\n')

p53_df <- tibble(pos = c(175,220,245,248,249,273,282), protein = 'p53',
  label = c('R175H','Y220C','G245S','R248W*','R249S','R273H','R282W'))
kras_df <- tibble(pos = c(12,13,61,146), protein = 'KRAS',
  label = c('G12V/C/D/R*','G13D*','Q61H/L/R','A146T*'))
p53_df$has_effect <- p53_df$label %in% c('R248W*','R249S')
kras_df$has_effect <- kras_df$label %in% c('G12V/C/D/R*','G13D*','A146T*')
map_df <- bind_rows(p53_df, kras_df)

p_map <- ggplot(map_df, aes(x = pos, y = protein, color = has_effect)) +
  geom_segment(data = tibble(protein = c('p53','KRAS'), start = c(1,1), end = c(323,188)),
               aes(x = start, xend = end, y = protein, yend = protein),
               color = '#DDDDDD', linewidth = 3, inherit.aes = FALSE) +
  geom_point(size = 5) +
  geom_text_repel(aes(label = label), size = 3.0, fontface = 'bold',
                  direction = 'both', nudge_y = 0.3, max.overlaps = 20,
                  segment.size = 0.3, segment.color = '#999999', min.segment.length = 0.1) +
  scale_color_manual(values = c('TRUE' = NATURE_COLORS['red'], 'FALSE' = '#BBBBBB'),
                     labels = c('TRUE' = 'Epitope-altering', 'FALSE' = 'No effect'), name = NULL) +
  scale_x_continuous(limits = c(-10, 400)) +
  scale_y_discrete(expand = expansion(mult = c(0.4, 0.4))) +
  labs(title = 'Mutation Hotspot Map — p53 & KRAS',
       subtitle = '* = epitope-altering | Red = alters MHC-I binding prediction',
       x = 'Amino Acid Position', y = NULL) +
  nature_theme

ggsave('figures/Figure5_mutation_map.pdf', p_map, width = 12, height = 3.5)
cat('  -> Figure5 saved\n')

############################################################
# FIGURE 6: Pan-Cancer Oncoplot (Nature style)
############################################################
cat('[6/8] Pan-cancer oncoplot...\n')

cancer_types <- c('PAAD\nPancreatic', 'COAD\nColorectal', 'LUAD\nLung Adeno',
                  'LUSC\nLung Squam', 'STAD\nStomach', 'UCEC\nEndometrial',
                  'OV\nOvarian', 'LIHC\nLiver', 'BRCA\nBreast')
kras_g12v <- c(23.0, 7.2, 5.0, 0.8, 1.8, 2.7, 0.3, 0.2, 0.1)
tp53_r248w <- c(1.0, 3.6, 1.0, 1.6, 2.8, 0.8, 1.0, 0.3, 0.5)

x <- 1:9; w <- 0.35
oncoplot_df <- bind_rows(
  tibble(cancer = cancer_types, freq = kras_g12v, mut = 'KRAS G12V', x = x - w/2),
  tibble(cancer = cancer_types, freq = tp53_r248w, mut = 'TP53 R248W', x = x + w/2))

p_onco <- ggplot(oncoplot_df, aes(x = x, y = freq, fill = mut)) +
  geom_col(position = position_identity(), width = w, alpha = 0.9) +
  geom_text(aes(label = ifelse(freq > 1, paste0(freq,'%'), ''), y = freq + 0.6),
            size = 2.8, fontface = 'bold', color = '#333333') +
  scale_fill_manual(values = c('KRAS G12V' = NATURE_COLORS['red'],
                                'TP53 R248W' = NATURE_COLORS['blue'])) +
  scale_x_continuous(breaks = x, labels = cancer_types) +
  labs(title = 'KRAS G12V and TP53 R248W Pan-Cancer Mutation Frequencies',
       subtitle = 'TCGA PanCancer Atlas | Data: cBioPortal / GENIE / COSMIC',
       x = NULL, y = 'Mutation Frequency (% of All Tumors)', fill = NULL) +
  ylim(0, 27) + nature_theme

ggsave('figures/Figure6_pan_cancer_oncoplot.pdf', p_onco, width = 12, height = 5.5)
cat('  -> Figure6 saved\n')

############################################################
# FIGURE 7: ESM-2 vs BLOSUM62 Comparison
############################################################
cat('[7/8] ESM-2 encoding comparison...\n')

esm_df <- tibble(
  Encoding = factor(c('BLOSUM62', 'ESM-2 t6\nmean pool', 'ESM-2 t12\nper-pos', 'ESM-2 t12\nper-pos reg', 'ESM-2 t6\nper-pos'),
                    levels = rev(c('BLOSUM62', 'ESM-2 t6\nmean pool', 'ESM-2 t12\nper-pos', 'ESM-2 t12\nper-pos reg', 'ESM-2 t6\nper-pos'))),
  Accuracy = c(91.9, 65.9, 90.9, 91.5, 93.3),
  Dims = c(180, 320, 4320, 4320, 2880),
  Type = c('Baseline', 'PLM failure', 'PLM overfit', 'PLM overfit reg', 'PLM best'))

p_esm <- ggplot(esm_df, aes(x = Accuracy, y = Encoding, fill = Type)) +
  geom_col(width = 0.6, alpha = 0.9) +
  geom_text(aes(label = sprintf('%.1f%% (%d-dim)', Accuracy, Dims)), hjust = -0.1, size = 3.2, fontface = 'bold') +
  geom_vline(xintercept = 91.9, linetype = 'dashed', color = '#999999', linewidth = 0.5) +
  scale_fill_manual(values = c('Baseline' = '#999999', 'PLM failure' = NATURE_COLORS['red'],
                                'PLM overfit' = NATURE_COLORS['orange'],
                                'PLM overfit reg' = '#E6A817',
                                'PLM best' = NATURE_COLORS['green'])) +
  labs(title = 'Protein Language Model vs BLOSUM62 Encoding',
       subtitle = 'ESM-2 per-position embeddings | Peptide-MHC Binding Prediction (HLA-A*02:01)',
       x = 'Test Accuracy (%)', y = NULL, fill = NULL) +
  xlim(60, 98) + nature_theme + theme(legend.position = 'right')

ggsave('figures/Figure7_esm2_comparison.pdf', p_esm, width = 10, height = 4.5)
cat('  -> Figure7 saved\n')

############################################################
# FIGURE 8: Real Negative Benchmark
############################################################
cat('[8/8] Real negative benchmark...\n')

neg_df <- tibble(
  Source = rep(c('Random\n9-mers', 'IEDB\nnon-binders', 'MHCflurry\nNB', 'Protein\nwindows', 'Homo-\npolymers'), 2),
  Specificity = c(0.94, 0.867, 0.96, 0.98, 0.80, 0.84, 0.867, 0.92, 0.94, 0.75),
  Classifier = rep(c('Binding Score >= 0.95', 'Argmax (baseline)'), each = 5))

p_neg <- ggplot(neg_df, aes(x = Source, y = Specificity, fill = Classifier)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6, alpha = 0.9) +
  geom_text(aes(label = sprintf('%.0f%%', Specificity * 100), group = Classifier),
            position = position_dodge(width = 0.7), vjust = -0.5, size = 3.5, fontface = 'bold') +
  geom_hline(yintercept = 0.80, linetype = 'dashed', color = '#999999', linewidth = 0.5) +
  scale_fill_manual(values = c('Binding Score >= 0.95' = NATURE_COLORS['green'],
                                'Argmax (baseline)' = NATURE_COLORS['gray'])) +
  labs(title = 'Per-Source Specificity: Real Negative Benchmark',
       subtitle = '4-source diverse negatives (n=185) + 20 homopolymer controls | Deep FFN',
       x = 'Negative Peptide Source', y = 'Specificity', fill = 'Classifier') +
  ylim(0, 1.1) + nature_theme

ggsave('figures/FigureS5_negative_benchmark.pdf', p_neg, width = 9, height = 5)
cat('  -> FigureS5 saved\n')

############################################################
# DONE
############################################################
cat('\n========================================\n')
cat('  ALL FIGURES REGENERATED (Nature style)\n')
cat('  Output: Analysis/figures/\n')
cat('========================================\n')
cat('\nFigure list:\n')
cat('  Figure1_model_comparison.pdf\n')
cat('  FigureS1_per_class_f1.pdf\n')
cat('  FigureS2_cv_folds.pdf\n')
cat('  Figure2_benchmark_roc.pdf\n')
cat('  FigureS3_confusion_matrix.pdf\n')
cat('  Figure3_blosum62_heatmap.pdf\n')
cat('  Figure4_mutation_delta.pdf\n')
cat('  FigureS4_mutation_paired.pdf\n')
cat('  Figure5_mutation_map.pdf\n')
cat('  Figure6_pan_cancer_oncoplot.pdf\n')
cat('  Figure7_esm2_comparison.pdf\n')
cat('  FigureS5_negative_benchmark.pdf\n')
