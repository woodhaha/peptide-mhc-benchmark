library(tidyverse)

# ---- Data ----
pssm <- tribble(
  ~Model,            ~Accuracy, ~Macro_F1, ~NB_F1, ~WB_F1, ~SB_F1,
  "Deep FFN",        94.8,      0.948,     0.976,  0.925,  0.944,
  "FFN (Jessen)",    93.2,      0.932,     0.970,  0.900,  0.926,
  "CNN",             93.2,      0.932,     0.976,  0.901,  0.920,
  "ResNet",          88.4,      0.881,     0.917,  0.813,  0.914,
  "LSTM",            87.8,      0.879,     0.954,  0.825,  0.858,
  "Random Forest",   87.3,      0.875,     0.954,  0.823,  0.848
) %>% mutate(Source = "PSSM")

mhcflurry <- tribble(
  ~Model,            ~Accuracy, ~Macro_F1, ~NB_F1, ~WB_F1, ~SB_F1,
  "Deep FFN",        91.9,      0.921,     0.969,  0.880,  0.913,
  "FFN (Jessen)",    90.9,      0.911,     0.969,  0.875,  0.891,
  "CNN",             90.0,      0.901,     0.963,  0.860,  0.882,
  "ResNet",          84.3,      0.847,     0.907,  0.792,  0.843,
  "LSTM",            83.3,      0.836,     0.935,  0.752,  0.822,
  "Random Forest",   81.1,      0.814,     0.927,  0.723,  0.793
) %>% mutate(Source = "MHCflurry")

cv_data <- tribble(
  ~Source,     ~Fold, ~Accuracy,
  "PSSM",      1,     90.42,
  "PSSM",      2,     91.71,
  "PSSM",      3,     91.72,
  "PSSM",      4,     89.90,
  "PSSM",      5,     91.82,
  "MHCflurry", 1,     90.86,
  "MHCflurry", 2,     89.38,
  "MHCflurry", 3,     89.88,
  "MHCflurry", 4,     89.29,
  "MHCflurry", 5,     88.79
)

both <- bind_rows(pssm, mhcflurry) %>%
  mutate(Model = factor(Model, levels = rev(c(
    "Deep FFN", "FFN (Jessen)", "CNN", "ResNet", "LSTM", "Random Forest"
  ))))

# ---- Plot 1: Accuracy comparison ----
p1 <- ggplot(both, aes(x = Accuracy, y = Model, fill = Source)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  geom_text(aes(label = sprintf("%.1f%%", Accuracy),
                group = Source),
            position = position_dodge(width = 0.7),
            hjust = -0.15, size = 3.8, fontface = "bold") +
  scale_fill_manual(values = c("PSSM" = "#4472C4", "MHCflurry" = "#ED7D31")) +
  labs(title = "Model Accuracy: PSSM vs MHCflurry",
       subtitle = "Peptide–MHC Class I Binding Prediction (HLA-A*02:01)",
       x = "Accuracy (%)", y = NULL, fill = "Label Source") +
  xlim(0, 100) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "bottom",
        plot.title = element_text(face = "bold"),
        panel.grid.major.y = element_blank())

# ---- Plot 2: Per-class F1 (Deep FFN) ----
deep_ffn <- both %>% filter(Model == "Deep FFN") %>%
  select(Source, NB_F1, WB_F1, SB_F1) %>%
  pivot_longer(-Source, names_to = "Class", values_to = "F1") %>%
  mutate(Class = factor(Class, levels = c("SB_F1", "WB_F1", "NB_F1"),
                        labels = c("SB\n(Strong)", "WB\n(Weak)", "NB\n(Non)")))

p2 <- ggplot(deep_ffn, aes(x = Class, y = F1, fill = Source)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  geom_text(aes(label = sprintf("%.3f", F1), group = Source),
            position = position_dodge(width = 0.7),
            vjust = -0.5, size = 4.5, fontface = "bold") +
  scale_fill_manual(values = c("PSSM" = "#4472C4", "MHCflurry" = "#ED7D31")) +
  labs(title = "Per-Class F1: Deep FFN (Best Model)",
       subtitle = "PSSM vs MHCflurry — Weak Binder class shows largest gap",
       y = "F1 Score", x = NULL, fill = "Label Source") +
  ylim(0, 1.05) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "bottom",
        plot.title = element_text(face = "bold"))

# ---- Plot 3: CV fold consistency ----
p3 <- ggplot(cv_data, aes(x = factor(Fold), y = Accuracy,
                          color = Source, group = Source)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3.5) +
  geom_hline(data = cv_data %>% group_by(Source) %>% summarise(m = mean(Accuracy)),
             aes(yintercept = m, color = Source), linetype = "dashed", alpha = 0.5) +
  geom_text(data = cv_data %>% group_by(Source) %>% summarise(m = mean(Accuracy)),
            aes(y = m, label = sprintf("μ=%.1f%%", m), x = 0.7, color = Source),
            hjust = 0, size = 4, fontface = "italic") +
  scale_color_manual(values = c("PSSM" = "#4472C4", "MHCflurry" = "#ED7D31")) +
  labs(title = "5-Fold Cross-Validation: Fold-by-Fold",
       subtitle = "PSSM tighter variance (±0.9%) vs MHCflurry (±0.8%)",
       x = "CV Fold", y = "Accuracy (%)", color = "Label Source") +
  ylim(87, 93) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "bottom",
        plot.title = element_text(face = "bold"))

# ---- Plot 4: Macro F1 gap across all models ----
both_long <- both %>%
  mutate(gap = Accuracy - lag(Accuracy, default = first(Accuracy)),
         .by = Model) %>%
  filter(Source == "MHCflurry") %>%
  mutate(f1_gap = NA_real_)

# Compute the PSSM - MHCflurry accuracy gap per model
gap_df <- both %>%
  select(Model, Source, Accuracy) %>%
  pivot_wider(names_from = Source, values_from = Accuracy) %>%
  mutate(Gap = PSSM - MHCflurry)

p4 <- ggplot(gap_df, aes(x = Gap, y = reorder(Model, Gap))) +
  geom_col(fill = "#C00000", alpha = 0.7, width = 0.55) +
  geom_text(aes(label = sprintf("+%.1f%%", Gap)), hjust = -0.15,
            size = 4.5, fontface = "bold") +
  labs(title = "Accuracy Gap (PSSM − MHCflurry)",
       subtitle = "How much harder is real ML-predicted binding data?",
       x = "Accuracy Gap (%)", y = NULL) +
  xlim(0, 7) +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold"),
        panel.grid.major.y = element_blank())

# ---- Save ----
dir.create("plots", showWarnings = FALSE)
ggsave("plots/comparison_accuracy.png",        p1, width = 10, height = 5,  dpi = 150)
ggsave("plots/comparison_per_class_f1.png",    p2, width = 8,  height = 5,  dpi = 150)
ggsave("plots/comparison_cv_folds.png",        p3, width = 8,  height = 5,  dpi = 150)
ggsave("plots/comparison_gap.png",             p4, width = 8,  height = 4.5, dpi = 150)

cat("Plots saved to plots/comparison_*.png\n")
