# Mutation delta score plot
library(tidyverse)
setwd('D:/Researching/Peptide epitope')

res <- read_csv("data/mutation_scan_results.csv", show_col_types = FALSE)

# Filter to mutation-affected positions only
delta <- res %>%
  filter(effect != "unchanged") %>%
  mutate(
    mut_clean = str_replace(mutation, " -- .*", ""),
    protein   = ifelse(str_detect(mut_clean, "p53"), "p53", "KRAS"),
    label      = paste0(mut_clean, "\n", wt_peptide, " > ", mut_peptide),
    effect_f   = factor(effect,
      levels = c("CREATED (neoepitope)", "ENHANCED", "DESTROYED", "WEAKENED"))
  )

# ---- Plot 1: Delta score dot plot ----
p1 <- ggplot(delta, aes(x = reorder(label, delta_score), y = delta_score)) +
  geom_col(aes(fill = effect_f), width = 0.55, alpha = 0.9) +
  geom_hline(yintercept = 0, linewidth = 0.8, color = "gray40") +
  geom_text(aes(
    label = sprintf("%+.3f", delta_score),
    y = ifelse(delta_score > 0, delta_score + 0.03, delta_score - 0.03),
    color = effect_f
  ), size = 4.5, fontface = "bold") +
  scale_fill_manual(values = c(
    "CREATED (neoepitope)" = "#228B22",
    "ENHANCED"             = "#4472C4",
    "DESTROYED"            = "#DC143C",
    "WEAKENED"             = "#FF8C00"
  )) +
  scale_color_manual(values = c(
    "CREATED (neoepitope)" = "#228B22",
    "ENHANCED"             = "#4472C4",
    "DESTROYED"            = "#DC143C",
    "WEAKENED"             = "#FF8C00"
  ), guide = "none") +
  coord_flip() +
  labs(
    title    = "Mutation Epitope Score Changes -- p53 & KRAS Hotspots",
    subtitle = paste0(
      "HLA-A*02:01 binding score delta (mutant - wild-type) | ",
      "Green = neoepitope created, Red = destroyed, Blue = enhanced"
    ),
    x = NULL,
    y = expression(Delta * " Binding Score (SB*1.0 + WB*0.5)"),
    fill = "Effect"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title        = element_text(face = "bold", size = 16),
    plot.subtitle     = element_text(size = 10, color = "gray40"),
    axis.text.y       = element_text(size = 9),
    legend.position   = "bottom",
    panel.grid.major.y = element_blank()
  )

ggsave("plots/mutation_delta_scores.png", p1,
       width = 12, height = 5.5, dpi = 180)

# ---- Plot 2: Paired WT vs Mutant score ----
paired <- delta %>%
  select(mut_clean, wt_peptide, mut_peptide, wt_score, mut_score, effect_f) %>%
  pivot_longer(c(wt_score, mut_score), names_to = "type", values_to = "score") %>%
  mutate(
    type = ifelse(type == "wt_score", "Wild-Type", "Mutant"),
    type = factor(type, levels = c("Wild-Type", "Mutant")),
    label = paste0(mut_clean, "\n", wt_peptide)
  )

p2 <- ggplot(paired, aes(x = type, y = score, group = label)) +
  geom_line(aes(color = effect_f), linewidth = 1.2, alpha = 0.7) +
  geom_point(aes(color = effect_f, shape = type), size = 4) +
  scale_color_manual(values = c(
    "CREATED (neoepitope)" = "#228B22",
    "ENHANCED"             = "#4472C4",
    "DESTROYED"            = "#DC143C",
    "WEAKENED"             = "#FF8C00"
  )) +
  scale_shape_manual(values = c("Wild-Type" = 1, "Mutant" = 16)) +
  scale_y_continuous(limits = c(0, 0.9)) +
  facet_wrap(~label, nrow = 2) +
  labs(
    title    = "WT vs Mutant Binding Scores -- p53 & KRAS Hotspots",
    subtitle = "HLA-A*02:01 | Green=Neoepitope, Red=Destroyed, Blue=Enhanced",
    x = NULL, y = "Binding Score",
    color = "Effect", shape = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title    = element_text(face = "bold", size = 16),
    strip.text    = element_text(size = 8.5, face = "bold"),
    legend.position = "bottom",
    panel.spacing = unit(1.2, "lines")
  )

ggsave("plots/mutation_wt_vs_mutant.png", p2,
       width = 14, height = 7, dpi = 180)

# ---- Plot 3: p53 and KRAS protein maps with mutation positions ----
p53_positions <- c(175, 220, 245, 248, 249, 273, 282)
kras_positions <- c(12, 13, 61, 146)

p53_df <- tibble(pos = p53_positions, protein = "p53",
  label = c("R175H", "Y220C", "G245S", "R248W*", "R249S", "R273H", "R282W"))
kras_df <- tibble(pos = kras_positions, protein = "KRAS",
  label = c("G12V/C/D/R*", "G13D*", "Q61H/L/R", "A146T*"))

# Mark which have epitope effects
p53_df$has_effect <- p53_df$label %in% c("R248W*", "R249S")
kras_df$has_effect <- kras_df$label %in% c("G12V/C/D/R*", "G13D*", "A146T*")

map_df <- bind_rows(p53_df, kras_df)

p3 <- ggplot(map_df, aes(x = pos, y = protein, color = has_effect)) +
  geom_segment(data = tibble(protein = c("p53","KRAS"), start = c(1,1), end = c(323,188)),
               aes(x = start, xend = end, y = protein, yend = protein),
               color = "gray70", linewidth = 3, inherit.aes = FALSE) +
  geom_point(size = 5) +
  geom_text(aes(label = label), hjust = -0.15, vjust = -0.8, size = 3.5,
            fontface = "bold") +
  scale_color_manual(values = c("TRUE" = "#DC143C", "FALSE" = "gray50"),
                     labels = c("TRUE" = "Epitope-altering", "FALSE" = "Silent"),
                     name = NULL) +
  scale_x_continuous(limits = c(-10, 370)) +
  labs(
    title = "Mutation Hotspot Map -- p53 & KRAS",
    subtitle = "* = epitope-altering mutation  |  Red = alters MHC-I binding",
    x = "Amino Acid Position", y = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "bottom"
  )

ggsave("plots/mutation_protein_map.png", p3,
       width = 12, height = 3.5, dpi = 180)

cat("Saved:\n")
cat("  plots/mutation_delta_scores.png\n")
cat("  plots/mutation_wt_vs_mutant.png\n")
cat("  plots/mutation_protein_map.png\n")
