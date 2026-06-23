#!/usr/bin/env Rscript
# Redraw BLOSUM feature maps with soft, Nature-appropriate palettes
options(warn = -1)
Sys.setenv(TF_CPP_MIN_LOG_LEVEL = "3")
suppressMessages({library(keras); library(tensorflow); library(tidyverse); library(viridis)})
setwd('D:/Researching/Peptide epitope')

# ---- BLOSUM62 matrix ----
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

# ---- Top epitopes ----
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
  "ALMDKSLHV",   "MART-1 56-64",       "Novel (100% SB)",
  "KIADYNYKL",   "Spike RBD 87-95",    "Novel (99.9% SB)",
  "RMPEAAPPV",   "p53 65-73",          "Novel (100% SB)",
  "LLTEVETYV",   "M1 3-11",            "Novel (100% SB)",
  "RLLQTGIHV",   "pp65 40-48",         "Novel (99.9% SB)",
  "SLGEQQYSV",   "WT1 187-195",        "Novel (100% SB)",
  "ILRGSVAHK",   "NP 265-273 (MISS)",  "False Negative",
  "QYDPVAALF",   "pp65 341-349 (MISS)","False Negative",
  "TLGIVCPIC",   "E7 86-94 (MISS)",    "False Negative",
  "AAAAAAAAA",   "Poly-A (Neg Ctrl)",  "Negative Control",
  "LLLLLLLLL",   "Poly-L (FP)",        "False Positive",
  "DDDDDDDDD",   "Poly-D (Neg Ctrl)",  "Negative Control"
)

cat(sprintf("Encoding %d epitopes...\n", nrow(epitopes)))
x <- encode_blosum62(epitopes$peptide)

# ---- Nature-soft theme ----
nature_soft <- theme_minimal(base_size = 9) +
  theme(
    plot.title = element_text(size = 10, face = "bold", family = "mono", color = "#333333"),
    plot.subtitle = element_text(size = 8, color = "#777777"),
    axis.text.x = element_text(size = 6, angle = 90, hjust = 1, vjust = 0.5, color = "#555555"),
    axis.text.y = element_text(size = 7, color = "#555555"),
    legend.position = "right",
    legend.key.height = unit(0.8, "cm"),
    legend.text = element_text(size = 7),
    panel.grid = element_blank(),
    plot.background = element_rect(fill = "white", color = NA)
  )

dir.create("plots", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

# ==============================================================
# VERSION A: Soft blue gradient ("Gentle Blue")
# ==============================================================
cat("\n[A] Soft blue palette...\n")
plot_list_blue <- list()
for (i in 1:nrow(epitopes)) {
  mat <- x[i,,]
  df <- expand.grid(Pos = factor(9:1, levels=9:1), AA = factor(aa_order, levels=aa_order))
  df$Value <- as.vector(t(mat[nrow(mat):1, ]))

  p <- ggplot(df, aes(x = AA, y = Pos, fill = Value)) +
    geom_tile(color = "white", linewidth = 0.3) +
    scale_fill_gradient(low = "#F4F8FD", high = "#5B8DB8",
                        limits = c(0, 1), name = NULL) +
    labs(title = sprintf("%s  |  %s", epitopes$peptide[i], epitopes$label[i]),
         subtitle = epitopes$type[i], x = NULL, y = NULL) +
    nature_soft
  plot_list_blue[[i]] <- p
}

pdf("plots/feature_blosum_images.pdf", width = 20, height = 20)
gridExtra::grid.arrange(
  grobs = plot_list_blue, ncol = 4,
  top = grid::textGrob("BLOSUM62 Encoding of Top Epitopes — HLA-A*02:01",
                       gp = grid::gpar(fontface = "bold", fontsize = 20, col = "#333333"))
)
dev.off()
cat("  -> plots/feature_blosum_images.pdf\n")

# ==============================================================
# VERSION B: Soft warm gradient ("Gentle Warm") replaces HOT
# ==============================================================
cat("[B] Soft warm palette...\n")
plot_list_warm <- list()
for (i in 1:nrow(epitopes)) {
  mat <- x[i,,]
  df <- expand.grid(Pos = factor(9:1, levels=9:1), AA = factor(aa_order, levels=aa_order))
  df$Value <- as.vector(t(mat[nrow(mat):1, ]))

  p <- ggplot(df, aes(x = AA, y = Pos, fill = Value)) +
    geom_tile(color = "white", linewidth = 0.3) +
    scale_fill_gradient(low = "#FFF8F4", high = "#D4856B",
                        limits = c(0, 1), name = NULL) +
    labs(title = sprintf("%s  |  %s", epitopes$peptide[i], epitopes$label[i]),
         subtitle = epitopes$type[i], x = NULL, y = NULL) +
    nature_soft
  plot_list_warm[[i]] <- p
}

pdf("plots/feature_blosum_images_warm.pdf", width = 20, height = 20)
gridExtra::grid.arrange(
  grobs = plot_list_warm, ncol = 4,
  top = grid::textGrob("BLOSUM62 Encoding of Top Epitopes — HLA-A*02:01",
                       gp = grid::gpar(fontface = "bold", fontsize = 20, col = "#333333"))
)
dev.off()
cat("  -> plots/feature_blosum_images_warm.pdf\n")

# ==============================================================
# VERSION C: Soft viridis ("Gentle Viridis")
# ==============================================================
cat("[C] Viridis palette...\n")
plot_list_viridis <- list()
for (i in 1:nrow(epitopes)) {
  mat <- x[i,,]
  df <- expand.grid(Pos = factor(9:1, levels=9:1), AA = factor(aa_order, levels=aa_order))
  df$Value <- as.vector(t(mat[nrow(mat):1, ]))

  p <- ggplot(df, aes(x = AA, y = Pos, fill = Value)) +
    geom_tile(color = "white", linewidth = 0.3) +
    scale_fill_viridis_c(option = "D", limits = c(0, 1), name = NULL,
                         begin = 0.05, end = 0.95) +
    labs(title = sprintf("%s  |  %s", epitopes$peptide[i], epitopes$label[i]),
         subtitle = epitopes$type[i], x = NULL, y = NULL) +
    nature_soft
  plot_list_viridis[[i]] <- p
}

pdf("plots/feature_blosum_images_viridis.pdf", width = 20, height = 20)
gridExtra::grid.arrange(
  grobs = plot_list_viridis, ncol = 4,
  top = grid::textGrob("BLOSUM62 Encoding of Top Epitopes — HLA-A*02:01",
                       gp = grid::gpar(fontface = "bold", fontsize = 20, col = "#333333"))
)
dev.off()
cat("  -> plots/feature_blosum_images_viridis.pdf\n")

# ==============================================================
# VERSION D: Soft cool-warm diverging ("Gentle Diverging")
# Best for showing both low and high values elegantly
# ==============================================================
cat("[D] Soft diverging palette...\n")
plot_list_div <- list()
for (i in 1:nrow(epitopes)) {
  mat <- x[i,,]
  df <- expand.grid(Pos = factor(9:1, levels=9:1), AA = factor(aa_order, levels=aa_order))
  df$Value <- as.vector(t(mat[nrow(mat):1, ]))

  p <- ggplot(df, aes(x = AA, y = Pos, fill = Value)) +
    geom_tile(color = "white", linewidth = 0.3) +
    scale_fill_gradient2(low = "#7FAFD4", mid = "#FFFDF9", high = "#E8916A",
                         midpoint = 0.5, limits = c(0, 1), name = NULL) +
    labs(title = sprintf("%s  |  %s", epitopes$peptide[i], epitopes$label[i]),
         subtitle = epitopes$type[i], x = NULL, y = NULL) +
    nature_soft
  plot_list_div[[i]] <- p
}

pdf("plots/feature_blosum_images_diverging.pdf", width = 20, height = 20)
gridExtra::grid.arrange(
  grobs = plot_list_div, ncol = 4,
  top = grid::textGrob("BLOSUM62 Encoding of Top Epitopes — HLA-A*02:01",
                       gp = grid::gpar(fontface = "bold", fontsize = 20, col = "#333333"))
)
dev.off()
cat("  -> plots/feature_blosum_images_diverging.pdf\n")

# ---- Copy to figures ----
file.copy("plots/feature_blosum_images.pdf", "figures/feature_blosum_images.pdf", overwrite = TRUE)
file.copy("plots/feature_blosum_images_warm.pdf", "figures/feature_blosum_images_warm.pdf", overwrite = TRUE)
file.copy("plots/feature_blosum_images_viridis.pdf", "figures/feature_blosum_images_viridis.pdf", overwrite = TRUE)
file.copy("plots/feature_blosum_images_diverging.pdf", "figures/feature_blosum_images_diverging.pdf", overwrite = TRUE)
file.copy("plots/feature_blosum_images.pdf", "Analysis/figures/feature_blosum_images.pdf", overwrite = TRUE)
file.copy("plots/feature_blosum_images_warm.pdf", "Analysis/figures/feature_blosum_images_warm.pdf", overwrite = TRUE)

cat("\n========================================\n")
cat("  DONE — 4 gentle palettes generated\n")
cat("========================================\n")
cat("  A: Soft blue       → feature_blosum_images.pdf\n")
cat("  B: Soft warm       → feature_blosum_images_warm.pdf\n")
cat("  C: Viridis D       → feature_blosum_images_viridis.pdf\n")
cat("  D: Soft diverging  → feature_blosum_images_diverging.pdf\n")
