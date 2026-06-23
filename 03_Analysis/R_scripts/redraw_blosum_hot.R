# Redraw BLOSUM feature images with "hot" color palette
library(tidyverse)
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
  "GILGFVFTL",   "M1 58-66 (Immunodominant)",           "Confirmed",
  "NLVPMVATV",   "pp65 495-503 (Immunodominant)",       "Confirmed",
  "RMFPNAPYL",   "WT1 126-134",                         "Confirmed",
  "LLGRNSFEV",   "p53 264-272",                         "Confirmed",
  "IMDQVPFSV",   "gp100 209-217",                       "Confirmed",
  "YMNGTMSQV",   "Tyrosinase 369-377",                  "Confirmed",
  "YLQPRTFLL",   "Spike 269-277 (Immunodominant)",      "Confirmed",
  "SLLMWITQC",   "NY-ESO-1 157-165",                    "Confirmed",
  "ALMDKSLHV",   "MART-1 56-64",                        "Novel (100% SB)",
  "KIADYNYKL",   "Spike RBD 87-95",                     "Novel (99.9% SB)",
  "RMPEAAPPV",   "p53 65-73",                           "Novel (100% SB)",
  "LLTEVETYV",   "M1 3-11",                             "Novel (100% SB)",
  "RLLQTGIHV",   "pp65 40-48",                          "Novel (99.9% SB)",
  "SLGEQQYSV",   "WT1 187-195",                         "Novel (100% SB)",
  "ILRGSVAHK",   "NP 265-273 (p9=K)",                   "False Negative",
  "QYDPVAALF",   "pp65 341-349 (p9=F)",                 "False Negative",
  "TLGIVCPIC",   "E7 86-94 (p9=C)",                     "False Negative",
  "AAAAAAAAA",   "Poly-A",                              "Negative Control",
  "LLLLLLLLL",   "Poly-L (False Positive)",              "False Positive",
  "DDDDDDDDD",   "Poly-D",                              "Negative Control"
)

# ---- Encode ----
x <- encode_blosum62(epitopes$peptide)

# ---- HOT color palette (MATLAB-style) ----
hot_colors <- colorRampPalette(
  c("#000000", "#8B0000", "#FF0000", "#FF4500", "#FF8C00", "#FFD700", "#FFFF00", "#FFFFFF")
)(100)

# ---- Plot each peptide as 9x20 heatmap with HOT colormap ----
plot_list <- list()
for (i in 1:nrow(epitopes)) {
  mat <- x[i,,]
  df <- expand.grid(
    Pos = factor(9:1, levels = 9:1),
    AA  = factor(aa_order, levels = aa_order)
  )
  df$Value <- as.vector(t(mat[nrow(mat):1, ]))

  p <- ggplot(df, aes(x = AA, y = Pos, fill = Value)) +
    geom_tile(color = "gray30", linewidth = 0.25) +
    scale_fill_gradientn(
      colors = hot_colors,
      limits = c(0, 1),
      name = "BLOSUM\nScore"
    ) +
    labs(
      title = sprintf("%s  |  %s", epitopes$peptide[i], epitopes$label[i]),
      subtitle = epitopes$type[i],
      x = "Amino Acid", y = "Position (p1-p9)"
    ) +
    theme_minimal(base_size = 9) +
    theme(
      plot.title        = element_text(size = 11, face = "bold", family = "mono"),
      plot.subtitle     = element_text(size = 8.5, color = "gray60"),
      axis.text.x       = element_text(size = 6.5, angle = 90, hjust = 1, vjust = 0.5),
      axis.text.y       = element_text(size = 8),
      axis.title.x      = element_text(size = 8, color = "gray50"),
      legend.position   = "right",
      legend.key.height = unit(0.7, "cm"),
      legend.title      = element_text(size = 7),
      legend.text       = element_text(size = 6),
      panel.grid        = element_blank(),
      plot.background   = element_rect(fill = "gray15", color = NA),
      plot.title.position = "panel"
    )
  plot_list[[i]] <- p
}

# ---- Arrange 5 rows x 4 cols ----
png("plots/feature_blosum_images_hot.png",
    width = 5000, height = 5200, res = 220)

gridExtra::grid.arrange(
  grobs = plot_list,
  ncol  = 4,
  top   = grid::textGrob(
    "BLOSUM62 Encoding of Top Epitopes -- HLA-A*02:01  |  HOT Colormap",
    gp = grid::gpar(fontface = "bold", fontsize = 22, col = "white")
  ),
  padding = unit(2, "line"),
  newpage = TRUE
)

# Dark background for the whole page
grid::grid.rect(gp = grid::gpar(fill = "gray15", col = NA))
grid::grid.draw(ggplot2:::ggplotGrob(
  cowplot::plot_grid(plotlist = plot_list, ncol = 4)
))

dev.off()

# Alternative: light-background version
plot_list_light <- list()
for (i in 1:nrow(epitopes)) {
  mat <- x[i,,]
  df <- expand.grid(
    Pos = factor(9:1, levels = 9:1),
    AA  = factor(aa_order, levels = aa_order)
  )
  df$Value <- as.vector(t(mat[nrow(mat):1, ]))

  p <- ggplot(df, aes(x = AA, y = Pos, fill = Value)) +
    geom_tile(color = "gray80", linewidth = 0.25) +
    scale_fill_gradientn(
      colors = hot_colors,
      limits = c(0, 1),
      name = "BLOSUM\nScore"
    ) +
    labs(
      title = sprintf("%s  |  %s", epitopes$peptide[i], epitopes$label[i]),
      subtitle = epitopes$type[i],
      x = "Amino Acid", y = "Position"
    ) +
    theme_minimal(base_size = 9) +
    theme(
      plot.title        = element_text(size = 11, face = "bold", family = "mono"),
      plot.subtitle     = element_text(size = 8.5, color = "gray40"),
      axis.text.x       = element_text(size = 6.5, angle = 90, hjust = 1, vjust = 0.5),
      axis.text.y       = element_text(size = 8),
      legend.position   = "right",
      legend.key.height = unit(0.8, "cm"),
      panel.grid        = element_blank()
    )
  plot_list_light[[i]] <- p
}

png("plots/feature_blosum_images.png",
    width = 5000, height = 5200, res = 220)

gridExtra::grid.arrange(
  grobs = plot_list_light,
  ncol  = 4,
  top   = grid::textGrob(
    "BLOSUM62 Encoding of Top Epitopes -- HLA-A*02:01  |  HOT Colormap",
    gp = grid::gpar(fontface = "bold", fontsize = 22)
  )
)

dev.off()

cat("Saved:\n")
cat("  plots/feature_blosum_images.png (light bg)\n")
cat("  plots/feature_blosum_images_hot.png (dark bg)\n")
