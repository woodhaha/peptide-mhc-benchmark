# Feature maps: BLOSUM encoding images + learned feature heatmap for top epitopes
options(warn = -1)
Sys.setenv(TF_CPP_MIN_LOG_LEVEL = "3")
suppressMessages({library(keras); library(tensorflow); library(tidyverse)})
setwd('D:/Researching/Peptide epitope')

# ---- 1. Load model ----
model <- load_model_hdf5("models/FFN_Deep.h5")

# ---- 2. BLOSUM62 matrix ----
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

# ---- 3. Top epitopes to visualize ----
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
x_4d <- array_reshape(x, c(dim(x)[1], 9, 20, 1))

# ---- 4. Plot 1: pep_plot_images() -- BLOSUM encoding as heatmaps ----
cat("Generating BLOSUM encoding images...\n")

plot_list <- list()
for (i in 1:nrow(epitopes)) {
  mat <- x[i,,]  # 9x20
  df <- expand.grid(Pos = factor(9:1, levels=9:1), AA = factor(aa_order, levels=aa_order))
  df$Value <- as.vector(t(mat[nrow(mat):1, ]))

  p <- ggplot(df, aes(x = AA, y = Pos, fill = Value)) +
    geom_tile(color = "white", linewidth = 0.3) +
    scale_fill_gradient(low = "white", high = "#4472C4",
                        limits = c(0, 1), name = "BLOSUM") +
    labs(
      title = sprintf("%s  |  %s", epitopes$peptide[i], epitopes$label[i]),
      subtitle = epitopes$type[i],
      x = NULL, y = "Position"
    ) +
    theme_minimal(base_size = 9) +
    theme(
      plot.title = element_text(size = 10, face = "bold", family = "mono"),
      plot.subtitle = element_text(size = 8, color = "gray40"),
      axis.text.x = element_text(size = 6, angle = 90, hjust = 1, vjust = 0.5),
      axis.text.y = element_text(size = 7),
      legend.position = "right",
      legend.key.height = unit(0.8, "cm"),
      panel.grid = element_blank()
    )
  plot_list[[i]] <- p
}

# Arrange 5 rows x 4 cols
png("plots/feature_blosum_images.png",
    width = 4800, height = 4800, res = 220)
gridExtra::grid.arrange(
  grobs = plot_list,
  ncol = 4,
  top = grid::textGrob(
    "BLOSUM62 Encoding of Top Epitopes -- HLA-A*02:01",
    gp = grid::gpar(fontface = "bold", fontsize = 20)
  )
)
dev.off()
cat("Saved: plots/feature_blosum_images.png\n")

# ---- 5. Extract learned features from penultimate layer ----
cat("Extracting learned features...\n")
layer_names <- sapply(model$layers, function(l) l$name)
dense_layers <- grep("dense_", layer_names, value = TRUE)
penultimate <- dense_layers[length(dense_layers) - 1]
cat(sprintf("Penultimate layer: %s\n", penultimate))

feature_model <- keras_model(
  inputs = model$input,
  outputs = get_layer(model, penultimate)$output
)

# Predict learned features
x_flat <- array_reshape(x_4d, c(dim(x_4d)[1], 180))
learned <- predict(feature_model, x_flat, verbose = 0)
cat(sprintf("Learned features shape: [%d, %d]\n", nrow(learned), ncol(learned)))

# Scale features for heatmap
learned_scaled <- scale(learned)
learned_scaled[learned_scaled > 2] <- 2
learned_scaled[learned_scaled < -2] <- -2

# ---- 6. Plot 2: Learned feature heatmap ----
cat("Generating learned feature heatmap...\n")
n_feats <- ncol(learned)

feat_df <- expand.grid(
  Peptide = epitopes$peptide,
  Feature = 1:n_feats
)
feat_df$Value <- as.vector(t(learned_scaled))

# Order peptides by type
feat_df$Peptide <- factor(feat_df$Peptide,
  levels = rev(epitopes$peptide))

# Annotate with type
type_colors <- c(
  "Confirmed" = "#228B22",
  "Novel (100% SB)" = "#4472C4",
  "Novel (99.9% SB)" = "#6495ED",
  "False Negative" = "#DC143C",
  "False Positive" = "#FF8C00",
  "Negative Control" = "#808080"
)

# Row annotation
row_ann <- epitopes %>%
  mutate(Peptide = factor(peptide, levels = rev(epitopes$peptide)))

png("plots/feature_learned_heatmap.png",
    width = 3600, height = 2200, res = 200)

p_heat <- ggplot(feat_df, aes(x = Feature, y = Peptide, fill = Value)) +
  geom_tile() +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       midpoint = 0, name = "Z-score") +
  # Add type annotation on left
  geom_text(data = row_ann,
            aes(x = -5, y = Peptide, label = type, color = type),
            inherit.aes = FALSE, hjust = 1, size = 3.5, fontface = "bold") +
  scale_color_manual(values = type_colors, guide = "none") +
  scale_x_continuous(expand = c(0, 0)) +
  labs(
    title = "Learned Feature Representations -- Deep FFN Penultimate Layer (90-dim)",
    subtitle = paste0("Top epitopes colored by validation status | ",
                      "Features from ", penultimate, " layer (", n_feats, "-dim)"),
    x = "Feature Index (1-90)", y = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(size = 10, color = "gray40"),
    axis.text.y = element_text(size = 9, family = "mono", face = "bold"),
    axis.text.x = element_text(size = 6),
    panel.grid = element_blank(),
    legend.position = "bottom"
  )

print(p_heat)
dev.off()
cat("Saved: plots/feature_learned_heatmap.png\n")

# ---- 7. Plot 3: Feature correlation with binding score ----
cat("Generating feature importance plot...\n")

preds <- predict(model, x_flat, verbose = 0)
binding_score <- preds[,3] + preds[,2] * 0.5

# Correlation of each feature with binding score
n_feats <- ncol(learned)
feature_cor <- sapply(1:n_feats, function(j) cor(learned[,j], binding_score))

cor_df <- tibble(
  Feature = 1:n_feats,
  Correlation = feature_cor
) %>% arrange(desc(abs(Correlation)))

p_cor <- ggplot(cor_df, aes(x = Feature, y = Correlation)) +
  geom_col(aes(fill = Correlation > 0), width = 0.7) +
  scale_fill_manual(values = c("TRUE" = "#B2182B", "FALSE" = "#2166AC"),
                    labels = c("TRUE" = "Positive", "FALSE" = "Negative"),
                    guide = "none") +
  labs(
    title = "Feature-Binding Score Correlation",
    subtitle = paste0("Top features: ", paste(head(cor_df$Feature, 5), collapse=", ")),
    x = "Feature Index", y = "Pearson r"
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))

ggsave("plots/feature_correlation.png", p_cor,
       width = 10, height = 4, dpi = 150)
cat("Saved: plots/feature_correlation.png\n")

# ---- 8. Plot 4: Clustered heatmap (skip if constant features break clustering) ----
cat("Generating clustered heatmap...\n")
tryCatch({
  dist_mat <- dist(learned_scaled)
  hc <- hclust(dist_mat, method = "ward.D2")
  pep_order <- hc$labels[hc$order]
  feat_df$Peptide_clust <- factor(feat_df$Peptide, levels = rev(pep_order))

  varying <- which(apply(learned_scaled, 2, sd) > 1e-6)
  feat_dist <- dist(t(learned_scaled[, varying, drop = FALSE]))
  feat_hc <- hclust(feat_dist, method = "ward.D2")
  feat_order <- colnames(learned_scaled[, varying])[feat_hc$order]
  if (is.null(feat_order) || length(feat_order) == 0) feat_order <- 1:n_feats
  feat_df$Feature_clust <- factor(feat_df$Feature, levels = feat_order)

  png("plots/feature_clustered_heatmap.png", width = 3600, height = 2400, res = 200)
  p_clust <- ggplot(feat_df, aes(x = Feature_clust, y = Peptide_clust, fill = Value)) +
    geom_tile() +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                         midpoint = 0, name = "Z-score") +
    labs(title = "Clustered Feature Heatmap -- Top Epitope Representations",
         subtitle = paste0("Rows: Ward.D2 clustering | Columns: ", n_feats, " features"),
         x = "Feature Index (clustered)", y = NULL) +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(face = "bold", size = 15),
          plot.subtitle = element_text(size = 10, color = "gray40"),
          axis.text.y = element_text(size = 9, family = "mono", face = "bold"),
          axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          panel.grid = element_blank(), legend.position = "bottom")
  print(p_clust)
  dev.off()
  cat("Saved: plots/feature_clustered_heatmap.png\n")
}, error = function(e) {
  cat(sprintf("Clustered heatmap skipped: %s\n", e$message))
})

cat("\nDone -- all feature maps generated.\n")
