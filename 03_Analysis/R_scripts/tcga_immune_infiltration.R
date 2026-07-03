#!/usr/bin/env Rscript
#
# TCGA Immune Infiltration Analysis for p53 R248W and KRAS G12V
# =============================================================
# Compares MHC-I (HLA-A/B/C) and CD8 T cell marker (CD8A/B) expression
# between mutant and wild-type tumors using cBioPortal REST API.
#
# Target tumor types:
#   p53 R248W  -> SKCM (skin melanoma), COADREAD (colorectal), ESCA (esophageal)
#   KRAS G12V  -> PAAD (pancreatic), LUAD (lung adeno), COADREAD (colorectal)
#
# Output: 02_Data/cleaned/tcga_immune_infiltration.csv
#         02_Data/cleaned/tcga_immune_test_results.csv
#         03_Analysis/figures/tcga_immune_*.png

options(warn = -1)
suppressMessages({
  library(httr)
  library(jsonlite)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggpubr)
})

CBIOPORTAL <- "https://www.cbioportal.org/api"
OUTPUT_DIR  <- normalizePath("D:/Researching/Peptide epitope")
DATA_DIR    <- file.path(OUTPUT_DIR, "02_Data", "cleaned")
FIGURE_DIR  <- file.path(OUTPUT_DIR, "03_Analysis", "figures")
dir.create(DATA_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(FIGURE_DIR, showWarnings = FALSE, recursive = TRUE)

`%||%` <- function(a, b) if (is.null(a)) b else a

# ---- Helper: GET with retry ----
cbioget <- function(url, max_retries = 3, timeout_sec = 30) {
  for (i in seq_len(max_retries)) {
    r <- tryCatch(GET(URLencode(url), timeout(timeout_sec * 1000)), error = function(e) NULL)
    if (!is.null(r) && status_code(r) == 200) return(content(r, as = "parsed"))
    if (i < max_retries) Sys.sleep(2)
  }
  NULL
}

# ---- Gene entrez IDs ----
GENE_ENTREZ <- list(
  "HLA-A" = 3105, "HLA-B" = 3106, "HLA-C" = 3107,
  "B2M"   = 567,  "CD8A"  = 925,  "CD8B"  = 926,
  "TAP1"  = 6890, "TAP2"  = 6891,
  "TP53"  = 7157, "KRAS"  = 3845
)

# ==== Query mutations ====
query_mutations <- function(study_id, entrez_id) {
  for (sl_suffix in c("_sequenced", "_all")) {
    sl <- paste0(study_id, sl_suffix)
    url <- paste0(CBIOPORTAL, "/molecular-profiles/", study_id, "_mutations",
                  "/mutations?sampleListId=", sl,
                  "&entrezGeneId=", entrez_id, "&pageSize=10000")
    muts <- cbioget(url)
    if (!is.null(muts) && length(muts) > 0) {
      return(list(
        samples = do.call(rbind, lapply(muts, function(m) {
          data.frame(sampleId = m$sampleId, proteinChange = m$proteinChange %||% "", stringsAsFactors = FALSE)
        })),
        sample_list = sl
      ))
    }
  }
  NULL
}

# ==== Find best expression profile for a study ====
find_best_profile <- function(study_id, mut_samples) {
  # Try RNA-seq first, then microarray; pick the one covering most mutant samples
  profiles_to_try <- list(
    list(id = paste0(study_id, "_rna_seq_v2_mrna"), sl = paste0(study_id, "_rna_seq_v2_mrna"), label = "RNA-seq"),
    list(id = paste0(study_id, "_mrna"),             sl = paste0(study_id, "_mrna"),             label = "microarray")
  )

  best <- NULL
  best_covered <- -1

  for (prof in profiles_to_try) {
    # Check profile exists
    r <- cbioget(paste0(CBIOPORTAL, "/molecular-profiles/", prof$id))
    if (is.null(r)) next

    # Get a single gene's data to identify which samples are present
    url <- paste0(CBIOPORTAL, "/molecular-profiles/", prof$id,
                  "/molecular-data?sampleListId=", prof$sl,
                  "&entrezGeneId=567&pageSize=10000")  # B2M as reference
    data <- cbioget(url)
    if (is.null(data) || length(data) == 0) next

    avail_samples <- unique(sapply(data, function(d) d$sampleId))
    covered <- sum(mut_samples %in% avail_samples)

    if (covered > best_covered || (covered == best_covered && prof$label == "RNA-seq")) {
      best <- list(
        id = prof$id,
        sl = prof$sl,
        avail_samples = avail_samples,
        label = prof$label,
        n_covered = covered,
        n_total = length(avail_samples)
      )
      best_covered <- covered
    }
  }
  best
}

# ==== Query expression via the best profile ====
query_expression <- function(profile_id, sample_list_id, entrez_ids, gene_names) {
  all_expr <- list()
  for (i in seq_along(entrez_ids)) {
    eid <- entrez_ids[i]; gnm <- gene_names[i]
    url <- paste0(CBIOPORTAL, "/molecular-profiles/", profile_id,
                  "/molecular-data?sampleListId=", sample_list_id,
                  "&entrezGeneId=", eid, "&pageSize=10000")
    data <- cbioget(url)
    if (is.null(data) || length(data) == 0) next

    df <- do.call(rbind, lapply(data, function(d) {
      data.frame(sampleId = d$sampleId, expression = as.numeric(d$value %||% d$expressionValue %||% NA),
                 stringsAsFactors = FALSE)
    }))
    if (nrow(df) > 0) df$gene <- gnm
    all_expr[[gnm]] <- df
  }
  if (length(all_expr) == 0) return(NULL)
  bind_rows(all_expr)
}

# ==== Sample type ====
sample_type_code <- function(sample_id) {
  code <- strsplit(sample_id, "-")[[1]]
  code <- code[length(code)]
  if (code == "01") return("Primary Solid Tumor")
  if (code %in% c("06", "07")) return("Metastatic")
  if (code == "11") return("Solid Tissue Normal")
  if (code == "02") return("Recurrent Solid Tumor")
  return("Other")
}

# ==== Run ====
run_analysis <- function() {
  cat("==============================================================\n")
  cat("  TCGA IMMUNE INFILTRATION ANALYSIS\n")
  cat("  p53 R248W / KRAS G12V -- MHC-I & CD8 expression\n")
  cat("==============================================================\n\n")

  scenario <- list(
    list(label = "p53.R248W", gene = "TP53", entrez = 7157, mut_pattern = "R248W",
         studies = c("skcm_tcga", "coadread_tcga", "esca_tcga"),
         cancer_names = c("SKCM (Melanoma)", "COADREAD (Colorectal)", "ESCA (Esophageal)")),
    list(label = "KRAS.G12V", gene = "KRAS", entrez = 3845, mut_pattern = "G12V",
         studies = c("paad_tcga", "luad_tcga", "coadread_tcga"),
         cancer_names = c("PAAD (Pancreatic)", "LUAD (Lung Adeno)", "COADREAD (Colorectal)"))
  )

  immune_genes <- c("HLA-A", "HLA-B", "HLA-C", "B2M", "CD8A", "CD8B", "TAP1", "TAP2")
  immune_entrez <- unname(unlist(GENE_ENTREZ[immune_genes]))
  all_data <- list()

  for (sc in scenario) {
    cat(sprintf("\n=== %s (%s) ===\n", sc$label, sc$gene))
    for (i in seq_along(sc$studies)) {
      sid <- sc$studies[i]; cnm <- sc$cancer_names[i]
      cat(sprintf("\n--- %s (%s) ---\n", cnm, sid))

      # 1. Mutations
      mut <- query_mutations(sid, sc$entrez)
      if (is.null(mut)) { cat("    No mutation data\n"); next }
      cat(sprintf("    Total %s mutations: %d\n", sc$gene, nrow(mut$samples)))

      carriers <- unique(mut$samples$sampleId[grepl(sc$mut_pattern, mut$samples$proteinChange)])
      cat(sprintf("    %s carriers: %d\n", sc$mut_pattern, length(carriers)))

      # 2. Find best expression profile
      prof <- find_best_profile(sid, carriers)
      if (is.null(prof)) { cat("    No expression profile\n"); next }
      cat(sprintf("    Best profile: %s (%s) -- %d samples, %d/%d mutants covered\n",
                  prof$id, prof$label, prof$n_total, prof$n_covered, length(carriers)))

      if (prof$n_covered < 2) {
        cat("    SKIP: < 2 mutant samples with expression data\n")
        next
      }

      # 3. Expression
      expr <- query_expression(prof$id, prof$sl, immune_entrez, immune_genes)
      if (is.null(expr)) { cat("    No expression data returned\n"); next }

      expr <- expr %>%
        mutate(
          cancer_type = cnm, study_id = sid,
          mutant = sampleId %in% carriers,
          mutation_status = if_else(mutant, sc$mut_pattern, "Wild-type"),
          mutation_label = sc$label,
          sample_type = sapply(sampleId, sample_type_code)
        )
      cat(sprintf("    Mutants: %d | WT: %d | Rows: %d\n",
                  n_distinct(expr$sampleId[expr$mutant]),
                  n_distinct(expr$sampleId[!expr$mutant]),
                  nrow(expr)))

      all_data[[paste(sc$label, sid, sep = "|")]] <- expr
    }
  }

  if (length(all_data) == 0) { cat("\nERROR: No data collected.\n"); return(invisible(NULL)) }

  comb <- bind_rows(all_data, .id = "batch")
  comb$mutation_group <- sub("\\|.*", "", comb$batch)

  # ---- Save ----
  fout <- file.path(DATA_DIR, "tcga_immune_infiltration.csv")
  write.csv(comb, fout, row.names = FALSE)
  cat(sprintf("\nRaw saved: %s (%d rows)\n", fout, nrow(comb)))

  # ---- Stats ----
  stats <- comb %>%
    group_by(mutation_group, cancer_type, gene, mutant) %>%
    summarise(n = n(), mean = mean(expression, na.rm = TRUE),
              median = median(expression, na.rm = TRUE), .groups = "drop")

  # ---- Wilcoxon ----
  tests <- comb %>%
    group_by(mutation_group, cancer_type, gene) %>%
    summarise(
      n_mut = sum(mutant), n_wt = sum(!mutant),
      mean_mut = mean(expression[mutant], na.rm = TRUE),
      mean_wt  = mean(expression[!mutant], na.rm = TRUE),
      p_value = tryCatch(wilcox.test(expression[mutant], expression[!mutant])$p.value,
                         error = function(e) NA_real_),
      .groups = "drop"
    ) %>%
    mutate(log2FC = log2(mean_mut / pmax(mean_wt, 1e-10)),
           p_adj = p.adjust(p_value, method = "BH"))

  tests_show <- tests %>% filter(!is.na(p_value)) %>% arrange(p_adj)
  fstat <- file.path(DATA_DIR, "tcga_immune_test_results.csv")
  write.csv(tests_show, fstat, row.names = FALSE)

  # ---- Print summary ----
  cat("\n\n=== KEY FINDINGS ===\n")
  for (mg in unique(tests_show$mutation_group)) {
    cat(sprintf("\n--- %s ---\n", mg))
    sub_t <- tests_show %>% filter(mutation_group == mg)
    sig <- sub_t %>% filter(p_adj < 0.05)
    if (nrow(sig) > 0) {
      for (j in 1:nrow(sig))
        cat(sprintf("  ** %s in %s: log2FC=%.3f, p_adj=%.4f (%s in mutant)\n",
                    sig$gene[j], sig$cancer_type[j], sig$log2FC[j], sig$p_adj[j],
                    ifelse(sig$log2FC[j] > 0, "up", "down")))
    } else {
      cat("  No significant differences after BH correction.\n")
    }
    for (j in 1:nrow(sub_t)) {
      if (is.na(sub_t$p_value[j]))
        cat(sprintf("    %s in %s: insufficient data\n", sub_t$gene[j], sub_t$cancer_type[j]))
      else
        cat(sprintf("    %s in %s: mut=%.0f wt=%.0f log2FC=%+.3f p=%.4f\n",
                    sub_t$gene[j], sub_t$cancer_type[j],
                    sub_t$mean_mut[j], sub_t$mean_wt[j],
                    sub_t$log2FC[j], sub_t$p_value[j]))
    }
  }

  # ---- Figures ----
  cat("\n\n--- Figures ---\n")

  for (mg in unique(comb$mutation_group)) {
    sub <- comb %>% filter(mutation_group == mg)
    p <- ggplot(sub, aes(x = cancer_type, y = log2(expression + 1), fill = mutation_status)) +
      geom_boxplot(outlier.size = 0.2, alpha = 0.7, coef = 1.5) +
      facet_wrap(~ gene, scales = "free_y", ncol = 4) +
      scale_fill_manual(values = c("Wild-type" = "grey70", "R248W" = "#E64B35", "G12V" = "#4DBBD5")) +
      labs(title = paste("TCGA:", mg, "- Immune gene expression"), x = NULL, y = "log2(Exp+1)") +
      theme_bw(base_size = 9) + theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom")
    fn <- file.path(FIGURE_DIR, paste0("tcga_immune_", gsub("\\.", "_", mg), ".png"))
    ggsave(fn, p, width = 10, height = 7, dpi = 150)
    cat(sprintf("  %s\n", fn))
  }

  fc_df <- tests_show %>% filter(is.finite(log2FC))
  if (nrow(fc_df) > 0) {
    p2 <- ggplot(fc_df, aes(x = cancer_type, y = gene, fill = log2FC)) +
      geom_tile(color = "white") +
      facet_wrap(~ mutation_group, ncol = 2, scales = "free_x") +
      scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                           midpoint = 0, name = "log2FC") +
      labs(title = "log2FC (Mutant vs Wild-type)") + theme_bw(base_size = 10) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    fn2 <- file.path(FIGURE_DIR, "tcga_immune_heatmap_log2FC.png")
    ggsave(fn2, p2, width = 8, height = 6, dpi = 150)
  }

  cd8 <- comb %>% filter(gene %in% c("CD8A", "CD8B")) %>%
    group_by(mutation_group, cancer_type, sampleId, mutant, mutation_status) %>%
    summarise(cd8 = mean(expression, na.rm = TRUE), .groups = "drop") %>% filter(is.finite(cd8))
  if (nrow(cd8) > 0) {
    p3 <- ggplot(cd8, aes(x = cancer_type, y = log2(cd8 + 1), fill = mutation_status)) +
      geom_boxplot(outlier.size = 0.2, alpha = 0.7) +
      facet_wrap(~ mutation_group, ncol = 2, scales = "free_x") +
      scale_fill_manual(values = c("Wild-type" = "grey70", "R248W" = "#E64B35", "G12V" = "#4DBBD5")) +
      labs(title = "CD8 expression (proxy for CD8+ T cell infiltration)", x = NULL, y = "log2(mean CD8 + 1)") +
      theme_bw(base_size = 10) + theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom")
    fn3 <- file.path(FIGURE_DIR, "tcga_immune_cd8_infiltration.png")
    ggsave(fn3, p3, width = 9, height = 5.5, dpi = 150)
    cat(sprintf("  %s\n", fn3))
  }

  mhc <- comb %>% filter(gene %in% c("HLA-A", "HLA-B", "HLA-C")) %>%
    group_by(mutation_group, cancer_type, sampleId, mutant, mutation_status) %>%
    summarise(mhc = mean(expression, na.rm = TRUE), .groups = "drop") %>% filter(is.finite(mhc))
  if (nrow(mhc) > 0) {
    p4 <- ggplot(mhc, aes(x = cancer_type, y = log2(mhc + 1), fill = mutation_status)) +
      geom_boxplot(outlier.size = 0.2, alpha = 0.7) +
      facet_wrap(~ mutation_group, ncol = 2, scales = "free_x") +
      scale_fill_manual(values = c("Wild-type" = "grey70", "R248W" = "#E64B35", "G12V" = "#4DBBD5")) +
      labs(title = "MHC-I (HLA-A/B/C mean) expression (antigen presentation proxy)", x = NULL, y = "log2(mean MHC-I + 1)") +
      theme_bw(base_size = 10) + theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom")
    fn4 <- file.path(FIGURE_DIR, "tcga_immune_mhc1_expression.png")
    ggsave(fn4, p4, width = 9, height = 5.5, dpi = 150)
    cat(sprintf("  %s\n", fn4))
  }

  cat("\n==============================================================\n")
  cat("  ANALYSIS COMPLETE\n")
  cat("==============================================================\n")
}

run_analysis()
