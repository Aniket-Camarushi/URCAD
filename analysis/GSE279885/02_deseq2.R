# ==============================================================================
# analysis/GSE279885/02_deseq2.R
#
# Purpose: Differential expression analysis using DESeq2.
#
# WHY DESEQ2:
#   DESeq2 uses a negative-binomial generalized linear model to test for
#   differential expression. It includes its own dispersion shrinkage (which
#   stabilizes variance estimates for low-replicate studies) and uses the
#   Wald test for pairwise comparisons. It is the established standard for
#   count-based bulk RNA-seq DE.
#
# REQUIRED COHORT (runs automatically):
#   Design: ~ treatment  (VCD vs Vehicle_Control)
#   Contrast: VCD / Vehicle_Control
#
# OPTIONAL COHORT (uncomment the "OPTIONAL DATASET" blocks):
#   Design: ~ cohort + treatment
#   WHY add cohort to the design: the optional samples come from a
#   different experiment (shRNA knockdown). Including "cohort" as a
#   covariate removes batch-like variation between the two experimental
#   setups, allowing treatment effect to be estimated cleanly.
#   Note: this combined analysis is exploratory. The two cohorts test
#   different biological interventions (estrogen depletion vs Esrra KD)
#   and results should be interpreted accordingly.
#
# Outputs saved to results/GSE279885/tables/:
#   - dds_required.rds           DESeq2 dataset (required)
#   - vsd_required.rds           Variance-stabilised data (required, for PCA)
#   - de_results_required.csv    Full DE table (required)
#   - de_sig_required.csv        Significant DEGs only (required)
# ==============================================================================

suppressPackageStartupMessages({
  suppressWarnings({
    library(here)
    library(DESeq2)
    library(tximport)
    library(dplyr)
    library(readr)
    library(ggplot2)
    library(ggrepel)
    library(pheatmap)
    library(tibble)
    library(writexl)
  })
})

source(here("analysis", "helpers", "utils.R"))
source(here("analysis", "helpers", "themes.R"))

log_section("GSE279885 | Step 2 | DESeq2 Differential Expression Analysis")

dirs <- create_result_dirs("GSE279885")

# Load outputs from Step 1
txi      <- readRDS(here("results", "GSE279885", "tables", "txi_required.rds"))
metadata <- readRDS(here("results", "GSE279885", "tables", "metadata_required.rds"))
sym_map  <- readRDS(here("results", "GSE279885", "tables", "ensembl_to_symbol.rds"))

# ── [OPTIONAL DATASET — UNCOMMENT TO INCLUDE] ────────────────────────────────
# txi      <- readRDS(here("results", "GSE279885", "tables", "txi_combined.rds"))
# metadata <- readRDS(here("results", "GSE279885", "tables", "metadata_combined.rds"))
# ─────────────────────────────────────────────────────────────────────────────

# Align metadata row order to count matrix column order
metadata <- metadata[match(colnames(txi$counts), metadata$srr_id), ]
stopifnot(all(colnames(txi$counts) == metadata$srr_id))

# Factor levels: reference = control group (first level = baseline for DE)
metadata$treatment <- factor(metadata$treatment,
                             levels = c("Vehicle_Control", "VCD"))

# ── [OPTIONAL DATASET — UNCOMMENT TO INCLUDE] ────────────────────────────────
# When cohorts are combined, add a cohort covariate.
# metadata$cohort    <- factor(metadata$cohort, levels = c("required", "optional"))
# metadata$treatment <- factor(metadata$treatment,
#                              levels = c("Vehicle_Control", "VCD",
#                                         "Control_shRNA",  "Esrra_shRNA"))
# ─────────────────────────────────────────────────────────────────────────────

# ── 1. Build DESeqDataSet ─────────────────────────────────────────────────────
#
# DESeqDataSetFromTximport: directly accepts tximport output.
# design = ~ treatment: estimate the effect of treatment on expression.
#   For the combined dataset, use ~ cohort + treatment to partial out
#   between-cohort variance.

dds <- load_or_compute(
  path = here("results", "GSE279885", "tables", "dds_required.rds"),
  expr = {
    log_msg("Building DESeqDataSet ...")
    dds_obj <- DESeqDataSetFromTximport(
      txi      = txi,
      colData  = metadata,
      design   = ~ treatment
      # ── [OPTIONAL DATASET] use this design instead: ──────────────────────
      # design = ~ cohort + treatment
      # ─────────────────────────────────────────────────────────────────────
    )

    # Low-count filter: remove genes with fewer than 10 counts across all samples.
    # WHY: very low counts are dominated by noise; removing them reduces the
    # multiple testing burden and speeds up computation.
    keep <- rowSums(counts(dds_obj)) >= PARAMS$deseq2$min_count
    dds_obj <- dds_obj[keep, ]
    log_msg("Genes after low-count filter: ", nrow(dds_obj))

    # Run DESeq2
    DESeq(dds_obj)
  }
)

log_msg("DESeq2 complete. Genes tested: ", nrow(dds))

# ── 2. Variance-stabilised transformation (VSD) for visualisation ────────────
#
# WHY VSD not raw counts for PCA/heatmaps:
#   Raw counts have mean-variance dependence (highly expressed genes dominate).
#   VST removes this dependence, making distances between samples meaningful.
#   We use blind = FALSE so the dispersion estimates use the fitted model
#   (more accurate for small n).

vsd <- load_or_compute(
  path = here("results", "GSE279885", "tables", "vsd_required.rds"),
  expr = { varianceStabilizingTransformation(dds, blind = FALSE) }
)

# ── 3. Extract DE results ─────────────────────────────────────────────────────
#
# contrast = c("treatment", "VCD", "Vehicle_Control"):
#   Log2FC > 0 means higher in VCD (estrogen-depleted) vs Vehicle (control).
# lfcShrink with type = "apeglm": shrinks noisy log2FC estimates for low-count
#   genes toward zero. Essential for volcano plots and ranked gene lists.

res_raw <- results(
  dds,
  contrast = c("treatment", "VCD", "Vehicle_Control"),
  alpha    = PARAMS$deseq2$alpha
)

res_shrunk <- lfcShrink(
  dds,
  coef = "treatment_VCD_vs_Vehicle_Control",
  type = "apeglm"
)

# ── [OPTIONAL DATASET — UNCOMMENT TO INCLUDE] ────────────────────────────────
# For the combined design, extract each contrast separately.
# VCD vs Vehicle from the required cohort samples:
# res_vcd <- results(dds, contrast = c("treatment", "VCD", "Vehicle_Control"), alpha = 0.05)
#
# Esrra KD vs Control shRNA from the optional cohort samples:
# res_esrra <- results(dds, contrast = c("treatment", "Esrra_shRNA", "Control_shRNA"), alpha = 0.05)
# ─────────────────────────────────────────────────────────────────────────────

# Annotate results
de_results <- as.data.frame(res_shrunk) %>%
  annotate_deseq2_results(
    alpha = PARAMS$deseq2$alpha,
    lfc   = PARAMS$deseq2$lfc_threshold
  ) %>%
  tibble::rownames_to_column("ENSEMBL") %>%
  dplyr::left_join(sym_map, by = "ENSEMBL") %>%
  dplyr::mutate(gene = dplyr::coalesce(SYMBOL, ENSEMBL)) %>%
  dplyr::select(-SYMBOL)

de_sig <- de_results %>% dplyr::filter(significant)

log_msg("Total genes tested: ",   nrow(de_results))
log_msg("Significant DEGs: ",     nrow(de_sig))
log_msg("  Upregulated in VCD: ", sum(de_sig$direction == "up"))
log_msg("  Downregulated in VCD:", sum(de_sig$direction == "down"))

# Save tables
write_csv(de_results, here("results", "GSE279885", "tables", "de_results_required.csv"))
write_csv(de_sig,     here("results", "GSE279885", "tables", "de_sig_required.csv"))
log_msg("DE tables saved.")

# ── 4. PCA plot ───────────────────────────────────────────────────────────────
#
# PCA on the top 500 most variable genes (default for plotPCA).
# This is the first QC check after DE: samples should separate by treatment.

pca_data <- plotPCA(vsd, intgroup = c("treatment", "srr_id"), returnData = TRUE)
pct_var  <- round(100 * attr(pca_data, "percentVar"), 1)

p_pca <- ggplot(pca_data, aes(PC1, PC2,
                               colour = treatment, label = srr_id)) +
  geom_point(size = 3, alpha = 0.9) +
  ggrepel::geom_text_repel(size = 3, show.legend = FALSE) +
  scale_color_manual(values = PALETTE_VCD) +
  labs(
    title    = "PCA — GSE279885 (required cohort)",
    subtitle = "VSD-transformed counts, top 500 variable genes",
    x        = paste0("PC1 (", pct_var[1], "%)"),
    y        = paste0("PC2 (", pct_var[2], "%)"),
    colour   = "Treatment"
  )

save_figure(p_pca, "01_pca", dirs$figures, width = 7, height = 5)

# ── 5. Volcano plot ───────────────────────────────────────────────────────────
#
# Volcano plot: x = shrunk log2FC, y = -log10(adjusted p-value).
# Labels the top 15 up and top 15 down genes by adjusted p-value.

top_labels <- bind_rows(
  de_results %>% filter(direction == "up")   %>% arrange(padj) %>% head(15),
  de_results %>% filter(direction == "down") %>% arrange(padj) %>% head(15)
)

p_volcano <- ggplot(de_results,
                    aes(log2FoldChange, -log10(padj), colour = direction)) +
  geom_point(alpha = 0.5, size = 1.2) +
  ggrepel::geom_text_repel(data = top_labels,
                            aes(label = gene),
                            size = 2.8, max.overlaps = 20,
                            show.legend = FALSE) +
  geom_vline(xintercept = c(-PARAMS$deseq2$lfc_threshold,
                              PARAMS$deseq2$lfc_threshold),
             linetype = "dashed", colour = "grey50") +
  geom_hline(yintercept = -log10(PARAMS$deseq2$alpha),
             linetype = "dashed", colour = "grey50") +
  scale_color_manual(values = PALETTE_DE,
                     labels = c("up" = "Up in VCD",
                                "down" = "Down in VCD",
                                "ns" = "Not significant")) +
  labs(
    title    = "Volcano plot - VCD vs Vehicle (GSE279885)",
    subtitle = sprintf("DEGs (|LFC| > %.1f, FDR < %.2f): %d up, %d down",
                       PARAMS$deseq2$lfc_threshold, PARAMS$deseq2$alpha,
                       sum(de_sig$direction == "up"),
                       sum(de_sig$direction == "down")),
    x        = "Shrunk log2 fold change (VCD / Vehicle)",
    y        = expression(-log[10](adjusted~p-value)),
    colour   = NULL
  )

save_figure(p_volcano, "02_volcano", dirs$figures, width = 7, height = 6)

# ── 6. Heatmap of top 50 DEGs ────────────────────────────────────────────────

top50 <- de_sig %>%
  arrange(padj) %>%
  head(50)

# Extract VSD matrix for top 50, scale by row (z-score per gene)
mat <- assay(vsd)[top50$ENSEMBL, ]
rownames(mat) <- top50$gene
mat_scaled <- t(scale(t(mat)))

ann_col <- data.frame(
  Treatment = metadata$treatment,
  row.names = metadata$srr_id
)

ann_colors <- list(Treatment = PALETTE_VCD)

pdf(here("results", "GSE279885", "figures", "03_heatmap_top50.pdf"),
    width = 8, height = 10)

pheatmap(
  mat_scaled,
  annotation_col  = ann_col,
  annotation_colors = ann_colors,
  show_rownames   = TRUE,
  show_colnames   = TRUE,
  cluster_cols    = TRUE,
  cluster_rows    = TRUE,
  color           = colorRampPalette(c("#2166AC", "white", "#D73027"))(100),
  fontsize_row    = 7,
  fontsize_col    = 9,
  main            = "Top 50 DEGs - GSE279885 (VCD vs Vehicle)"
)

dev.off()

log_msg("Saved: 03_heatmap_top50.pdf")

log_section("GSE279885 | Step 2 | COMPLETE")
message("Next: run analysis/GSE279885/03_enrichment.R")
