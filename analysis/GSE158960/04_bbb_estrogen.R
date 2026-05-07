# ==============================================================================
# analysis/GSE158960/04_bbb_estrogen.R
#
# Purpose: Characterize estrogen receptor expression and BBB integrity markers
#          in the BBB-relevant cell types (endothelial, pericyte, astrocyte).
#
# SCIENTIFIC RATIONALE:
#   Estrogen (via ERα / Esr1) maintains BBB integrity by upregulating tight
#   junction proteins (Cldn5, Tjp1, Ocln) and suppressing inflammatory
#   mediators (Il1b, Tnf). In this dataset:
#     - WT vs Esr1cKO compares cells with intact vs ablated ERα signaling.
#     - Control vs Tamoxifen compares untreated vs tamoxifen-treated animals.
#   We test whether BBB gene expression differs across these four groups,
#   specifically in the endothelial and astrocyte populations.
#
# Approach:
#   1. Module scoring: summarize pathway activity per cell as a single score.
#   2. Differential expression: Wilcoxon rank-sum test via FindMarkers.
#   3. Violin / Feature plots: visualize score distributions and key genes.
# ==============================================================================

suppressPackageStartupMessages({
  suppressWarnings({
    library(here)
    library(Seurat)
    library(DESeq2)
    library(Matrix)
    library(dplyr)
    library(ggplot2)
    library(patchwork)
    library(writexl)
  })
})

source(here("analysis", "helpers", "utils.R"))
source(here("analysis", "helpers", "themes.R"))

log_section("GSE158960 | Step 4 | BBB & Estrogen Analysis")
dirs <- create_result_dirs("GSE158960")

seurat_bbb <- readRDS(
  here("results", "GSE158960", "tables", "seurat_bbb.rds")
)

# ── Gene sets for module scoring ──────────────────────────────────────────────

GENE_SETS <- list(
  # Estrogen receptors and downstream response genes
  estrogen_receptors = c("Esr1", "Esr2", "Gper1", "Esrra", "Esrrb", "Esrrg"),

  # Classical estrogen-response genes (ERE-containing promoters)
  estrogen_response = c("Pgr", "Prlr", "Ltf", "Muc1", "Tnfsf11",
                         "Igfbp3", "Tff1", "Tff3", "Cited2"),

  # BBB structural integrity (tight junctions, adherens junctions, transporters)
  bbb_integrity = c("Cldn5", "Tjp1", "Tjp2", "Ocln", "Cdh5",
                     "Esam", "Pecam1", "Abcb1a", "Slc2a1", "Mfsd2a"),

  # BBB-supportive factors from astrocyte endfeet
  astrocyte_bbb = c("Aqp4", "Gfap", "Aldh1l1", "Lrp1", "Apoe",
                     "Gja1", "S1pr3"),

  # Pericyte / mural cell markers and BBB regulators
  pericyte_bbb = c("Pdgfrb", "Rgs5", "Notch3", "Acta2",
                    "Angpt1", "Tek", "Des"),

  # Neuroinflammation (disrupts BBB in injury / estrogen deficiency)
  neuroinflammation = c("Il1b", "Il6", "Tnf", "Ccl2", "Cxcl10",
                         "Nos2", "Ptgs2", "Icam1", "Vcam1", "Sele")
)

# Intersect with genes present in the dataset
genes_in_data <- lapply(GENE_SETS, function(g) {
  intersect(g, rownames(seurat_bbb))
})

lapply(names(genes_in_data), function(n) {
  log_msg(sprintf("  %s: %d/%d genes present",
                  n, length(genes_in_data[[n]]), length(GENE_SETS[[n]])))
})

# ── 1. Module scoring ─────────────────────────────────────────────────────────
#
# AddModuleScore calculates a per-cell score for a gene set by averaging
# the expression of set genes and subtracting the average of randomly
# selected control genes with similar expression levels.
# WHY: reduces the effect of overall transcriptional activity on the score.

seurat_bbb <- AddModuleScore(seurat_bbb,
  features = genes_in_data,
  name     = names(genes_in_data),
  ctrl     = 20,
  seed     = 42
)

# Rename auto-generated columns (AddModuleScore appends numbers)
score_cols_old <- paste0(names(genes_in_data), seq_along(genes_in_data))
score_cols_new <- paste0("score_", names(genes_in_data))
for (i in seq_along(score_cols_old)) {
  seurat_bbb[[score_cols_new[i]]] <- seurat_bbb[[score_cols_old[i]]]
  seurat_bbb[[score_cols_old[i]]] <- NULL
}

# ── 2. Module score violin plots ─────────────────────────────────────────────

# score_plots <- lapply(score_cols_new, function(sc) {
#   VlnPlot(seurat_bbb,
#           features  = sc,
#           group.by  = "group",
#           pt.size   = 0,
#           cols      = PALETTE_GSE158960,
#           split.by  = "cell_type") +
#     theme_publication(base_size = 9) +
#     labs(title = gsub("score_", "", sc), x = NULL) +
#     theme(legend.position = "none",
#           axis.text.x = element_text(angle = 45, hjust = 1))
# })
# 
# p_scores <- patchwork::wrap_plots(score_plots, ncol = 2) +
#   patchwork::plot_annotation(
#     title    = "Module scores by group and cell type — GSE158960 BBB subset",
#     subtitle = "WT vs Esr1cKO × Control vs Tamoxifen",
#     theme    = theme_publication()
#   )

seurat_wt <- subset(seurat_bbb, subset = group %in% c("WT_Control", "WT_Tamoxifen"))

seurat_wt$group <- factor(
  seurat_wt$group,
  levels = c("WT_Control", "WT_Tamoxifen"),
  labels = c("Control",    "Tamoxifen")
)

# ── BBB integrity: endothelial cells, WT_Control vs WT_Tamoxifen ──
p_bbb_integrity <- VlnPlot(
  subset(seurat_wt, subset = cell_type == "Endothelial"),
  features = "score_bbb_integrity",
  group.by = "group",
  pt.size  = 0,
  cols     = c("#4393C3", "#D6604D")
) +
  theme_publication(base_size = 9) +
  labs(
    title = "BBB Integrity (Endothelial)",
    x     = NULL,
    y     = "Module score",
    fill  = "Group"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right")

p_bbb_integrity

save_figure(
  p_bbb_integrity,
  "07a_module_score_bbb_integrity_endothelial_WT",
  dirs$figures,
  width = 5, height = 4
)

# ── Astrocyte BBB: astrocytes, WT_Control vs WT_Tamoxifen ──
p_astro_bbb <- VlnPlot(
  subset(seurat_wt, subset = cell_type == "Astrocyte"),
  features = "score_astrocyte_bbb",
  group.by = "group",
  pt.size  = 0,
  cols     = c("#4393C3", "#D6604D")
) +
  theme_publication(base_size = 9) +
  labs(
    title = "Astrocyte BBB module (Astrocytes)",
    x     = NULL,
    y     = "Module score",
    fill  = "Group"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right")

p_astro_bbb

save_figure(
  p_astro_bbb,
  "07b_module_score_astrocyte_bbb_astro_WT",
  dirs$figures,
  width = 4, height = 4
)

# ── Neuroinflammation: astrocytes, WT_Control vs WT_Tamoxifen ──
p_neuroinf <- VlnPlot(
  subset(seurat_wt, subset = cell_type == "Astrocyte"),
  features = "score_neuroinflammation",
  group.by = "group",
  pt.size  = 0,
  cols     = c("#4393C3", "#D6604D")
) +
  theme_publication(base_size = 9) +
  labs(
    title = "Neuroinflammation (Astrocytes)",
    x     = NULL,
    y     = "Module score",
    fill  = "Group"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right")

p_neuroinf

save_figure(
  p_neuroinf,
  "07c_module_score_neuroinflammation_astro_WT",
  dirs$figures,
  width = 4, height = 4
)


# score_plots <- lapply(score_cols_new, function(sc) {
#   VlnPlot(seurat_bbb,
#           features  = sc,
#           group.by  = "group",
#           pt.size   = 0,
#           cols      = PALETTE_GSE158960,
#           split.by  = "cell_type") +
#     theme_publication(base_size = 9) +
#     labs(title = gsub("score_", "", sc), x = NULL) +
#     theme(legend.position = "none",
#           axis.text.x = element_text(angle = 45, hjust = 1))
# })
# 
# p_scores <- patchwork::wrap_plots(score_plots, ncol = 2) +
#   patchwork::plot_annotation(
#     title    = "Module scores by group and cell type - GSE158960 BBB subset",
#     subtitle = "Control vs Tamoxifen",
#     theme    = theme_publication()
#   )
# 
# save_figure(p_scores, "07_module_scores_bbb",
#             dirs$figures, width = 14, height = 16)

# ── 3. UMAP of module scores ──────────────────────────────────────────────────

umap_score_plots <- lapply(score_cols_new[1:4], function(sc) {
  FeaturePlot(seurat_bbb,
              features  = sc,
              reduction = "umap",
              cols      = c("grey90", "#D73027"),
              order     = TRUE) +
    theme_publication(base_size = 9) +
    labs(title = gsub("score_", "", sc))
})

p_umap_scores <- patchwork::wrap_plots(umap_score_plots, ncol = 2)

save_figure(p_umap_scores, "08_umap_module_scores",
            dirs$figures, width = 12, height = 10)

# ── 4. Differential expression in BBB cell types ─────────────────────────────
#
# FindMarkers with Wilcoxon rank-sum test (non-parametric; appropriate for
# count-like single-cell data without distributional assumptions).
#
# We compare:
#   a) WT vs Esr1cKO (genotype effect = proxy for ER-alpha signaling)
#   b) Control vs Tamoxifen (treatment effect), within each genotype

seurat_bbb[["RNA"]] <- CreateAssayObject(
  counts = GetAssayData(seurat_bbb, assay = "SCT", layer = "counts")
)

DefaultAssay(seurat_bbb) <- "RNA"

seurat_bbb <- SCTransform(
  seurat_bbb,
  vst.flavor      = "v2",
  vars.to.regress = "percent.mt",
  verbose         = FALSE
)

seurat_bbb <- PrepSCTFindMarkers(seurat_bbb, verbose = FALSE)

Idents(seurat_bbb) <- "group"

de_list <- list()
comparisons <- list(
  WT_Tamox_vs_Ctrl  = c("WT_Tamoxifen",       "WT_Control"),
  CKO_Tamox_vs_Ctrl = c("Esr1cKO_Tamoxifen", "Esr1cKO_Control"),
  WT_vs_CKO_ctrl    = c("WT_Control",          "Esr1cKO_Control")
)

for (comp_name in names(comparisons)) {
  grp <- comparisons[[comp_name]]
  log_msg("DE: ", comp_name)
  de <- tryCatch(
    FindMarkers(seurat_bbb,
                ident.1 = grp[1], 
                ident.2 = grp[2],
                logfc.threshold = 0.01,
                min.pct = 0.1,
                test.use = "wilcox",
                verbose = FALSE),
    error = function(e) {
      warning("FindMarkers failed for ", comp_name, ": ", e$message)
      NULL
    }
  )
  
  if (!is.null(de)) {
    de$gene      <- rownames(de)
    de$comparison <- comp_name
    de_list[[comp_name]] <- de
  }
}

# Save all DE results to Excel (one sheet per comparison)
if (length(de_list) > 0) {
  writexl::write_xlsx(de_list,
    here("results", "GSE158960", "tables", "de_bbb_comparisons.xlsx"))
  log_msg("Saved: de_bbb_comparisons.xlsx")
}

# ── 5. Key BBB gene violin plot ───────────────────────────────────────────────

# Show the top BBB tight-junction and estrogen-receptor genes
key_genes <- intersect(
  c("Cldn5", "Tjp1", "Ocln", "Esr1", "Gper1", "Aqp4", "Pdgfrb"),
  rownames(seurat_bbb)
)

p_key <- VlnPlot(seurat_bbb,
                  features = key_genes,
                  group.by = "group",
                  pt.size  = 0,
                  cols     = PALETTE_GSE158960,
                  ncol     = 4) &
  theme_publication(base_size = 9) &
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

save_figure(p_key, "09_key_genes_violin",
            dirs$figures, width = 14, height = 8)

saveRDS(seurat_bbb,
        here("results", "GSE158960", "tables", "seurat_bbb_scored.rds"))

# ── 5. Pseudobulk DESeq2 (avoids pseudoreplication) ──────────────────────────

# WHY PSEUDOBULK:
# FindMarkers treats every cell as an independent replicate, but cells from
# the same mouse share genetics, hormone levels, and handling effects.
# Pseudobulk collapses all cells from one sample × cell_type into a single
# count vector, giving ~4 real biological replicates per group.
# DESeq2 then tests those 4 vs 4 observations — statistically correct.
# Reference: Squair et al. (2021) Nature Communications.

pseudobulk_counts <- function(seurat_obj) {
  DefaultAssay(seurat_obj) <- "RNA"
  meta   <- seurat_obj@meta.data
  groups <- paste(meta$sample_id, meta$cell_type, sep = "__")
  mat    <- GetAssayData(seurat_obj, assay = "RNA", layer = "counts")
  
  agg <- sapply(unique(groups), function(grp) {
    idx <- which(groups == grp)
    if (length(idx) == 1) mat[, idx, drop = FALSE]
    else Matrix::rowSums(mat[, idx, drop = FALSE])
  })
  
  col_meta <- as.data.frame(do.call(rbind, strsplit(colnames(agg), "__")))
  colnames(col_meta) <- c("sample_id", "cell_type")
  col_meta$group    <- meta$group[match(col_meta$sample_id, meta$sample_id)]
  col_meta$genotype <- meta$genotype[match(col_meta$sample_id, meta$sample_id)]
  col_meta$treatment <- meta$treatment[match(col_meta$sample_id, meta$sample_id)]
  rownames(col_meta) <- colnames(agg)
  
  list(counts = as.matrix(agg), meta = col_meta)
}

pb <- pseudobulk_counts(seurat_bbb)
log_msg("Pseudobulk: ", ncol(pb$counts), " sample×cell_type pseudo-samples")

pb_de_list <- list()

for (ct in c("Endothelial", "Pericyte", "Astrocyte")) {
  idx <- pb$meta$cell_type == ct
  if (sum(idx) < 4) {
    log_msg("Skipping ", ct, ": fewer than 4 pseudobulk samples"); next
  }
  
  counts_ct <- pb$counts[, idx, drop = FALSE]
  meta_ct   <- pb$meta[idx, , drop = FALSE]
  meta_ct$genotype  <- factor(meta_ct$genotype,
                              levels = c("Wildtype", "Esr1cKO"))
  meta_ct$treatment <- factor(meta_ct$treatment,
                              levels = c("Control", "Tamoxifen"))
  
  # Low-count gene filter: >= 10 counts in at least 2 samples
  keep      <- rowSums(counts_ct >= PARAMS$deseq2$min_count) >= 2
  counts_ct <- counts_ct[keep, ]
  
  dds <- DESeqDataSetFromMatrix(
    countData = counts_ct,
    colData   = meta_ct,
    design    = ~ genotype + treatment   # main effects only (too few reps for interaction)
  )
  dds <- DESeq(dds, quiet = TRUE)
  
  # Genotype contrast: WT vs Esr1cKO (controlling for treatment)
  res_geno <- as.data.frame(results(dds,
                                    contrast = c("genotype", "Wildtype", "Esr1cKO"),
                                    alpha    = PARAMS$deseq2$alpha))
  res_geno <- annotate_deseq2_results(res_geno,
                                      alpha = PARAMS$deseq2$alpha, lfc = 0)
  res_geno$cell_type  <- ct
  res_geno$comparison <- "WT_vs_CKO_pseudobulk"
  
  # Treatment contrast: Tamoxifen vs Control (controlling for genotype)
  res_treat <- as.data.frame(results(dds,
                                     contrast = c("treatment", "Tamoxifen", "Control"),
                                     alpha    = PARAMS$deseq2$alpha))
  res_treat <- annotate_deseq2_results(res_treat,
                                       alpha = PARAMS$deseq2$alpha, lfc = 0)
  res_treat$cell_type  <- ct
  res_treat$comparison <- "Tamoxifen_vs_Control_pseudobulk"
  
  pb_de_list[[paste0(ct, "_geno")]]  <- res_geno
  pb_de_list[[paste0(ct, "_treat")]] <- res_treat
  
  log_msg(sprintf("  %s: %d sig genes (genotype), %d sig genes (treatment)",
                  ct,
                  sum(!is.na(res_geno$padj)  & res_geno$padj  < PARAMS$deseq2$alpha),
                  sum(!is.na(res_treat$padj) & res_treat$padj < PARAMS$deseq2$alpha)))
}

if (length(pb_de_list) > 0) {
  writexl::write_xlsx(pb_de_list,
                      here("results", "GSE158960", "tables", "pseudobulk_deseq2_bbb.xlsx"))
  log_msg("Saved: pseudobulk_deseq2_bbb.xlsx")
}


log_section("GSE158960 | Step 4 | COMPLETE")
message("Next: run analysis/GSE158960/05_cellchat.R")
