# ==============================================================================
# analysis/integration/01_build_master.R
#
# Purpose: Combine all four datasets into two integration layers:
#   Layer A — Single-cell/nucleus object: scRNA (GSE158960) + snRNA (GSE299005)
#             integrated via Harmony in a single Seurat v5 multi-layer object.
#   Layer B — Results master table: DE results, pathway enrichment, and
#             module scores from all four datasets in one tidy data frame.
#
# WHY TWO LAYERS:
#   The sc/sn object is a proper molecular integration — cells and nuclei
#   from both datasets live in a shared embedding space, enabling cross-
#   dataset cell type comparison and combined ligand-receptor analysis.
#   The results table is a metadata-level integration — it allows you to
#   ask "which genes are differentially expressed in ALL datasets?" without
#   forcing incompatible assay types (RNA-seq counts vs NanoString counts)
#   into the same matrix.
#
# WHY HARMONY AGAIN (not RPCA or CCA):
#   Harmony is best when the biological variation across datasets is larger
#   than the technical variation. scRNA (hypothalamus, estrogen model) vs
#   snRNA (brain, FPI model) differ in cell type composition AND biological
#   context; Harmony's iterative correction handles this better than Seurat
#   CCA, which can over-integrate when groups have genuinely different biology.
#
# Outputs:
#   - seurat_integrated.rds      Seurat v5 sc+sn combined object
#   - results_master_table.csv   Cross-dataset DE / pathway summary
# ==============================================================================

suppressPackageStartupMessages({
  suppressWarnings({
    library(here)
    library(Seurat)
    library(harmony)
    library(dplyr)
    library(readr)
    library(ggplot2)
    library(patchwork)
    library(tibble)
    library(ggrepel)
  })
})

source(here("analysis", "helpers", "utils.R"))
source(here("analysis", "helpers", "themes.R"))

log_section("Integration | Step 1 | Build Master Object")
dirs <- create_result_dirs("integration")

# ── Layer A: sc + sn Seurat integration ──────────────────────────────────────

seurat_sc <- readRDS(
  here("results", "GSE158960", "tables", "seurat_annotated.rds")
)
seurat_sn <- readRDS(
  here("results", "GSE299005", "tables", "seurat_sn_annotated.rds")
)

# Before merging, stamp each object clearly with its source so cells are
# always traceable back to their original dataset and assay type.
seurat_sc$dataset    <- "GSE158960"
seurat_sc$data_type  <- "scRNA"
seurat_sc$model      <- "estrogen_model"   # WT/Esr1cKO tamoxifen

seurat_sn$dataset    <- "GSE299005"
seurat_sn$data_type  <- "snRNA"
seurat_sn$model      <- "FPI"              # fluid percussion injury

# Harmonise cell type labels across both objects so they map to the same
# canonical names and PALETTE_CELL_TYPES colours.
harmonise_cell_types <- function(obj) {
  obj$cell_type <- dplyr::recode(
    obj$cell_type,
    "Endothelial cells" = "Endothelial",
    "Pericytes"         = "Pericyte",
    "Astrocytes"        = "Astrocyte",
    "Oligodendrocytes"  = "Oligodendrocyte",
    "Microglia"         = "Microglia",
    "Neurons"           = "Neuron"
  )
  obj$cell_type[is.na(obj$cell_type)] <- "Other"
  obj
}

seurat_sc <- harmonise_cell_types(seurat_sc)
seurat_sn <- harmonise_cell_types(seurat_sn)

# ── Merge: Seurat v5 with per-sample layers ───────────────────────────────────
#
# In Seurat v5, merge() creates a multi-layer object where each sample's
# counts live in its own layer (counts.WT_Control_s01, counts.GSM9030357, ...).
# This is more memory-efficient than earlier versions and preserves the
# original counts for pseudo-bulk analysis.

seurat_integrated <- load_or_compute(
  path = here("results", "integration", "seurat_integrated.rds"),
  expr = {
    log_msg("Merging sc + sn objects ...")
    merged <- merge(
      seurat_sc,
      y           = seurat_sn,
      add.cell.ids = c("sc", "sn"),
      project      = "URCAD_integrated",
      merge.data   = TRUE
    )
    log_msg("Total cells after merge: ", ncol(merged))

    # Normalize merged object
    log_msg("SCTransform on merged object ...")
    merged <- SCTransform(merged,
                          vst.flavor      = "v2",
                          vars.to.regress = "percent.mt",
                          verbose         = FALSE)

    # PCA
    merged <- RunPCA(merged, npcs = PARAMS$seurat$n_pcs, verbose = FALSE)

    # Harmony: correct for dataset (scRNA vs snRNA) AND sample_id
    # WHY two group variables: dataset captures the modality-level effect
    # (sc vs sn); sample_id captures individual sample variation within each.
    log_msg("Harmony integration (dataset + sample_id) ...")
    merged <- RunHarmony(
      merged,
      group.by.vars  = c("dataset", "sample_id"),
      reduction      = "pca",
      reduction.save = "harmony",
      verbose        = FALSE
    )

    # UMAP + clustering
    merged <- RunUMAP(merged, reduction = "harmony",
                      dims = 1:PARAMS$seurat$n_pcs, verbose = FALSE)
    merged <- FindNeighbors(merged, reduction = "harmony",
                            dims = 1:PARAMS$seurat$n_pcs, verbose = FALSE)
    merged <- FindClusters(merged, resolution = PARAMS$seurat$resolution,
                           verbose = FALSE)
    merged
  }
)

log_msg("Integrated object cells: ", ncol(seurat_integrated))
log_msg("Datasets: ", paste(unique(seurat_integrated$dataset), collapse = ", "))

# ── UMAP figures ──────────────────────────────────────────────────────────────

p_int_dataset <- DimPlot(seurat_integrated, reduction = "umap",
                          group.by = "dataset",
                          cols = c("GSE158960" = "#4393C3",
                                   "GSE299005" = "#D6604D"),
                          alpha = 0.5) +
  theme_publication() +
  labs(title = "Integrated UMAP — dataset of origin")

p_int_celltype <- DimPlot(seurat_integrated, reduction = "umap",
                           group.by = "cell_type",
                           label = TRUE, repel = TRUE,
                           label.size = 3,
                           cols = PALETTE_CELL_TYPES) +
  theme_publication() +
  labs(title = "Integrated UMAP — cell type")

p_int_datatype <- DimPlot(seurat_integrated, reduction = "umap",
                           group.by = "data_type",
                           cols = c("scRNA" = "#4393C3", "snRNA" = "#D6604D"),
                           alpha = 0.5) +
  theme_publication() +
  labs(title = "Integrated UMAP — assay type (sc vs sn)")

p_int_combined <- (p_int_dataset | p_int_celltype) / p_int_datatype
save_figure(p_int_combined, "01_integrated_umap",
            dirs$figures, width = 14, height = 12)

# ── Layer B: Results master table ─────────────────────────────────────────────
#
# Pull the DE results from each dataset into a unified tidy data frame.
# Columns that all datasets share: gene, log2FC (or equivalent), adj_p, direction.
# Dataset-specific columns are preserved with a prefix.

build_master_table <- function() {

  # GSE279885 bulk DE
  de_bulk <- readr::read_csv(
    here("results", "GSE279885", "tables", "de_results_required.csv"),
    show_col_types = FALSE
  ) %>%
    dplyr::select(gene,
                  bulk_lfc   = log2FoldChange,
                  bulk_padj  = padj,
                  bulk_dir   = direction)

  # GSE160651 NanoString DE
  de_nano <- readr::read_csv(
    here("results", "GSE160651", "tables", "de_deseq2_results_tmp.csv"),
    show_col_types = FALSE
  ) %>%
    dplyr::select(gene,
                  nano_lfc   = logFC,
                  nano_padj  = adj.P.Val,
                  nano_dir   = direction)

  # Join on gene symbol
  master <- dplyr::full_join(de_bulk, de_nano, by = "gene")

  # Count in how many datasets each gene is significant
  master$n_datasets_sig <- rowSums(
    cbind(
      bulk = !is.na(master$bulk_padj) & master$bulk_padj < 0.05,
      nano = !is.na(master$nano_padj) & master$nano_padj < 0.05
    ),
    na.rm = TRUE
  )

  # Flag genes significant in BOTH datasets (high-confidence cross-dataset DEGs)
  master$cross_dataset_sig <- master$n_datasets_sig >= 2

  master %>% dplyr::arrange(desc(n_datasets_sig), bulk_padj)
}

master_table <- build_master_table()
readr::write_csv(master_table,
                 here("results", "integration", "tables", "results_master_table.csv"))

log_msg("Master table rows (unique genes): ", nrow(master_table))
log_msg("Cross-dataset significant genes: ",
        sum(master_table$cross_dataset_sig, na.rm = TRUE))

# Show top cross-dataset DEGs
log_msg("Top cross-dataset DEGs:")
master_table %>%
  dplyr::filter(cross_dataset_sig) %>%
  dplyr::arrange(bulk_padj) %>%
  dplyr::select(gene, bulk_lfc, bulk_padj, nano_lfc, nano_padj) %>%
  head(20) %>%
  print()

log_section("Integration | Step 1 | COMPLETE")
message("Next: run analysis/integration/02_cross_dataset_pathways.R")


PADJ_CUTOFF <- 0.05
LFC_CUTOFF  <- 0.5

# ── 1. Load GSE279885 bulk DE ──────────────────────────────────────────────
de_bulk <- read_csv(
  here("results", "GSE279885", "tables", "de_results_required.csv"),
  show_col_types = FALSE
) %>%
  select(gene,
         bulk_lfc  = log2FoldChange,
         bulk_padj = padj,
         bulk_dir  = direction) %>%
  mutate(sig_bulk = !is.na(bulk_padj) & bulk_padj < PADJ_CUTOFF)

de_bulk_unique <- de_bulk %>%
  group_by(gene) %>%
  arrange(bulk_padj, .by_group = TRUE) %>%
  slice(1) %>%
  ungroup()

log_msg("Bulk DEGs loaded: ", nrow(de_bulk),
        " | significant: ", sum(de_bulk$sig_bulk, na.rm = TRUE))

# ── 2. Load GSE160651 NanoString/DESeq2 — LONG format (3 contrasts) ───────
de_nano_long <- read_csv(
  here("results", "GSE160651", "tables", "de_deseq2_results_tmp.csv"),
  show_col_types = FALSE
) %>%
  select(gene, contrast,
         nano_lfc  = log2FoldChange,
         nano_padj = padj,
         nano_dir  = direction) %>%
  mutate(sig_nano = !is.na(nano_padj) &
           nano_padj < PADJ_CUTOFF &
           abs(nano_lfc) > LFC_CUTOFF)

for (cname in unique(de_nano_long$contrast)) {
  n <- sum(de_nano_long$sig_nano[de_nano_long$contrast == cname], na.rm = TRUE)
  log_msg("  Nano sig [", cname, "]: ", n)
}

# ── 3. Pivot NanoString to wide: one row per gene ──────────────────────────
de_nano_wide <- de_nano_long %>%
  tidyr::pivot_wider(
    id_cols    = gene,
    names_from = contrast,
    values_from = c(nano_lfc, nano_padj, nano_dir, sig_nano),
    names_sep  = "_"
  )

lfc_cols  <- grep("^nano_lfc_",  names(de_nano_wide), value = TRUE)
padj_cols <- grep("^nano_padj_", names(de_nano_wide), value = TRUE)

de_nano_wide <- de_nano_wide %>%
  mutate(across(all_of(lfc_cols),  as.numeric),
         across(all_of(padj_cols), as.numeric))

# ── 4. Full join (now 1-to-1) ─────────────────────────────────────────────
# master <- full_join(de_bulk, de_nano_wide, by = "gene")
master <- full_join(de_bulk_unique, de_nano_wide, by = "gene")

# ── 5. Temporal persistence flags ─────────────────────────────────────────
contrasts <- c("TBI_1d_vs_Control", "TBI_7d_vs_Control", "TBI_30d_vs_Control")

# for (ct in contrasts) {
#   master[[paste0("sig_", ct)]] <- {
#     col_p <- paste0("nano_padj_", ct)
#     col_l <- paste0("nano_lfc_", ct)
#     !is.na(master[[col_p]]) &
#       master[[col_p]] < PADJ_CUTOFF &
#       abs(master[[col_l]]) > LFC_CUTOFF
#   }
# }

for (ct in contrasts) {
  col_p <- paste0("nano_padj_", ct)
  col_l <- paste0("nano_lfc_",  ct)
  
  if (!col_p %in% names(master) || !col_l %in% names(master)) {
    warning("Missing columns for contrast ", ct, ": ",
            col_p, ", ", col_l)
    master[[paste0("sig_", ct)]] <- FALSE
  } else {
    master[[paste0("sig_", ct)]] <-
      !is.na(master[[col_p]]) &
      master[[col_p]] < PADJ_CUTOFF &
      abs(master[[col_l]]) > LFC_CUTOFF
  }
}

# master$sig_any_tbi    <- Reduce(`|`, lapply(contrasts, function(ct) master[[paste0("sig_", ct)]]))
# master$n_tbi_contrasts <- Reduce(`+`, lapply(contrasts, function(ct) as.integer(master[[paste0("sig_", ct)]])))

master$sig_any_tbi      <- Reduce(`|`, lapply(contrasts, function(ct) master[[paste0("sig_", ct)]]))
master$n_tbi_contrasts  <- Reduce(`+`, lapply(contrasts, function(ct) as.integer(master[[paste0("sig_", ct)]])))
master$cross_dataset_sig <- master$sig_bulk & master$sig_any_tbi

log_msg("Cross-dataset significant genes: ",
        sum(master$cross_dataset_sig, na.rm = TRUE))

# ── 6. Cross-dataset flag: sig in bulk VCD AND any TBI timepoint ──────────
master$cross_dataset_sig <- master$sig_bulk & master$sig_any_tbi

log_msg("Total master table rows: ", nrow(master))
log_msg("Cross-dataset significant genes: ", sum(master$cross_dataset_sig, na.rm = TRUE))

# ── 7. Save master table ──────────────────────────────────────────────────
write_csv(master, here("results", "integration", "tables", "results_master_table.csv"))

cross <- master %>%
  filter(cross_dataset_sig) %>%
  arrange(desc(n_tbi_contrasts), bulk_padj)

write_csv(cross, here("results", "integration", "tables", "cross_dataset_sig_genes.csv"))

# Excel with two tabs
write_xlsx(
  list(master_table = master, cross_dataset_hits = cross),
  here("results", "integration", "tables", "integration_results.xlsx")
)



plot_df <- master %>%
  filter(!is.na(bulk_lfc), !is.na(nano_lfc_TBI_7d_vs_Control)) %>%
  mutate(
    direction_match = sign(bulk_lfc) == sign(nano_lfc_TBI_7d_vs_Control),
    category = case_when(
      sig_bulk & direction_match ~ "Bulk sig + same direction",
      sig_bulk & !direction_match ~ "Bulk sig + opposite",
      sig_any_tbi ~ "Nano sig only",
      TRUE ~ "Background"
    )
  )

plot_df

p_scatter <- ggplot(plot_df, aes(bulk_lfc, nano_lfc_TBI_7d_vs_Control)) +
  geom_point(aes(color = category), alpha = 0.6, size = 1.8) +
  geom_smooth(method = "lm", color = "black", linewidth = 1) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  
  scale_color_manual(values = c(
    "Bulk sig + same direction" = "#20b2aa",
    "Bulk sig + opposite" = "#ff6b6b",
    "Nano sig only" = "#f59e0b",
    "Background" = "grey70"
  )) +
  
  labs(
    title = "Cross-Model Concordance of Gene Regulation",
    subtitle = "Shared directional shifts despite platform and model differences",
    x = "log₂FC VCD (Bulk RNA-seq)",
    y = "log₂FC TBI 7d (NanoString)"
  ) +
  theme_publication()

cor.test(plot_df$bulk_lfc, plot_df$nano_lfc_TBI_7d_vs_Control)


rank_bulk <- master %>%
  filter(!is.na(bulk_lfc)) %>%
  arrange(desc(bulk_lfc)) %>%
  pull(gene)

rank_nano <- master %>%
  filter(!is.na(nano_lfc_TBI_7d_vs_Control)) %>%
  arrange(desc(nano_lfc_TBI_7d_vs_Control)) %>%
  pull(gene)

# Step 2

# Step 3
pathway_df <- full_join(gsea_bulk, gsea_nano, by = "pathway")

ggplot(pathway_df, aes(NES_bulk, NES_nano)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(
    title = "Shared Pathway Dysregulation",
    subtitle = "Strong concordance at pathway level despite gene-level differences",
    x = "NES (VCD bulk)",
    y = "NES (TBI 7d)"
  ) +
  theme_publication()



bulk_sig_genes <- master %>%
  filter(sig_bulk) %>%
  pull(gene)

nano_ranked <- master %>%
  filter(!is.na(nano_lfc_TBI_7d_vs_Control)) %>%
  arrange(desc(nano_lfc_TBI_7d_vs_Control)) %>%
  pull(gene)

# ── 8. FIGURE: LFC Scatter — Bulk VCD vs TBI 7d ───────────────────────────
plot_df <- master %>%
  filter(!is.na(bulk_lfc), !is.na(nano_lfc_TBI7dvsControl)) %>%
  mutate(
    category = case_when(
      sig_bulk & sig_TBI7dvsControl ~ "Cross-dataset hit",
      sig_bulk                      ~ "Bulk VCD only",
      sig_TBI7dvsControl            ~ "TBI 7d only",
      TRUE                          ~ "Not significant"
    ),
    label = if_else(cross_dataset_sig, gene, NA_character_)
  )

p_scatter <- ggplot(plot_df, aes(bulk_lfc, nano_lfc_TBI7dvsControl, colour = category)) +
  geom_point(data = filter(plot_df, category == "Not significant"),
             alpha = 0.2, size = 0.8) +
  geom_point(data = filter(plot_df, category != "Not significant"),
             alpha = 0.8, size = 2.5) +
  geom_text_repel(aes(label = label), size = 2.8, max.overlaps = 20,
                  fontface = "italic", show.legend = FALSE) +
  geom_vline(xintercept = c(-1, 0, 1), linetype = c("dotted","dashed","dotted"), colour = "grey60") +
  geom_hline(yintercept = c(-1, 0, 1), linetype = c("dotted","dashed","dotted"), colour = "grey60") +
  scale_colour_manual(values = c("Cross-dataset hit" = "#20b2aa",
                                 "Bulk VCD only"     = "#a78bfa",
                                 "TBI 7d only"       = "#ff6b6b",
                                 "Not significant"   = "grey50")) +
  labs(title    = "Convergent Dysregulation: VCD Estrogen Loss × TBI 7d",
       subtitle = "Bottom-left quadrant = genes down in both models (high-confidence BBB targets)",
       x = "log₂FC  VCD vs Vehicle (Bulk, GSE279885)",
       y = "log₂FC  TBI 7d vs Sham (NanoString, GSE160651)",
       colour = NULL) +
  theme_publication()

save_figure(p_scatter, "01_scatter_bulk_vs_tbi7d", dirs$figures, width = 9, height = 7)

# ── 9. FIGURE: LFC Heatmap all cross-dataset hits ─────────────────────────
heat_df <- cross %>%
  select(gene, n_tbi_contrasts,
         `VCD Bulk`  = bulk_lfc,
         `TBI 1d`   = nano_lfc_TBI1dvsControl,
         `TBI 7d`   = nano_lfc_TBI7dvsControl,
         `TBI 30d`  = nano_lfc_TBI30dvsControl) %>%
  tidyr::pivot_longer(-c(gene, n_tbi_contrasts), names_to = "condition", values_to = "lfc") %>%
  mutate(condition = factor(condition, levels = c("VCD Bulk","TBI 1d","TBI 7d","TBI 30d")),
         gene      = factor(gene, levels = rev(cross$gene)))

p_heat <- ggplot(heat_df, aes(condition, gene, fill = lfc)) +
  geom_tile(colour = "grey20", linewidth = 0.3) +
  scale_fill_gradient2(low = "#4575b4", mid = "white", high = "#d73027",
                       midpoint = 0, name = "log₂FC") +
  labs(title = "Cross-Dataset Significant Genes — log₂FC",
       subtitle = "VCD estrogen depletion (left) and TBI at three timepoints (right)",
       x = NULL, y = NULL) +
  theme_publication(base_size = 9) +
  theme(axis.text.y = element_text(face = "italic"))

save_figure(p_heat, "02_heatmap_cross_genes", dirs$figures, width = 6, height = 8)


