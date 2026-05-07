# ==============================================================================
# analysis/GSE158960/03_annotate.R
#
# Purpose: Automated cell type annotation using SingleR, manual refinement,
#          and extraction of BBB-relevant cell types for downstream analysis.
#
# WHY SINGLER FOR INITIAL ANNOTATION:
#   SingleR correlates each cell's expression profile against curated reference
#   transcriptomes (e.g. MouseBrainAtlas) and assigns the most similar label.
#   This removes subjectivity from the first-pass annotation. We then manually
#   inspect marker genes to confirm or correct assignments for BBB-relevant
#   populations, which are the scientific priority.
#
# BBB-relevant cell types extracted for all downstream analysis:
#   - Endothelial cells: form the physical BBB (tight junctions, Cldn5, Tjp1)
#   - Pericytes/Mural cells: regulate BBB permeability (Pdgfrb, Rgs5)
#   - Astrocytes: endfeet contact endothelial cells; produce BBB support factors
#
# Outputs:
#   - seurat_annotated.rds     Full object with cell_type column
#   - seurat_bbb.rds           BBB cell subset (endothelial + pericyte + astrocyte)
# ==============================================================================

suppressPackageStartupMessages({
  suppressWarnings({
    library(here)
    library(Seurat)
    library(SingleR)
    library(celldex)
    library(dplyr)
    library(ggplot2)
    library(patchwork)
  })
})

source(here("analysis", "helpers", "utils.R"))
source(here("analysis", "helpers", "themes.R"))

log_section("GSE158960 | Step 3 | Cell Type Annotation")
dirs <- create_result_dirs("GSE158960")

seurat_merged <- readRDS(
  here("results", "GSE158960", "tables", "seurat_clustered.rds")
)

# ── 1. Marker gene sets for manual validation ─────────────────────────────────

MARKER_GENES <- list(
  Endothelial  = c("Cldn5", "Pecam1", "Cdh5", "Tjp1", "Ocln", "Esam"),
  Pericyte     = c("Pdgfrb", "Rgs5", "Acta2", "Notch3", "Des"),
  Astrocyte    = c("Gfap", "Aqp4", "Aldh1l1", "S100b", "Slc1a2"),
  Microglia    = c("Cx3cr1", "P2ry12", "Tmem119", "Iba1", "Csf1r"),
  Neuron       = c("Rbfox3", "Map2", "Tubb3", "Snap25", "Syn1"),
  Oligodendrocyte = c("Mog", "Plp1", "Mbp", "Cnp", "Mag"),
  OPC          = c("Pdgfra", "Cspg4", "Sox10", "Olig2"),
  # Estrogen receptor expression — key for this project
  Estrogen_receptors = c("Esr1", "Esr2", "Gper1", "Esrra", "Esrrb", "Esrrg")
)

# ── 2. SingleR automated annotation ──────────────────────────────────────────
#
# Reference: ImmGen + MouseBrainAtlas from celldex.
# WHY two references: ImmGen captures immune cells (microglia);
# MouseBrainAtlas captures neural/glial cell types.

singler_results <- load_or_compute(
  path = here("results", "GSE158960", "tables", "singler_results.rds"),
  expr = {
    log_msg("Loading SingleR reference ...")
    ref_brain <- celldex::MouseRNAseqData()

    # Get normalised count matrix for SingleR
    log_msg("Running SingleR ...")
    sce_counts <- GetAssayData(seurat_merged,
                                assay = "SCT", layer = "data")
    sce_broad <- SingleR(
      test      = sce_counts,
      ref       = ref_brain,
      labels    = ref_brain$label.main
    )
    
    sce_fine <- SingleR(
      test      = sce_counts,
      ref       = ref_brain,
      labels    = ref_brain$label.fine
    )
  }
)

# Add SingleR labels to metadata
seurat_merged$singler_label <- singler_results$labels
seurat_merged$singler_score <- apply(singler_results$scores, 1, max)

# ── 3. Cluster-level cell type assignment ────────────────────────────────────
#
# WHY assign per cluster not per cell:
#   SingleR per-cell labels are noisy. We aggregate by cluster majority vote,
#   then manually inspect the top markers for ambiguous clusters.

cluster_labels <- seurat_merged@meta.data %>%
  group_by(seurat_clusters, singler_label) %>%
  tally() %>%
  slice_max(n, n = 1) %>%
  ungroup() %>%
  dplyr::select(seurat_clusters, singler_label)

log_msg("Cluster -> SingleR label majority vote:")
print(cluster_labels)

plotScoreHeatmap(sce_broad, max.labels = 26, 
                 clusters = seurat_merged$seurat_clusters, 
                 order.by = "clusters", show_colnames = F)

plotScoreHeatmap(sce_fine,
                 clusters = seurat_merged$seurat_clusters, 
                 order.by = "clusters", show_colnames = F)

# ── 4. Manual cell type map ───────────────────────────────────────────────────
#
# After inspecting the cluster markers (see FeaturePlot below), apply your
# manual assignments here. The vector maps seurat_cluster number to cell type.
# BBB-relevant types use the exact names expected by subset_bbb_cells().
#
# IMPORTANT: Run this script once to get the UMAP and FeaturePlots, inspect
# them visually, then fill in this map before re-running to generate
# seurat_annotated.rds.
#

cluster_cell_type_map <- c(
  "0" = "Astrocyte",
  "1" = "Oligodendrocyte",
  "2" = "Endothelial",
  "3" = "Neuron",
  "4" = "Microglia",
  "5" = "NPCs",         # Finer SingleR analysis showed it is aNPCs
  "6" = "Pericyte",
  "7" = "Neuron",
  "8" = "Astrocyte",
  "9" = "Epithelial",
  "10" = "Neuron",
  "11" = "Pericyte",
  "12" = "Neuron",
  "13" = "Microglia",
  "14" = "Neuron",
  "15" = "Oligodendrocyte",
  "16" = "Pericyte",
  "17" = "Oligodendrocyte",
  "18" = "Neuron"
  
  # setNames(cluster_labels$singler_label,
  #          as.character(cluster_labels$seurat_clusters))
)

# Standardise label names to match PALETTE_CELL_TYPES
# Map any SingleR-specific labels to canonical names used across this project
# label_harmonise <- c(
#   "Endothelial cells" = "Endothelial",
#   "Pericytes"         = "Pericyte",
#   "Astrocytes"        = "Astrocyte",
#   "Oligodendrocytes"  = "Oligodendrocyte",
#   "Microglia"         = "Microglia",
#   "Neurons"           = "Neuron"
# )

seurat_merged$cell_type <- dplyr::recode(
  cluster_cell_type_map[as.character(seurat_merged$seurat_clusters)],
  !!!label_harmonise
)

seurat_merged$cell_type[is.na(seurat_merged$cell_type)] <- "Other"

# ── 5. Marker gene FeaturePlots for validation ────────────────────────────────
#
# WHY: these figures are how you verify the annotation. Open them in RStudio
# Viewer or the PDF and confirm that Cldn5+ cells match the "Endothelial"
# cluster, Gfap+ cells match "Astrocyte", etc.

for (cell_type in names(MARKER_GENES)) {
  genes_present <- intersect(MARKER_GENES[[cell_type]], rownames(seurat_merged))
  if (length(genes_present) == 0) next

  p <- FeaturePlot(seurat_merged,
                    features  = genes_present[1:min(6, length(genes_present))],
                    reduction = "umap",
                    ncol      = 3,
                    cols      = c("grey90", "#D73027"),
                    order     = TRUE) &
    theme_publication(base_size = 9)

  safe_name <- gsub("[^A-Za-z0-9_]", "_", cell_type)
  save_figure(p,
              paste0("04_markers_", safe_name),
              dirs$figures, width = 12, height = 8)
}

# ── 6. Annotated UMAP ────────────────────────────────────────────────────────

p_annotated <- DimPlot(seurat_merged, reduction = "umap",
                        group.by = "cell_type",
                        label = TRUE, repel = TRUE,
                        label.size = 3.5) +
  theme_publication() +
  labs(title   = "GSE158960 - Cell type annotation")

save_figure(p_annotated, "05_umap_annotated",
            dirs$figures, width = 10, height = 8)

saveRDS(seurat_merged,
        here("results", "GSE158960", "tables", "seurat_annotated_tmp.rds"))
log_msg("Saved: seurat_annotated.rds")

# ── 7. Extract BBB cell types ─────────────────────────────────────────────────
#
# This is the core subset for BBB analysis: endothelial cells, pericytes,
# and astrocytes. All downstream ligand-receptor and estrogen analyses
# will operate primarily on this subset.

seurat_merged@graphs <- list()

seurat_bbb <- subset_bbb_cells(
  seurat_merged,
  cell_types = c("Endothelial", "Pericyte", "Astrocyte")
)

log_msg("BBB cells (Endothelial + Pericyte + Astrocyte): ", ncol(seurat_bbb))
log_msg("BBB cell type breakdown:")
print(table(seurat_bbb$cell_type))

# Re-cluster BBB subset at higher resolution to resolve sub-populations
seurat_bbb <- FindNeighbors(seurat_bbb, reduction = "harmony",
                             dims = 1:PARAMS$seurat$n_pcs, verbose = FALSE)
seurat_bbb <- FindClusters(seurat_bbb, resolution = 0.8, verbose = FALSE)
seurat_bbb <- RunUMAP(seurat_bbb, reduction = "harmony",
                       dims = 1:PARAMS$seurat$n_pcs, verbose = FALSE)

p_bbb_umap <- DimPlot(seurat_bbb, reduction = "umap",
                       group.by = "cell_type",
                       cols = PALETTE_CELL_TYPES) +
  DimPlot(seurat_bbb, reduction = "umap",
          group.by = "group",
          cols = PALETTE_GSE158960, alpha = 0.6) +
  patchwork::plot_layout(ncol = 2) +
  patchwork::plot_annotation(
    title = "BBB cell populations - GSE158960",
    theme = theme_publication()
  )

save_figure(p_bbb_umap, "06_umap_bbb_subset",
            dirs$figures, width = 14, height = 6)

saveRDS(seurat_bbb,
        here("results", "GSE158960", "tables", "seurat_bbb.rds"))
log_msg("Saved: seurat_bbb.rds")

log_section("GSE158960 | Step 3 | COMPLETE")
message("Next: run analysis/GSE158960/04_bbb_estrogen.R")
