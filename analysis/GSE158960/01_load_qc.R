# ==============================================================================
# analysis/GSE158960/01_load_qc.R
#
# Purpose: Load all 16 GSE158960 scRNA-seq samples, calculate per-cell QC
#          metrics, visualize them, and apply threshold-based cell filtering.
#
# WHY PER-SAMPLE QC BEFORE MERGING:
#   Merging first and then filtering loses the per-sample context needed to
#   diagnose problems. A sample with 5% mitochondrial reads in one group and
#   20% in another suggests a tissue-dissociation quality difference — that
#   information disappears once you pool all cells. We filter per sample,
#   then merge the cleaned objects.
#
# QC metrics and their biology:
#   nFeature_RNA (genes detected):  200 lower bound removes empty droplets;
#                                   4000 upper bound removes doublets.
#   percent.mt (mito reads):        >15% in hypothalamus suggests dying cells
#                                   (mito transcripts leak out of damaged cells
#                                   before nuclear RNA does).
#   percent.ribo (ribo reads):      tracked for QC but not used as a hard filter
#                                   (ribosomal content varies naturally by cell type).
#   percent.hb (hemoglobin reads):  >5% indicates red blood cell contamination.
#
# Outputs saved to results/GSE158960/tables/:
#   - seurat_list_filtered.rds   List of 16 per-sample filtered Seurat objects
# ==============================================================================

suppressPackageStartupMessages({
  suppressWarnings({
    library(here)
    library(Seurat)
    library(dplyr)
    library(ggplot2)
    library(patchwork)
  })
})

source(here("analysis", "helpers", "utils.R"))
source(here("analysis", "helpers", "themes.R"))

log_section("GSE158960 | Step 1 | Load & QC")
dirs <- create_result_dirs("GSE158960")

# ── Sample map ────────────────────────────────────────────────────────────────
# Each row defines one sample. The data_dir column points to the folder
# containing barcodes.tsv.gz / features.tsv.gz / matrix.mtx.gz.

sample_map <- tibble::tribble(
  ~sample_id,             ~genotype, ~treatment, ~group,
  "WT_Control_s01",       "Wildtype", "Control",  "WT_Control",
  "WT_Control_s05",       "Wildtype", "Control",  "WT_Control",
  "WT_Control_s10",       "Wildtype", "Control",  "WT_Control",
  "WT_Tamoxifen_s03",     "Wildtype", "Tamoxifen","WT_Tamoxifen",
  "WT_Tamoxifen_s08",     "Wildtype", "Tamoxifen","WT_Tamoxifen",
  "WT_Tamoxifen_s09",     "Wildtype", "Tamoxifen","WT_Tamoxifen",
  "WT_Tamoxifen_s13",     "Wildtype", "Tamoxifen","WT_Tamoxifen",
  "WT_Tamoxifen_s16",     "Wildtype", "Tamoxifen","WT_Tamoxifen",
  "Esr1cKO_Control_s04",  "Esr1cKO", "Control",  "Esr1cKO_Control",
  "Esr1cKO_Control_s07",  "Esr1cKO", "Control",  "Esr1cKO_Control",
  "Esr1cKO_Control_s12",  "Esr1cKO", "Control",  "Esr1cKO_Control",
  "Esr1cKO_Control_s14",  "Esr1cKO", "Control",  "Esr1cKO_Control",
  "Esr1cKO_Tamoxifen_s02","Esr1cKO", "Tamoxifen","Esr1cKO_Tamoxifen",
  "Esr1cKO_Tamoxifen_s06","Esr1cKO", "Tamoxifen","Esr1cKO_Tamoxifen",
  "Esr1cKO_Tamoxifen_s11","Esr1cKO", "Tamoxifen","Esr1cKO_Tamoxifen",
  "Esr1cKO_Tamoxifen_s15","Esr1cKO", "Tamoxifen","Esr1cKO_Tamoxifen"
) %>%
  dplyr::mutate(
    data_dir = here("data", "raw", "GSE158960", sample_id)
  )

for (i in seq_len(nrow(sample_map))) {
  sid  <- sample_map$sample_id[i]
  ddir <- sample_map$data_dir[i]

  if (dir.exists(ddir)) {
    message("Data directory found for sample: ", sid, " at path: ", ddir)
  }
}

log_msg("Samples defined: ", nrow(sample_map))

# ── QC thresholds from config ─────────────────────────────────────────────────
qc <- PARAMS$qc$scrna

# ── 1. Load and filter each sample ───────────────────────────────────────────

seurat_list_filtered <- load_or_compute(
  path  = here("results", "GSE158960", "tables", "seurat_list_filtered.rds"),
  expr  = {
    filtered <- list()

    for (i in seq_len(nrow(sample_map))) {
      sid  <- sample_map$sample_id[i]
      ddir <- sample_map$data_dir[i]

      if (!dir.exists(ddir)) {
        warning("Data directory not found: ", ddir, ". Skipping.")
        next
      }

      log_msg("Loading: ", sid)

      # Read 10X CellRanger output (barcodes / features / matrix)
      counts <- Read10X(data.dir = ddir)

      # Create Seurat object. min.cells = 3 drops genes seen in fewer than
      # 3 cells — removes very sparse genes that add noise without information.
      obj <- CreateSeuratObject(
        counts    = counts,
        project   = sid,
        min.cells = 3,
        min.features = qc$min_features
      )

      # ── QC metrics ──
      # WHY PercentageFeatureSet: Seurat calculates the fraction of UMIs
      # attributable to a feature pattern. Mouse mito genes start with "mt-".
      obj$percent.mt   <- PercentageFeatureSet(obj, pattern = "^mt-")
      obj$percent.ribo <- PercentageFeatureSet(obj, pattern = "^Rp[sl]")
      obj$percent.hb   <- PercentageFeatureSet(obj, pattern = "^Hb[ab]")
      obj$sample_id    <- sid
      obj$genotype     <- sample_map$genotype[i]
      obj$treatment    <- sample_map$treatment[i]
      obj$group        <- sample_map$group[i]
      obj$dataset      <- "GSE158960"
      obj$data_type    <- "scRNA"

      cells_before <- ncol(obj)

      # Apply thresholds
      obj <- subset(
        obj,
        subset = nFeature_RNA >= qc$min_features &
                 nFeature_RNA <= qc$max_features &
                 percent.mt  <= qc$max_pct_mt    &
                 percent.hb  <= qc$max_pct_hb
      )

      cells_after <- ncol(obj)
      log_msg(sprintf("  %s: %d -> %d cells (%.1f%% retained)",
                      sid, cells_before, cells_after,
                      100 * cells_after / cells_before))

      filtered[[sid]] <- obj
    }
    filtered
  }
)

log_msg("Samples loaded and filtered: ", length(seurat_list_filtered))
total_cells <- sum(sapply(seurat_list_filtered, ncol))
log_msg("Total cells after filtering: ", total_cells)

# ── 2. QC summary figure ─────────────────────────────────────────────────────
#
# WHY a combined QC plot before moving on:
#   This is the one figure that tells you if any sample is an outlier.
#   A sample with dramatically different nFeature distribution or percent.mt
#   should be investigated before it contaminates downstream clustering.

qc_summary <- purrr::map_dfr(seurat_list_filtered, function(obj) {
  obj@meta.data
})

p_vln_features <- ggplot(qc_summary,
                          aes(x = sample_id, y = nFeature_RNA, fill = group)) +
  geom_violin(scale = "width", alpha = 0.8) +
  geom_boxplot(width = 0.1, outlier.size = 0.5, fill = "white") +
  scale_fill_manual(values = PALETTE_GSE158960) +
  labs(title = "Genes detected per cell — post-filter",
       x = NULL, y = "nFeature_RNA") +
  theme_publication() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        legend.position = "none")

p_vln_mt <- ggplot(qc_summary,
                    aes(x = sample_id, y = percent.mt, fill = group)) +
  geom_violin(scale = "width", alpha = 0.8) +
  geom_boxplot(width = 0.1, outlier.size = 0.5, fill = "white") +
  scale_fill_manual(values = PALETTE_GSE158960) +
  geom_hline(yintercept = qc$max_pct_mt, linetype = "dashed", colour = "grey40") +
  labs(title = "Mitochondrial % — post-filter",
       x = NULL, y = "% mitochondrial reads") +
  theme_publication() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        legend.position = "none")

p_qc_combined <- p_vln_features / p_vln_mt +
  plot_annotation(
    title = "GSE158960 — QC metrics after filtering",
    subtitle = "Violin + boxplots of nFeature_RNA and percent.mt per sample"
  ) &
  theme_publication(base_size = 10)

save_figure(p_qc_combined, "01_qc_postfilter_tmp",
            dirs$figures, width = 14, height = 8)

log_section("GSE158960 | Step 1 | COMPLETE")
message("Next: run analysis/GSE158960/02_doublets_merge.R")
