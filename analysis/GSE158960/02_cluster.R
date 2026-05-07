# ==============================================================================
# analysis/GSE158960/02_cluster.R
#
# Purpose: Doublet removal, merge all 16 samples, normalize, find variable
#          features, scale, reduce dimensions, and cluster.
#
# WHY DOUBLET DETECTION PER SAMPLE (not on merged object):
#   DoubletFinder simulates doublets by averaging randomly selected cell pairs
#   and asks which real cells look like those simulations. This only works
#   correctly when run on one sample at a time because doublets are specific
#   to a sequencing run — a "doublet" in the merged object could just be a
#   transition state between two real cell types from different samples.
#
# WHY SCTRANSFORM v2 FOR NORMALIZATION:
#   Standard log-normalization assumes constant library size across cells.
#   SCTransform v2 fits a regularized negative-binomial regression to the UMI
#   counts, removing technical variation (sequencing depth) while preserving
#   biological signal. It is particularly well-suited to scRNA-seq where
#   library size varies >10-fold between cells.
#
# WHY HARMONY FOR BATCH CORRECTION:
#   We merge 16 samples from 4 experimental groups; each sample is a separate
#   10X capture, so there is inevitable sample-level batch variation.
#   Harmony corrects for this in the PCA embedding space (not in gene
#   expression space), which is faster and avoids over-correction.
#
# Outputs:
#   - seurat_merged.rds    Merged, normalized, clustered Seurat v5 object
# ==============================================================================

suppressPackageStartupMessages({
  suppressWarnings({
    library(here)
    library(Seurat)
    library(harmony)
    library(dplyr)
    library(ggplot2)
    library(patchwork)
    library(future)
  })
})

source(here("analysis", "helpers", "utils.R"))
source(here("analysis", "helpers", "themes.R"))

log_section("GSE158960 | Step 2 | Doublets, Merge, Cluster")
dirs <- create_result_dirs("GSE158960")

seurat_list <- readRDS(here("results", "GSE158960", "tables",
                             "seurat_list_filtered.rds"))

plan("multisession", workers = 2)        # use 4 CPU cores
options(future.globals.maxSize = 4000 * 1024^2)  # 8 GB RAM limit per worker

# ── 1. Doublet detection with DoubletFinder ───────────────────────────────────
#
# The workflow per sample is:
#   a) Quick normalise + PCA + UMAP (DoubletFinder needs a UMAP embedding)
#   b) Estimate the expected doublet rate (0.8% per 1000 cells; 10X guideline)
#   c) Run DoubletFinder; mark cells as Singlet / Doublet
#   d) Remove doublets before merging
#
# This is optional but strongly recommended. If DoubletFinder is not installed,
# the loop is skipped and a warning is printed.

run_doublet_finder <- requireNamespace("DoubletFinder", quietly = TRUE)

if (!run_doublet_finder) {
  warning("DoubletFinder not installed. Skipping doublet removal.\n",
          "Install with: remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')")
}

seurat_list_clean <- load_or_compute(
  path = here("results", "GSE158960", "tables", "seurat_list_nodbl.rds"),
  expr = {
    clean <- list()
    for (sid in names(seurat_list)) {
      obj <- seurat_list[[sid]]

      if (run_doublet_finder) {
        log_msg("DoubletFinder: ", sid)

        # Quick pre-processing needed by DoubletFinder
        obj <- NormalizeData(obj, verbose = FALSE)
        obj <- FindVariableFeatures(obj, nfeatures = 2000, verbose = FALSE)
        obj <- ScaleData(obj, verbose = FALSE)
        obj <- RunPCA(obj, npcs = 20, verbose = FALSE)
        obj <- RunUMAP(obj, dims = 1:20, verbose = FALSE)

        # Expected doublet rate: ~0.8% per 1000 cells captured
        n_exp_dbl <- round(ncol(obj) * 0.008 * (ncol(obj) / 1000))
        n_exp_dbl <- max(n_exp_dbl, 1)

        # pK = neighbourhood size; 0.09 is a commonly used default
        # For a fully rigorous analysis, run paramSweep to find optimal pK
        obj <- DoubletFinder::doubletFinder(
          obj, PCs = 1:20, pK = 0.09, nExp = n_exp_dbl
        )

        # The classification column is named dynamically; find it
        dbl_col <- grep("DF.classifications", colnames(obj@meta.data),
                        value = TRUE)[1]
        n_dbl <- sum(obj@meta.data[[dbl_col]] == "Doublet")
        log_msg(sprintf("  Removed %d doublets (%.1f%%)",
                        n_dbl, 100 * n_dbl / ncol(obj)))
        obj <- subset(obj, subset = !!rlang::sym(dbl_col) == "Singlet")

        # Clean up extra metadata columns added by DoubletFinder
        drop_cols <- grep("pANN|DF.classifications", colnames(obj@meta.data),
                          value = TRUE)
        obj@meta.data[drop_cols] <- NULL
      }

      clean[[sid]] <- obj
    }
    clean
  }
)

# ── 2. Merge all samples into one Seurat v5 object ───────────────────────────
#
# merge() with add.cell.ids prepends the sample name to each cell barcode,
# preventing barcode collisions across samples (each 10X run uses the same
# ~6700 barcodes). The result is a Seurat v5 object with one layer per sample.

seurat_merged <- load_or_compute(
  path = here("results", "GSE158960", "tables", "seurat_merged.rds"),
  expr = {
    log_msg("Merging ", length(seurat_list_clean), " samples ...")
    obj_list <- unname(seurat_list_clean)
    merged <- merge(
      obj_list[[1]],
      y           = obj_list[-1],
      add.cell.ids = names(seurat_list_clean),
      project      = "GSE158960"
    )
    log_msg("Total cells: ", ncol(merged))
    merged
  }
)

# ── 3. Normalisation with SCTransform v2 ─────────────────────────────────────
#
# vars.to.regress: we regress out percent.mt to remove the confounding effect
#   of mitochondrial read fraction on gene expression. We do NOT regress out
#   cell cycle by default because cell cycle differences between WT and
#   Esr1cKO are biologically relevant to estrogen signaling.
#
# vst.flavor = "v2": uses the improved theta estimation from Choudhary & Satija 2022.

plan("sequential")

seurat_merged <- load_or_compute(
  path = here("results", "GSE158960", "tables", "seurat_sct.rds"),
  expr = {
    log_msg("Running SCTransform v2 ...")
    SCTransform(
      seurat_merged,
      vst.flavor     = "v2",
      vars.to.regress = "percent.mt",
      verbose        = FALSE
    )
  }
)

# ── 4. PCA, Harmony batch correction, UMAP, clustering ───────────────────────

seurat_merged <- load_or_compute(
  path = here("results", "GSE158960", "tables", "seurat_clustered.rds"),
  force = F,
  expr = {
    n_pcs <- PARAMS$seurat$n_pcs
    
    obj <- seurat_merged
    DefaultAssay(obj) <- "SCT"

    log_msg("PCA ...")
    obj <- RunPCA(obj,
                  npcs    = n_pcs,
                  verbose = FALSE)

    # Elbow plot to guide PC selection — save for inspection
    p_elbow <- ElbowPlot(obj, ndims = n_pcs) +
      theme_publication() +
      labs(title = "PCA elbow plot — GSE158960",
           subtitle = "Choose PCs where curve flattens")
    save_figure(p_elbow, "02_elbow_plot", dirs$figures, width = 6, height = 4)

    # Harmony: corrects the PCA embedding for sample-level batch effects.
    # group.by.vars = "sample_id": treat each 10X capture as one batch.
    log_msg("Harmony batch correction ...")
    obj <- RunHarmony(
      obj,
      group.by.vars = "sample_id",
      reduction     = "pca",
      reduction.save = "harmony",
      verbose       = FALSE
    )

    # UMAP on Harmony-corrected embedding
    log_msg("UMAP ...")
    obj <- RunUMAP(obj, reduction = "harmony",
                   dims = 1:n_pcs, verbose = FALSE)

    # Nearest-neighbour graph + Louvain clustering
    # WHY two resolutions: 0.5 gives broad clusters for cell type annotation;
    # fine-tune later if needed.
    log_msg("Clustering ...")
    obj <- FindNeighbors(obj, reduction = "harmony",
                         dims = 1:n_pcs, verbose = FALSE)
    obj <- FindClusters(obj, resolution = PARAMS$seurat$resolution,
                        verbose = FALSE)
    obj
  }
)

log_msg("Clusters: ",  nlevels(seurat_merged$seurat_clusters))
log_msg("Total cells: ", ncol(seurat_merged))

# ── 5. UMAP overview figures ──────────────────────────────────────────────────

p_umap_cluster <- DimPlot(seurat_merged, reduction = "umap",
                            group.by = "seurat_clusters",
                            label = TRUE, repel = TRUE, label.size = 3.5) +
  theme_publication() +
  labs(title = "GSE158960 — UMAP by cluster")

p_umap_group <- DimPlot(seurat_merged, reduction = "umap",
                         group.by = "group",
                         cols = PALETTE_GSE158960,
                         alpha = 0.5) +
  theme_publication() +
  labs(title = "GSE158960 — UMAP by experimental group")

save_figure(p_umap_cluster / p_umap_group, "03_umap_overview_tmp",
            dirs$figures, width = 10, height = 12)

log_section("GSE158960 | Step 2 | COMPLETE")
message("Next: run analysis/GSE158960/03_annotate.R")
