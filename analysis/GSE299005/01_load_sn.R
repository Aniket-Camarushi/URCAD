# ==============================================================================
# analysis/GSE299005/01_load_sn.R
#
# Purpose: Load GSE299005 snRNA-seq data (FPI male mice), run QC, normalize,
#          integrate two DGE matrices, cluster, and annotate cell types.
#
# DATA FORMAT NOTE:
#   GSE299005 provides data as sparse DGE matrices (.mtx.gz) + cell metadata
#   (.csv.gz), NOT as standard 10X CellRanger output. This requires a
#   different loading strategy than Read10X().
#   Files:
#     GSM9030357_DGE1.mtx.gz + GSM9030357_cell_metadata1.csv.gz  (sample 1)
#     GSM9030358_DGE2.mtx.gz + GSM9030358_cell_metadata2.csv.gz  (sample 2)
#
# SNRNA-SEQ VS SCRNA-SEQ — KEY DIFFERENCES:
#   1. Nuclei not cells: nuclear RNA only, so cytoplasmic markers are absent.
#      Mitochondrial reads are very low (<5%) because mitochondria stay in
#      the cytoplasm; cells with >5% mt are likely cytoplasmic contamination.
#   2. Higher proportion of pre-mRNA (intronic reads): splicing is incomplete
#      in nuclei. This slightly increases nFeature_RNA for the same cell type.
#   3. Cell type composition may differ: neuron-dense tissues yield more
#      oligodendrocytes and fewer endothelial cells in snRNA-seq because
#      nuclei from small cells (capillary endothelium) are harder to capture.
#
# FPI MODEL CONTEXT:
#   Fluid Percussion Injury causes diffuse TBI. Male mice only (no estrogen).
#   This dataset provides the injury context: how does the BBB respond to TBI
#   in the absence of estrogen? Cross-referencing with GSE158960 (estrogen
#   present/absent via ERα KO) enables us to ask whether estrogen loss
#   phenocopies aspects of TBI-induced BBB dysfunction.
#
# Outputs:
#   - seurat_sn_annotated.rds   Annotated snRNA object
#   - seurat_sn_bbb.rds         BBB subset (endothelial + pericyte + astrocyte)
# ==============================================================================

suppressPackageStartupMessages({
  suppressMessages({
    library(here)
    library(Seurat)
    library(Matrix)
    library(dplyr)
    library(readr)
    library(ggplot2)
    library(patchwork)
    library(purrr)
    library(knitr)
    library(SingleR)
    library(celldex)
    library(DoubletFinder)
  })
})

source(here("analysis", "helpers", "utils.R"))
source(here("analysis", "helpers", "themes.R"))

options(future.globals.maxSize = 8 * 1024^3)

log_section("GSE299005 | Step 1 | Load snRNA-seq")
dirs <- create_result_dirs("GSE299005")

data_dir <- here("data", "raw", "GSE299005")

# ── Sample metadata ───────────────────────────────────────────────────────────
# Both samples are from male FPI mice (no female, no sham controls in this GEO
# submission). GSM9030357 and GSM9030358 are two biological replicates.

sample_info <- tibble::tribble(
  ~gsm_id,        ~dge_file,                          ~meta_file,                              ~condition,
  "GSM9030357",   "GSM9030357_DGE1.mtx.gz",           "GSM9030357_cell_metadata1.csv.gz",       "FPI",
  "GSM9030358",   "GSM9030358_DGE2.mtx.gz",           "GSM9030358_cell_metadata2.csv.gz",       "FPI"
)

# ── 1. Load DGE matrices ───────────────────────────────────────────────────────
#
# The DGE (Digital Gene Expression) format stores counts as a sparse matrix
# where rows = genes and columns = cell barcodes.
# readMM reads the MatrixMarket format (.mtx); we need to attach gene/cell
# names from the metadata CSV.

load_dge_sample <- function(gsm_id, dge_file, meta_file, condition) {
  mtx_path  <- file.path(data_dir, dge_file)
  meta_path <- file.path(data_dir, meta_file)

  if (!file.exists(mtx_path)) {
    warning("DGE matrix not found: ", mtx_path)
    return(NULL)
  }

  log_msg("Loading: ", gsm_id)

  # Read sparse matrix
  mat  <- Matrix::readMM(gzcon(file(mtx_path, "rb")))

  # Read cell metadata (contains barcodes and any pre-computed annotations)
  meta <- readr::read_csv(meta_path, show_col_types = FALSE)
  
  mat <- Matrix::t(mat)

  # The first column of meta should be cell barcodes
  barcodes <- meta[[1]]

  # Dimensions: rows = genes, cols = cells
  # If meta has gene names, attach them. Otherwise rows are gene indices.
  # if ("gene" %in% colnames(meta) || nrow(mat) == nrow(meta)) {
  #   # Gene names in metadata rows
  #   gene_names <- if ("gene" %in% colnames(meta)) meta$gene else
  #                 paste0("Gene_", seq_len(nrow(mat)))
  #   rownames(mat) <- gene_names
  #   colnames(mat) <- barcodes
  # } else {
  #   # Standard orientation: rows = genes, cols = barcodes
  #   colnames(mat) <- barcodes
  # }
 
  # Attaching gene names
  genes_path <- file.path(data_dir, "all_genes.csv")
  genes      <- readr::read_csv(genes_path, show_col_types = FALSE)
  if (nrow(genes) != nrow(mat))
    stop(gsm_id, ": gene count mismatch — ", nrow(genes), " vs ", nrow(mat))
  rownames(mat) <- make.unique(genes$gene_name)
  
  if (ncol(mat) != length(barcodes)) {
    stop(gsm_id, ": ncol(mat) ", ncol(mat),
         " != nrow(meta) ", length(barcodes),
         " — check matrix orientation")
  }
  colnames(mat) <- paste0(gsm_id, "_", barcodes)

  # Create Seurat object
  obj <- CreateSeuratObject(
    counts       = mat,
    project      = gsm_id,
    min.cells    = 3,
    min.features = PARAMS$qc$snrna$min_features
  )
  
  rm(mat); gc()

  # Attach provided metadata columns
  shared_cols <- intersect(colnames(meta), c("cell_type", "cluster",
                                              "umap_1", "umap_2"))
  for (col in shared_cols) {
    matched <- meta[[col]][match(
      sub(paste0("^", gsm_id, "_"), "", colnames(obj)),  # strip prefix to match
      barcodes
    )]
    
    obj[[col]] <- matched
  }

  # Add project-level metadata
  obj$sample_id  <- gsm_id
  obj$condition  <- condition
  obj$dataset    <- "GSE299005"
  obj$data_type  <- "snRNA"
  obj$sex        <- "male"
  obj$model      <- "FPI"

  obj
}

sn_list_raw <- mapply(
  load_dge_sample,
  gsm_id    = sample_info$gsm_id,
  dge_file  = sample_info$dge_file,
  meta_file = sample_info$meta_file,
  condition = sample_info$condition,
  SIMPLIFY  = FALSE
)

sn_list_raw <- Filter(Negate(is.null), sn_list_raw)

# ── 2. QC metrics and filtering ───────────────────────────────────────────────

qc <- PARAMS$qc$snrna

sn_list_filtered <- load_or_compute(
  path = here("results", "GSE299005", "tables", "sn_list_filtered.rds"),
  force = F,
  expr = {
    lapply(sn_list_raw, function(obj) {
      # mt < 5% for nuclei (cytoplasmic contamination threshold)
      obj$percent.mt   <- PercentageFeatureSet(obj, pattern = "^mt-")
      obj$percent.ribo <- PercentageFeatureSet(obj, pattern = "^Rp[sl]")

      cells_before <- ncol(obj)
      obj <- subset(obj,
                    subset = nFeature_RNA >= qc$min_features &
                             nFeature_RNA <= qc$max_features &
                             percent.mt   <= qc$max_pct_mt)
      log_msg(sprintf("  %s: %d -> %d nuclei",
                      obj$sample_id[1], cells_before, ncol(obj)))
      obj
    })
  }
)

rm(sn_list_raw); gc()

qc_summary <- purrr::map_dfr(sn_list_filtered, function(obj) {
  df <- obj@meta.data %>%
    dplyr::select(sample_id, nFeature_RNA, nCount_RNA, percent.mt, percent.ribo) %>%
    dplyr::mutate(group = "FPI") %>%  # Single condition
    tibble::as_tibble()
  return(df)
})

# Violin + boxplot QC panels
p_vln_features <- ggplot(qc_summary, aes(x = sample_id, y = nFeature_RNA, fill = group)) +
  geom_violin(scale = "width", alpha = 0.8, draw_quantiles = c(0.25, 0.5, 0.75)) +
  geom_boxplot(width = 0.1, outlier.size = 0.5, fill = "white", alpha = 0.9) +
  scale_fill_manual(values = "#E31A1C") +  # FPI red (adapt PALETTE_GSE158960)
  labs(title = "Genes detected per nucleus — post-filter",
       x = NULL, y = "nFeature_RNA") +
  theme_publication() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
        legend.position = "none")

p_vln_count <- ggplot(qc_summary, aes(x = sample_id, y = nCount_RNA, fill = group)) +
  geom_violin(scale = "width", alpha = 0.8, draw_quantiles = c(0.25, 0.5, 0.75)) +
  geom_boxplot(width = 0.1, outlier.size = 0.5, fill = "white", alpha = 0.9) +
  scale_fill_manual(values = "#E31A1C") +
  labs(title = "Reads per nucleus — post-filter",
       x = NULL, y = "nCount_RNA") +
  theme_publication() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
        legend.position = "none")

p_vln_mt <- ggplot(qc_summary, aes(x = sample_id, y = percent.mt, fill = group)) +
  geom_violin(scale = "width", alpha = 0.8, draw_quantiles = c(0.25, 0.5, 0.75)) +
  geom_boxplot(width = 0.1, outlier.size = 0.5, fill = "white", alpha = 0.9) +
  scale_fill_manual(values = "#E31A1C") +
  geom_hline(yintercept = qc$max_pct_mt, linetype = "dashed", colour = "grey40", linewidth = 0.8) +
  labs(title = "Mitochondrial % — post-filter",
       x = NULL, y = "% mitochondrial reads") +
  theme_publication() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
        legend.position = "none")

# Combined QC figure
p_qc_combined <- (p_vln_features | p_vln_count) / 
  p_vln_mt +
  plot_annotation(
    title = "GSE299005 — snRNA-seq QC (FPI male mice)",
    subtitle = "Post-filtering metrics by sample (nFeature_RNA, nCount_RNA, percent.mt)"
  ) &
  theme(plot.title = element_text(size = 14, face = "bold"),
        plot.subtitle = element_text(size = 11))

save_figure(p_qc_combined, "01_qc_violins_postfilter_tmp", 
            dirs$figures, width = 12, height = 8)

qc_summary %>%
  group_by(sample_id) %>%
  summarise(
    n_nuclei = n(),
    nFeature_median = median(nFeature_RNA),
    nFeature_q95 = quantile(nFeature_RNA, 0.95),
    nCount_median = median(nCount_RNA),
    nCount_q95 = quantile(nCount_RNA, 0.95),
    pct_mt_median = median(percent.mt),
    .groups = "drop"
  ) %>%
  knitr::kable(digits = 0) %>%
  print()

# ── 3. Merge and normalize ─────────────────────────────────────────────────────

sn_list_sct <- load_or_compute(
  path  = here("results", "GSE299005", "tables", "sn_list_sct.rds"),
  force = F,
  expr  = {
    log_msg("SCTransforming ", length(sn_list_filtered), " samples individually ...")
    lapply(sn_list_filtered, function(obj) {
      log_msg("  SCTransform: ", unique(obj$sample_id),
              "  (", ncol(obj), " nuclei x ", nrow(obj), " genes)")
      SCTransform(
        obj,
        vst.flavor      = "v2",
        vars.to.regress = "percent.mt",
        verbose         = TRUE,          # shows per-gene progress
        conserve.memory = TRUE,          # set TRUE if RAM < 32 GB
        return.only.var.genes = TRUE
      )
    })
  }
)

rm(sn_list_filtered); gc()

seurat_sn <- load_or_compute(
  path  = here("results", "GSE299005", "tables", "seurat_sn_merged.rds"),
  force = F,
  expr  = {
    
    log_msg("Stripping RNA assay to reduce memory footprint ...")
    sn_list_sct <- lapply(sn_list_sct, function(obj) {
      obj[["RNA"]] <- NULL
      obj
    })
    gc(); gc()
    
    # Compute var features BEFORE merge — avoids holding list + merged simultaneously
    var_features <- SelectIntegrationFeatures(
      object.list = sn_list_sct,
      nfeatures   = 3000
    )
    
    # Unlist, then immediately free named list
    log_msg("Merging ", length(sn_list_sct), " SCTransformed samples ...")
    obj_list <- unname(sn_list_sct)
    rm(sn_list_sct); gc(); gc()
    
    # Merge, then immediately free source objects
    merged <- merge(obj_list[[1]], y = obj_list[-1], project = "GSE299005")
    rm(obj_list); gc(); gc()
    
    DefaultAssay(merged)     <- "SCT"
    VariableFeatures(merged) <- var_features
    merged
    
    # log_msg("Merging ", length(sn_list_sct), " SCTransformed samples ...")
    # obj_list <- unname(sn_list_sct)
    # merged   <- merge(
    #   obj_list[[1]],
    #   y            = obj_list[-1],
    #   add.cell.ids = NULL,
    #   project      = "GSE299005"
    # )
    # DefaultAssay(merged) <- "SCT"
    # 
    # VariableFeatures(merged) <- SelectIntegrationFeatures(
    #   object.list = sn_list_sct,
    #   nfeatures   = 3000
    # )
    # 
    # rm(sn_list_sct, obj_list); gc()
    # merged
  }
)

if (exists("sn_list_sct")) { rm(sn_list_sct); gc() }

# ── Block 2: Dim reduction + clustering ──────────────────────────────────────


# seurat_sn <- load_or_compute(
#   path  = here("results", "GSE299005", "tables", "seurat_sn_merged.rds"),
#   force = FALSE,
#   expr  = {
#     log_msg("Merging ", length(sn_list_sct), " SCTransformed samples ...")
#     # obj_list <- unname(sn_list_sct)
#     # merged   <- merge(
#     #   obj_list[[1]],
#     #   y            = obj_list[-1],
#     #   add.cell.ids = NULL,
#     #   project      = "GSE299005"
#     # )
#     # DefaultAssay(merged) <- "SCT"
#     # 
#     # # Select variable features that are consistently variable across both
#     # # samples — required after per-sample SCT before RunPCA
#     # VariableFeatures(merged) <- SelectIntegrationFeatures(
#     #   object.list = sn_list_sct,
#     #   nfeatures   = 3000
#     # )
#     # 
#     # rm(sn_list_sct); gc()
#     # 
#     # merged
#     
#     n_pcs <- PARAMS$seurat$n_pcs
#     obj   <- seurat_
#     
#     obj <- JoinLayers(obj)   # ← ADD THIS before RunPCA
#     gc()
#     
#     obj <- RunPCA(obj, npcs = n_pcs, verbose = FALSE)
#     obj <- harmony::RunHarmony(obj, group.by.vars = "sample_id",
#                                reduction = "pca", verbose = FALSE)
#     obj <- RunUMAP(obj,       reduction = "harmony", dims = 1:n_pcs, verbose = FALSE)
#     obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:n_pcs, verbose = FALSE)
#     obj <- FindClusters(obj,  resolution = 0.5,      verbose = FALSE)
#     
#     obj
#   }
# )
# 
# if (exists("sn_list_sct")) { rm(sn_list_sct); gc() }
# 
# seurat_sn <- load_or_compute(
#   path = here("results", "GSE299005", "tables", "seurat_sn_merged.rds"),
#   expr = {
#     log_msg("Merging ", length(sn_list_filtered), " snRNA samples ...")
#     obj_list <- unname(sn_list_filtered)
#     merged <- merge(obj_list[[1]],
#                     y           = obj_list[-1],
#                     add.cell.ids = NULL,
#                     project      = "GSE299005")
#     # SCTransform v2
#     SCTransform(merged, vst.flavor = "v2",
#                 vars.to.regress = "percent.mt", verbose = FALSE)
#   }
# )

# ── 4. Dimensionality reduction and clustering ────────────────────────────────

seurat_sn <- load_or_compute(
  path = here("results", "GSE299005", "tables", "seurat_sn_clustered.rds"),
  force = F,
  expr = {
    n_pcs <- PARAMS$seurat$n_pcs
    obj   <- seurat_sn
    gc()
    # obj <- JoinLayers(obj)
    obj <- RunPCA(obj, npcs = n_pcs, verbose = FALSE)
    obj <- harmony::RunHarmony(
      obj,
      group.by.vars = "sample_id",
      reduction.use = "pca",        # ← was "reduction", now "reduction.use"
      verbose       = FALSE
    )
    obj <- RunUMAP(obj, reduction = "harmony", dims = 1:n_pcs, verbose = FALSE)
    obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:n_pcs, verbose = FALSE)
    obj <- FindClusters(obj, resolution = 0.5,      verbose = FALSE)
    
    obj
    
    # n_pcs <- PARAMS$seurat$n_pcs
    # obj   <- seurat_sn          # ← now seurat_sn exists from Block 1
    # 
    # # obj <- JoinLayers(obj)
    # gc()
    # 
    # obj <- RunPCA(obj,        npcs = n_pcs,        verbose = FALSE)
    # obj <- harmony::RunHarmony(obj, group.by.vars = "sample_id",
    #                            reduction = "pca",  verbose = FALSE)
    # obj <- RunUMAP(obj,       reduction = "harmony", dims = 1:n_pcs, verbose = FALSE)
    # obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:n_pcs, verbose = FALSE)
    # obj <- FindClusters(obj,  resolution = 0.5,      verbose = FALSE)
    # obj
  }
)

# seurat_sn <- load_or_compute(
#   path = here("results", "GSE299005", "tables", "seurat_sn_clustered.rds"),
#   expr = {
#     n_pcs <- PARAMS$seurat$n_pcs
#     obj <- RunPCA(seurat_sn, npcs = n_pcs, verbose = FALSE)
#     # Harmony batch correction across the two DGE samples
#     obj <- harmony::RunHarmony(obj, group.by.vars = "sample_id",
#                                reduction = "pca", verbose = FALSE)
#     obj <- RunUMAP(obj, reduction = "harmony", dims = 1:n_pcs, verbose = FALSE)
#     obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:n_pcs, verbose = FALSE)
#     obj <- FindClusters(obj, resolution = 0.5, verbose = FALSE)
#     obj
#   }
# )

# ── 5. Cell type annotation ───────────────────────────────────────────────────
# If the original metadata included cluster labels, use them; otherwise use SingleR.

# if ("cell_type" %in% colnames(seurat_sn@meta.data)) {
#   log_msg("Using cell_type labels from original metadata.")
# } else {
#   log_msg("No pre-existing labels; running SingleR ...")
#   ref <- celldex::MouseRNAseqData()
#   
#   sce_pseudo <- AggregateExpression(
#     seurat_sn,
#     assays        = "SCT",
#     return.seurat = FALSE,
#     group.by      = "seurat_clusters"
#   )$SCT
#   
#   singleR_res <- SingleR::SingleR(
#     test   = sce_pseudo,
#     ref    = ref,
#     labels = ref$label.main
#   )
#   
#   # Diagnostic — confirm rownames match cluster IDs
#   log_msg("SingleR cluster IDs: ", paste(rownames(singleR_res), collapse = ", "))
#   log_msg("seurat_clusters levels: ", 
#           paste(levels(seurat_sn$seurat_clusters), collapse = ", "))
#   
#   cluster_map <- singleR_res$labels
#   names(cluster_map) <- rownames(singleR_res)
#   
#   # unname() prevents Seurat's name-matching check from hitting NAs
#   seurat_sn$cell_type <- unname(
#     cluster_map[as.character(seurat_sn$seurat_clusters)])
#   
#   # sce  <- GetAssayData(seurat_sn, assay = "SCT", layer = "data")
#   # singleR_res <- SingleR::SingleR(test = sce, ref = ref,
#   #                                  labels = ref$label.main)
#   # seurat_sn$cell_type <- singleR_res$labels
# }

# lapply(c("dplyr", "Seurat", "HGNChelper", "openxlsx"), library, character.only = TRUE)
# 
# source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
# source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")
# 
# # ScType cell marker DB — use "Brain" tissue type
# gs_list <- gene_sets_prepare(
#   "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx",
#   cell_type = "Brain"
# )
# 
# # Score using SCT scaled data
# es.max <- sctype_score(
#   scRNAseqData = seurat_sn[["SCT"]]@scale.data,
#   scaled       = TRUE,
#   gs           = gs_list$gs_positive,
#   gs2          = gs_list$gs_negative   # negative markers reduce false positives
# )
# 
# # Assign best label per cluster
# sctype_results <- do.call("rbind",
#                           lapply(unique(seurat_sn$seurat_clusters), function(cl) {
#                             cells <- which(seurat_sn$seurat_clusters == cl)
#                             scores <- sort(rowSums(es.max[, cells, drop = FALSE]), decreasing = TRUE)
#                             data.frame(
#                               cluster   = cl,
#                               cell_type = names(scores)[1],
#                               score     = scores[1],
#                               ncells    = length(cells)
#                             )
#                           })
# )
# 
# sctype_results$cell_type[sctype_results$score / sctype_results$ncells < 0] <- "Unknown"
# 
# print(sctype_results)
# 
# # Attach to Seurat object
# label_map <- setNames(sctype_results$cell_type, sctype_results$cluster)
# seurat_sn$cell_type_sctype <- unname(label_map[as.character(seurat_sn$seurat_clusters)])
# table(seurat_sn$cell_type_sctype, useNA = "always")


singler_results <- load_or_compute(
  path  = here("results", "GSE299005", "tables", "singler_results.rds"),
  force = T,
  expr  = {
    log_msg("Loading SingleR reference: MouseRNAseqData ...")
    ref_brain <- celldex::MouseRNAseqData()
    
    log_msg("Running SingleR on SCT data layer (",
            ncol(seurat_sn), " nuclei) ...")
    sce_counts <- GetAssayData(seurat_sn, assay = "SCT", layer = "data")
    
    sce_broad <- SingleR(
      test   = sce_counts,
      ref    = ref_brain,
      labels = ref_brain$label.main
    )
    
    sce_fine <- SingleR(
      test   = sce_counts,
      ref    = ref_brain,
      labels = ref_brain$label.fine
    )
    
    
  }
)

seurat_sn$singler_label <- singler_results$labels
seurat_sn$singler_score <- apply(singler_results$scores, 1, max)

plotScoreHeatmap(sce_broad, max.labels = 26, 
                 clusters = seurat_sn$seurat_clusters, 
                 order.by = "clusters", show_colnames = F)

plotScoreHeatmap(sce_fine,
                 clusters = seurat_sn$seurat_clusters, 
                 order.by = "clusters", show_colnames = F)

# ref <- celldex::MouseRNAseqData()   # or ImmGenData() for immune subtypes
# 
# # Pseudo-bulk per cluster (faster than per-cell, avoids noise)
# pb <- AggregateExpression(seurat_sn, assays = "SCT",
#                           group.by = "seurat_clusters",
#                           return.seurat = FALSE)$SCT
# 
# sr_res <- SingleR(test   = pb,
#                   ref    = ref,
#                   labels = ref$label.main)
# 
# print(data.frame(cluster = rownames(sr_res), label = sr_res$labels,
#                  pruned  = sr_res$pruned.labels))


marker_sets <- list(
  Endothelial     = c("Flt1", "Slco1a4", "Egfl7", "Pltp", "Abcg2"),
  Astrocyte       = c("Aqp4", "Gfap", "Vim"),
  Oligodendrocyte = c("Mbp", "Cldn11", "Mal", "Ermn", "Opalin"),
  Microglia       = c("Tmem119", "P2ry12", "Aif1"),
  Neuron          = c("Syt1", "Snap25"),
  Pericyte        = c("Rgs5", "P2ry14", "Ptn", "Cox4i2", "Rgs4", "Nbl1")
  # "Epithelial" = c("Epcam", "Cldn6", "Cdh1", "Ttr", "Folr1", "Prlr"),
  # "Astrocytes" = c("Gfap", "Aldh1l1", "Aqp4", "S100b"),
  # "Oligodendrocytes" = c("Olig1", "Olig2", "Mbp", "Cnp"),
  # "Endothelial_cells" = c("Cldn5", "Flt1", "Esam", "Cdh5", "Pecam1"),
  # "Microglia" = c("Iba1", "Aif1", "Tmem119", "P2ry12", "Cx3cr1"),
  # "Neurons" = c("Slc17a6", "Gad1", "Snap25", "Thy1", "Map2"),
  # "Mural_cells" = c("Cspg4", "Pdgfra", "Acta2", "Tagln"),
  # "Ependymal_cells" = c("Foxj1", "Ccdc67", "Dnah5", "Pifo"),
  # "Fibroblasts" = c("Dcn", "Col1a1", "Fbn1"),
  # "Pericytes" = c("Abcc9", "Kcnj8", "Rgs5", "Cspg4")
)

marker_vln <- function(srt, marker_list, marker_now) {
  vln <- VlnPlot(srt, features = marker_list[[marker_now]], pt.size = 0.1) +
    patchwork::plot_annotation(title = marker_now,
                               theme = theme(title = element_text(size = 20)))
  return(vln)
}

vlnList_p1 <- lapply(names(marker_sets[1:6]), marker_vln, 
                     srt = seurat_sn, marker_list = marker_sets)

ga_p1 <- ggarrange(plotlist = vlnList_p1)
ga_p1

ggsave(ga_p1, filename = "Violin Plots of Manual Annotation.png", path = ".", width = 36, height = 22)

# vlnList_p2 <- lapply(names(cell_type_markers[6:10]), marker_vln, 
#                      srt = seurat_merged, marker_list = cell_type_markers)

# seurat_sn <- PrepSCTFindMarkers(seurat_sn, verbose = T)
avg_exp <- AverageExpression(
  seurat_sn,
  features = unlist(marker_sets),
  assays   = "SCT",
  group.by = "seurat_clusters",
  layer    = "data"
)$SCT

score_mat <- sapply(marker_sets, function(genes) {
  genes <- genes[genes %in% rownames(avg_exp)]
  colMeans(avg_exp[genes, , drop = FALSE])
})

score_scaled <- apply(score_mat, 2, function(x) x / max(x + 1e-9))
cluster_labels <- colnames(score_scaled)[apply(score_scaled, 1, which.max)]
names(cluster_labels) <- rownames(score_mat)

print(cluster_labels)

names(cluster_labels) <- sub("^g", "", names(cluster_labels))

seurat_sn$cell_type <- unname(
  cluster_labels[as.character(seurat_sn$seurat_clusters)]
)

table(seurat_sn$cell_type, useNA = "always")

marker_sets <- list(
  Endothelial     = c("Egfl7"),                          # Only truly specific one
  Astrocyte       = c("Aqp4", "Gfap", "Vim"),
  Oligodendrocyte = c("Mbp", "Cldn11", "Mal", "Ermn", "Opalin"),
  Microglia       = c("Tmem119", "P2ry12", "Aif1"),
  Neuron          = c("Satb2","Bcl11b","Sox5","Rorb",    # Nuclear-enriched TFs
                      "Nrgn","Camk2a","Grin1","Neurod6","Tbr1"),
  Pericyte        = c("Rgs5", "P2ry14", "Ptn", "Cox4i2", "Rgs4", "Nbl1")
)

avg_exp <- AverageExpression(
  seurat_sn,
  features = unlist(marker_sets),
  assays   = "SCT",
  group.by = "seurat_clusters",
  layer    = "data"
)$SCT

score_mat <- sapply(marker_sets, function(genes) {
  genes <- genes[genes %in% rownames(avg_exp)]
  colMeans(avg_exp[genes, , drop = FALSE])
})

# Row-wise z-score
score_z <- t(apply(score_mat, 1, function(x) {
  sd_x <- sd(x)
  if (sd_x == 0) return(rep(0, length(x)))
  (x - mean(x)) / sd_x
}))
colnames(score_z) <- colnames(score_mat)

cluster_labels <- colnames(score_z)[apply(score_z, 1, which.max)]
names(cluster_labels) <- sub("^g", "", rownames(score_z))

print(cluster_labels)
print(round(score_mat, 4))   # Keep this to verify no single gene dominates again

seurat_sn$cell_type <- unname(
  cluster_labels[as.character(seurat_sn$seurat_clusters)]
)
table(seurat_sn$cell_type, useNA = "always")

# Standardise labels
# seurat_sn$cell_type <- dplyr::recode(
#   seurat_sn$cell_type,
#   "Endothelial cells" = "Endothelial",
#   "Pericytes"         = "Pericyte",
#   "Astrocytes"        = "Astrocyte",
#   "Oligodendrocytes"  = "Oligodendrocyte",
#   "Microglia"         = "Microglia",
#   "Neurons"           = "Neuron"
# )

seurat_sn$cell_type <- dplyr::case_match(
  seurat_sn$cell_type,
  "Endothelial cells" ~ "Endothelial",
  "Pericytes"         ~ "Pericyte",
  "Astrocytes"        ~ "Astrocyte",
  "Oligodendrocytes"  ~ "Oligodendrocyte",
  "Microglia"         ~ "Microglia",
  "Neurons"           ~ "Neuron",
  .default = seurat_sn$cell_type   # keep unmatched labels as-is
)

# ── 6. Save annotated object and BBB subset ───────────────────────────────────

saveRDS(seurat_sn,
        here("results", "GSE299005", "tables", "seurat_sn_annotated.rds"))

seurat_sn_bbb <- subset_bbb_cells(seurat_sn)
seurat_sn_bbb$dataset   <- "GSE299005"
seurat_sn_bbb$data_type <- "snRNA"

saveRDS(seurat_sn_bbb,
        here("results", "GSE299005", "tables", "seurat_sn_bbb.rds"))

log_msg("snRNA BBB cells: ", ncol(seurat_sn_bbb))

# UMAP figure
p_sn_umap <- DimPlot(seurat_sn, group.by = "cell_type",
                      label = TRUE, repel = TRUE,
                      cols = PALETTE_CELL_TYPES) +
  theme_publication() +
  labs(title = "GSE299005 — snRNA-seq (FPI male mice)",
       subtitle = "Cell type annotation")

save_figure(p_sn_umap, "01_umap_sn_annotated",
            dirs$figures, width = 9, height = 7)

log_section("GSE299005 | Step 1 | COMPLETE")
message("Next: run analysis/GSE160651/01_load_nanostring.R")
