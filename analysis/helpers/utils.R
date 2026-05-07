# ==============================================================================
# analysis/helpers/utils.R
#
# Shared utility functions for the URCAD project.
# Source this at the top of every analysis script:
#   source(here::here("analysis", "helpers", "utils.R"))
#
# WHY A HELPERS FILE:
#   Prevents copy-pasting the same RDS caching logic, directory creation,
#   and logging into every script. Changes in one place propagate everywhere.
# ==============================================================================

suppressPackageStartupMessages({
  library(here)
  library(yaml)
  library(dplyr)
})

# ── Project-wide configuration ───────────────────────────────────────────────

#' Load the pipeline params YAML once and cache it.
#' Returns a named list matching the structure of config/pipeline_params.yaml.
load_params <- function(path = here("config", "pipeline_params.yaml")) {
  if (!file.exists(path))
    stop("pipeline_params.yaml not found at: ", path)
  yaml::read_yaml(path)
}

PARAMS <- load_params()

# ── Logging ──────────────────────────────────────────────────────────────────

#' Print a timestamped section header to console.
log_section <- function(title) {
  bar <- strrep("═", 60)
  message("\n", bar)
  message("  ", title)
  message("  ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
  message(bar)
}

#' Print a timestamped message.
log_msg <- function(...) {
  message(format(Sys.time(), "[%H:%M:%S]"), " ", ...)
}

# ── Directory management ─────────────────────────────────────────────────────

#' Create a directory (and any parent directories) if it does not exist.
#' Silently returns the path so it can be used in pipelines.
ensure_dir <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
    log_msg("Created directory: ", path)
  }
  invisible(path)
}

#' Create the standard results subdirectories for a given dataset ID.
#' Returns a named list of paths.
create_result_dirs <- function(dataset_id) {
  base   <- here("results", dataset_id)
  dirs   <- list(
    base    = base,
    figures = file.path(base, "figures"),
    tables  = file.path(base, "tables")
  )
  lapply(dirs, ensure_dir)
  invisible(dirs)
}

# ── RDS caching (load_or_compute) ────────────────────────────────────────────

#' Load an RDS file if it exists, otherwise compute with `expr`, save, and return.
#'
#' This is the single most important function for interactive development:
#' it makes each step idempotent so you never re-run a 2-hour Seurat integration
#' just because you re-opened RStudio.
#'
#' Usage:
#'   seurat_obj <- load_or_compute(
#'     path = here("results", "GSE158960", "seurat_merged.rds"),
#'     expr = {
#'       # ... expensive computation ...
#'       my_seurat_obj
#'     }
#'   )
#'
#' @param path  Character. Full path to the RDS file.
#' @param expr  Expression to evaluate if the cache is missing.
#' @param force Logical. If TRUE, always recompute (ignores cache).
load_or_compute <- function(path, expr, force = FALSE) {
  if (!force && file.exists(path)) {
    log_msg("Loading cached: ", basename(path))
    return(readRDS(path))
  }
  log_msg("Computing: ", basename(path))
  result <- expr
  ensure_dir(dirname(path))
  saveRDS(result, path)
  log_msg("Saved: ", path)
  result
}

# ── Gene annotation helpers ──────────────────────────────────────────────────

#' Convert Ensembl IDs to gene symbols using org.Mm.eg.db.
#' Returns a data.frame with columns ensembl_id and gene_symbol.
#' Drops unmapped IDs silently (set drop = FALSE to keep NAs).
ensembl_to_symbol <- function(ensembl_ids, drop = TRUE) {
  if (!requireNamespace("AnnotationDbi", quietly = TRUE) ||
      !requireNamespace("org.Mm.eg.db", quietly = TRUE))
    stop("Install AnnotationDbi and org.Mm.eg.db")

  res <- AnnotationDbi::select(
    org.Mm.eg.db::org.Mm.eg.db,
    keys    = ensembl_ids,
    columns = c("ENSEMBL", "SYMBOL"),
    keytype = "ENSEMBL"
  )
  names(res) <- c("ensembl_id", "gene_symbol")
  if (drop) res <- res[!is.na(res$gene_symbol), ]
  res
}

#' Convert gene symbols to Entrez IDs for clusterProfiler.
symbol_to_entrez <- function(symbols) {
  if (!requireNamespace("AnnotationDbi", quietly = TRUE) ||
      !requireNamespace("org.Mm.eg.db", quietly = TRUE))
    stop("Install AnnotationDbi and org.Mm.eg.db")

  AnnotationDbi::mapIds(
    org.Mm.eg.db::org.Mm.eg.db,
    keys    = symbols,
    column  = "ENTREZID",
    keytype = "SYMBOL",
    multiVals = "first"
  )
}

# ── Sample metadata loaders ──────────────────────────────────────────────────

#' Load GSE279885 sample metadata.
#' @param cohort One of "required", "optional", or "both".
load_gse279885_metadata <- function(cohort = "required") {
  req_path <- here("config", "samples_GSE279885_required.csv")
  opt_path <- here("config", "samples_GSE279885_optional.csv")

  switch(cohort,
    required = readr::read_csv(req_path, show_col_types = FALSE),
    optional = readr::read_csv(opt_path, show_col_types = FALSE),
    both     = dplyr::bind_rows(
                 readr::read_csv(req_path, show_col_types = FALSE),
                 readr::read_csv(opt_path, show_col_types = FALSE)
               ),
    stop("cohort must be 'required', 'optional', or 'both'. Got: ", cohort)
  )
}

# ── DESeq2 result helpers ────────────────────────────────────────────────────

#' Annotate a DESeq2 results data frame with gene symbols and significance flag.
#' @param res   data.frame from as.data.frame(DESeq2::results(...))
#' @param alpha Adjusted p-value cutoff (default: 0.05)
#' @param lfc   |log2FC| cutoff (default: 1.0)
annotate_deseq2_results <- function(res, alpha = 0.05, lfc = 1.0) {
  res$gene        <- rownames(res)
  res$significant <- !is.na(res$padj) &
                     res$padj < alpha &
                     abs(res$log2FoldChange) >= lfc
  res$direction   <- dplyr::case_when(
    res$significant & res$log2FoldChange > 0 ~ "up",
    res$significant & res$log2FoldChange < 0 ~ "down",
    TRUE                                     ~ "ns"
  )
  res[order(res$padj, na.last = TRUE), ]
}

# ── Seurat helpers ───────────────────────────────────────────────────────────

#' Add a dataset-level metadata column to a Seurat object.
#' Useful before merging sc and sn objects so you can always trace each cell.
add_dataset_label <- function(seurat_obj, dataset_id) {
  seurat_obj$dataset <- dataset_id
  seurat_obj
}

#' Extract BBB-relevant cell clusters from a Seurat object.
#' Filters by cell_type annotation column.
#' @param seurat_obj  Seurat object with `cell_type` metadata column.
#' @param cell_types  Character vector of cell types to keep.
subset_bbb_cells <- function(seurat_obj,
                             cell_types = c("Endothelial",
                                            "Pericyte",
                                            "Astrocyte")) {
  if (!"cell_type" %in% colnames(seurat_obj@meta.data))
    stop("'cell_type' column not found. Run annotation first.")
  
  seurat_obj@graphs <- list()

  present <- intersect(cell_types, unique(seurat_obj$cell_type))
  if (length(present) == 0)
    warning("None of the requested cell types found. Check cell_type labels.")

  subset(seurat_obj, subset = cell_type %in% present)
}

message("[utils.R] Loaded URCAD helper functions.")
