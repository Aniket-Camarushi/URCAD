# ==============================================================================
# analysis/00_install_packages.R
#
# Install all R packages required across the URCAD project.
#
# Run this ONCE after cloning the repository or when setting up renv.
# If you are using renv (recommended), run renv::restore() instead —
# that will install the exact versions recorded in renv.lock.
#
# When renv is NOT active (e.g. clean Docker build or first-time setup):
#   source("analysis/00_install_packages.R")
#
# WHY THIS FILE EXISTS:
#   renv::restore() handles exact reproducibility, but the first time you
#   set up the project you need to install and snapshot packages. This file
#   serves as the authoritative list of what the project needs so you or a
#   collaborator can recreate the environment from scratch.
# ==============================================================================

# ── Helper: install if missing ──────────────────────────────────────────────

install_if_missing <- function(pkg, source = "CRAN", repo = NULL) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message(sprintf("[install] %s (%s)", pkg, source))
    if (source == "CRAN") {
      install.packages(pkg, repos = "https://cloud.r-project.org", quiet = TRUE)
    } else if (source == "Bioconductor") {
      if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager", repos = "https://cloud.r-project.org")
      BiocManager::install(pkg, update = FALSE, ask = FALSE)
    } else if (source == "GitHub" && !is.null(repo)) {
      if (!requireNamespace("remotes", quietly = TRUE))
        install.packages("remotes", repos = "https://cloud.r-project.org")
      remotes::install_github(repo, upgrade = "never")
    }
  } else {
    message(sprintf("[ok]      %s", pkg))
  }
}

# ── 1. Core infrastructure ──────────────────────────────────────────────────

# renv: R environment reproducibility (already installed if using renv workflow)
install_if_missing("renv")

# here: find project root reliably; solves the "where is my working directory?" problem
install_if_missing("here")

# yaml: read pipeline_params.yaml and dataset_manifest.yaml from R
install_if_missing("yaml")

# ── 2. Shared data wrangling and visualization ──────────────────────────────

for (pkg in c(
  "dplyr",          # data manipulation (filter, mutate, group_by, summarise)
  "tidyr",          # reshaping data (pivot_wider, pivot_longer)
  "readr",          # fast CSV/TSV reading
  "tibble",         # modern data frames
  "stringr",        # string operations
  "forcats",        # factor reordering
  "purrr",          # functional programming (map, reduce)
  "ggplot2",        # base plotting (all figures)
  "ggrepel",        # non-overlapping labels on scatter plots
  "patchwork",      # combining multiple ggplot objects
  "pheatmap",       # heatmaps with hierarchical clustering
  "RColorBrewer",   # color palettes
  "viridis",        # perceptually uniform color scales (colorblind-friendly)
  "cowplot",        # publication-ready ggplot themes
  "scales",         # axis formatting (percent, comma)
  "writexl",        # export results to Excel
  "gt",             # beautiful publication-quality tables
  "openxlsx"        # Excel I/O
)) {
  install_if_missing(pkg)
}

# ── 3. GSE279885 — Bulk RNA-seq (tximport + DESeq2) ─────────────────────────

# tximport: import transcript-level Kallisto estimates into gene-level
# WHY: Kallisto outputs transcript abundances; DESeq2 expects gene-level counts.
#      tximport performs the aggregation and accounts for transcript-length bias.
install_if_missing("tximport",      "Bioconductor")
install_if_missing("tximeta",       "Bioconductor")  # metadata-aware tximport

# DESeq2: negative-binomial model for RNA-seq differential expression
# WHY: the industry standard for bulk RNA-seq DE; handles overdispersion properly
install_if_missing("DESeq2",        "Bioconductor")

# Annotation packages for mouse
install_if_missing("AnnotationDbi", "Bioconductor")  # gene ID conversion
install_if_missing("org.Mm.eg.db",  "Bioconductor")  # mouse gene annotations
install_if_missing("biomaRt",       "Bioconductor")  # Ensembl ID lookup

# Enrichment analysis
install_if_missing("clusterProfiler","Bioconductor") # ORA and GSEA
install_if_missing("fgsea",          "Bioconductor") # fast pre-ranked GSEA
install_if_missing("msigdbr",        "CRAN")         # MSigDB gene sets in R
install_if_missing("enrichplot",     "Bioconductor") # clusterProfiler visualizations
install_if_missing("DOSE",           "Bioconductor") # disease ontology enrichment

# ── 4. GSE158960 + GSE299005 — scRNA-seq / snRNA-seq (Seurat v5) ──────────

# Seurat v5: the primary framework for single-cell analysis
# WHY: Seurat v5 has native multi-assay layers, enabling sc+sn integration.
install_if_missing("Seurat")
install_if_missing("SeuratObject")

# Harmony: fast cross-dataset integration by iterative PCA correction
# WHY: faster and more stable than Seurat CCA for datasets of this size;
#      works well when cell-type composition differs between datasets
install_if_missing("harmony")

# BPCells: on-disk sparse matrix operations for very large Seurat objects
# WHY: GSE158960 (16 samples) can exceed Windows RAM; BPCells keeps matrices
#      on disk and operates on chunks, avoiding R memory limits.
# NOTE: requires HDF5 system library (libhdf5-dev)
install_if_missing("BPCells", "GitHub", "bnprks/BPCells")

# DoubletFinder: probabilistic doublet detection per sample
# WHY: scRNA-seq captures two cells in one droplet ~5-10% of the time;
#      doublets appear as intermediate clusters and confound cell type calling
install_if_missing("DoubletFinder", "GitHub", "chris-mcginnis-ucsf/DoubletFinder")

# SingleR: automated cell type annotation using reference transcriptomes
install_if_missing("SingleR",       "Bioconductor")
install_if_missing("celldex",       "Bioconductor")  # reference datasets for SingleR

# CellChat v2: ligand-receptor interaction analysis
# WHY: the primary tool for inferring cell-cell communication from scRNA-seq;
#      v2 added support for spatial data and improved the interaction database
install_if_missing("CellChat", "GitHub", "jinworks/CellChat")

# NMF: non-negative matrix factorization (used internally by CellChat)
install_if_missing("NMF")

# Matrix utilities
install_if_missing("Matrix")      # sparse matrix operations
install_if_missing("MatrixExtra") # additional sparse matrix utilities

# ── 5. GSE160651 — NanoString bulk RNA-seq ──────────────────────────────────

# NACHO: NanoString QC, normalization, and downstream analysis
# WHY: NACHO is the purpose-built Bioconductor package for NanoString nCounter
#      data; it implements all the manufacturer-recommended QC steps
install_if_missing("NACHO",  "Bioconductor")

# limma: linear models for microarray/RNA-seq; standard for NanoString DE
# WHY: NanoString data after normalization behaves like log-normal;
#      limma's voom/lmFit framework handles this appropriately
install_if_missing("limma",  "Bioconductor")
install_if_missing("edgeR",  "Bioconductor")  # for TMM normalization

# ── 6. Integration ───────────────────────────────────────────────────────────

# NicheNet: ligand-receptor inference with prior interaction knowledge
install_if_missing("nichenetr", "GitHub", "saeyslab/nichenetr")

# ggvenn / VennDiagram: for DEG overlap figures
install_if_missing("ggvenn")
install_if_missing("VennDiagram")

# ComplexHeatmap: highly customizable heatmaps for integration figures
install_if_missing("ComplexHeatmap", "Bioconductor")
install_if_missing("circlize",       "CRAN")  # required by ComplexHeatmap

# ── 7. Snapshot into renv ────────────────────────────────────────────────────

message("\n========================================")
message("All packages installed.")
message("Run renv::snapshot() to record versions.")
message("========================================\n")

# Uncomment to automatically snapshot after installation:
# renv::snapshot()
# rm(list = ls())
