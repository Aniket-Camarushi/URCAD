# ==============================================================================
# analysis/integration/02_cross_dataset_pathways.R
#
# Purpose: Identify estrogen-BBB pathway signatures conserved across all four
#          datasets and perform integrated ligand-receptor analysis.
#
# SCIENTIFIC GOAL:
#   We want to identify genes and pathways that are consistently associated
#   with estrogen signaling and BBB regulation across:
#     - Female VCD mice (bulk RNA-seq, estrogen depletion)
#     - WT vs Esr1cKO neurons/BBB cells (scRNA-seq, receptor ablation)
#     - Male FPI mice (snRNA-seq, injury-induced BBB disruption)
#     - Male FPI mice (NanoString, injury panel)
#   The intersection of these datasets reveals the core estrogen-BBB
#   regulatory program that is robust across sex, model, and assay type.
#
# Analysis steps:
#   1. Venn diagram of DEG overlaps (requires cross_dataset_sig genes)
#   2. Shared pathway enrichment (ORA on cross-dataset DEGs)
#   3. NicheNet: which ligands are predicted to regulate BBB target genes?
# ==============================================================================

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(ggvenn)
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(patchwork)
  library(pheatmap)
})

source(here("analysis", "helpers", "utils.R"))
source(here("analysis", "helpers", "themes.R"))

log_section("Integration | Step 2 | Cross-Dataset Pathways")
dirs <- create_result_dirs("integration")

# ── Load DE results from all datasets ────────────────────────────────────────

de_bulk <- readr::read_csv(
  here("results", "GSE279885", "tables", "de_results_required.csv"),
  show_col_types = FALSE)

de_nano <- readr::read_csv(
  here("results", "GSE160651", "tables", "de_limma_results.csv"),
  show_col_types = FALSE)

master_table <- readr::read_csv(
  here("results", "integration", "tables", "results_master_table.csv"),
  show_col_types = FALSE)

# scRNA-seq BBB DE (from GSE158960)
de_scrna_path <- here("results", "GSE158960", "tables",
                       "de_bbb_comparisons.xlsx")
if (file.exists(de_scrna_path)) {
  de_scrna <- openxlsx::read.xlsx(de_scrna_path, sheet = "WT_vs_CKO_ctrl")
} else {
  de_scrna <- tibble(gene = character(), avg_log2FC = numeric(),
                     p_val_adj = numeric(), direction = character())
  warning("de_bbb_comparisons.xlsx not found. Run GSE158960/04_bbb_estrogen.R first.")
}

# ── 1. Venn diagram: DEG overlap across datasets ──────────────────────────────
#
# For each dataset, get the set of significant upregulated genes.
# The overlap is the cross-dataset conserved estrogen-BBB signature.

sig_sets <- list(
  "Bulk VCD\n(GSE279885)" = de_bulk %>%
    filter(direction == "up") %>% pull(gene),

  "NanoString FPI\n(GSE160651)" = de_nano %>%
    filter(direction == "up") %>% pull(gene),

  "scRNA BBB\n(GSE158960)" = de_scrna %>%
    filter(direction == "up") %>% pull(gene)
)
sig_sets <- Filter(function(x) length(x) > 0, sig_sets)

if (length(sig_sets) >= 2) {
  p_venn <- ggvenn(sig_sets,
                   fill_color = c("#4393C3", "#D6604D", "#74C476"),
                   stroke_size = 0.5,
                   set_name_size = 4) +
    labs(title    = "Upregulated DEG overlap across datasets",
         subtitle = "Intersection = cross-dataset estrogen-BBB signature") +
    theme_publication()

  save_figure(p_venn, "02_venn_up_DEGs",
              dirs$figures, width = 7, height = 6)
}

# ── 2. Pathway enrichment on cross-dataset conserved genes ───────────────────
#
# Find genes that are upregulated in at least 2 of the 3 datasets.

conserved_up <- Reduce(union, lapply(sig_sets, identity))  # all sig in any
# Genes in 2+ datasets
if (length(sig_sets) >= 2) {
  gene_counts <- table(unlist(sig_sets))
  conserved_up_shared <- names(gene_counts)[gene_counts >= 2]
} else {
  conserved_up_shared <- conserved_up
}
log_msg("Genes upregulated in 2+ datasets: ", length(conserved_up_shared))

if (length(conserved_up_shared) >= 5) {
  entrez_shared <- symbol_to_entrez(conserved_up_shared)
  entrez_shared <- entrez_shared[!is.na(entrez_shared)]

  ora_shared <- tryCatch(
    enrichGO(
      gene          = entrez_shared,
      OrgDb         = org.Mm.eg.db,
      ont           = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff  = 0.05,
      readable      = TRUE
    ),
    error = function(e) NULL
  )

  if (!is.null(ora_shared) && nrow(as.data.frame(ora_shared)) > 0) {
    write_csv(as.data.frame(ora_shared),
              here("results", "integration", "tables",
                   "conserved_ora_go.csv"))

    p_ora_shared <- dotplot(ora_shared, showCategory = 20) +
      theme_publication() +
      labs(title    = "GO enrichment — cross-dataset conserved upregulated genes",
           subtitle = sprintf("%d genes upregulated in 2+ datasets",
                              length(conserved_up_shared)))

    save_figure(p_ora_shared, "03_conserved_ora_go",
                dirs$figures, width = 9, height = 8)
  }
}

# ── 3. BBB gene heatmap across datasets ───────────────────────────────────────
#
# Show the fold change (or log2FC) of key BBB and estrogen genes across
# all datasets in one figure. This is the summary figure for a paper.

bbb_estrogen_genes <- c(
  # Estrogen receptors
  "Esr1", "Esr2", "Gper1", "Esrra",
  # BBB tight junctions
  "Cldn5", "Tjp1", "Tjp2", "Ocln",
  # BBB adhesion / transport
  "Cdh5", "Esam", "Abcb1a", "Slc2a1", "Mfsd2a",
  # Pericyte
  "Pdgfrb", "Angpt1", "Notch3",
  # Astrocyte
  "Aqp4", "Gja1",
  # Neuroinflammation
  "Il1b", "Tnf", "Icam1", "Vcam1", "Ccl2"
)

# Build a matrix: rows = genes, columns = datasets
build_lfc_matrix <- function(gene_set) {
  cols <- list(
    `Bulk VCD\n(GSE279885)` = de_bulk %>%
      dplyr::select(gene, lfc = log2FoldChange),
    `NanoString FPI\n(GSE160651)` = de_nano %>%
      dplyr::select(gene, lfc = logFC)
  )
  if (nrow(de_scrna) > 0) {
    cols[["scRNA BBB\n(GSE158960)"]] <- de_scrna %>%
      dplyr::select(gene, lfc = avg_log2FC)
  }

  # Merge all into a matrix
  mat_df <- purrr::reduce(
    lapply(names(cols), function(n) {
      cols[[n]] %>% dplyr::filter(gene %in% gene_set) %>%
        dplyr::rename(!!n := lfc)
    }),
    dplyr::full_join, by = "gene"
  ) %>%
    tibble::column_to_rownames("gene")

  # Keep only genes present in at least one dataset
  mat_df <- mat_df[rowSums(!is.na(mat_df)) > 0, , drop = FALSE]
  as.matrix(mat_df)
}

lfc_mat <- build_lfc_matrix(bbb_estrogen_genes)
lfc_mat[is.na(lfc_mat)] <- 0  # set NA to 0 (not tested in that dataset)

if (nrow(lfc_mat) > 0) {
  pdf(here("results", "integration", "figures",
           "04_bbb_estrogen_cross_dataset_heatmap.pdf"),
      width = 8, height = 10)
  pheatmap(
    lfc_mat,
    color         = colorRampPalette(c("#2166AC", "white", "#D73027"))(100),
    breaks        = seq(-2, 2, length.out = 101),
    cluster_cols  = FALSE,
    cluster_rows  = TRUE,
    show_rownames = TRUE,
    show_colnames = TRUE,
    fontsize_row  = 9,
    fontsize_col  = 10,
    main          = "BBB & Estrogen genes — log2FC across datasets\n(blue = down, red = up)",
    na_col        = "grey90",
    border_color  = NA
  )
  dev.off()
  log_msg("Saved: 04_bbb_estrogen_cross_dataset_heatmap.pdf")
}

# ── 4. Proteomics scaffold ────────────────────────────────────────────────────
#
# Proteomics data is not yet available. This section is scaffolded so that
# when mass spectrometry data (e.g. MaxQuant output) is added, it integrates
# naturally into the pipeline.
#
# To activate: place MaxQuant output (proteinGroups.txt or similar) in
#   data/raw/proteomics/
# Then uncomment and complete the block below.

# ── [PROTEOMICS — UNCOMMENT WHEN DATA IS AVAILABLE] ─────────────────────────
#
# suppressPackageStartupMessages({
#   library(limma)    # for differential abundance
#   library(MSnbase)  # for proteomics data containers (Bioconductor)
#   library(DEP)      # for DESeq2-style proteomics DE
# })
#
# prot_path <- here("data", "raw", "proteomics", "proteinGroups.txt")
# if (file.exists(prot_path)) {
#   # Load MaxQuant output
#   # Columns "LFQ intensity *" contain label-free quantification values
#   prot_raw <- read_tsv(prot_path)
#
#   # Filter: remove reverse hits, contaminants, only-identified-by-site
#   prot_clean <- prot_raw %>%
#     filter(!`Reverse` %in% "+",
#            !`Potential contaminant` %in% "+",
#            !`Only identified by site` %in% "+")
#
#   # Proceed with DEP / limma for differential abundance testing
#   # Match significant proteins to the mRNA DEG lists for concordance analysis
# }
# ─────────────────────────────────────────────────────────────────────────────

log_section("Integration | Step 2 | COMPLETE")
message("Analysis pipeline complete. Review figures in results/integration/figures/")
