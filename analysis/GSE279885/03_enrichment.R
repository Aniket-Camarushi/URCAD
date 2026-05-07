# ==============================================================================
# analysis/GSE279885/03_enrichment.R
#
# Purpose: Pathway and gene set enrichment analysis on GSE279885 DEGs.
#
# WHY TWO APPROACHES (ORA + GSEA):
#   ORA (Over-Representation Analysis): asks "are DEGs enriched in pathway X?"
#     - Uses only the significant DEG list as input.
#     - Fast and easy to interpret; but misses genes just below the threshold.
#   GSEA (Gene Set Enrichment Analysis): uses the FULL ranked gene list.
#     - Ranks all genes by shrunk log2FC.
#     - Detects subtle, coordinated shifts across a pathway even when no
#       single gene clears the significance threshold.
#     - The preferred method when you want to report on a specific pathway
#       (e.g. estrogen signaling, BBB integrity).
#
# Gene set sources:
#   - GO Biological Process (org.Mm.eg.db): broad functional annotation
#   - KEGG: metabolism and signaling pathways
#   - MSigDB Hallmarks (msigdbr): curated, low-redundancy gene sets
#     Specifically: HALLMARK_ESTROGEN_RESPONSE_EARLY/LATE
#
# Outputs saved to results/GSE279885/tables/ and figures/:
#   - enrichment_ora_go.rds / .csv
#   - enrichment_gsea_hallmark.rds / .csv
#   - 04_ora_dotplot.pdf
#   - 05_gsea_estrogen.pdf
# ==============================================================================

suppressPackageStartupMessages({
  suppressWarnings({
    library(here)
    library(dplyr)
    library(readr)
    library(ggplot2)
    library(clusterProfiler)
    library(org.Mm.eg.db)
    library(msigdbr)
    library(enrichplot)    
  })
})

source(here("analysis", "helpers", "utils.R"))
source(here("analysis", "helpers", "themes.R"))

log_section("GSE279885 | Step 3 | Enrichment Analysis")

dirs <- create_result_dirs("GSE279885")

de_results <- read_csv(
  here("results", "GSE279885", "tables", "de_results_required.csv"),
  show_col_types = FALSE
)

# ── [OPTIONAL DATASET — UNCOMMENT TO INCLUDE] ────────────────────────────────
# de_results <- read_csv(
#   here("results", "GSE279885", "tables", "de_results_combined.csv"),
#   show_col_types = FALSE
# )
# ─────────────────────────────────────────────────────────────────────────────

# ── 1. Prepare gene lists ─────────────────────────────────────────────────────

# ORA input: significant upregulated genes as Entrez IDs
sig_up <- de_results %>%
  filter(direction == "up") %>%
  pull(gene)

entrez_up <- symbol_to_entrez(sig_up)
entrez_up <- entrez_up[!is.na(entrez_up)]

# GSEA input: ALL genes ranked by shrunk log2FC (descending)
# WHY shrunk LFC not raw: apeglm shrinkage gives more reliable rankings for
# low-count genes, which would otherwise dominate the ranked list.
ranked_list <- de_results %>%
  filter(!is.na(log2FoldChange)) %>%
  arrange(desc(log2FoldChange)) %>%
  pull(log2FoldChange, name = gene)

# Convert gene symbols to Entrez for clusterProfiler
ranked_entrez <- symbol_to_entrez(names(ranked_list))
ranked_list_entrez <- ranked_list[!is.na(ranked_entrez)]
names(ranked_list_entrez) <- ranked_entrez[!is.na(ranked_entrez)]
ranked_list_entrez <- ranked_list_entrez[!duplicated(names(ranked_list_entrez))]
ranked_list_entrez <- sort(ranked_list_entrez, decreasing = TRUE)

# Universe = all tested genes (for ORA background)
universe_entrez <- symbol_to_entrez(de_results$gene)
universe_entrez <- universe_entrez[!is.na(universe_entrez)]

log_msg("ORA input (sig. up): ",    length(entrez_up))
log_msg("GSEA ranked list size: ",  length(ranked_list_entrez))
log_msg("Background universe: ",    length(universe_entrez))

# ── 2. ORA — GO Biological Process ───────────────────────────────────────────

ora_go <- load_or_compute(
  path = here("results", "GSE279885", "tables", "enrichment_ora_go.rds"),
  expr = {
    enrichGO(
      gene          = entrez_up,
      universe      = universe_entrez,
      OrgDb         = org.Mm.eg.db,
      ont           = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff  = PARAMS$enrichment$pvalue_cutoff,
      qvalueCutoff  = PARAMS$enrichment$qvalue_cutoff,
      minGSSize     = PARAMS$enrichment$min_gene_set_size,
      maxGSSize     = PARAMS$enrichment$max_gene_set_size,
      readable      = TRUE   # convert Entrez IDs back to symbols in results
    )
  }
)

if (nrow(as.data.frame(ora_go)) > 0) {
  write_csv(as.data.frame(ora_go),
            here("results", "GSE279885", "tables", "enrichment_ora_go.csv"))
  log_msg("ORA GO terms found: ", nrow(as.data.frame(ora_go)))

  p_ora <- dotplot(ora_go, showCategory = 20, font.size = 9) +
    labs(title = "ORA — GO Biological Process",
         subtitle = "Upregulated in VCD vs Vehicle (GSE279885)") +
    theme_publication()
  save_figure(p_ora, "04_ora_go_dotplot", dirs$figures, width = 8, height = 8)
} else {
  log_msg("ORA: no significant GO terms found.")
}

# ── 3. GSEA — MSigDB Hallmark gene sets ──────────────────────────────────────
#
# WHY Hallmarks: the 50 Hallmark gene sets are manually curated, coherent
# biological states with minimal redundancy. This makes them ideal for
# hypothesis-driven analysis. We specifically look for:
#   HALLMARK_ESTROGEN_RESPONSE_EARLY
#   HALLMARK_ESTROGEN_RESPONSE_LATE
#   HALLMARK_ANDROGEN_RESPONSE    (sex hormone cross-talk)
#   HALLMARK_INFLAMMATORY_RESPONSE

hallmarks_mm <- msigdbr(species = "Mus musculus", collection = "H") %>%
  dplyr::select(gs_name, ncbi_gene) %>%
  dplyr::rename(entrez_gene = ncbi_gene) %>%
  dplyr::mutate(entrez_gene = as.character(entrez_gene))

gsea_hallmark <- load_or_compute(
  path = here("results", "GSE279885", "tables", "enrichment_gsea_hallmark.rds"),
  expr = {
    GSEA(
      geneList      = ranked_list_entrez,
      TERM2GENE     = hallmarks_mm,
      pAdjustMethod = "BH",
      pvalueCutoff  = 1,     # keep all for visualisation; filter downstream
      nPermSimple   = PARAMS$enrichment$n_permutations,
      seed          = 42
    )
  }
)

write_csv(as.data.frame(gsea_hallmark),
          here("results", "GSE279885", "tables", "enrichment_gsea_hallmark.csv"))

# Highlight estrogen-relevant gene sets
estrogen_sets <- grep("ESTROGEN|ANDROGEN|INFLAMMATORY",
                      gsea_hallmark@result$ID,
                      value = TRUE)

if (length(estrogen_sets) > 0) {
  log_msg("Estrogen/inflammatory Hallmark sets:")
  gsea_hallmark@result %>%
    filter(ID %in% estrogen_sets) %>%
    dplyr::select(ID, NES, pvalue, p.adjust) %>%
    print()

  # GSEA enrichment plot for estrogen response
  for (set_name in estrogen_sets[1:min(3, length(estrogen_sets))]) {
    # safe_name <- gsub("[^A-Za-z0-9_]", "_", set_name)
    # p_gsea <- gseaplot2(gsea_hallmark, geneSetID = set_name,
    #                      title = set_name, color = "#D6604D", base_size = 11)
    # p_gsea <- p_gsea &
    #   theme(
    #     panel.background  = element_rect(fill = "transparent", colour = NA),
    #     plot.background   = element_rect(fill = "transparent", colour = NA),
    #     legend.background = element_rect(fill = "transparent", colour = NA),
    #     legend.key        = element_rect(fill = "transparent", colour = NA)
    #   )
    
    # p_gsea <- p_gsea +
    #   theme(
    #     panel.background = element_rect(fill = NA, colour = NA),
    #     plot.background  = element_rect(fill = NA, colour = NA)
    #   )
    # 
    # png_path <- file.path(dirs$figures, paste0("05_tmp_gsea_", safe_name, ".png"))
    # 
    # ggplot2::ggsave(
    #   filename = png_path,
    #   plot     = p_gsea,
    #   width    = 8,
    #   height   = 5,
    #   dpi      = 300,
    #   bg       = "transparent",
    #   device   = ragg::agg_png
    # )
    
    # save_figure(p_gsea,
    #             paste0("05_gsea_", safe_name),
    #             dirs$figures, width = 8, height = 5, bg = "transparent")
    
    
    # library(grid)
    # 
    # p_gsea <- gseaplot2(
    #   gsea_hallmark,
    #   geneSetID = set_name,
    #   title = set_name,
    #   color = "#D6604D",
    #   base_size = 11
    # )
    # 
    # safe_name <- gsub("[^A-Za-z0-9_]", "_", set_name)
    png_path <- file.path(dirs$figures, paste0("05_tmp_gsea_", safe_name, ".png"))
    # 
    # png(
    #   filename = png_path,
    #   width = 8,
    #   height = 5,
    #   units = "in",
    #   res = 300,
    #   bg = "transparent"
    # )
    # 
    # grid::grid.draw(p_gsea)
    # 
    # dev.off()
    
    safe_name <- gsub("[^A-Za-z0-9_]", "_", set_name)
    
    p_gsea <- gseaplot2(
      gsea_hallmark,
      geneSetID = set_name,
      title     = set_name,
      color     = "#D6604D",
      base_size = 11
    )
    
    # Apply transparent theme to all patchwork panels with &
    p_gsea <- p_gsea &
      theme(
        panel.background  = element_rect(fill = "transparent", colour = NA),
        plot.background   = element_rect(fill = "transparent", colour = NA),
        legend.background = element_rect(fill = "transparent", colour = NA),
        legend.key        = element_rect(fill = "transparent", colour = NA)
      )
    
    save_figure(
      p_gsea,
      paste0("05_gsea_", safe_name),
      dirs$figures,
      width  = 8,
      height = 5
    )
    
    png(filename = png_path, width = 8, height = 5,
        units = "in", res = 300, bg = "transparent")
    print(p_gsea)
    dev.off() 
    
  }
}

# ── 4. GSEA — KEGG pathways ──────────────────────────────────────────────────

kegg_mm <- msigdbr(species = "Mus musculus", collection = "C2", subcategory = "CP:KEGG_LEGACY") %>%
  dplyr::select(gs_name, ncbi_gene) %>%
  dplyr::rename(entrez_gene = ncbi_gene) %>%
  dplyr::mutate(entrez_gene = as.character(entrez_gene))

gsea_kegg <- load_or_compute(
  path = here("results", "GSE279885", "tables", "enrichment_gsea_kegg.rds"),
  expr = {
    GSEA(
      geneList      = ranked_list_entrez,
      TERM2GENE     = kegg_mm,
      pAdjustMethod = "BH",
      pvalueCutoff  = 1,
      nPermSimple   = PARAMS$enrichment$n_permutations,
      seed          = 42
    )
  }
)

sig_kegg <- as.data.frame(gsea_kegg) %>%
  filter(p.adjust < PARAMS$enrichment$pvalue_cutoff) %>%
  arrange(p.adjust)

write_csv(sig_kegg,
          here("results", "GSE279885", "tables", "enrichment_gsea_kegg_sig.csv"))
log_msg("Significant KEGG pathways: ", nrow(sig_kegg))

if (nrow(sig_kegg) > 0) {
  p_kegg <- dotplot(gsea_kegg, showCategory = 20, font.size = 9) +
    labs(title = "GSEA — KEGG Pathways",
         subtitle = "Ranked by shrunk log2FC: VCD vs Vehicle") +
    theme_publication()
  save_figure(p_kegg, "06_gsea_kegg_dotplot", dirs$figures, width = 9, height = 8)
}

log_section("GSE279885 | Step 3 | COMPLETE")
message("Next: run analysis/GSE158960/01_load_qc.R")
