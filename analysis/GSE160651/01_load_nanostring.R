# ==============================================================================
# analysis/GSE160651/01_load_nanostring.R
#
# Purpose: Load NanoString nCounter RCC files, run QC with NACHO,
#          normalize, and perform differential expression with limma.
#
# WHY NANOSTRING IS DIFFERENT FROM RNA-SEQ:
#   NanoString nCounter uses molecular barcodes to directly count mRNA
#   molecules without amplification. This eliminates PCR bias but produces
#   count data that is NOT directly comparable to RNA-seq counts:
#     - Counts are raw hybridization fluorescent signals, not sequence reads.
#     - Technical variation is captured by positive/negative control probes
#       included in every cartridge lane.
#     - Normalization uses housekeeping genes (not TMM or DESeq2 size factors).
#
# WHY NACHO:
#   NACHO implements the manufacturer-recommended QC pipeline:
#     1. Positive control scaling: corrects for lane-to-lane efficiency variation
#     2. Negative control subtraction: removes background hybridization signal
#     3. Housekeeping normalization: corrects for sample input amount
#   It also produces standardized QC plots (imaging QC, binding density,
#   positive control linearity) that reviewers expect.
#
# WHY LIMMA NOT DESEQ2:
#   After NACHO normalization, the data is approximately log-normal rather
#   than negative-binomial. limma's lmFit/eBayes framework is appropriate
#   here and is the standard for NanoString DE in the literature.
#
# Outputs:
#   - nacho_normalized.rds     NACHO normalized data object
#   - de_limma_results.csv     Full DE table
#   - de_limma_sig.csv         Significant DEGs
# ==============================================================================

suppressPackageStartupMessages({
  suppressWarnings({
    library(here)
    library(NACHO)
    library(limma)
    library(dplyr)
    library(readr)
    library(ggplot2)
    library(patchwork)
    library(pheatmap)
    library(ggrepel)
    library(DESeq2)
  })
})

source(here("analysis", "helpers", "utils.R"))
source(here("analysis", "helpers", "themes.R"))

log_section("GSE160651 | Step 1 | NanoString QC & DE")
dirs <- create_result_dirs("GSE160651")

rcc_dir <- here("data", "raw", "GSE160651")
rcc_files <- list.files(rcc_dir, pattern = "\\.RCC$", full.names = TRUE)

# ── Complete GEO metadata from GSE160651 ──────────────────────────────────────
# Vehicle-only samples extracted from GEO (PLX5622 samples excluded)
geo_meta <- tibble::tribble(
  ~gsm_id,        ~animal_id, ~treatment,  ~injury,    ~timepoint,
  "GSM4876458",   "684",      "Vehicle",   "TBI",      "30d",
  "GSM4876460",   "686",      "Vehicle",   "TBI",      "7d",
  "GSM4876464",   "691",      "Vehicle",   "Control",  "0d",
  "GSM4876467",   "695",      "Vehicle",   "Control",  "7d",
  "GSM4876469",   "697",      "Vehicle",   "TBI",      "1d",
  "GSM4876473",   "702",      "Vehicle",   "TBI",      "7d",
  "GSM4876475",   "704",      "Vehicle",   "TBI",      "1d",
  "GSM4876477",   "706",      "Vehicle",   "Control",  "0d",
  "GSM4876478",   "707",      "Vehicle",   "TBI",      "30d",
  "GSM4876481",   "710",      "Vehicle",   "TBI",      "7d",
  "GSM4876483",   "712",      "Vehicle",   "TBI",      "1d",
  "GSM4876485",   "714",      "Vehicle",   "Control",  "0d",
  "GSM4876487",   "716",      "Vehicle",   "TBI",      "30d",
  "GSM4876489",   "718",      "Vehicle",   "TBI",      "7d",
  "GSM4876491",   "720",      "Vehicle",   "TBI",      "1d",
  "GSM4876493",   "722",      "Vehicle",   "TBI",      "30d",
  "GSM4876495",   "724",      "Vehicle",   "TBI",      "7d",
  "GSM4876497",   "726",      "Vehicle",   "TBI",      "1d",
  "GSM4876499",   "728",      "Vehicle",   "TBI",      "30d",
  "GSM4876501",   "730",      "Vehicle",   "TBI",      "7d",
  "GSM4876503",   "732",      "Vehicle",   "Control",  "0d",
  "GSM4876506",   "735",      "Vehicle",   "TBI",      "7d",
  "GSM4876508",   "737",      "Vehicle",   "Control",  "0d",
  "GSM4876511",   "741",      "Vehicle",   "TBI",      "1d",
  "GSM4876512",   "742",      "Vehicle",   "Control",  "0d"
)

# Build composite group column
geo_meta <- geo_meta %>%
  mutate(
    group = case_when(
      injury == "Control"                  ~ "Vehicle_Control",
      injury == "TBI" & timepoint == "1d"  ~ "Vehicle_TBI_1d",
      injury == "TBI" & timepoint == "7d"  ~ "Vehicle_TBI_7d",
      injury == "TBI" & timepoint == "30d" ~ "Vehicle_TBI_30d"
    )
  )

# Match GSM IDs to actual RCC file paths on disk
rcc_files <- list.files(rcc_dir, pattern = "\\.RCC$", full.names = TRUE)

# RCC filenames start with the GSM ID: GSM4876458_...RCC
geo_meta <- geo_meta %>%
  rowwise() %>%
  mutate(
    RCC_FILE = {
      all_rcc <- list.files(rcc_dir, pattern = "\\.RCC$", full.names = TRUE)
      m <- rcc_files[startsWith(basename(rcc_files), gsm_id)]
      if (length(m) == 1) m else NA_character_
    }
  ) %>%
  ungroup()

# Extract cartridge (batch) from filename for downstream QC
geo_meta <- geo_meta %>%
  mutate(
    batch = case_when(
      grepl("1st", RCC_FILE) ~ "batch1",
      grepl("2nd", RCC_FILE) ~ "batch2",
      grepl("3rd", RCC_FILE) ~ "batch3",
      grepl("4th", RCC_FILE) ~ "batch4",
      grepl("5th", RCC_FILE) ~ "batch5",
      TRUE                   ~ "unknown"
    )
  )

# Verify all files found
missing <- geo_meta %>% filter(is.na(RCC_FILE))
if (nrow(missing) > 0) {
  warning("Could not match RCC files for: ",
          paste(missing$gsm_id, collapse = ", "))
}

message(nrow(geo_meta) - nrow(missing), " / ", nrow(geo_meta),
        " samples matched to RCC files")
print(table(geo_meta$group))

# Save
out_path <- here("data", "raw", "GSE160651", "GSE160651_design_manual.csv")
readr::write_csv(geo_meta, out_path)
message("Saved: ", out_path)

metadata <- geo_meta %>% filter(!is.na(RCC_FILE) & file.exists(RCC_FILE))
log_msg("Samples loaded: ", nrow(metadata))
log_msg("Groups: ", paste(sort(unique(metadata$group)), collapse = ", "))

if (nrow(metadata) == 0) {
  stop("No valid RCC files found. Check rcc_dir: ", rcc_dir)
}

# ── 2. NACHO QC and normalization ─────────────────────────────────────────────

nacho_raw <- load_or_compute(
  path = here("results", "GSE160651", "tables", "nacho_raw.rds"),
  force = T,
  expr = {
    log_msg("Running NACHO load_rcc ...")
    
    # NACHO needs the filename only, not the full path, in the id column
    ssheet <- metadata %>%
      dplyr::mutate(RCC_FILE = basename(RCC_FILE))
    
    load_rcc(
      data_directory      = rcc_dir,
      ssheet_csv          = ssheet,
      id_colname          = "RCC_FILE",
      housekeeping_genes  = NULL,       # use defaults from the panel
      housekeeping_predict = TRUE,      # let NACHO find stable HKGs
      normalisation_method = "GEO",    # geometric mean (manufacturer default)
      n_comp              = 5
    )
  }
)

hkg <- c("Aars", "Asb7", "Ccdc127", "Cnot10", "Csnk2a2", "Lars", "Mto1")

nacho_obj <- load_or_compute(
  path = here("results", "GSE160651", "tables", "nacho_normalized.rds"),
  force = T,
  expr = {
    log_msg("Running NACHO normalise ...")
    normalise(
      nacho_object         = nacho_raw,
      housekeeping_genes   = hkg,
      housekeeping_predict = FALSE,
      housekeeping_norm    = TRUE,
      normalisation_method = "GEO",
      remove_outliers      = TRUE
    )
  }
)

log_msg("HKGs used: ", paste(nacho_obj$housekeeping_genes, collapse = ", "))

log_msg("NACHO normalised. Genes: ",
        length(unique(nacho_obj$nacho$Name)),
        "  Samples: ", length(unique(nacho_obj$nacho$RCC_FILE)))


# ── 3. QC figures ────────────────────────────────────────────────────────────
#
# NACHO has built-in plot methods. We call them and save as PDFs.
# These are the figures reviewers will look for in a NanoString paper.

pdf(here("results", "GSE160651", "figures", "01_nacho_qc.pdf"),
    width = 12, height = 8)
autoplot(nacho_obj, x = "FoV") + theme_publication() + labs(title = "NACHO QC: Imaging")
dev.off()
log_msg("Saved: 01_nacho_qc.pdf")

qc_dir <- here("results", "GSE160651", "figures", "qc")
dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)

label_map <- setNames(
  metadata$group,
  basename(metadata$RCC_FILE)
)

qc_metrics <- c(
  "BD",           # Binding Density — must be 0.1–2.25 (FLEX/MAX)
  "FoV",          # Field of View / Imaging quality — must be > 75%
  "PCL",          # Positive Control Linearity — must be > 0.95
  "LoD",          # Limit of Detection
  "Positive",     # Positive control counts
  "Negative",     # Negative control counts (background)
  "Housekeeping", # Housekeeping gene expression
  "PN",           # Positive vs. Negative controls
  "ACBD",         # Average Counts vs. Binding Density
  "ACMC",         # Average Counts vs. Median Counts
  "PCA12",        # PC1 vs. PC2
  "PCAi",         # PCA scree plot
  "PFNF",         # Positive Factor vs. Negative Factor
  "HF",           # Housekeeping Factor
  "NORM"          # Normalisation Factor
)

pdf(here("results", "GSE160651", "figures", "01_nacho_qc.pdf"),
    width = 12, height = 8)
for (metric in qc_metrics) {
  p <- tryCatch({
    autoplot(
      object      = nacho_obj,
      x           = metric,
      colour      = "group",       
      size        = 1.8,
      show_legend = TRUE
    ) + 
    theme_publication(base_size = 10) +
    labs(title = paste("NACHO QC:", metric))
    },
    error = function(e) {
      message("Skipping metric '", metric, "': ", conditionMessage(e))
      NULL
    }
  )
  if (!is.null(p)) print(p)
  
  p <- tryCatch(
    p +
      scale_x_discrete(labels = label_map) +
      theme(axis.text.x  = element_text(angle = 45, hjust = 1, size = 7),
            axis.title.x = element_blank()),
    error = function(e) p
  )
  
  save_figure(p,
              name   = paste0("qc_", tolower(metric)),
              outdir = qc_dir,
              width  = 10,
              height = 6)
}

# dev.off()

log_msg("Saved: 01_nacho_qc.pdf (", length(qc_metrics), " QC panels)")


# ── 4. Prepare expression matrix for limma ───────────────────────────────────
#
# Extract the normalized log2 count matrix from NACHO.
# Rows = genes, Columns = samples.

# norm_counts <- nacho_obj$nacho %>%
#   dplyr::filter(CodeClass == "Endogenous") %>%   # genes only, not controls
#   dplyr::select(Name, RCC_FILE, Count_Norm) %>%
#   tidyr::pivot_wider(names_from  = RCC_FILE,
#                      values_from = Count_Norm) %>%
#   tibble::column_to_rownames("Name") %>%
#   as.matrix()
# 
# # log2-transform (NACHO normalised counts are on raw scale)
# norm_counts <- log2(norm_counts + 1)
# 
# log_msg("Normalised count matrix: ", nrow(norm_counts), " x ", ncol(norm_counts))
# 
# # Re-align metadata to matrix column order
# common_samples <- intersect(colnames(norm_counts), basename(metadata$RCC_FILE))
# norm_counts    <- norm_counts[, common_samples]
# metadata_match <- metadata %>%
#   dplyr::mutate(RCC_FILE_base = basename(RCC_FILE)) %>%
#   dplyr::filter(RCC_FILE_base %in% common_samples) %>%
#   dplyr::arrange(match(RCC_FILE_base, common_samples))



# Pull RAW counts (Count, not Count_Norm) from NACHO
raw_counts <- nacho_obj$nacho %>%
  dplyr::filter(CodeClass == "Endogenous") %>%
  dplyr::select(Name, RCC_FILE, Count) %>%         # RAW counts
  tidyr::pivot_wider(names_from  = RCC_FILE,
                     values_from = Count) %>%
  tibble::column_to_rownames("Name") %>%
  as.matrix() %>%
  round()                                           # DESeq2 requires integers

# Align metadata
common_samples <- intersect(colnames(raw_counts), basename(metadata$RCC_FILE))
raw_counts     <- raw_counts[, common_samples]
metadata_deseq <- metadata %>%
  mutate(RCC_FILE_base = basename(RCC_FILE)) %>%
  filter(RCC_FILE_base %in% common_samples) %>%
  arrange(match(RCC_FILE_base, common_samples)) %>%
  mutate(
    group = factor(group,
                   levels = c("Vehicle_Control", "Vehicle_TBI_1d",
                              "Vehicle_TBI_7d",  "Vehicle_TBI_30d")),
    batch = factor(batch)
  ) %>%
  tibble::column_to_rownames("RCC_FILE_base")

qc_samples <- unique(nacho_obj$nacho$RCC_FILE)

metadata_match <- metadata %>%
  dplyr::mutate(RCC_FILE_base = basename(RCC_FILE)) %>%
  dplyr::filter(RCC_FILE_base %in% qc_samples) %>%
  dplyr::arrange(match(RCC_FILE_base, qc_samples))

pivot_matrix <- function(value_col) {
  nacho_obj$nacho %>%
    dplyr::filter(CodeClass == "Endogenous") %>%
    dplyr::select(Name, RCC_FILE, all_of(value_col)) %>%
    tidyr::pivot_wider(names_from = RCC_FILE, values_from = all_of(value_col)) %>%
    tibble::column_to_rownames("Name") %>%
    as.matrix()
}

raw_counts  <- round(pivot_matrix("Count"))
norm_counts <- log2(pivot_matrix("Count_Norm") + 1)

raw_counts  <- raw_counts[,  metadata_match$RCC_FILE_base]
norm_counts <- norm_counts[, metadata_match$RCC_FILE_base]

log_msg("Raw count matrix:  ", nrow(raw_counts),  " x ", ncol(raw_counts))
log_msg("Norm count matrix: ", nrow(norm_counts), " x ", ncol(norm_counts))

# 
# 
# dds <- DESeqDataSetFromMatrix(
#   countData = raw_counts,
#   colData   = metadata_deseq,
#   design    = ~ batch + group          # batch covariate + group
# )
# 
# dds <- DESeq(dds)
# 
# # Extract each contrast — paper used p < 0.05 (NOT adjusted)
# contrasts_list <- list(
#   TBI_1d_vs_Control  = c("group", "Vehicle_TBI_1d",  "Vehicle_Control"),
#   TBI_7d_vs_Control  = c("group", "Vehicle_TBI_7d",  "Vehicle_Control"),
#   TBI_30d_vs_Control = c("group", "Vehicle_TBI_30d", "Vehicle_Control")
# )
# 
# de_results_deseq <- lapply(names(contrasts_list), function(cname) {
#   res <- results(dds,
#                  contrast      = contrasts_list[[cname]],
#                  independentFiltering = TRUE) %>%
#     as.data.frame() %>%
#     tibble::rownames_to_column("gene") %>%
#     mutate(
#       contrast    = cname,
#       # Paper threshold: p < 0.05 unadjusted
#       significant = pvalue < 0.05 & abs(log2FoldChange) >= 1,
#       direction   = case_when(
#         significant & log2FoldChange > 0 ~ "up",
#         significant & log2FoldChange < 0 ~ "down",
#         TRUE ~ "ns"
#       )
#     )
# }) %>% bind_rows()

# ── 5. Design matrix and contrasts ────────────────────────────────────────────
#
# Experimental design: male FPI mice vs sham controls (if present).
# Identify the group column — adjust if your metadata uses a different name.

# group_factor <- factor(metadata_match[[group_col]])
# log_msg("Groups: ", paste(levels(group_factor), collapse = " vs "))

# group_factor <- factor(
#   metadata_match$group,
#   levels = c("Vehicle_Control", "Vehicle_TBI_1d",
#              "Vehicle_TBI_7d",  "Vehicle_TBI_30d")
# )
# 
# batch_factor <- factor(metadata_match$batch)
# 
# log_msg("Groups (n per group):")
# print(table(group_factor))
# 
# log_msg("Batch (n per batch):")
# print(table(batch_factor))
# 
# design <- model.matrix(~ 0 + group_factor + batch_factor)
# # colnames(design) <- levels(group_factor)
# colnames(design) <- gsub("^group_factor", "", colnames(design))
# colnames(design) <- gsub("^batch_factor",  "batch_", colnames(design))
# 
# log_msg("Design: ", nrow(design), " samples x ", ncol(design), " terms")
# 
# # Contrast: FPI vs Sham (adjust to your actual group names)
# group_levels <- levels(group_factor)
# if (length(group_levels) >= 2) {
#   contrast_str <- paste0(group_levels[2], " - ", group_levels[1])
#   contrast_matrix <- makeContrasts(
#     contrasts = contrast_str,
#     levels    = design
#   )
# }
# 
# contrast_matrix <- makeContrasts(
#   TBI_1d_vs_Control  = Vehicle_TBI_1d  - Vehicle_Control,
#   TBI_7d_vs_Control  = Vehicle_TBI_7d  - Vehicle_Control,
#   TBI_30d_vs_Control = Vehicle_TBI_30d - Vehicle_Control,
#   levels = design
# )

col_data <- data.frame(
  group = factor(metadata_match$group,
                 levels = c("Vehicle_Control", "Vehicle_TBI_1d",
                            "Vehicle_TBI_7d",  "Vehicle_TBI_30d")),
  batch = factor(metadata_match$batch),
  row.names = metadata_match$RCC_FILE_base
)

log_msg("Groups (n per group):"); print(table(col_data$group))
log_msg("Batch (n per batch):");  print(table(col_data$batch))

dds <- DESeqDataSetFromMatrix(
  countData = raw_counts,
  colData   = col_data,
  design    = ~ batch + group
)
dds <- DESeq(dds)
log_msg("DESeq2 complete. Coefficients: ",
        paste(resultsNames(dds), collapse = ", "))


# ── 6. DESeq2 differential expression ─────────────────────────────────────────

# fit1 <- lmFit(norm_counts, design)
# fit2 <- contrasts.fit(fit1, contrast_matrix)
# fit2 <- eBayes(fit2)
# 
# de_results_ns <- topTable(fit2, number = Inf, adjust.method = "BH",
#                            sort.by = "P")
# de_results_ns$gene        <- rownames(de_results_ns)
# de_results_ns$significant <- de_results_ns$adj.P.Val < PARAMS$deseq2$alpha &
#                               abs(de_results_ns$logFC) >= PARAMS$deseq2$lfc_threshold
# de_results_ns$direction   <- dplyr::case_when(
#   de_results_ns$significant & de_results_ns$logFC > 0 ~ "up",
#   de_results_ns$significant & de_results_ns$logFC < 0 ~ "down",
#   TRUE ~ "ns"
# )
# 
# de_results_ns <- lapply(colnames(contrast_matrix), function(cname) {
#   res <- topTable(fit2, coef = cname, number = Inf,
#                   adjust.method = "BH", sort.by = "P")
#   res$gene       <- rownames(res)
#   res$contrast   <- cname
#   res$significant <- res$adj.P.Val < PARAMS$deseq2$alpha &
#                      abs(res$logFC) >= PARAMS$deseq2$lfc_threshold
#   res$direction  <- dplyr::case_when(
#     res$significant & res$logFC > 0 ~ "up",
#     res$significant & res$logFC < 0 ~ "down",
#     TRUE ~ "ns"
#   )
#   res
# }) %>% dplyr::bind_rows()
# 
# write_csv(de_results_ns,
#           here("results", "GSE160651", "tables", "de_limma_results.csv"))
# write_csv(de_results_ns %>% filter(significant),
#           here("results", "GSE160651", "tables", "de_limma_sig.csv"))
# 
# log_msg("Total gene-contrast tests: ", nrow(de_results_ns))
# log_msg("Significant DEGs per contrast:")
# de_results_ns %>%
#   filter(significant) %>%
#   dplyr::count(contrast, direction) %>%
#   print()

contrasts_list <- list(
  TBI_1d_vs_Control  = c("group", "Vehicle_TBI_1d",  "Vehicle_Control"),
  TBI_7d_vs_Control  = c("group", "Vehicle_TBI_7d",  "Vehicle_Control"),
  TBI_30d_vs_Control = c("group", "Vehicle_TBI_30d", "Vehicle_Control")
)

de_results_ns <- lapply(names(contrasts_list), function(cname) {
  results(dds, contrast = contrasts_list[[cname]],
          independentFiltering = TRUE) %>%
    as.data.frame() %>%
    tibble::rownames_to_column("gene") %>%
    dplyr::mutate(
      contrast    = cname,
      significant = !is.na(pvalue) & pvalue < 0.05 &
        abs(log2FoldChange) >= PARAMS$deseq2$lfc_threshold,
      direction   = dplyr::case_when(
        significant & log2FoldChange > 0 ~ "up",
        significant & log2FoldChange < 0 ~ "down",
        TRUE ~ "ns"
      )
    )
}) %>% dplyr::bind_rows()

write_csv(de_results_ns,
          here("results", "GSE160651", "tables", "de_deseq2_results_tmp.csv"))
write_csv(de_results_ns %>% filter(significant),
          here("results", "GSE160651", "tables", "de_deseq2_sig_tmp.csv"))

log_msg("Total gene-contrast tests: ", nrow(de_results_ns))
log_msg("Significant DEGs (p<0.05, |LFC|>=1) per contrast:")
de_results_ns %>%
  filter(significant) %>%
  dplyr::count(contrast, direction) %>%
  print()

# ── 7. Volcano plot ───────────────────────────────────────────────────────────

de_results_ns <- load_or_compute(
  path = here("results", "GSE160651", "tables", "de_deseq2_results_tmp.csv"),
  force = T,
  expr = read_csv(here("results", "GSE160651", "tables", "de_deseq2_results_tmp.csv"))
)

volcano_plots <- lapply(names(contrasts_list), function(cname) {
  res <- de_results_ns %>% dplyr::filter(contrast == cname)
  
  top_labels <- dplyr::bind_rows(
    res %>% filter(direction == "up")   %>% arrange(pvalue) %>% head(10),
    res %>% filter(direction == "down") %>% arrange(pvalue) %>% head(10)
  )
  
  # ggplot(res, aes(logFC, -log10(adj.P.Val), colour = direction)) +
  #   geom_point(alpha = 0.5, size = 1.2) +
  #   ggrepel::geom_text_repel(data = top_labels, aes(label = gene),
  #                            size = 2.8, max.overlaps = 20,
  #                            show.legend = FALSE) +
  #   geom_vline(xintercept = c(-PARAMS$deseq2$lfc_threshold,
  #                             PARAMS$deseq2$lfc_threshold),
  #              linetype = "dashed", colour = "grey50") +
  #   geom_hline(yintercept = -log10(PARAMS$deseq2$alpha),
  #              linetype = "dashed", colour = "grey50") +
  #   scale_color_manual(values = PALETTE_DE) +
  #   theme_publication(base_size = 10) +
  #   labs(
  #     title    = gsub("_", " ", cname),
  #     subtitle = paste0("Up: ", sum(res$direction == "up"),
  #                       "  Down: ", sum(res$direction == "down")),
  #     x        = "log2 fold change",
  #     y        = expression(-log[10](adjusted~p-value)),
  #     colour   = NULL
  #   )
  
  ggplot(res, aes(log2FoldChange, -log10(pvalue), colour = direction)) +
    geom_point(alpha = 0.5, size = 1.2) +
    ggrepel::geom_text_repel(data = top_labels, aes(label = gene),
                             size = 2.8, max.overlaps = 20,
                             show.legend = FALSE) +
    geom_vline(xintercept = c(-PARAMS$deseq2$lfc_threshold,
                              PARAMS$deseq2$lfc_threshold),
               linetype = "dashed", colour = "grey50") +
    geom_hline(yintercept = -log10(PARAMS$deseq2$alpha),          # PARAMS$deseq2$alpha = 0.05
               linetype = "dashed", colour = "grey50") +
    scale_color_manual(values = PALETTE_DE) +
    theme_publication(base_size = 10) +
    labs(
      title    = gsub("_", " ", cname),
      subtitle = paste0("Up: ", sum(res$direction == "up"),
                        "  Down: ", sum(res$direction == "down")),
      x        = "log2 fold change",
      y        = expression(-log[10](p-value)),
      colour   = NULL
    )
  
})

names(volcano_plots) <- names(contrasts_list)

p_volcanoes <- patchwork::wrap_plots(volcano_plots, ncol = 3) +
  patchwork::plot_annotation(
    title    = "GSE160651 - TBI vs Control (DESeq2)",
    subtitle = "NanoString Neuropathology Panel | p < 0.05, |LFC| ≥ 1"
  )

save_figure(p_volcanoes, "02_nanostring_volcanoes_tmp",
            dirs$figures, width = 21, height = 7)

# Quick check
de_results_ns %>%
  filter(contrast == "TBI_1d_vs_Control", significant) %>%
  arrange(pvalue) %>%
  select(gene, log2FoldChange, pvalue, direction) %>%
  head(20) %>%
  print()

# p_ns_volcano <- ggplot(de_results_ns,
#                         aes(logFC, -log10(adj.P.Val), colour = direction)) +
#   geom_point(alpha = 0.5, size = 1.2) +
#   ggrepel::geom_text_repel(data = top_ns_labels,
#                             aes(label = gene),
#                             size = 2.8, max.overlaps = 20,
#                             show.legend = FALSE) +
#   geom_vline(xintercept = c(-PARAMS$deseq2$lfc_threshold,
#                               PARAMS$deseq2$lfc_threshold),
#              linetype = "dashed", colour = "grey50") +
#   geom_hline(yintercept = -log10(PARAMS$deseq2$alpha),
#              linetype = "dashed", colour = "grey50") +
#   scale_color_manual(values = PALETTE_DE) +
#   theme_publication(base_size = 10) +
#   labs(
#     title    = sprintf("Volcano plot — NanoString FPI (GSE160651)"),
#     subtitle = sprintf("Comparison: %s", contrast_str),
#     x        = "log2 fold change",
#     y        = expression(-log[10](adjusted~p-value)),
#     colour   = NULL
#   )
# 
# save_figure(p_ns_volcano, "02_nanostring_volcano",
#             dirs$figures, width = 7, height = 6)

# ── 8. Save normalized matrix for integration ────────────────────────────────

saveRDS(list(
  norm_counts = norm_counts,
  metadata    = metadata_match,
  de_results  = de_results_ns
), here("results", "GSE160651", "tables", "nanostring_processed_tmp.rds"))

log_section("GSE160651 | Step 1 | COMPLETE")
message("Next: run analysis/integration/01_build_master.R")
