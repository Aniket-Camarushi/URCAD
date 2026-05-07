# ==============================================================================
# analysis/GSE279885/01_setup_tximport.R
#
# Purpose: Load Kallisto pseudo-alignment outputs into R as gene-level
#          count and TPM matrices ready for DESeq2.
#
# WHY TXIMPORT:
#   Kallisto outputs transcript-level estimated counts and TPMs. DESeq2 needs
#   gene-level counts. tximport aggregates transcripts to genes and applies the
#   "scaled" method that corrects for differences in transcript length across
#   samples — a bias that naive summing would ignore. This step is the bridge
#   between the HPC pipeline and the R analysis.
#
# REQUIRED COHORT (runs automatically):
#   7 samples: 3 Vehicle_Control + 4 VCD
#
# OPTIONAL COHORT (uncomment the "OPTIONAL DATASET" blocks below):
#   7 additional samples: 4 Control_shRNA + 3 Esrra_shRNA
#   When combined, DESeq2 includes a "cohort" term in the design to account
#   for the different experimental batches.
#
# Outputs saved to results/GSE279885/tables/:
#   - txi_required.rds      tximport object for required samples
#   - txi_combined.rds      tximport object for required + optional
#   - metadata_required.rds sample metadata data frame (required)
#   - metadata_combined.rds sample metadata data frame (combined)
# ==============================================================================

suppressPackageStartupMessages({
  suppressWarnings({
    library(here)
    library(tximport)
    library(dplyr)
    library(readr)
    library(AnnotationDbi)
    library(org.Mm.eg.db)
    library(GenomicFeatures)
    library(txdbmaker)
  })
})

source(here("analysis", "helpers", "utils.R"))
source(here("analysis", "helpers", "themes.R"))

log_section("GSE279885 | Step 1 | Setup & tximport")

# ── Output directories ────────────────────────────────────────────────────────

dirs <- create_result_dirs("GSE279885")

# ── 1. Build transcript-to-gene mapping from GTF ──────────────────────────────
#
# WHY GTF not org.Mm.eg.db: org.Mm.eg.db only covers RefSeq-registered
# transcripts (22% of Kallisto index). The GTF used to build the Kallisto
# index gives 100% coverage across all transcripts.
#
# This mapping is stored once as an RDS so it is not re-queried every run.

gtf_path <- here("data", "reference", "Mus_musculus.GRCm39.108.gtf")

tx2gene_ensembl <- load_or_compute(
  path = here("results", "GSE279885", "tables", "tx2gene.rds"),
  expr = {
    log_msg("Building tx2gene from GTF (txdbmaker) ...")
    txdb <- txdbmaker::makeTxDbFromGFF(gtf_path, format = "gtf")
    AnnotationDbi::select(
      txdb,
      keys    = keys(txdb, keytype = "TXNAME"),
      columns = "GENEID",
      keytype = "TXNAME"
    ) %>%
      dplyr::rename(target_id = TXNAME, gene_name = GENEID) %>%
      dplyr::filter(!is.na(gene_name)) %>%
      dplyr::distinct()
  }
)

log_msg("tx2gene rows: ", nrow(tx2gene_ensembl))

# ── 2. Load sample metadata ───────────────────────────────────────────────────
#
# REQUIRED cohort: VCD vs Vehicle_Control
meta_required <- load_gse279885_metadata("required")

# ── [OPTIONAL DATASET — UNCOMMENT TO INCLUDE] ────────────────────────────────
# meta_optional <- load_gse279885_metadata("optional")
# meta_combined <- bind_rows(meta_required, meta_optional)
# ─────────────────────────────────────────────────────────────────────────────

log_msg("Required samples: ", nrow(meta_required))

# ── 3. Build file paths to Kallisto abundance.h5 files ───────────────────────
#
# WHY .h5 not .tsv: the HDF5 binary file stores bootstrap replicates needed
# for uncertainty estimation. tximport can read both formats; h5 is preferred.

kallisto_base <- here("data", "processed", "GSE279885", "kallisto")

make_file_paths <- function(metadata) {
  paths <- file.path(kallisto_base, metadata$srr_id, "abundance.h5")
  names(paths) <- metadata$srr_id
  missing <- paths[!file.exists(paths)]
  if (length(missing) > 0) {
    warning("Missing Kallisto output files:\n",
            paste(" ", missing, collapse = "\n"),
            "\nRun the SLURM pipeline first.")
  }
  paths
}

files_required <- make_file_paths(meta_required)

# ── [OPTIONAL DATASET — UNCOMMENT TO INCLUDE] ────────────────────────────────
# files_combined <- make_file_paths(meta_combined)
# ─────────────────────────────────────────────────────────────────────────────

# ── 4. Run tximport ──────────────────────────────────────────────────────────
#
# type = "kallisto": tells tximport the input format
# txOut = FALSE:     aggregate to gene level (not transcript level)
# countsFromAbundance = "lengthScaledTPM": the recommended method for DESeq2.
#   It creates count-scale values that account for transcript-length bias.
#   See: Soneson et al. 2015 (F1000Research) for the statistical justification.

txi_required <- load_or_compute(
  path = here("results", "GSE279885", "tables", "txi_required.rds"),
  expr = {
    log_msg("Running tximport for required cohort ...")
    tximport(
      files_required,
      type                 = "kallisto",
      tx2gene              = tx2gene_ensembl,
      txOut                = FALSE,
      countsFromAbundance  = "lengthScaledTPM",
      ignoreTxVersion      = TRUE   # strip version suffix from Ensembl IDs
    )
  }
)

log_msg("txi_required: ",
        nrow(txi_required$counts), " genes x ",
        ncol(txi_required$counts), " samples")

# ── [OPTIONAL DATASET — UNCOMMENT TO INCLUDE] ────────────────────────────────
# txi_combined <- load_or_compute(
#   path = here("results", "GSE279885", "tables", "txi_combined.rds"),
#   expr = {
#     log_msg("Running tximport for combined cohort ...")
#     tximport(
#       files_combined,
#       type                = "kallisto",
#       tx2gene             = tx2gene,
#       txOut               = FALSE,
#       countsFromAbundance = "lengthScaledTPM",
#       ignoreTxVersion     = TRUE
#     )
#   }
# )
# ─────────────────────────────────────────────────────────────────────────────

# ── 5. Save metadata RDS objects ─────────────────────────────────────────────

saveRDS(meta_required, here("results", "GSE279885", "tables", "metadata_required.rds"))
log_msg("Saved: metadata_required.rds")

# ── [OPTIONAL DATASET — UNCOMMENT TO INCLUDE] ────────────────────────────────
# saveRDS(meta_combined, here("results", "GSE279885", "tables", "metadata_combined.rds"))
# ─────────────────────────────────────────────────────────────────────────────

# ── 6. Quick sanity check ─────────────────────────────────────────────────────

log_msg("Sample names in txi match metadata: ",
        all(colnames(txi_required$counts) == meta_required$srr_id))

# TPM summary: per-sample total (should all be ~1,000,000 after TPM scaling)
tpm_totals <- colSums(txi_required$abundance)
log_msg("Per-sample TPM totals (should be ~1e6):")
print(round(tpm_totals, 0))

# ── 7. Build Ensembl → symbol lookup for DESeq2 output ───────────────────────

ensembl_to_symbol <- AnnotationDbi::select(
  org.Mm.eg.db,
  keys    = rownames(txi_required$counts),
  columns = "SYMBOL",
  keytype = "ENSEMBL"
) %>%
  dplyr::distinct(ENSEMBL, .keep_all = TRUE)

saveRDS(ensembl_to_symbol,
        here("results", "GSE279885", "tables", "ensembl_to_symbol.rds"))
log_msg("Symbols mapped: ", sum(!is.na(ensembl_to_symbol$SYMBOL)),
        " / ", nrow(ensembl_to_symbol))

log_section("GSE279885 | Step 1 | COMPLETE")
message("Next: run analysis/GSE279885/02_deseq2.R")
