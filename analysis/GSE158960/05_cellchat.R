# ==============================================================================
# analysis/GSE158960/05_cellchat.R
#
# Purpose: Infer cell-cell communication via ligand-receptor (LR) pairs
#          in the BBB microenvironment using CellChat v2.
#
# WHY CELLCHAT:
#   CellChat infers signaling networks by combining expression evidence
#   (is the ligand expressed in sender cells? is the receptor expressed in
#   receiver cells?) with a curated database of experimentally validated
#   LR interactions. It then computes interaction probability using the
#   law of mass action.
#   For estrogen-BBB research, we specifically look for:
#     - PDGF signaling (pericyte maintenance)
#     - VEGF signaling (permeability regulation)
#     - Notch signaling (tight junction regulation)
#     - Complement signaling (neuroinflammation)
#     - Estrogen-induced pathways (ANGPT, WNT, TGFb in endothelial)
#
# Strategy:
#   Run CellChat separately on WT_Control and Esr1cKO_Control (no tamoxifen)
#   to isolate the estrogen receptor effect. Then use compareCellChat to
#   identify pathways that are gained or lost when ERα is ablated.
#
# Outputs:
#   - cellchat_WT.rds        CellChat object: WT_Control
#   - cellchat_CKO.rds       CellChat object: Esr1cKO_Control
#   - 10_cellchat_comparison.pdf
# ==============================================================================

suppressPackageStartupMessages({
  library(here)
  library(Seurat)
  library(CellChat)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

source(here("analysis", "helpers", "utils.R"))
source(here("analysis", "helpers", "themes.R"))

log_section("GSE158960 | Step 5 | CellChat LR Analysis")
dirs <- create_result_dirs("GSE158960")

seurat_bbb <- readRDS(
  here("results", "GSE158960", "tables", "seurat_bbb_scored.rds")
)

# CellChat requires cell type labels in the "idents" slot
Idents(seurat_bbb) <- "cell_type"

# ── Helper: build and run a CellChat object ───────────────────────────────────
#
# This function encapsulates the full CellChat workflow for one group.
# Running it as a function avoids code duplication between WT and CKO.

run_cellchat <- function(seurat_subset, group_label) {
  log_msg("Building CellChat for: ", group_label)

  # Extract normalized expression matrix and metadata
  data_input  <- GetAssayData(seurat_subset,
                               assay = "SCT", layer = "data")
  meta_data   <- seurat_subset@meta.data[, "cell_type", drop = FALSE]
  colnames(meta_data) <- "labels"

  # Create CellChat object
  cc <- createCellChat(object   = data_input,
                        meta     = meta_data,
                        group.by = "labels")

  # Load mouse CellChat database (secreted + ECM + contact interactions)
  # WHY CellChatDB.mouse: the interaction probabilities are calibrated for
  # mouse proteins; do not use the human database for mouse data.
  CellChatDB <- CellChatDB.mouse

  # Subset to "Secreted Signaling" and "Cell-Cell Contact" categories.
  # We exclude "ECM-Receptor" in this first pass to keep the analysis focused
  # on paracrine and juxtacrine signals relevant to the BBB microenvironment.
  cc@DB <- subsetDB(CellChatDB,
                    search = c("Secreted Signaling", "Cell-Cell Contact"))

  # Pre-process: identify over-expressed genes/interactions
  cc <- subsetData(cc)
  cc <- identifyOverExpressedGenes(cc)
  cc <- identifyOverExpressedInteractions(cc)

  # Compute communication probabilities
  # population.size = TRUE: corrects for differences in cell population size
  cc <- computeCommunProb(cc, type = "triMean",
                           population.size = TRUE)

  # Filter out interactions supported by very few cells
  cc <- filterCommunication(cc, min.cells = 10)

  # Summarise at pathway level
  cc <- computeCommunProbPathway(cc)

  # Aggregate network (counts + weights)
  cc <- aggregateNet(cc)

  saveRDS(cc, here("results", "GSE158960", "tables",
                    paste0("cellchat_", group_label, ".rds")))
  log_msg("Saved: cellchat_", group_label, ".rds")
  cc
}

# ── 1. Run CellChat on WT_Control and Esr1cKO_Control ────────────────────────

cc_wt <- load_or_compute(
  path = here("results", "GSE158960", "tables", "cellchat_WT_Control.rds"),
  expr = {
    wt_cells <- subset(seurat_bbb, subset = group == "WT_Control")
    run_cellchat(wt_cells, "WT_Control")
  }
)

cc_cko <- load_or_compute(
  path = here("results", "GSE158960", "tables", "cellchat_CKO_Control.rds"),
  expr = {
    cko_cells <- subset(seurat_bbb, subset = group == "Esr1cKO_Control")
    run_cellchat(cko_cells, "Esr1cKO_Control")
  }
)

# ── 2. Interaction chord diagram ─────────────────────────────────────────────

pdf(here("results", "GSE158960", "figures", "10_cellchat_chord_WT.pdf"),
    width = 8, height = 8)
netVisual_circle(cc_wt@net$count,
                 vertex.weight   = as.numeric(table(cc_wt@idents)),
                 weight.scale    = TRUE,
                 title.name      = "Interaction count — WT_Control BBB cells")
dev.off()

pdf(here("results", "GSE158960", "figures", "10_cellchat_chord_CKO.pdf"),
    width = 8, height = 8)
netVisual_circle(cc_cko@net$count,
                 vertex.weight   = as.numeric(table(cc_cko@idents)),
                 weight.scale    = TRUE,
                 title.name      = "Interaction count — Esr1cKO_Control BBB cells")
dev.off()

# ── 3. Compare WT vs CKO: which pathways are gained/lost? ───────────────────
#
# WHY compare not just describe:
#   The biological question is specifically about what ERα ablation (Esr1cKO)
#   does to BBB cell-cell communication. compareCellChat quantifies the
#   differential information flow through each signaling pathway.

cc_list <- list(WT_Control = cc_wt, Esr1cKO_Control = cc_cko)

# Merge object sizes to same cell type set
all_idents <- union(levels(cc_wt@idents), levels(cc_cko@idents))
cc_wt  <- liftCellChat(cc_wt,  group.new = all_idents)
cc_cko <- liftCellChat(cc_cko, group.new = all_idents)

# cc_wt  <- liftCellChat(cc_wt,  group.new = levels(cc_cko@idents))
# cc_cko <- liftCellChat(cc_cko, group.new = levels(cc_wt@idents))

cc_combined <- mergeCellChat(list(cc_wt, cc_cko),
                              add.names = names(cc_list))

# Information flow comparison
pdf(here("results", "GSE158960", "figures", "11_cellchat_infoflow.pdf"),
    width = 9, height = 6)
p_infoflow <- rankNet(cc_combined,
                       mode        = "comparison",
                       stacked     = TRUE,
                       do.stat     = TRUE,
                       color.use   = c("#4393C3", "#74C476"))
print(p_infoflow)
dev.off()

# ── 4. Bubble plot: LR pairs between BBB cell types ──────────────────────────
#
# Shows which specific LR pairs are active in WT but not CKO (or vice versa).
# Focus on endothelial and astrocyte as receiver cells (key BBB cell types).

pdf(here("results", "GSE158960", "figures", "12_cellchat_bubble_endo_recvr.pdf"),
    width = 10, height = 8)
p_bubble <- netVisual_bubble(cc_wt,
                              sources.use = c("Pericyte", "Astrocyte"),
                              targets.use = "Endothelial",
                              remove.isolate = TRUE)
print(p_bubble)
dev.off()

log_msg("CellChat analysis complete.")
log_msg("Key outputs: cellchat_WT_Control.rds, 10_cellchat_chord_WT.pdf, 11_cellchat_infoflow.pdf")

log_section("GSE158960 | Step 5 | COMPLETE")
message("Next: run analysis/GSE299005/01_load_sn.R")
