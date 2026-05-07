# ==============================================================================
# analysis/helpers/themes.R
#
# Publication-quality ggplot2 theme and color palettes for the URCAD project.
# All figures are white-background, clean panel, suitable for journal submission.
#
# Source this at the top of every analysis script:
#   source(here::here("analysis", "helpers", "themes.R"))
#
# WHY A SEPARATE THEMES FILE:
#   Changing one line here updates every figure in the project.
#   This is how you maintain a consistent visual identity across
#   a multi-dataset paper.
# ==============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(RColorBrewer)
  library(viridis)
  library(here)
})

# ── Core publication theme ────────────────────────────────────────────────────

#' A clean, publication-ready ggplot2 theme.
#'
#' Design choices:
#'   - White background (no grey panel) for journals that require it.
#'   - Black axes and ticks for contrast in print.
#'   - Minimal gridlines (major only, very light) for readability.
#'   - No right/top axis border (classic academic style).
#'   - Bold facet strip labels for multi-panel figures.
#'
#' @param base_size   Base font size in pt (default 12).
#' @param base_family Font family (default "sans" which resolves to Helvetica/Arial).
theme_publication <- function(base_size = 12, base_family = "sans") {
  theme_bw(base_size = base_size, base_family = base_family) %+replace%
    theme(
      # Panel
      panel.background = element_rect(fill = "transparent",   colour = NA),
      plot.background  = element_rect(fill = "transparent",   colour = NA),
      panel.border     = element_rect(fill = NA, colour = "black", linewidth = 0.5),
      panel.grid.major = element_line(colour = "grey92", linewidth = 0.3),
      panel.grid.minor = element_blank(),

      # Axes
      axis.line        = element_line(colour = "black", linewidth = 0.4),
      axis.ticks       = element_line(colour = "black", linewidth = 0.3),
      axis.ticks.length = unit(3, "pt"),
      axis.text        = element_text(size  = base_size * 0.85, colour = "black"),
      axis.title       = element_text(size  = base_size,        colour = "black"),
      axis.title.x     = element_text(margin = margin(t = 8)),
      axis.title.y     = element_text(margin = margin(r = 8), angle = 90),

      # Legend
      legend.background = element_rect(fill = "transparent", colour = NA),
      legend.key        = element_rect(fill = "transparent", colour = NA),
      legend.text       = element_text(size = base_size * 0.80),
      legend.title      = element_text(size = base_size * 0.90, face = "bold"),
      legend.key.size   = unit(0.4, "cm"),

      # Facets
      strip.background  = element_rect(fill = "grey95", colour = "black",
                                        linewidth = 0.5),
      strip.text        = element_text(size = base_size * 0.90, face = "bold",
                                        colour = "black"),

      # Titles
      plot.title        = element_text(size = base_size * 1.10, face = "bold",
                                        hjust = 0, colour = "black"),
      plot.subtitle     = element_text(size = base_size * 0.95, hjust = 0,
                                        colour = "grey40"),
      plot.caption      = element_text(size = base_size * 0.75, hjust = 1,
                                        colour = "grey50"),
      plot.margin       = margin(10, 10, 10, 10)
    )
}

# Set as default theme for the session
theme_set(theme_publication())

# ── Color palettes ────────────────────────────────────────────────────────────

# Experimental groups for GSE279885
PALETTE_VCD <- c(
  "Vehicle_Control" = "#2166AC",   # blue
  "VCD"             = "#D6604D"    # red-orange
)

PALETTE_SHRNA <- c(
  "Control_shRNA" = "#1A9850",     # green
  "Esrra_shRNA"   = "#D73027"      # red
)

# Genotype / treatment for GSE158960
PALETTE_GSE158960 <- c(
  "WT_Control"         = "#4393C3",
  "WT_Tamoxifen"       = "#D6604D",
  "Esr1cKO_Control"   = "#74C476",
  "Esr1cKO_Tamoxifen" = "#9E9AC8"
)

# DE direction
PALETTE_DE <- c(
  "up"   = "#D6604D",
  "down" = "#4393C3",
  "ns"   = "grey70"
)

# Cell types (BBB-relevant highlighted)
PALETTE_CELL_TYPES <- c(
  "Endothelial"  = "#D6604D",
  "Pericyte"     = "#FDAE61",
  "Astrocyte"    = "#74ADD1",
  "Microglia"    = "#A6D96A",
  "Neuron"       = "#ABD9E9",
  "Oligodendrocyte" = "#9E9AC8",
  "OPC"          = "#C2A5CF",
  "Other"        = "grey75"
)

# Continuous scale: expression / module scores
scale_color_expression <- function(...) {
  scale_color_gradient2(
    low  = "#2166AC",
    mid  = "white",
    high = "#D73027",
    midpoint = 0,
    ...
  )
}

# Continuous scale for UMAP density / count
scale_color_umap <- function(...) {
  viridis::scale_color_viridis(option = "magma", ...)
}

# ── Figure saving ─────────────────────────────────────────────────────────────

#' Save a ggplot to both PDF and PNG at publication quality.
#'
#' WHY PDF + PNG:
#'   PDF is vector (infinite resolution, required by most journals).
#'   PNG is raster (required by Word/PowerPoint and for previewing in GitHub).
#'
#' @param plot   A ggplot object (defaults to last plot).
#' @param name   Base filename without extension (e.g. "01_pca").
#' @param outdir Output directory (created if absent).
#' @param width  Width in inches.
#' @param height Height in inches.
#' @param dpi    Resolution for PNG (default 300).
#' @param bg     Background colour (default "transparent").
save_figure <- function(plot    = ggplot2::last_plot(),
                        name,
                        outdir,
                        width   = 6,
                        height  = 5,
                        dpi     = 300,
                        bg      = "transparent") {
  if (!dir.exists(outdir))
    dir.create(outdir, recursive = TRUE)

  pdf_path <- file.path(outdir, paste0(name, ".pdf"))
  png_path <- file.path(outdir, paste0(name, ".png"))

  ggplot2::ggsave(pdf_path, plot = plot, width = width, height = height,
                  bg = bg, device = cairo_pdf)
  ggplot2::ggsave(png_path, plot = plot, width = width, height = height,
                  dpi = dpi, bg = bg)

  message(sprintf("[save_figure] %s.pdf / .png -> %s", name, outdir))
  invisible(list(pdf = pdf_path, png = png_path))
}

message("[themes.R] Publication theme set. Palettes loaded.")
