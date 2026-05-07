# ==============================================================================
# Dockerfile — URCAD R Analysis Environment
#
# Base image: Bioconductor Docker (R 4.5 + Bioconductor 3.21 pre-installed).
# This image ships with most Bioconductor system dependencies already present,
# saving ~2 hours of compilation compared to building from scratch.
#
# Build:
#   docker build -t urcad:latest .
#
# Run RStudio Server interactively (maps to http://localhost:8787):
#   docker run -d -p 8787:8787 \
#     -v "D:/Bioinfomatics Scripts/HPC/URCAD:/home/rstudio/URCAD" \
#     -e PASSWORD=urcad2025 \
#     urcad:latest
#
# Run a script non-interactively:
#   docker run --rm \
#     -v "D:/Bioinfomatics Scripts/HPC/URCAD:/URCAD" \
#     urcad:latest \
#     Rscript /URCAD/analysis/GSE279885/01_setup.R
#
# ==============================================================================

FROM bioconductor/bioconductor_docker:devel

LABEL maintainer="pi_carmunoz"
LABEL description="URCAD: Estrogen effects on BBB - scRNA/snRNA/bulk/NanoString analysis"
LABEL r_version="4.5"
LABEL bioconductor_version="3.21"

# ── System libraries required by R packages ──────────────────────────────────

RUN apt-get update && apt-get install -y --no-install-recommends \
    libglpk-dev       \
    libhdf5-dev       \
    libigraph-dev     \
    libgsl-dev        \
    libmagick++-dev   \
    libharfbuzz-dev   \
    libfribidi-dev    \
    libgit2-dev       \
    libssl-dev        \
    libcurl4-openssl-dev \
    libxml2-dev       \
    libudunits2-dev   \
    libgdal-dev       \
    libproj-dev       \
    libjpeg-dev       \
    libpng-dev        \
    libtiff-dev       \
    libfreetype6-dev  \
    pandoc            \
    && rm -rf /var/lib/apt/lists/*

# ── renv bootstrap ────────────────────────────────────────────────────────────
# We copy renv.lock first so Docker can cache the layer when only
# source files change (no need to re-install all packages).

WORKDIR /URCAD
COPY renv.lock .
COPY .Rprofile .

RUN R -e "install.packages('renv', repos = 'https://cloud.r-project.org')"
RUN R -e "renv::restore(prompt = FALSE)"

# ── Copy project files ────────────────────────────────────────────────────────

COPY analysis/ analysis/
COPY config/   config/

# ── Default command: start RStudio Server ────────────────────────────────────

CMD ["/init"]
