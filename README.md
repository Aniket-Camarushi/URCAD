# URCAD: Estrogen–TBI Blood-Brain Barrier Analysis Pipeline

Multi-dataset R pipeline investigating how estrogen depletion interacts with
traumatic brain injury (TBI) at the blood-brain barrier (BBB). Integrates
single-cell RNA-seq, NanoString neuropathology panel, and bulk RNA-seq data
across four public GEO datasets.

Developed as an undergraduate honors thesis at the University of
Maryland, Baltimore County (UMBC) in the Munoz-Ballester Lab.

---

## Biological Question

Estrogen is thought to be neuroprotective, but the molecular overlap between
estrogen loss and TBI-driven BBB disruption is not well characterized. This
pipeline asks: do estrogen-deficient models and TBI models dysregulate the
same BBB genes, and which cell types drive that overlap?

---

## Datasets

| GEO ID | Data type | Model | Comparisons |
|---|---|---|---|
| GSE158960 | scRNA-seq (10x Genomics) | Mouse BBB — WT vs *Esr1*cKO × Control vs Tamoxifen | 16 samples, 4 groups |
| GSE160651 | NanoString Neuropathology Panel | Mouse TBI | TBI 1d / 7d / 30d vs Sham |
| GSE279885 | Bulk RNA-seq | Mouse VCD estrogen depletion | VCD vs Vehicle Control |
| GSE299005 | Single-nucleus RNA-seq | Mouse (in progress) | — |

Raw data is not included. See [Data Access](#data-access) below.

---

## Pipeline Structure
analysis/
├── GSE158960/
│ ├── 01_load_qc.R Per-sample QC, filtering, merging
│ ├── 02_cluster.R SCTransform, PCA, UMAP, Leiden clustering
│ ├── 03_annotate.R Cell type annotation (SingleR + manual)
│ ├── 04_bbb_estrogen.R BBB module scoring, pseudobulk DESeq2
│ └── 05_cellchat.R CellChat ligand-receptor analysis
├── GSE160651/
│ └── 01_load_nanostring.R NACHO QC, normalization, DESeq2, GSEA
├── GSE279885/
│ ├── 01_setup_tximport.R Kallisto → tximport, tx2gene mapping
│ ├── 02_deseq2.R DESeq2 differential expression
│ └── 03_enrichment.R GSEA (MSigDB Hallmark), KEGG dotplots
├── GSE299005/
│ └── 01_load_sn.R snRNA-seq QC (in progress)
├── integration/
│ ├── 01_build_master.R Cross-dataset gene overlap table
│ └── 02_cross_dataset_pathways.R Pathway-level integration
└── helpers/
├── utils.R Shared helper functions
└── themes.R ggplot2 publication theme


Key genes tracked across datasets: *Cldn5*, *Tjp1*, *Ocln* (tight junctions);
*Esr1*, *Esrra*, *Gper1* (estrogen receptors); *Aqp4*, *Gfap* (astrocyte);
*Pdgfrb*, *Notch3* (pericyte); *Il6*, *Il1b*, *Tnf* (neuroinflammation).

---

## Requirements

- **R** ≥ 4.5.0
- **renv** for package management (all versions pinned in `renv.lock`)
- **Kallisto** for bulk RNA-seq pseudo-alignment (GSE279885)
- **STAR** alignment outputs can substitute for Kallisto (see `scripts/`)
- Optional: Docker / Snakemake (see below)

---

## Getting Started

### 1. Clone the repo

```bash
git clone https://github.com/Aniket-Camarushi/URCAD.git
cd URCAD
```

### 2. Restore the R environment

```r
# In R, from the project root:
install.packages("renv")
renv::restore()
```

This installs all packages at the exact versions recorded in `renv.lock`.
Expect ~15–20 min on first run due to Bioconductor packages.

### 3. Download raw data

```bash
# GSE279885 (bulk RNA-seq) — uses SRA Toolkit
bash scripts/download_sra_fastq.sh config/accessions_GSE279885_required.txt

# GSE158960, GSE160651, GSE299005 — download manually from GEO
# and place in data/raw/<GSEID>/
```

### 4. Run the pipeline

**Option A — Run scripts sequentially in R/RStudio**

Open `URCAD.Rproj` and source each script in numbered order within each
dataset folder. Outputs go to `results/<GSEID>/`.

**Option B — Snakemake (bulk processing)**

```bash
snakemake --snakefile workflow/Snakefile --cores 8
```

> The Snakemake workflow covers GSE279885 alignment and quantification.
> scRNA-seq scripts are designed to be run interactively in RStudio.

**Option C — Docker**

```bash
docker-compose up
```

> Docker image includes R 4.5, all renv packages, and Kallisto.
> Image has not been pushed to a registry; build locally with `docker build`.

---

## Configuration

All tunable parameters are in `config/pipeline_params.yaml`:

```yaml
# Key defaults
qc_scrna:
  min_features: 200
  max_features: 4000
  max_pct_mt: 15
deseq2:
  alpha: 0.05
  min_count: 10
gsea:
  n_perm: 1000
  min_gs_size: 15
```

Sample metadata and SRA accession lists are in `config/`.

---

## Data Access

Raw data is excluded from this repository (too large for GitHub; see
`.gitignore`). All datasets are publicly available on NCBI GEO:

- [GSE158960](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE158960)
- [GSE160651](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE160651)
- [GSE279885](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE279885)
- [GSE299005](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE299005)

Reference genome: Mus musculus GRCm39 (Ensembl release 108).

---

## Dependencies (key packages)

| Package | Purpose |
|---|---|
| Seurat ≥ 5 | scRNA-seq processing |
| DESeq2 | Differential expression |
| CellChat | Ligand-receptor communication |
| fgsea / clusterProfiler | GSEA and pathway analysis |
| tximport | Transcript-level import |
| NACHO | NanoString QC and normalization |
| here, dplyr, ggplot2, patchwork | Infrastructure |

Full version list: `renv.lock`

---

## License

MIT License — see `LICENSE`.

---

## Acknowledgements

Developed as a departmental honors thesis in the Munoz-Ballester Lab at UMBC.
Recipient of the Biological Sciences Departmental Honors in Research Award
from the UMBC Department of Biological Sciences.

Datasets from publicly available GEO submissions by the original authors;
this repository contains only analysis code, not the underlying data.
