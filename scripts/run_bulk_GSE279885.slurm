#!/bin/bash
#
# scripts/run_bulk_GSE279885.slurm
#
# SLURM batch script to run the bulk RNA-seq pipeline for GSE279885
# on the UMBC chip-cpu cluster.
#
# Changes from original:
#   - Renamed to run_bulk_GSE279885.slurm (descriptive name)
#   - Added COHORT variable so you can download required, optional, or both
#     before running the pipeline in one submission.
#   - Added --skip-gatk flag (optional; remove if you need variant calling)
#   - Uses updated path conventions (data/raw/ not data_raw/)
#
# Usage:
#   # Download + process required cohort only (default):
#   sbatch scripts/run_bulk_GSE279885.slurm
#
#   # To run with optional cohort included, change COHORT=both below.
#
# ─────────────────────────────────────────────────────────────────────────────

#SBATCH --job-name=bulk_GSE279885
#SBATCH --cluster=chip-cpu
#SBATCH --partition=general
#SBATCH --account=pi_carmunoz
#SBATCH --qos=long
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=30
#SBATCH --mem=200G
#SBATCH --time=1-00:00:00
#SBATCH --output=logs/bulk_GSE279885_%j.out
#SBATCH --error=logs/bulk_GSE279885_%j.err

# ─────────────────────────────────────────────────────────────────────────────
# Configuration
# Edit COHORT here before submitting.
#   required  -> 7 samples (VCD vs Vehicle); already downloaded
#   optional  -> 7 samples (Esrra_shRNA vs Control_shRNA); needs download
#   both      -> all 14 samples
# ─────────────────────────────────────────────────────────────────────────────

PAPER_ID="sun_jk_GSE279885"
COHORT="both"           # <-- change to "optional" or "both" as needed
THREADS=30

# ─────────────────────────────────────────────────────────────────────────────
# Environment
# ─────────────────────────────────────────────────────────────────────────────

module load anaconda
conda init
source ~/.bashrc
conda activate ngs_pipeline

echo "Active conda environment:"
conda list | grep -E 'fastqc|kallisto|star |samtools|picard|gatk4|multiqc|sra-tools'

# ─────────────────────────────────────────────────────────────────────────────
# Step A: Download FASTQs for the requested cohort
#
# The download script is idempotent: it skips any SRR whose FASTQs already
# exist in data/raw/sun_jk_GSE279885/fastq/. Safe to re-run.
# ─────────────────────────────────────────────────────────────────────────────

echo "=========================================="
echo "Step A: Downloading FASTQ files"
echo "Paper ID : ${PAPER_ID}"
echo "Cohort   : ${COHORT}"
echo "Started  : $(date)"
echo "=========================================="

bash ./download_sra_fastq.sh \
    --paper-id "${PAPER_ID}" \
    --cohort   "${COHORT}"  \
    --threads  "${THREADS}"

echo "Download complete: $(date)"

# ─────────────────────────────────────────────────────────────────────────────
# Step B: Run the bulk RNA-seq pipeline
#
# FastQC -> Kallisto -> STAR -> Picard -> (GATK) -> MultiQC
# Remove --skip-gatk below if you need variant calls.
# ─────────────────────────────────────────────────────────────────────────────

echo "=========================================="
echo "Step B: Running bulk RNA-seq pipeline"
echo "Started: $(date)"
echo "=========================================="

bash ./run_pipeline.sh \
    --paper-id "${PAPER_ID}" \
    --threads  "${THREADS}"              # remove if you need GATK variant calls

echo "Pipeline complete: $(date)"
echo "All steps finished for ${PAPER_ID} (${COHORT} cohort)."
