#!/usr/bin/env bash
#
# scripts/download_sra_fastq.sh
#
# Download FASTQ files from NCBI SRA using prefetch + fasterq-dump.
#
# Key change from original: adds --cohort flag so the caller can download
# just the required samples, just the optional samples, or both in one call.
# This is critical for the GSE279885 project which has two experimental
# cohorts (VCD/Vehicle = required; Esrra_shRNA/Control_shRNA = optional).
#
# Usage:
#   bash scripts/download_sra_fastq.sh \
#       --paper-id sun_jk_GSE279885 \
#       --cohort required \
#       --threads 10
#
#   bash scripts/download_sra_fastq.sh \
#       --paper-id sun_jk_GSE279885 \
#       --cohort both \
#       --threads 10
#
# Exit codes:
#   0  success
#   1  configuration / argument error
#   2  missing tools
#   3  download or conversion error
#

set -euo pipefail

# ─────────────────────────────────────────────────────────────────────────────
# Defaults
# ─────────────────────────────────────────────────────────────────────────────

PROJECT_ROOT="$(pwd)"
PAPER_ID=""
COHORT="both"          # options: required | optional | both | custom
ACCESSION_FILE=""          # if --cohort custom, provide your own file

DATA_RAW_DIR="${PROJECT_ROOT}/data/raw"
CACHEDIR="${DATA_RAW_DIR}/sra_cache"
TMPDIR="${PROJECT_ROOT}/tmp/sra_tmp"
LOGDIR="${PROJECT_ROOT}/logs"
THREADS=10

# ─────────────────────────────────────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────────────────────────────────────

usage() {
    cat <<EOF
Usage: $(basename "$0") [options]

Options:
  --paper-id NAME       Short identifier used to locate data/raw/NAME/fastq
                        (e.g. sun_jk_GSE279885).
  --cohort COHORT       Which accession list to download.
                        One of: required | optional | both | custom
                        (default: required)
  --accession-file FILE Only used when --cohort custom. Path to a plain-text
                        file with one SRA accession per non-comment line.
  --threads N           Threads for fasterq-dump (default: ${THREADS}).
  -h, --help            Show this message.

Cohort files are resolved automatically from config/:
  required  -> config/accessions_<PAPER_ID>_required.txt
               Falls back to config/accessions_<PAPER_ID>.txt
  optional  -> config/accessions_<PAPER_ID>_optional.txt
  both      -> required + optional merged

Examples:
  bash scripts/download_sra_fastq.sh --paper-id sun_jk_GSE279885 --cohort required
  bash scripts/download_sra_fastq.sh --paper-id sun_jk_GSE279885 --cohort both --threads 8
EOF
}

log() { printf "[%s] %s\n" "$(date '+%Y-%m-%d %H:%M:%S')" "$*" >&2; }
die() { log "ERROR: $*"; exit 1; }

# ─────────────────────────────────────────────────────────────────────────────
# Argument parsing
# ─────────────────────────────────────────────────────────────────────────────

while [[ $# -gt 0 ]]; do
    case "$1" in
        --paper-id)        PAPER_ID="${2:-}";        shift 2 ;;
        --cohort)          COHORT="${2:-}";           shift 2 ;;
        --accession-file)  ACCESSION_FILE="${2:-}";  shift 2 ;;
        --threads)         THREADS="${2:-}";          shift 2 ;;
        -h|--help)         usage; exit 0 ;;
        *) die "Unknown argument: $1 (use --help)" ;;
    esac
done

[[ -z "${PAPER_ID}" ]] && die "--paper-id is required."

OUTDIR="${DATA_RAW_DIR}/${PAPER_ID}/fastq"

# ─────────────────────────────────────────────────────────────────────────────
# Tool checks
# ─────────────────────────────────────────────────────────────────────────────

command -v prefetch     >/dev/null 2>&1 || die "prefetch not found. Activate SRA-tools conda env."
command -v fasterq-dump >/dev/null 2>&1 || die "fasterq-dump not found. Activate SRA-tools conda env."

# ─────────────────────────────────────────────────────────────────────────────
# Build accession list from cohort
# ─────────────────────────────────────────────────────────────────────────────

RUNS=()

load_accessions_from_file() {
    local file="$1"
    [[ -f "$file" ]] || die "Accession file not found: $file"
    local arr
    mapfile -t arr < <(grep -E '^[^#[:space:]]' "$file")
    [[ ${#arr[@]} -gt 0 ]] || die "No accessions found in $file"
    log "  Loaded ${#arr[@]} accessions from $(basename "$file")"
    RUNS+=("${arr[@]}")
}

# Derive canonical config filenames from PAPER_ID.
# Strip any GEO prefix (e.g. sun_jk_GSE279885 -> GSE279885) for the config.
GEO_ID="${PAPER_ID##*_}"   # last segment after final underscore

REQ_FILE="${PROJECT_ROOT}/config/accessions_${GEO_ID}_required.txt"
OPT_FILE="${PROJECT_ROOT}/config/accessions_${GEO_ID}_optional.txt"
# Legacy fallback: if paper used a single accessions file
LEGACY_FILE="${PROJECT_ROOT}/config/accessions_${GEO_ID}.txt"

case "${COHORT}" in
    required)
        if [[ -f "$REQ_FILE" ]]; then
            load_accessions_from_file "$REQ_FILE"
        elif [[ -f "$LEGACY_FILE" ]]; then
            log "No *_required.txt found; falling back to $LEGACY_FILE"
            load_accessions_from_file "$LEGACY_FILE"
        else
            die "Required accession file not found: $REQ_FILE"
        fi
        ;;
    optional)
        load_accessions_from_file "$OPT_FILE"
        ;;
    both)
        load_accessions_from_file "$REQ_FILE"
        load_accessions_from_file "$OPT_FILE"
        ;;
    custom)
        [[ -n "$ACCESSION_FILE" ]] || die "--accession-file is required when --cohort custom."
        load_accessions_from_file "$ACCESSION_FILE"
        ;;
    *)
        die "Unknown --cohort value: ${COHORT}. Use required | optional | both | custom"
        ;;
esac

[[ ${#RUNS[@]} -gt 0 ]] || die "No accessions loaded."

# ─────────────────────────────────────────────────────────────────────────────
# Prepare directories and logging
# ─────────────────────────────────────────────────────────────────────────────

mkdir -p "${OUTDIR}" "${CACHEDIR}" "${TMPDIR}" "${LOGDIR}"

LOGFILE="${LOGDIR}/download_${PAPER_ID}_${COHORT}.log"
log "Logging to: ${LOGFILE}"
exec > >(tee -a "${LOGFILE}") 2>&1

log "Paper ID   : ${PAPER_ID}"
log "Cohort     : ${COHORT}"
log "Output dir : ${OUTDIR}"
log "Accessions : ${#RUNS[@]}"
log ""
for acc in "${RUNS[@]}"; do log "  - ${acc}"; done
log ""

# ─────────────────────────────────────────────────────────────────────────────
# Main download loop
# ─────────────────────────────────────────────────────────────────────────────

for ACC in "${RUNS[@]}"; do
    log "────────────────────────────────────────"
    log "Processing: ${ACC}"

    # Skip if FASTQ already exists (idempotent re-runs)
    existing=$(find "${OUTDIR}" -maxdepth 1 -name "${ACC}*.fastq" 2>/dev/null | wc -l)
    if [[ "${existing}" -gt 0 ]]; then
        log "  FASTQs already present for ${ACC}. Skipping."
        continue
    fi

    # Step 1: prefetch (downloads the SRA object to cache)
    log "  prefetch ${ACC} ..."
    prefetch "${ACC}" -O "${CACHEDIR}" \
        || die "prefetch failed for ${ACC}"

    # Locate the downloaded .sra file
    SRA_PATH=""
    [[ -d "${CACHEDIR}/${ACC}" ]] && \
        SRA_PATH=$(find "${CACHEDIR}/${ACC}" -maxdepth 1 -name "*.sra" | head -n1 || true)
    [[ -z "$SRA_PATH" && -f "${CACHEDIR}/${ACC}.sra" ]] && \
        SRA_PATH="${CACHEDIR}/${ACC}.sra"
    [[ -z "$SRA_PATH" ]] && \
        die "Could not locate .sra for ${ACC} under ${CACHEDIR}"

    # Step 2: fasterq-dump (converts SRA to FASTQ, outputs paired-end _1/_2)
    log "  fasterq-dump ${ACC} ..."
    (
        cd "${OUTDIR}"
        fasterq-dump "${SRA_PATH}" \
            -t "${TMPDIR}" \
            -e "${THREADS}" \
            -p \
            --split-files
    ) || die "fasterq-dump failed for ${ACC}"

    log "  Done: ${ACC}"
done

log ""
log "All accessions for cohort '${COHORT}' processed successfully."
log "FASTQs are in: ${OUTDIR}"
