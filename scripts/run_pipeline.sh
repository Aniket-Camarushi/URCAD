#!/usr/bin/env bash
#
# scripts/run_pipeline.sh
#
# Bulk RNA-seq pipeline: FASTQ -> FastQC -> Kallisto -> STAR -> Picard -> GATK -> MultiQC
#
# Changes from original:
#   - data_raw/  renamed to  data/raw/   (unified URCAD convention)
#   - processed outputs go to  data/processed/<paper-id>/
#   - All path logic is relative to PROJECT_ROOT (pwd when called)
#   - Added --skip-gatk flag (GATK variant calling is optional for RNA-seq DE)
#
# Run from project root:
#   bash scripts/run_pipeline.sh --paper-id sun_jk_GSE279885 --threads 20
#   bash scripts/run_pipeline.sh --paper-id sun_jk_GSE279885 --threads 20 --skip-gatk

set -euo pipefail

# ─────────────────────────────────────────────────────────────────────────────
# Defaults
# ─────────────────────────────────────────────────────────────────────────────

PAPER_ID=""
THREADS=20
SKIP_GATK=false

PROJECT_ROOT="$(pwd)"
DATA_DIR="${PROJECT_ROOT}/data"
REF_DIR="${PROJECT_ROOT}/data/reference"
LOGDIR="${PROJECT_ROOT}/logs"

BAM_RAM_LIMIT=30000000000   # 30 GB for STAR sort
MAX_JAVA_MEM="32g"

# ─────────────────────────────────────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────────────────────────────────────

usage() {
    cat <<EOF
Usage: $(basename "$0") --paper-id ID [--threads N] [--skip-gatk]

Options:
  --paper-id ID   Dataset identifier (FASTQs expected at data/raw/ID/fastq/).
  --threads N     CPU threads (default: ${THREADS}).
  --skip-gatk     Skip GATK variant-calling steps (faster; not needed for DE).
  -h, --help      Show this message.

Example:
  bash scripts/$(basename "$0") --paper-id sun_jk_GSE279885 --threads 20
EOF
}

log() { printf "[%s] %s\n" "$(date '+%Y-%m-%d %H:%M:%S')" "$*" >&2; }
die() { log "ERROR: $*"; exit 1; }

# ─────────────────────────────────────────────────────────────────────────────
# Argument parsing
# ─────────────────────────────────────────────────────────────────────────────

while [[ $# -gt 0 ]]; do
    case "$1" in
        --paper-id)    PAPER_ID="${2:-}"; shift 2 ;;
        --threads)     THREADS="${2:-}";  shift 2 ;;
        --skip-gatk)   SKIP_GATK=true;   shift   ;;
        -h|--help)     usage; exit 0 ;;
        *) die "Unknown argument: $1" ;;
    esac
done

[[ -z "${PAPER_ID}" ]] && die "--paper-id is required."

# Derive working paths
FASTQ_DIR="${DATA_DIR}/raw/${PAPER_ID}/fastq"
PROC_DIR="${DATA_DIR}/processed/${PAPER_ID}"
FASTQC_DIR="${PROC_DIR}/fastqc"
KALLISTO_DIR="${PROC_DIR}/kallisto"
STAR_DIR="${PROC_DIR}/star"
GATK_DIR="${PROC_DIR}/gatk"
MULTIQC_DIR="${PROC_DIR}/multiqc"
LOGFILE="${LOGDIR}/pipeline_${PAPER_ID}.log"

# ─────────────────────────────────────────────────────────────────────────────
# Pre-flight checks
# ─────────────────────────────────────────────────────────────────────────────

[[ -d "${FASTQ_DIR}" ]] || die "FASTQ directory not found: ${FASTQ_DIR}"

for tool in fastqc multiqc kallisto STAR samtools picard gatk; do
    command -v "${tool}" >/dev/null 2>&1 || die "${tool} not in PATH."
done

CDNA_FA="${REF_DIR}/Mus_musculus.GRCm39.cdna.all.fa"
GENOME_FA="${REF_DIR}/Mus_musculus.GRCm39.dna.primary_assembly.fa"
GTF_FILE="${REF_DIR}/Mus_musculus.GRCm39.108.gtf"
[[ -f "${CDNA_FA}" ]]   || die "Missing cDNA FASTA: ${CDNA_FA}"
[[ -f "${GENOME_FA}" ]] || die "Missing genome FASTA: ${GENOME_FA}"
[[ -f "${GTF_FILE}" ]]  || die "Missing GTF: ${GTF_FILE}"

# ─────────────────────────────────────────────────────────────────────────────
# Set up dirs and logging
# ─────────────────────────────────────────────────────────────────────────────

mkdir -p "${FASTQC_DIR}" "${KALLISTO_DIR}" "${STAR_DIR}" \
         "${GATK_DIR}" "${MULTIQC_DIR}" "${LOGDIR}"

log "Logging to ${LOGFILE}"
exec > >(tee -a "${LOGFILE}") 2>&1

log "Paper ID   : ${PAPER_ID}"
log "FASTQ dir  : ${FASTQ_DIR}"
log "Threads    : ${THREADS}"
log "Skip GATK  : ${SKIP_GATK}"

# ─────────────────────────────────────────────────────────────────────────────
# Step 1: FastQC
# ─────────────────────────────────────────────────────────────────────────────

log "━━━ Step 1: FastQC ━━━"
fastqc "${FASTQ_DIR}"/*.fastq -t "${THREADS}" -o "${FASTQC_DIR}"

# ─────────────────────────────────────────────────────────────────────────────
# Step 2: Kallisto quantification
# ─────────────────────────────────────────────────────────────────────────────

log "━━━ Step 2: Kallisto ━━━"

KALLISTO_INDEX="${KALLISTO_DIR}/mmusculus_GRCm39.cdna.index"
if [[ ! -f "${KALLISTO_INDEX}" ]]; then
    log "  Building Kallisto index..."
    kallisto index -i "${KALLISTO_INDEX}" "${CDNA_FA}"
fi

for fq1 in "${FASTQ_DIR}"/*_1.fastq; do
    [[ -e "${fq1}" ]] || { log "  No *_1.fastq found in ${FASTQ_DIR}"; break; }
    sample=$(basename "${fq1%_1.fastq}")
    fq2="${FASTQ_DIR}/${sample}_2.fastq"
    [[ -f "${fq2}" ]] || die "  Missing mate: ${fq2}"
    outdir="${KALLISTO_DIR}/${sample}"
    mkdir -p "${outdir}"
    log "  Quantifying ${sample}..."
    kallisto quant -i "${KALLISTO_INDEX}" -o "${outdir}" -t "${THREADS}" "${fq1}" "${fq2}"
done

# ─────────────────────────────────────────────────────────────────────────────
# Step 3: STAR genome index
# ─────────────────────────────────────────────────────────────────────────────

log "━━━ Step 3: STAR index ━━━"

STAR_INDEX_DIR="${STAR_DIR}/genome_index"
if [[ ! -d "${STAR_INDEX_DIR}" ]] || [[ -z "$(ls -A "${STAR_INDEX_DIR}" 2>/dev/null)" ]]; then
    mkdir -p "${STAR_INDEX_DIR}"
    STAR \
        --runThreadN "${THREADS}" \
        --runMode genomeGenerate \
        --genomeDir "${STAR_INDEX_DIR}" \
        --genomeFastaFiles "${GENOME_FA}" \
        --sjdbGTFfile "${GTF_FILE}" \
        --sjdbOverhang 100
else
    log "  STAR index already exists; skipping."
fi

# ─────────────────────────────────────────────────────────────────────────────
# Step 4: STAR alignment
# ─────────────────────────────────────────────────────────────────────────────

log "━━━ Step 4: STAR alignment ━━━"

for fq1 in "${FASTQ_DIR}"/*_1.fastq; do
    [[ -e "${fq1}" ]] || { log "  No *_1.fastq found"; break; }
    sample=$(basename "${fq1%_1.fastq}")
    fq2="${FASTQ_DIR}/${sample}_2.fastq"
    [[ -f "${fq2}" ]] || die "  Missing mate: ${fq2}"
    log "  Aligning ${sample}..."
    STAR \
        --genomeDir "${STAR_INDEX_DIR}" \
        --readFilesIn "${fq1}" "${fq2}" \
        --runThreadN "${THREADS}" \
        --limitBAMsortRAM "${BAM_RAM_LIMIT}" \
        --outSAMtype BAM SortedByCoordinate \
        --quantMode TranscriptomeSAM GeneCounts \
        --outFileNamePrefix "${STAR_DIR}/${sample}_" \
        --outSAMattrRGline ID:"${sample}" SM:"${sample}" LB:lib1 PL:ILLUMINA PU:unit1 \
        --outFilterType BySJout \
        --outFilterMultimapNmax 20 \
        --alignSJoverhangMin 8 \
        --alignSJDBoverhangMin 1 \
        --outFilterMismatchNmax 999 \
        --outFilterMismatchNoverLmax 0.04 \
        --alignIntronMin 20 \
        --alignIntronMax 1000000 \
        --alignMatesGapMax 1000000 \
        --twopassMode Basic \
        --outSAMattributes NH HI AS nM NM MD XS
done

# ─────────────────────────────────────────────────────────────────────────────
# Step 5: Picard MarkDuplicates
# ─────────────────────────────────────────────────────────────────────────────

log "━━━ Step 5: Picard MarkDuplicates ━━━"

mkdir -p "${STAR_DIR}/picard_tmp"
for bam in "${STAR_DIR}"/*Aligned.sortedByCoord.out.bam; do
    [[ -e "${bam}" ]] || { log "  No aligned BAMs found"; break; }
    base="${bam%Aligned.sortedByCoord.out.bam}"
    log "  MarkDuplicates for $(basename "${base}")..."
    picard -Xmx"${MAX_JAVA_MEM}" MarkDuplicates \
        I="${bam}" \
        O="${base}_dedup.bam" \
        M="${base}_dedup_metrics.txt" \
        VALIDATION_STRINGENCY=SILENT \
        TMP_DIR="${STAR_DIR}/picard_tmp"
    samtools index "${base}_dedup.bam"
done

# ─────────────────────────────────────────────────────────────────────────────
# Step 6: GATK variant calling (optional; skip with --skip-gatk)
# ─────────────────────────────────────────────────────────────────────────────

if [[ "${SKIP_GATK}" == "true" ]]; then
    log "━━━ Step 6: GATK skipped (--skip-gatk) ━━━"
else
    log "━━━ Step 6: GATK variant calling ━━━"
    mkdir -p "${GATK_DIR}"

    [[ -f "${GENOME_FA%.fa}.dict" ]] || \
        picard CreateSequenceDictionary R="${GENOME_FA}" O="${GENOME_FA%.fa}.dict"
    [[ -f "${GENOME_FA}.fai" ]] || samtools faidx "${GENOME_FA}"

    for bam in "${STAR_DIR}"/*_dedup.bam; do
        [[ -e "${bam}" ]] || { log "  No dedup BAMs found"; break; }
        base=$(basename "${bam}" _dedup.bam)
        log "  HaplotypeCaller: ${base}"
        gatk --java-options "-Xmx${MAX_JAVA_MEM}" HaplotypeCaller \
            -R "${GENOME_FA}" -I "${bam}" \
            -O "${GATK_DIR}/${base}_.g.vcf.gz" -ERC GVCF
    done

    log "  CombineGVCFs..."
    GVCFs=$(ls "${GATK_DIR}"/*.g.vcf.gz | tr '\n' ' ')
    [[ -n "${GVCFs}" ]] || die "No GVCFs in ${GATK_DIR}"
    gatk --java-options "-Xmx${MAX_JAVA_MEM}" CombineGVCFs \
        -R "${GENOME_FA}" \
        $(printf -- "--variant %s " ${GVCFs}) \
        -O "${GATK_DIR}/combined.g.vcf.gz"

    log "  GenotypeGVCFs..."
    gatk --java-options "-Xmx${MAX_JAVA_MEM}" GenotypeGVCFs \
        -R "${GENOME_FA}" \
        -V "${GATK_DIR}/combined.g.vcf.gz" \
        -O "${GATK_DIR}/joint_variants.vcf.gz"

    log "  VariantFiltration..."
    gatk --java-options "-Xmx${MAX_JAVA_MEM}" VariantFiltration \
        -R "${GENOME_FA}" \
        -V "${GATK_DIR}/joint_variants.vcf.gz" \
        --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || SOR > 3.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
        --filter-name "GATK_RNAseq_filter" \
        -O "${GATK_DIR}/filtered_variants.vcf.gz"
fi

# ─────────────────────────────────────────────────────────────────────────────
# Step 7: MultiQC
# ─────────────────────────────────────────────────────────────────────────────

log "━━━ Step 7: MultiQC ━━━"
multiqc "${PROC_DIR}" -o "${MULTIQC_DIR}"

log "Pipeline complete for ${PAPER_ID}."
log "Results in: ${PROC_DIR}"
