#!/usr/bin/env bash
set -euo pipefail

THREADS=4
BASE_DIR=${1:-$(pwd)}

OUTDIR="${OUTPUT_DIR:-fastp_50bp}"
CLEANDIR="${OUTDIR}/clean"
FINALDIR="${OUTDIR}/final_50bp"
LOGDIR="${OUTDIR}/logs"
HTMLDIR="${OUTDIR}/qc_html"

mkdir -p "$CLEANDIR" "$FINALDIR" "$LOGDIR" "$HTMLDIR"

echo "======================================"
echo " fastp + hard 50bp truncation pipeline"
echo "======================================"

clean_reads() {

    R1="$1"
    R2="$2"

    SAMPLE=$(basename "$R1")
    SAMPLE="${SAMPLE%_R1_001.fastq.gz}"
    SAMPLE="${SAMPLE%_R1.fastq.gz}"

    CLEAN_R1="${CLEANDIR}/${SAMPLE}_R1.clean.fq.gz"
    CLEAN_R2="${CLEANDIR}/${SAMPLE}_R2.clean.fq.gz"

    FINAL_R1="${FINALDIR}/${SAMPLE}_R1_50bp.fq.gz"
    FINAL_R2="${FINALDIR}/${SAMPLE}_R2_50bp.fq.gz"

    HTML="${HTMLDIR}/${SAMPLE}.html"
    JSON="${HTMLDIR}/${SAMPLE}.json"
    LOG="${LOGDIR}/${SAMPLE}.log"

    echo "[$SAMPLE] fastp cleaning..."

    fastp \
        -i "$R1" -I "$R2" \
        -o "$CLEAN_R1" -O "$CLEAN_R2" \
        --detect_adapter_for_pe \
        --trim_poly_g \
        --cut_right \
        --cut_right_window_size 4 \
        --cut_right_mean_quality 20 \
        --qualified_quality_phred 20 \
        --length_required 1 \
        --thread "$THREADS" \
        --html "$HTML" \
        --json "$JSON" \
        > "$LOG" 2>&1

    

  echo "[$SAMPLE] extracting bases 36-85 (50 bp) after fastp..."

seqkit subseq -r 36:85 "$CLEAN_R1" | gzip > "$FINAL_R1"
seqkit subseq -r 36:85 "$CLEAN_R2" | gzip > "$FINAL_R2"

    echo "[$SAMPLE] DONE → $FINALDIR"
}

export -f clean_reads
export CLEANDIR FINALDIR LOGDIR HTMLDIR THREADS

while IFS= read -r R1; do
    if [[ "$R1" == *_R1_001.fastq.gz ]]; then
        R2="${R1%_R1_001.fastq.gz}_R2_001.fastq.gz"
    else
        R2="${R1%_R1.fastq.gz}_R2.fastq.gz"
    fi
    [[ -f "$R2" ]] && clean_reads "$R1" "$R2"
done < <(find "$BASE_DIR" -type f \( -name "*_R1_001.fastq.gz" -o -name "*_R1.fastq.gz" \) | sort)

echo "Running MultiQC..."
multiqc "$OUTDIR" -o "$OUTDIR"

echo "✅ ALL DONE"
