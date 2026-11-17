#!/bin/bash -l
#$ -N single_cell_alleles
#$ -cwd
#$ -pe smp 1
#$ -l mem=12G
#$ -l h_rt=24:00:00
set -euo pipefail
: "${ASCAT_SAMPLE_ID:?Set ASCAT_SAMPLE_ID}"
: "${ASCAT_CELL_ID:?Set ASCAT_CELL_ID}"
: "${ALLELECOUNT_BAM:?Set ALLELECOUNT_BAM to the single-cell BAM}"
: "${ALLELECOUNT_LOCI_DIR:?Set ALLELECOUNT_LOCI_DIR}"
: "${ASCAT_INPUT_ROOT:?Set ASCAT_INPUT_ROOT}"
Rscript "${ASCAT_REPO_ROOT:?Source config/example_environment.sh before submission}/02_allele_counts/count_alleles.R" \
  --bam "$ALLELECOUNT_BAM" \
  --loci-dir "$ALLELECOUNT_LOCI_DIR" \
  --output-dir "${SINGLE_CELL_AC_OUTPUT_DIR:-$ASCAT_INPUT_ROOT/$ASCAT_SAMPLE_ID/bam_to_beagle_output/$ASCAT_CELL_ID}" \
  --sample-id "$ASCAT_CELL_ID" \
  --allele-counter "${ALLELECOUNTER_EXECUTABLE:-alleleCounter}"
