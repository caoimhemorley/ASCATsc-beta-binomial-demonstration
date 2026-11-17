#!/bin/bash -l
#$ -N allele_counts
#$ -cwd
#$ -pe smp 1
#$ -l mem=12G
#$ -l h_rt=24:00:00
set -euo pipefail
: "${ASCAT_SAMPLE_ID:?Set ASCAT_SAMPLE_ID}"
: "${ALLELECOUNT_BAM:?Set ALLELECOUNT_BAM}"
: "${ALLELECOUNT_LOCI_DIR:?Set ALLELECOUNT_LOCI_DIR}"
: "${BEAGLE_INPUT_ROOT:?Set BEAGLE_INPUT_ROOT}"
Rscript "${ASCAT_REPO_ROOT:?Source config/example_environment.sh before submission}/02_allele_counts/count_alleles.R" \
  --bam "$ALLELECOUNT_BAM" \
  --loci-dir "$ALLELECOUNT_LOCI_DIR" \
  --output-dir "${BULK_AC_OUTPUT_DIR:-$BEAGLE_INPUT_ROOT/$ASCAT_SAMPLE_ID/bam_to_beagle_output}" \
  --sample-id "$ASCAT_SAMPLE_ID" \
  --sex "${ASCAT_SEX:-female}" \
  --allele-counter "${ALLELECOUNTER_EXECUTABLE:-alleleCounter}"
