#!/usr/bin/env bash
# Edit the paths to point to existing local files and directories.
export ASCAT_REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
# Source this file before running the wrappers or ASCAT.sc.
export ASCAT_SAMPLE_IDS="sample_01,sample_02"
export ASCAT_SAMPLE_ID="sample_01"
export ASCAT_CELL_ID="cell_01"
export ASCAT_CELL_IDS="cell_01,cell_02"
export ASCAT_INPUT_ROOT="/path/to/input"
export ASCAT_OUTPUT_ROOT="/path/to/ascat_output"

# Independent ASCAT input roots. {sample} is replaced with each sample ID.
# For a single sample, these can be literal paths without {sample}.
export ASCAT_BAM_DIR="/path/to/{sample}/single_cell_bams"
export ASCAT_CELL_AC_ROOT="/path/to/{sample}/single_cell_allele_counts"
export ASCAT_PHASE_DIR="/path/to/{sample}/bulk_beagle_phasing"
# Set these patterns to match the existing validated BAM names.
export ASCAT_BAM_PATTERN='\.hg38\.sorted\.dedup\.bam$'
export ASCAT_BAM_CELL_SUFFIX='\.dedup\.bam$'
# If your files end in .markdup.bam, instead use:
# export ASCAT_BAM_PATTERN='\.hg38\.sorted\.markdup\.bam$'
# export ASCAT_BAM_CELL_SUFFIX='\.markdup\.bam$'

export ASCAT_SEX="female"
export ASCAT_REGION_CHROM="chr4"
export ASCAT_REGION_START_MB="10"
export ASCAT_REGION_END_MB="30"

export ALIGN_FASTQ_DIR="/path/to/trimmed_fastq"
export ALIGN_REFERENCE="/path/to/reference/GRCh38.fa"
export ALIGN_OUTPUT_DIR="/path/to/aligned_bams"
export MARKDUP_INPUT_DIR="/path/to/aligned_bams"
export MARKDUP_OUTPUT_DIR="/path/to/duplicate_marked_bams"

export ALLELECOUNT_BAM="/path/to/a.bam"
export ALLELECOUNT_LOCI_DIR="/path/to/GRCh38/loci"
export ALLELECOUNTER_EXECUTABLE="alleleCounter"
# Optional: explicit independent outputs, containing the allele_count/ directory.
export BULK_AC_OUTPUT_DIR="/path/to/bulk_allele_counts"
export SINGLE_CELL_AC_OUTPUT_DIR="/path/to/cell_01_allele_counts"

export BEAGLE_CONTAINER="/path/to/beagle_container.sif"
export BEAGLE_INPUT_ROOT="/path/to/bulk_data"
# These two phasing paths may be completely independent.
export BULK_AC_DIR="/path/to/{sample}/bulk_allele_counts"
export BEAGLE_PHASE_OUTPUT_DIR="/path/to/{sample}/bulk_beagle_phasing"
# Regional figure scripts can use their existing directory arrangement, or override root.
export ASCAT_PLOT_AC_ROOT="/path/to/allele_counts_single_cells"
export ASCAT_PLOT_PHASE_DIR="/path/to/bulk_beagle_phasing"
