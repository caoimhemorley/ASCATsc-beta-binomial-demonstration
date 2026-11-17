#!/bin/bash -l
#$ -N confident_beagle
#$ -cwd
#$ -tc 4
#$ -pe smp 4
#$ -l mem=16G
#$ -l h_rt=24:00:00
set -euo pipefail
: "${BEAGLE_CONTAINER:?Set BEAGLE_CONTAINER}"
: "${BEAGLE_INPUT_ROOT:?Set BEAGLE_INPUT_ROOT}"
IFS="," read -r -a DONORS <<< "${ASCAT_SAMPLE_IDS:?Set ASCAT_SAMPLE_IDS}"
CHROMS=({1..22} X)
task=$(( ${SGE_TASK_ID:?Submit as an SGE array} - 1 ))
if (( task < 0 || task >= ${#DONORS[@]} * 23 )); then
  echo "Invalid SGE_TASK_ID: ${SGE_TASK_ID}" >&2
  exit 1
fi
donor=${DONORS[$((task / 23))]}
chrom=${CHROMS[$((task % 23))]}
bulk_ac_root="${BULK_AC_DIR:-$BEAGLE_INPUT_ROOT/{sample}/bam_to_beagle_output}"
phase_root="${BEAGLE_PHASE_OUTPUT_DIR:-$BEAGLE_INPUT_ROOT/{sample}/bam_to_beagle_output/beagle_confident_v1}"
bulk_ac_root="${bulk_ac_root//\{sample\}/$donor}"
phase_root="${phase_root//\{sample\}/$donor}"
sex_var="ASCAT_SEX_${donor}"
sex="${!sex_var:-${ASCAT_SEX:-male}}"
script_dir="${ASCAT_REPO_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}/02_phasing"
python3 "$script_dir/phase_bulk_confident_SNPs.py" \
  --donor "$donor" --chrom "$chrom" \
  --sex "$sex" \
  --threads "${NSLOTS:-4}" \
  --input-root "$bulk_ac_root" \
  --out-root "$phase_root" \
  --loci-dir "$ALLELECOUNT_LOCI_DIR" \
  --beagle-container "$BEAGLE_CONTAINER"
