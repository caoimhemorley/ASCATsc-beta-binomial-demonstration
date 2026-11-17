#!/bin/bash -l
#$ -N ASCATsc_cohort
#$ -cwd
#$ -tc 2
#$ -pe smp 8
#$ -l mem=8G
#$ -l h_rt=48:00:00
set -euo pipefail
Rscript "${ASCAT_REPO_ROOT:?Source config/example_environment.sh before submission}/03_ascat/run_ASCATsc_cohort.R"
