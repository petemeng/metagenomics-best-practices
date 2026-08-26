#!/usr/bin/env bash
set -euo pipefail
export MICROBECENSUS_ROOT=${MICROBECENSUS_ROOT}
export BIOBAKERY_ENV_PREFIX=${BIOBAKERY_ENV_PREFIX}
PROJECT_ROOT=${PROJECT_ROOT} \
  ARTICLE21_INPUT_R1=${PROJECT_ROOT}/data/raw/article13/ERR9765746_clean_R1.fastq.gz \
  ARTICLE21_INPUT_R2=${PROJECT_ROOT}/data/raw/article13/ERR9765746_clean_R2.fastq.gz \
  bash ${PROJECT_ROOT}/scripts/run_article21_microbecensus.sh
