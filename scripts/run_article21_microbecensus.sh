#!/usr/bin/env bash
set -euo pipefail

project_root="${PROJECT_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
raw_dir="${ARTICLE21_RAW_DIR:-${project_root}/data/raw/article21}"
frozen_dir="${ARTICLE21_FROZEN_DIR:-${project_root}/data/small/21-table-semantics-frozen}"
source_root="${MICROBECENSUS_ROOT:?Set MICROBECENSUS_ROOT to the official v1.1.1 checkout}"
python_bin="${BIOBAKERY_ENV_PREFIX:-${HOME}/miniconda3/envs/metagenome-biobakery-2026.07}/bin/python"
input_r1="${ARTICLE21_INPUT_R1:-${project_root}/data/raw/article13/ERR9765746_clean_R1.fastq.gz}"
input_r2="${ARTICLE21_INPUT_R2:-${project_root}/data/raw/article13/ERR9765746_clean_R2.fastq.gz}"

mkdir -p "${raw_dir}" "${frozen_dir}" "${raw_dir}/tmp"

commit="$(git -C "${source_root}" rev-parse HEAD)"
test "${commit}" = "dfc42d356bfd7943633cde6c0fbfc0b116f29ae2"
test -f "${input_r1}"
test -f "${input_r2}"

run_branch() {
  local read_length="$1"
  local branch="read-length-${read_length}"
  TMPDIR="${raw_dir}/tmp" PYTHONHASHSEED=0 \
    /usr/bin/time -v \
    -o "${raw_dir}/${branch}.resource.txt" \
    "${python_bin}" "${project_root}/scripts/run_article21_microbecensus.py" \
      --source-root "${source_root}" \
      --input-r1 "${input_r1}" \
      --input-r2 "${input_r2}" \
      --output "${raw_dir}/${branch}.output.tsv" \
      --metadata-output "${raw_dir}/${branch}.metadata.json" \
      --read-length "${read_length}" \
      --threads 8 \
      --nreads 100000000 \
      --verbose \
      > "${raw_dir}/${branch}.log" 2>&1
}

run_branch 150
run_branch 100

"${python_bin}" "${project_root}/scripts/freeze_article21_microbecensus.py" \
  --project-root "${project_root}" \
  --raw-dir "${raw_dir}" \
  --source-root "${source_root}" \
  --frozen-dir "${frozen_dir}"
