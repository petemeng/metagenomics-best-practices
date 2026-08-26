#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat >&2 <<'USAGE'
Usage:
  run_article15_metaphlan4.sh \
    --project-root DIR \
    --environment-prefix DIR \
    --raw-dir DIR \
    --metadata-archive FILE \
    --bowtie2-archive FILE \
    --database-dir DIR \
    --frozen-dir DIR

The cleaned FASTQ, database archives, extracted database, paired subsamples,
and read-to-marker mapout files remain in Git-ignored storage. Only profiles,
aggregate evidence, normalized logs, and checksum manifests are frozen.
USAGE
}

project_root=""
environment_prefix=""
raw_dir=""
metadata_archive=""
bowtie2_archive=""
database_dir=""
frozen_dir=""

while [[ "$#" -gt 0 ]]; do
  case "$1" in
    --project-root)
      project_root="$2"
      shift 2
      ;;
    --environment-prefix)
      environment_prefix="$2"
      shift 2
      ;;
    --raw-dir)
      raw_dir="$2"
      shift 2
      ;;
    --metadata-archive)
      metadata_archive="$2"
      shift 2
      ;;
    --bowtie2-archive)
      bowtie2_archive="$2"
      shift 2
      ;;
    --database-dir)
      database_dir="$2"
      shift 2
      ;;
    --frozen-dir)
      frozen_dir="$2"
      shift 2
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "Unknown argument: $1" >&2
      usage
      exit 2
      ;;
  esac
done

if [[ -z "${project_root}" ||
      -z "${environment_prefix}" ||
      -z "${raw_dir}" ||
      -z "${metadata_archive}" ||
      -z "${bowtie2_archive}" ||
      -z "${database_dir}" ||
      -z "${frozen_dir}" ]]; then
  usage
  exit 2
fi

project_root="$(cd "${project_root}" && pwd)"
environment_prefix="$(cd "${environment_prefix}" && pwd)"
metadata_archive="$(realpath "${metadata_archive}")"
bowtie2_archive="$(realpath "${bowtie2_archive}")"
mkdir -p "${raw_dir}" "${database_dir}"
raw_dir="$(cd "${raw_dir}" && pwd)"
database_dir="$(cd "${database_dir}" && pwd)"

if [[ -e "${frozen_dir}" ]] &&
   [[ -n "$(find "${frozen_dir}" -mindepth 1 -maxdepth 1 -print -quit)" ]]; then
  echo "Refusing to overwrite non-empty frozen directory: ${frozen_dir}" >&2
  exit 1
fi
mkdir -p "${frozen_dir}"
frozen_dir="$(cd "${frozen_dir}" && pwd)"

work_dir="${raw_dir}/work"
if [[ -e "${work_dir}" ]] &&
   [[ -n "$(find "${work_dir}" -mindepth 1 -maxdepth 1 -print -quit)" ]]; then
  echo "Refusing to overwrite non-empty Article 15 work directory: ${work_dir}" >&2
  exit 1
fi
mkdir -p "${work_dir}" "${work_dir}/tmp" "${frozen_dir}/logs"

export LC_ALL=C
export TZ=UTC
export MPLCONFIGDIR="${work_dir}/.matplotlib"
export XDG_CACHE_HOME="${work_dir}/.cache"
export PATH="${environment_prefix}/bin:${PATH}"
mkdir -p "${MPLCONFIGDIR}" "${XDG_CACHE_HOME}"

metaphlan="${environment_prefix}/bin/metaphlan"
bowtie2="${environment_prefix}/bin/bowtie2"
python="${environment_prefix}/bin/python"

for executable in "${metaphlan}" "${bowtie2}" "${python}"; do
  if [[ ! -x "${executable}" ]]; then
    echo "Missing executable: ${executable}" >&2
    exit 1
  fi
done

index_name="mpa_vJan26_CHOCOPhlAnSGB_202605"
clean_r1="${project_root}/data/raw/article13/ERR9765746_clean_R1.fastq.gz"
clean_r2="${project_root}/data/raw/article13/ERR9765746_clean_R2.fastq.gz"

for input_file in "${clean_r1}" "${clean_r2}"; do
  if [[ ! -s "${input_file}" ]]; then
    echo "Missing Article 13 clean FASTQ: ${input_file}" >&2
    exit 1
  fi
done

verify_archive() {
  local archive="$1"
  local expected_bytes="$2"
  local expected_md5="$3"
  local observed_bytes observed_md5

  observed_bytes="$(stat -c '%s' "${archive}")"
  if [[ "${observed_bytes}" != "${expected_bytes}" ]]; then
    echo "Archive byte-count mismatch: ${archive}" >&2
    exit 1
  fi
  observed_md5="$(md5sum "${archive}" | cut -d ' ' -f 1)"
  if [[ "${observed_md5}" != "${expected_md5}" ]]; then
    echo "Archive MD5 mismatch: ${archive}" >&2
    exit 1
  fi
}

verify_archive \
  "${metadata_archive}" \
  6014095360 \
  7162b0c3493663dce9abef08ccc06aea
verify_archive \
  "${bowtie2_archive}" \
  41742510080 \
  ac93e1e9c0829629266f5b6ab19c318d

if [[ ! -s "${database_dir}/${index_name}.pkl" ]]; then
  /usr/bin/time -v \
    -o "${frozen_dir}/logs/extract-metadata.resources.txt" \
    tar -xf "${metadata_archive}" -C "${database_dir}" \
    > "${frozen_dir}/logs/extract-metadata.log" 2>&1
fi

index_file_count="$(
  find "${database_dir}" \
    -maxdepth 1 \
    -type f \
    -name "${index_name}*.bt2*" \
    -printf '.' \
    | wc -c
)"
if (( index_file_count == 0 )); then
  /usr/bin/time -v \
    -o "${frozen_dir}/logs/extract-bowtie2.resources.txt" \
    tar -xf "${bowtie2_archive}" -C "${database_dir}" \
    > "${frozen_dir}/logs/extract-bowtie2.log" 2>&1
elif (( index_file_count != 6 )); then
  echo "Incomplete Bowtie2 index in ${database_dir}: ${index_file_count}/6 files" >&2
  exit 1
fi

index_file_count="$(
  find "${database_dir}" \
    -maxdepth 1 \
    -type f \
    -name "${index_name}*.bt2*" \
    -size +0c \
    -printf '.' \
    | wc -c
)"
if (( index_file_count != 6 )); then
  echo "Bowtie2 extraction did not produce six non-empty files" >&2
  exit 1
fi

"${metaphlan}" \
  --install \
  --db_dir "${database_dir}" \
  --index "${index_name}" \
  --offline \
  > "${frozen_dir}/logs/database-offline-check.log" 2>&1

cp "$0" "${frozen_dir}/commands.sh"

{
  printf 'Tool\tVersion\tExecutable\n'
  printf 'MetaPhlAn\t%s\t%s\n' \
    "$("${metaphlan}" --version 2>&1 | sed -n '1s/^MetaPhlAn version \([^ ]*\).*/\1/p')" \
    '${BIOBAKERY_ENV_PREFIX}/bin/metaphlan'
  printf 'Bowtie2\t%s\t%s\n' \
    "$("${bowtie2}" --version 2>&1 | sed -n 's/.* version \([0-9.]*\)$/\1/p' | head -n 1)" \
    '${BIOBAKERY_ENV_PREFIX}/bin/bowtie2'
  printf 'Python\t%s\t%s\n' \
    "$("${python}" --version | sed 's/^Python //')" \
    '${BIOBAKERY_ENV_PREFIX}/bin/python'
} > "${frozen_dir}/tool-versions.tsv"

"${python}" "${project_root}/scripts/validate_article15_metaphlan4.py" \
  --project-root "${project_root}" \
  --environment-prefix "${environment_prefix}" \
  --frozen-dir "${frozen_dir}" \
  --raw-dir "${raw_dir}" \
  --work-dir "${work_dir}" \
  --metadata-archive "${metadata_archive}" \
  --bowtie2-archive "${bowtie2_archive}" \
  --database-dir "${database_dir}" \
  --preflight-only

main_mapout="${work_dir}/ERR9765746-full.mapout.bz2"

/usr/bin/time -v \
  -o "${frozen_dir}/logs/metaphlan-full.resources.txt" \
  "${metaphlan}" \
  "${clean_r1},${clean_r2}" \
  --input_type fastq \
  --db_dir "${database_dir}" \
  --index "${index_name}" \
  --offline \
  --bowtie2_exe "${bowtie2}" \
  --mapout "${main_mapout}" \
  --tmp_dir "${work_dir}/tmp" \
  --nproc 8 \
  --read_min_len 70 \
  --perc_nonzero 0.33 \
  --stat tavg_g \
  --stat_q 0.2 \
  -t rel_ab_w_read_stats \
  --tax_lev a \
  --sample_id ERR9765746_MOCK1 \
  --output_file "${frozen_dir}/profile-all.tsv" \
  > "${frozen_dir}/logs/metaphlan-full.log" 2>&1

# MetaPhlAn 4.2.5 validates that --mapout is present even for mapout input.
# The non-existing path below satisfies that validation but is not written.
classified_dummy_mapout="${work_dir}/classified-reprofile-unused.mapout"
/usr/bin/time -v \
  -o "${frozen_dir}/logs/metaphlan-classified-only.resources.txt" \
  "${metaphlan}" \
  "${main_mapout}" \
  --input_type mapout \
  --db_dir "${database_dir}" \
  --index "${index_name}" \
  --offline \
  --mapout "${classified_dummy_mapout}" \
  --nproc 8 \
  --perc_nonzero 0.33 \
  --stat tavg_g \
  --stat_q 0.2 \
  -t rel_ab_w_read_stats \
  --tax_lev a \
  --skip_unclassified_estimation \
  --sample_id ERR9765746_MOCK1_CLASSIFIED_ONLY \
  --output_file "${frozen_dir}/profile-classified-only.tsv" \
  > "${frozen_dir}/logs/metaphlan-classified-only.log" 2>&1

for pair_depth in 20000 50000; do
  paired_read_argument="$((pair_depth * 2))"
  /usr/bin/time -v \
    -o "${frozen_dir}/logs/metaphlan-depth-${pair_depth}.resources.txt" \
    "${metaphlan}" \
    --input_type fastq \
    -1 "${clean_r1}" \
    -2 "${clean_r2}" \
    --subsampling_paired "${paired_read_argument}" \
    --subsampling_seed 20260721 \
    --subsampling_output "${work_dir}/subsample-${pair_depth}.fastq.gz" \
    --db_dir "${database_dir}" \
    --index "${index_name}" \
    --offline \
    --bowtie2_exe "${bowtie2}" \
    --mapout "${work_dir}/ERR9765746-depth-${pair_depth}.mapout.bz2" \
    --tmp_dir "${work_dir}/tmp" \
    --nproc 8 \
    --read_min_len 70 \
    --perc_nonzero 0.33 \
    --stat tavg_g \
    --stat_q 0.2 \
    -t rel_ab_w_read_stats \
    --tax_lev a \
    --sample_id "ERR9765746_MOCK1_${pair_depth}_PAIRS" \
    --output_file "${work_dir}/profile-depth-${pair_depth}.tsv" \
    > "${frozen_dir}/logs/metaphlan-depth-${pair_depth}.log" 2>&1
done

"${python}" "${project_root}/scripts/validate_article15_metaphlan4.py" \
  --project-root "${project_root}" \
  --environment-prefix "${environment_prefix}" \
  --frozen-dir "${frozen_dir}" \
  --raw-dir "${raw_dir}" \
  --work-dir "${work_dir}" \
  --metadata-archive "${metadata_archive}" \
  --bowtie2-archive "${bowtie2_archive}" \
  --database-dir "${database_dir}" \
  --initialize-frozen

echo "Article 15 one-time MetaPhlAn analysis completed: ${frozen_dir}"
