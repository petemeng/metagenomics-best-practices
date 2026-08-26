#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat >&2 <<'USAGE'
Usage:
  run_article18_profiler_benchmark.sh \
    --project-root DIR \
    --environment-prefix DIR \
    --work-dir DIR \
    --motus-database-root DIR \
    --motus-database-archive FILE \
    --metaphlan-database-pkl FILE \
    --frozen-dir DIR

The one-time builder maps the Article 13 MOCK1 clean FASTQ once, calculates a
single MGC table, derives mOTUs profiles at g=1/3/6, and asks the validator to
build the strict NCBI-species/SGB/mOTU crosswalk. FASTQ, BAM, MGC, database
archives and indexes stay outside the frozen Git payload.
USAGE
}

project_root=""
environment_prefix=""
work_dir=""
motus_database_root=""
motus_database_archive=""
metaphlan_database_pkl=""
frozen_dir=""

while [[ "$#" -gt 0 ]]; do
  case "$1" in
    --project-root) project_root="$2"; shift 2 ;;
    --environment-prefix) environment_prefix="$2"; shift 2 ;;
    --work-dir) work_dir="$2"; shift 2 ;;
    --motus-database-root) motus_database_root="$2"; shift 2 ;;
    --motus-database-archive) motus_database_archive="$2"; shift 2 ;;
    --metaphlan-database-pkl) metaphlan_database_pkl="$2"; shift 2 ;;
    --frozen-dir) frozen_dir="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown argument: $1" >&2; usage; exit 2 ;;
  esac
done

if [[ -z "${project_root}" || -z "${environment_prefix}" ||
      -z "${work_dir}" || -z "${motus_database_root}" ||
      -z "${motus_database_archive}" || -z "${metaphlan_database_pkl}" ||
      -z "${frozen_dir}" ]]; then
  usage
  exit 2
fi

project_root="$(cd "${project_root}" && pwd)"
environment_prefix="$(cd "${environment_prefix}" && pwd)"
motus_database_root="$(cd "${motus_database_root}" && pwd)"
motus_database_archive="$(realpath "${motus_database_archive}")"
metaphlan_database_pkl="$(realpath "${metaphlan_database_pkl}")"
mkdir -p "${work_dir}" "${frozen_dir}"
work_dir="$(cd "${work_dir}" && pwd)"
frozen_dir="$(cd "${frozen_dir}" && pwd)"

if [[ -s "${frozen_dir}/run-summary.json" ]]; then
  echo "Refusing to overwrite completed frozen evidence: ${frozen_dir}" >&2
  exit 1
fi

mkdir -p \
  "${work_dir}/logs" \
  "${work_dir}/resources" \
  "${work_dir}/.matplotlib" \
  "${work_dir}/.cache" \
  "${frozen_dir}/logs"

export LC_ALL=C
export TZ=UTC
export PYTHONHASHSEED=0
export PYTHONNOUSERSITE=1
export MPLCONFIGDIR="${work_dir}/.matplotlib"
export XDG_CACHE_HOME="${work_dir}/.cache"
export PATH="${environment_prefix}/bin:${PATH}"

motus="${environment_prefix}/bin/motus"
python="${environment_prefix}/bin/python"
for executable in "${motus}" "${python}" "${environment_prefix}/bin/bwa"; do
  [[ -x "${executable}" ]] || {
    echo "Missing executable: ${executable}" >&2
    exit 1
  }
done

r1="${project_root}/data/raw/article13/ERR9765746_clean_R1.fastq.gz"
r2="${project_root}/data/raw/article13/ERR9765746_clean_R2.fastq.gz"
database_dir="${motus_database_root}/db_mOTU"
for input_file in \
  "${r1}" "${r2}" "${motus_database_archive}" \
  "${metaphlan_database_pkl}" "${database_dir}/mOTUsv4.1.db"; do
  [[ -s "${input_file}" ]] || {
    echo "Missing Article 18 input: ${input_file}" >&2
    exit 1
  }
done

archive_bytes="$(stat -c '%s' "${motus_database_archive}")"
archive_md5="$(md5sum "${motus_database_archive}" | awk '{print $1}')"
[[ "${archive_bytes}" == "5552983255" ]] || {
  echo "mOTUs archive byte-count mismatch: ${archive_bytes}" >&2
  exit 1
}
[[ "${archive_md5}" == "471ea128f0c0839f5c4629b949ea5f8a" ]] || {
  echo "mOTUs archive MD5 mismatch: ${archive_md5}" >&2
  exit 1
}

: > "${work_dir}/logs/database-md5-check.log"
while IFS=: read -r filename expected_md5; do
  filename="${filename%$'\t'}"
  expected_md5="${expected_md5//$'\t'/}"
  [[ "${filename}" == *.gz ]] || continue
  observed_md5="$(md5sum "${database_dir}/${filename}" | awk '{print $1}')"
  if [[ "${observed_md5}" != "${expected_md5}" ]]; then
    echo "FAIL ${filename} ${observed_md5}" >> "${work_dir}/logs/database-md5-check.log"
    exit 1
  fi
  echo "PASS ${filename} ${observed_md5}" >> "${work_dir}/logs/database-md5-check.log"
done < "${database_dir}/mOTUsv4.1.db"

bam="${work_dir}/ERR9765746.motus.bam"
mgc="${work_dir}/ERR9765746.mgc"

/usr/bin/time -v \
  -o "${work_dir}/resources/motus-map-tax.resources.txt" \
  "${motus}" map_tax \
    -f "${r1}" \
    -r "${r2}" \
    -db "${motus_database_root}" \
    -l 75 \
    -t 8 \
    -o "${bam}" \
  > "${work_dir}/logs/motus-map-tax.log" 2>&1

/usr/bin/time -v \
  -o "${work_dir}/resources/motus-calc-mgc.resources.txt" \
  "${motus}" calc_mgc \
    -i "${bam}" \
    -db "${motus_database_root}" \
    -l 75 \
    -o "${mgc}" \
  > "${work_dir}/logs/motus-calc-mgc.log" 2>&1

for minimum_markers in 1 3 6; do
  /usr/bin/time -v \
    -o "${work_dir}/resources/motus-calc-motu-g${minimum_markers}.resources.txt" \
    "${motus}" calc_motu \
      -i "${mgc}" \
      -db "${motus_database_root}" \
      -n ERR9765746_MOCK1 \
      -g "${minimum_markers}" \
      -y INSERT_SCALED \
      -o "${work_dir}/motus-profile-g${minimum_markers}.tsv" \
    > "${work_dir}/logs/motus-calc-motu-g${minimum_markers}.log" 2>&1
  cp \
    "${work_dir}/motus-profile-g${minimum_markers}.tsv" \
    "${frozen_dir}/motus-profile-g${minimum_markers}.tsv"
done

cp "$0" "${frozen_dir}/commands.sh"
cp "${work_dir}/logs/"*.log "${frozen_dir}/logs/"
cp "${work_dir}/resources/"*.txt "${frozen_dir}/logs/"

{
  printf 'Tool\tVersion\tPackageVersion\tExecutable\n'
  printf 'mOTUs\t4.1.0\t4.1.0\t%s\n' '${MOTUS_ENV_PREFIX}/bin/motus'
  printf 'BWA\t0.7.19-r1273\t0.7.19\t%s\n' '${MOTUS_ENV_PREFIX}/bin/bwa'
  printf 'VSEARCH\t2.31.0\t2.31.0\t%s\n' '${MOTUS_ENV_PREFIX}/bin/vsearch'
  printf 'Python\t3.12.13\t3.12.13\t%s\n' '${MOTUS_ENV_PREFIX}/bin/python'
  printf 'MetaPhlAn\t4.2.5\treused Article 15\tchecksum-locked frozen profile\n'
  printf 'Kraken2\t2.17.1\treused Articles 16-17\tchecksum-locked frozen profile\n'
  printf 'Bracken\t3.0.1\t3.1p1\tchecksum-locked frozen profile\n'
} > "${frozen_dir}/tool-versions.tsv"

"${python}" "${project_root}/scripts/validate_article18_profiler_benchmark.py" \
  --project-root "${project_root}" \
  --environment-prefix "${environment_prefix}" \
  --frozen-dir "${frozen_dir}" \
  --work-dir "${work_dir}" \
  --motus-database-root "${motus_database_root}" \
  --motus-database-archive "${motus_database_archive}" \
  --metaphlan-database-pkl "${metaphlan_database_pkl}" \
  --initialize-frozen

echo "PASS Article 18 frozen evidence: ${frozen_dir}"
