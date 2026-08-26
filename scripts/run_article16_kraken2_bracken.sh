#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat >&2 <<'USAGE'
Usage:
  run_article16_kraken2_bracken.sh \
    --project-root DIR \
    --environment-prefix DIR \
    --raw-dir DIR \
    --database-archive FILE \
    --database-dir DIR \
    --frozen-dir DIR

The clean FASTQ, 5.54 GiB database archive, extracted index, and per-fragment
Kraken output remain outside Git. Only aggregate reports, normalized logs,
resource measurements, and checksum-locked audit tables are frozen.
USAGE
}

project_root=""
environment_prefix=""
raw_dir=""
database_archive=""
database_dir=""
frozen_dir=""

while [[ "$#" -gt 0 ]]; do
  case "$1" in
    --project-root) project_root="$2"; shift 2 ;;
    --environment-prefix) environment_prefix="$2"; shift 2 ;;
    --raw-dir) raw_dir="$2"; shift 2 ;;
    --database-archive) database_archive="$2"; shift 2 ;;
    --database-dir) database_dir="$2"; shift 2 ;;
    --frozen-dir) frozen_dir="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown argument: $1" >&2; usage; exit 2 ;;
  esac
done

if [[ -z "${project_root}" || -z "${environment_prefix}" ||
      -z "${raw_dir}" || -z "${database_archive}" ||
      -z "${database_dir}" || -z "${frozen_dir}" ]]; then
  usage
  exit 2
fi

project_root="$(cd "${project_root}" && pwd)"
environment_prefix="$(cd "${environment_prefix}" && pwd)"
database_archive="$(realpath "${database_archive}")"
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
  echo "Refusing to overwrite non-empty Article 16 work directory: ${work_dir}" >&2
  exit 1
fi
mkdir -p "${work_dir}" "${frozen_dir}/logs"

export LC_ALL=C
export TZ=UTC
export PYTHONHASHSEED=0
export PYTHONNOUSERSITE=1
export MPLCONFIGDIR="${work_dir}/.matplotlib"
export XDG_CACHE_HOME="${work_dir}/.cache"
export PATH="${environment_prefix}/bin:${PATH}"
mkdir -p "${MPLCONFIGDIR}" "${XDG_CACHE_HOME}"

kraken2="${environment_prefix}/bin/kraken2"
kraken2_inspect="${environment_prefix}/bin/kraken2-inspect"
bracken="${environment_prefix}/bin/bracken"
python="${environment_prefix}/bin/python"
for executable in "${kraken2}" "${kraken2_inspect}" "${bracken}" "${python}"; do
  [[ -x "${executable}" ]] || { echo "Missing executable: ${executable}" >&2; exit 1; }
done

clean_r1="${project_root}/data/raw/article13/ERR9765746_clean_R1.fastq.gz"
clean_r2="${project_root}/data/raw/article13/ERR9765746_clean_R2.fastq.gz"
for input_file in "${clean_r1}" "${clean_r2}"; do
  [[ -s "${input_file}" ]] || { echo "Missing Article 13 clean FASTQ: ${input_file}" >&2; exit 1; }
done

observed_bytes="$(stat -c '%s' "${database_archive}")"
[[ "${observed_bytes}" == "5946578575" ]] || {
  echo "Database archive byte-count mismatch: ${observed_bytes}" >&2
  exit 1
}
observed_md5="$(md5sum "${database_archive}" | cut -d ' ' -f 1)"
[[ "${observed_md5}" == "7685f43cce057c2ca18511c925399b72" ]] || {
  echo "Database archive MD5 mismatch: ${observed_md5}" >&2
  exit 1
}

if [[ ! -s "${database_dir}/hash.k2d" ]]; then
  /usr/bin/time -v \
    -o "${frozen_dir}/logs/extract-database.resources.txt" \
    tar -xzf "${database_archive}" -C "${database_dir}" \
    > "${frozen_dir}/logs/extract-database.log" 2>&1
fi

(
  cd "${database_dir}"
  md5sum --check "${project_root}/data/small/16-standard8-files.md5"
) > "${frozen_dir}/logs/database-md5-check.log" 2>&1

cp "${project_root}/data/small/16-standard8-files.md5" \
  "${frozen_dir}/database-files.md5"
cp "$0" "${frozen_dir}/commands.sh"

{
  printf 'Tool\tVersion\tPackageVersion\tExecutable\n'
  printf 'Kraken2\t%s\t2.17.1\t%s\n' \
    "$("${kraken2}" --version | sed -n '1s/^Kraken version //p')" \
    '${KRAKEN_ENV_PREFIX}/bin/kraken2'
  printf 'Bracken\t%s\t3.1p1\t%s\n' \
    "$("${bracken}" -v 2>&1 | sed -n 's/^Bracken v//p')" \
    '${KRAKEN_ENV_PREFIX}/bin/bracken'
  printf 'Python\t%s\t3.12.13\t%s\n' \
    "$("${python}" --version | sed 's/^Python //')" \
    '${KRAKEN_ENV_PREFIX}/bin/python'
} > "${frozen_dir}/tool-versions.tsv"

"${kraken2_inspect}" \
  --db "${database_dir}" \
  --threads 8 \
  --skip-counts \
  > "${frozen_dir}/database-inspect.txt" \
  2> "${frozen_dir}/logs/kraken2-inspect.log"

/usr/bin/time -v \
  -o "${frozen_dir}/logs/kraken2-full.resources.txt" \
  "${kraken2}" \
    --db "${database_dir}" \
    --threads 8 \
    --paired \
    --gzip-compressed \
    --confidence 0.0 \
    --minimum-hit-groups 2 \
    --report "${frozen_dir}/kraken-report.tsv" \
    --output "${work_dir}/ERR9765746.kraken.output" \
    "${clean_r1}" "${clean_r2}" \
  > "${frozen_dir}/logs/kraken2-full.log" 2>&1

run_bracken() {
  local rank="$1"
  local read_length="$2"
  local threshold="$3"
  local label="$4"
  /usr/bin/time -v \
    -o "${frozen_dir}/logs/bracken-${label}.resources.txt" \
    "${bracken}" \
      -d "${database_dir}" \
      -i "${frozen_dir}/kraken-report.tsv" \
      -o "${frozen_dir}/bracken-${label}.tsv" \
      -w "${frozen_dir}/bracken-${label}.kreport.tsv" \
      -r "${read_length}" \
      -l "${rank}" \
      -t "${threshold}" \
    > "${frozen_dir}/logs/bracken-${label}.log" 2>&1
}

run_bracken S 150 10 species-r150-t10
run_bracken S 150 0 species-r150-t0
run_bracken S 100 10 species-r100-t10
run_bracken G 150 10 genus-r150-t10

"${python}" "${project_root}/scripts/validate_article16_kraken2_bracken.py" \
  --project-root "${project_root}" \
  --environment-prefix "${environment_prefix}" \
  --frozen-dir "${frozen_dir}" \
  --raw-dir "${raw_dir}" \
  --work-dir "${work_dir}" \
  --database-archive "${database_archive}" \
  --database-dir "${database_dir}" \
  --initialize-frozen

echo "Article 16 one-time Kraken2/Bracken analysis completed: ${frozen_dir}"
