#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat >&2 <<'USAGE'
Usage:
  run_article17_kraken2_database_confidence.sh \
    --project-root DIR \
    --environment-prefix DIR \
    --raw-dir DIR \
    --database-root DIR \
    --frozen-dir DIR

The script is restart-safe at the branch level. FASTQ, database archives and
indexes, and per-fragment Kraken output remain outside Git. Standard reports,
Bracken tables, aggregate audits, normalized commands and checksums are frozen.
USAGE
}

project_root=""
environment_prefix=""
raw_dir=""
database_root=""
frozen_dir=""

while [[ "$#" -gt 0 ]]; do
  case "$1" in
    --project-root) project_root="$2"; shift 2 ;;
    --environment-prefix) environment_prefix="$2"; shift 2 ;;
    --raw-dir) raw_dir="$2"; shift 2 ;;
    --database-root) database_root="$2"; shift 2 ;;
    --frozen-dir) frozen_dir="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown argument: $1" >&2; usage; exit 2 ;;
  esac
done

if [[ -z "${project_root}" || -z "${environment_prefix}" ||
      -z "${raw_dir}" || -z "${database_root}" || -z "${frozen_dir}" ]]; then
  usage
  exit 2
fi

project_root="$(cd "${project_root}" && pwd)"
environment_prefix="$(cd "${environment_prefix}" && pwd)"
mkdir -p "${raw_dir}" "${database_root}" "${frozen_dir}"
raw_dir="$(cd "${raw_dir}" && pwd)"
database_root="$(cd "${database_root}" && pwd)"
frozen_dir="$(cd "${frozen_dir}" && pwd)"

if [[ -s "${frozen_dir}/run-summary.json" ]]; then
  echo "Refusing to overwrite completed frozen evidence: ${frozen_dir}" >&2
  exit 1
fi

work_dir="${raw_dir}/work"
mkdir -p \
  "${work_dir}/per-fragment" \
  "${work_dir}/logs" \
  "${work_dir}/resources" \
  "${work_dir}/sentinels" \
  "${work_dir}/.matplotlib" \
  "${work_dir}/.cache" \
  "${frozen_dir}/reports" \
  "${frozen_dir}/bracken"

export LC_ALL=C
export TZ=UTC
export PYTHONHASHSEED=0
export PYTHONNOUSERSITE=1
export MPLCONFIGDIR="${work_dir}/.matplotlib"
export XDG_CACHE_HOME="${work_dir}/.cache"
export PATH="${environment_prefix}/bin:${PATH}"

kraken2="${environment_prefix}/bin/kraken2"
bracken="${environment_prefix}/bin/bracken"
python="${environment_prefix}/bin/python"
for executable in "${kraken2}" "${bracken}" "${python}"; do
  [[ -x "${executable}" ]] || { echo "Missing executable: ${executable}" >&2; exit 1; }
done

mock_r1="${project_root}/data/raw/article13/ERR9765746_clean_R1.fastq.gz"
mock_r2="${project_root}/data/raw/article13/ERR9765746_clean_R2.fastq.gz"
human_r1="${project_root}/data/raw/article14/ERR194147_prefix20k_R1.fastq.gz"
human_r2="${project_root}/data/raw/article14/ERR194147_prefix20k_R2.fastq.gz"
for input_file in "${mock_r1}" "${mock_r2}" "${human_r1}" "${human_r2}"; do
  [[ -s "${input_file}" ]] || { echo "Missing Article 17 input: ${input_file}" >&2; exit 1; }
done

declare -A db_archive=(
  [standard8]="${database_root}/archives/kraken2-standard8-20260626/k2_standard_08_GB_20260626.tar.gz"
  [standard16]="${database_root}/archives/kraken2-standard16-20260626/k2_standard_16_GB_20260626.tar.gz"
  [pluspf8]="${database_root}/archives/kraken2-pluspf8-20260626/k2_pluspf_08_GB_20260626.tar.gz"
)
declare -A db_dir=(
  [standard8]="${database_root}/standard-8-20260626"
  [standard16]="${database_root}/standard-16-20260626"
  [pluspf8]="${database_root}/pluspf-8-20260626"
)
declare -A db_bytes=(
  [standard8]="5946578575"
  [standard16]="11995707291"
  [pluspf8]="5933654083"
)
declare -A db_md5=(
  [standard8]="7685f43cce057c2ca18511c925399b72"
  [standard16]="f130daa49fd0befa688330b288a623de"
  [pluspf8]="79a153b99f045bc2ae95e6d57c17a02d"
)
declare -A db_files_md5=(
  [standard8]="${project_root}/data/small/16-standard8-files.md5"
  [standard16]="${project_root}/data/small/17-standard16-files.md5"
  [pluspf8]="${project_root}/data/small/17-pluspf8-files.md5"
)

for database in standard8 standard16 pluspf8; do
  archive="${db_archive[${database}]}"
  target_dir="${db_dir[${database}]}"
  [[ -s "${archive}" ]] || { echo "Missing database archive: ${archive}" >&2; exit 1; }
  observed_bytes="$(stat -c '%s' "${archive}")"
  observed_md5="$(md5sum "${archive}" | awk '{print $1}')"
  [[ "${observed_bytes}" == "${db_bytes[${database}]}" ]] || {
    echo "Archive byte-count mismatch for ${database}: ${observed_bytes}" >&2
    exit 1
  }
  [[ "${observed_md5}" == "${db_md5[${database}]}" ]] || {
    echo "Archive MD5 mismatch for ${database}: ${observed_md5}" >&2
    exit 1
  }
  if [[ ! -s "${target_dir}/hash.k2d" ]]; then
    if [[ -e "${target_dir}" ]] &&
       [[ -n "$(find "${target_dir}" -mindepth 1 -maxdepth 1 -print -quit)" ]]; then
      echo "Refusing to extract into a non-empty incomplete directory: ${target_dir}" >&2
      exit 1
    fi
    mkdir -p "${target_dir}"
    /usr/bin/time -v \
      -o "${work_dir}/resources/extract-${database}.txt" \
      tar -xzf "${archive}" -C "${target_dir}" \
      > "${work_dir}/logs/extract-${database}.log" 2>&1
  fi
  (
    cd "${target_dir}"
    md5sum --check "${db_files_md5[${database}]}"
  ) > "${work_dir}/logs/database-md5-${database}.log" 2>&1
done

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
  printf 'NCBI Datasets\t18.33.1\tone-time metadata snapshot\tnot required for routine QA\n'
} > "${frozen_dir}/tool-versions.tsv"

branch_design="${work_dir}/branch-design.tsv"
printf 'BranchID\tDatabaseID\tControlID\tConfidence\tMinimumHitGroups\tBrackenReadLength\tPrimaryMatrix\tHitGroupMatrix\n' > "${branch_design}"

confidence_label() {
  case "$1" in
    0.0) printf 'c000\n' ;;
    0.05) printf 'c005\n' ;;
    0.10) printf 'c010\n' ;;
    0.20) printf 'c020\n' ;;
    0.50) printf 'c050\n' ;;
    *) echo "Unsupported confidence: $1" >&2; return 2 ;;
  esac
}

verify_sentinel() {
  local sentinel="$1"
  local report="$2"
  local output="$3"
  local expected_lines="$4"
  [[ -s "${sentinel}" && -s "${report}" && -s "${output}" ]] || return 1
  [[ "$(wc -l < "${output}")" == "${expected_lines}" ]] || return 1
  local expected_report expected_output
  expected_report="$(awk -F '\t' '$1 == "report_sha256" {print $2}' "${sentinel}")"
  expected_output="$(awk -F '\t' '$1 == "output_sha256" {print $2}' "${sentinel}")"
  [[ "$(sha256sum "${report}" | awk '{print $1}')" == "${expected_report}" ]] || return 1
  [[ "$(sha256sum "${output}" | awk '{print $1}')" == "${expected_output}" ]] || return 1
}

run_branch() {
  local database="$1"
  local control="$2"
  local confidence="$3"
  local hit_groups="$4"
  local primary="$5"
  local hit_matrix="$6"
  local r1 r2 read_length expected_lines conf_tag branch

  if [[ "${control}" == "mock" ]]; then
    r1="${mock_r1}"
    r2="${mock_r2}"
    read_length="150"
    expected_lines="99991"
  else
    r1="${human_r1}"
    r2="${human_r2}"
    read_length="100"
    expected_lines="20000"
  fi
  conf_tag="$(confidence_label "${confidence}")"
  branch="${database}-${control}-${conf_tag}-h${hit_groups}"
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "${branch}" "${database}" "${control}" "${confidence}" "${hit_groups}" \
    "${read_length}" "${primary}" "${hit_matrix}" >> "${branch_design}"

  local report output sentinel tmp_report tmp_output tmp_resource tmp_log
  report="${frozen_dir}/reports/${branch}.kreport.tsv"
  output="${work_dir}/per-fragment/${branch}.output"
  sentinel="${work_dir}/sentinels/${branch}.kraken.done.tsv"
  if verify_sentinel "${sentinel}" "${report}" "${output}" "${expected_lines}"; then
    echo "SKIP verified Kraken branch ${branch}"
  else
    tmp_report="${report}.tmp.$$"
    tmp_output="${output}.tmp.$$"
    tmp_resource="${work_dir}/resources/${branch}.kraken.txt.tmp.$$"
    tmp_log="${work_dir}/logs/${branch}.kraken.log.tmp.$$"
    /usr/bin/time -v \
      -o "${tmp_resource}" \
      "${kraken2}" \
        --db "${db_dir[${database}]}" \
        --threads 8 \
        --paired \
        --gzip-compressed \
        --confidence "${confidence}" \
        --minimum-hit-groups "${hit_groups}" \
        --report "${tmp_report}" \
        --output "${tmp_output}" \
        "${r1}" "${r2}" \
      > "${tmp_log}" 2>&1
    [[ "$(wc -l < "${tmp_output}")" == "${expected_lines}" ]] || {
      echo "Per-fragment row-count mismatch for ${branch}" >&2
      exit 1
    }
    mv "${tmp_report}" "${report}"
    mv "${tmp_output}" "${output}"
    mv "${tmp_resource}" "${work_dir}/resources/${branch}.kraken.txt"
    mv "${tmp_log}" "${work_dir}/logs/${branch}.kraken.log"
    {
      printf 'branch_id\t%s\n' "${branch}"
      printf 'report_sha256\t%s\n' "$(sha256sum "${report}" | awk '{print $1}')"
      printf 'output_sha256\t%s\n' "$(sha256sum "${output}" | awk '{print $1}')"
      printf 'output_rows\t%s\n' "${expected_lines}"
    } > "${sentinel}.tmp.$$"
    mv "${sentinel}.tmp.$$" "${sentinel}"
  fi

  local bracken_table bracken_report bracken_sentinel
  bracken_table="${frozen_dir}/bracken/${branch}.tsv"
  bracken_report="${frozen_dir}/bracken/${branch}.kreport.tsv"
  bracken_sentinel="${work_dir}/sentinels/${branch}.bracken.done.tsv"
  local recorded_outcome recorded_report_sha
  recorded_outcome="$(awk -F '\t' '$1 == "outcome" {print $2}' "${bracken_sentinel}" 2>/dev/null || true)"
  [[ -n "${recorded_outcome}" ]] || recorded_outcome="estimated"
  recorded_report_sha="$(awk -F '\t' '$1 == "report_sha256" {print $2}' "${bracken_sentinel}" 2>/dev/null || true)"
  if [[ -s "${bracken_sentinel}" && -s "${bracken_table}" ]] &&
     [[ "$(sha256sum "${bracken_table}" | awk '{print $1}')" == "$(awk -F '\t' '$1 == "table_sha256" {print $2}' "${bracken_sentinel}")" ]] &&
     { [[ "${recorded_outcome}" == "no_eligible_species_above_threshold" &&
          "${recorded_report_sha}" == "not_created" && ! -e "${bracken_report}" ]] ||
       [[ "${recorded_outcome}" == "estimated" && -s "${bracken_report}" &&
          "$(sha256sum "${bracken_report}" | awk '{print $1}')" == "${recorded_report_sha}" ]]; }; then
    echo "SKIP verified Bracken branch ${branch}"
  else
    local tmp_table tmp_breport tmp_bresource tmp_blog bracken_exit outcome report_sha
    tmp_table="${bracken_table}.tmp.$$"
    tmp_breport="${bracken_report}.tmp.$$"
    tmp_bresource="${work_dir}/resources/${branch}.bracken.txt.tmp.$$"
    tmp_blog="${work_dir}/logs/${branch}.bracken.log.tmp.$$"
    set +e
    /usr/bin/time -v \
      -o "${tmp_bresource}" \
      "${bracken}" \
        -d "${db_dir[${database}]}" \
        -i "${report}" \
        -o "${tmp_table}" \
        -w "${tmp_breport}" \
        -r "${read_length}" \
        -l S \
        -t 10 \
      > "${tmp_blog}" 2>&1
    bracken_exit="$?"
    set -e
    if [[ "${bracken_exit}" -eq 0 ]]; then
      outcome="estimated"
      mv "${tmp_breport}" "${bracken_report}"
      report_sha="$(sha256sum "${bracken_report}" | awk '{print $1}')"
    elif [[ "${bracken_exit}" -eq 1 ]] &&
         grep -Fq "Error: no reads found. Please check your Kraken report" "${tmp_blog}"; then
      outcome="no_eligible_species_above_threshold"
      printf 'name\ttaxonomy_id\ttaxonomy_lvl\tkraken_assigned_reads\tadded_reads\tnew_est_reads\tfraction_total_reads\n' \
        > "${tmp_table}"
      report_sha="not_created"
    else
      echo "Bracken failed for ${branch} with exit status ${bracken_exit}" >&2
      sed -n '1,120p' "${tmp_blog}" >&2
      exit "${bracken_exit}"
    fi
    mv "${tmp_table}" "${bracken_table}"
    mv "${tmp_bresource}" "${work_dir}/resources/${branch}.bracken.txt"
    mv "${tmp_blog}" "${work_dir}/logs/${branch}.bracken.log"
    {
      printf 'branch_id\t%s\n' "${branch}"
      printf 'outcome\t%s\n' "${outcome}"
      printf 'tool_exit_status\t%s\n' "${bracken_exit}"
      printf 'table_sha256\t%s\n' "$(sha256sum "${bracken_table}" | awk '{print $1}')"
      printf 'report_sha256\t%s\n' "${report_sha}"
    } > "${bracken_sentinel}.tmp.$$"
    mv "${bracken_sentinel}.tmp.$$" "${bracken_sentinel}"
  fi
}

for database in standard8 standard16 pluspf8; do
  for confidence in 0.0 0.05 0.10 0.20 0.50; do
    for control in mock human; do
      hit_matrix="No"
      if [[ "${database}" == "standard16" && "${confidence}" == "0.10" ]]; then
        hit_matrix="Yes"
      fi
      run_branch "${database}" "${control}" "${confidence}" 2 Yes "${hit_matrix}"
    done
  done
done

for hit_groups in 1 3 4; do
  for control in mock human; do
    run_branch standard16 "${control}" 0.10 "${hit_groups}" No Yes
  done
done

cp "${branch_design}" "${frozen_dir}/branch-design.tsv"

"${python}" "${project_root}/scripts/validate_article17_kraken2_database_confidence.py" \
  --project-root "${project_root}" \
  --environment-prefix "${environment_prefix}" \
  --frozen-dir "${frozen_dir}" \
  --raw-dir "${raw_dir}" \
  --work-dir "${work_dir}" \
  --database-root "${database_root}" \
  --initialize-frozen

echo "Article 17 one-time database/confidence analysis completed: ${frozen_dir}"
