#!/usr/bin/env bash
set -euo pipefail

usage() {
  echo "Usage: download_article30_assembly_reads.sh <project-root> [raw-dir]" >&2
}

if [[ "$#" -lt 1 || "$#" -gt 2 ]]; then
  usage
  exit 2
fi

project_root="$(cd "$1" && pwd)"
raw_dir="${2:-${project_root}/data/raw/article30}"
mkdir -p "${raw_dir}/full" "${raw_dir}/selected"
raw_dir="$(cd "${raw_dir}" && pwd)"
manifest="${project_root}/data/small/30-source-manifest.tsv"
selector="${project_root}/scripts/select_article30_read_pairs.py"

download_one() {
  local url="$1"
  local target="$2"
  local expected_md5="$3"
  local expected_bytes="$4"
  local observed_md5 observed_bytes

  if [[ ! -f "${target}" ]]; then
    if command -v aria2c >/dev/null 2>&1; then
      local attempt downloaded="false"
      for attempt in 1 2 3 4 5; do
        if aria2c \
          --continue=true \
          --allow-overwrite=true \
          --auto-file-renaming=false \
          --file-allocation=none \
          --max-connection-per-server=4 \
          --split=4 \
          --min-split-size=20M \
          --max-tries=8 \
          --retry-wait=10 \
          --timeout=60 \
          --user-agent="metagenomics-best-practices/Article30 ENA downloader" \
          --dir "$(dirname "${target}")" \
          --out "$(basename "${target}").part" \
          "${url}"; then
          downloaded="true"
          break
        fi
        echo "aria2 attempt ${attempt}/5 failed; backing off before resume" >&2
        sleep "$((attempt * 10))"
      done
      if [[ "${downloaded}" != "true" ]]; then
        echo "aria2 could not complete after five attempts: ${url}" >&2
        exit 1
      fi
    else
      curl --fail --location --continue-at - \
        --retry 8 --retry-all-errors --retry-delay 5 \
        --connect-timeout 30 --speed-time 120 --speed-limit 1024 \
        --output "${target}.part" "${url}"
    fi
    mv "${target}.part" "${target}"
  fi
  observed_bytes="$(stat -c '%s' "${target}")"
  observed_md5="$(md5sum "${target}" | cut -d ' ' -f 1)"
  if [[ "${observed_bytes}" != "${expected_bytes}" || "${observed_md5}" != "${expected_md5}" ]]; then
    echo "ENA archive identity mismatch: ${target}" >&2
    exit 1
  fi
  printf 'verified\t%s\t%s\t%s\n' "${target}" "${observed_bytes}" "${observed_md5}"
}

while IFS=$'\t' read -r mock run sample mate url expected_md5 expected_bytes read_count base_count instrument layout checked; do
  if [[ "${mock}" == "Mock" ]]; then
    continue
  fi
  target="${raw_dir}/full/${run}_${mate}.fastq.gz"
  download_one "${url}" "${target}" "${expected_md5}" "${expected_bytes}"
done < "${manifest}"

select_mock() {
  local mock="$1"
  local run="$2"
  local total_pairs="$3"
  local out1="${raw_dir}/selected/${run}_selected2m_R1.fastq.gz"
  local out2="${raw_dir}/selected/${run}_selected2m_R2.fastq.gz"
  local summary="${raw_dir}/selected/${run}_selection-summary.json"

  if [[ -f "${summary}" && -s "${out1}" && -s "${out2}" ]]; then
    echo "selection already exists: ${run}"
    return
  fi
  if [[ -e "${summary}" || -e "${out1}" || -e "${out2}" ]]; then
    echo "Refusing incomplete Article 30 selection outputs for ${run}" >&2
    exit 1
  fi
  python3 "${selector}" \
    --r1 "${raw_dir}/full/${run}_R1.fastq.gz" \
    --r2 "${raw_dir}/full/${run}_R2.fastq.gz" \
    --output-r1 "${out1}" \
    --output-r2 "${out2}" \
    --summary "${summary}" \
    --sample "${mock}" \
    --run-accession "${run}" \
    --total-pairs "${total_pairs}" \
    --target-pairs 2000000 \
    --seed 20260730
}

select_mock MOCK1 ERR9765746 20597525
select_mock MOCK2 ERR9765747 23173964

echo "Article 30 ENA archives and deterministic paired subsets are ready in ${raw_dir}."
