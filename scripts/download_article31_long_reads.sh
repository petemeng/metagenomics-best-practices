#!/usr/bin/env bash
set -euo pipefail

usage() {
  echo "Usage: download_article31_long_reads.sh <project-root> [raw-dir]" >&2
}

if [[ "$#" -lt 1 || "$#" -gt 2 ]]; then
  usage
  exit 2
fi

project_root="$(cd "$1" && pwd)"
raw_dir="${2:-${project_root}/data/raw/article31}"
mkdir -p "${raw_dir}/full"
raw_dir="$(cd "${raw_dir}" && pwd)"
manifest="${project_root}/data/small/31-source-manifest.tsv"

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
          --user-agent="metagenomics-best-practices/Article31 ENA downloader" \
          --dir "$(dirname "${target}")" \
          --out "$(basename "${target}").part" \
          "${url}"; then
          downloaded="true"
          break
        fi
        echo "aria2 attempt ${attempt}/5 failed; the .part file will be resumed" >&2
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
    echo "expected bytes=${expected_bytes} md5=${expected_md5}" >&2
    echo "observed bytes=${observed_bytes} md5=${observed_md5}" >&2
    exit 1
  fi
  printf 'verified\t%s\t%s\t%s\n' "${target}" "${observed_bytes}" "${observed_md5}"
}

while IFS=$'\t' read -r platform label run sample layout url expected_md5 expected_bytes read_count base_count mean_length max_length identity boundary checked; do
  if [[ "${platform}" == "PlatformKey" ]]; then
    continue
  fi
  target="${raw_dir}/full/${run}.fastq.gz"
  download_one "${url}" "${target}" "${expected_md5}" "${expected_bytes}"
done < "${manifest}"

echo "Article 31 checksum-identified long-read archives are ready in ${raw_dir}/full."
