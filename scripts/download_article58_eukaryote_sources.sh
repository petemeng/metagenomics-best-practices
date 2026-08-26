#!/usr/bin/env bash
set -euo pipefail

root="${1:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
raw_rel="${2:-data/raw/article58}"
raw_dir="${root}/${raw_rel}"
mkdir -p "${raw_dir}"

download_verified() {
  local url="$1"
  local target="$2"
  local expected_bytes="$3"
  local expected_md5="$4"
  local partial="${target}.part"

  if [[ ! -f "${target}" ]]; then
    if command -v aria2c >/dev/null 2>&1; then
      aria2c \
        --allow-overwrite=true \
        --auto-file-renaming=false \
        --continue=true \
        --file-allocation=none \
        --max-connection-per-server=8 \
        --max-tries=0 \
        --min-split-size=10M \
        --retry-wait=5 \
        --split=8 \
        --connect-timeout=30 \
        --timeout=90 \
        --summary-interval=30 \
        --dir "$(dirname "${partial}")" \
        --out "$(basename "${partial}")" \
        "${url}"
    else
      curl --fail --location \
        --retry 50 --retry-all-errors --retry-delay 5 \
        --connect-timeout 30 --speed-time 180 --speed-limit 1024 \
        --continue-at - --output "${partial}" "${url}"
    fi

    local partial_bytes
    partial_bytes="$(stat -c '%s' "${partial}")"
    if [[ "${partial_bytes}" != "${expected_bytes}" ]]; then
      echo "Byte-count mismatch for ${partial}: ${partial_bytes} != ${expected_bytes}" >&2
      exit 1
    fi
    local partial_md5
    partial_md5="$(md5sum "${partial}" | awk '{print $1}')"
    if [[ "${partial_md5}" != "${expected_md5}" ]]; then
      echo "MD5 mismatch for ${partial}: ${partial_md5} != ${expected_md5}" >&2
      exit 1
    fi
    mv "${partial}" "${target}"
  fi

  local observed_bytes observed_md5
  observed_bytes="$(stat -c '%s' "${target}")"
  observed_md5="$(md5sum "${target}" | awk '{print $1}')"
  if [[ "${observed_bytes}" != "${expected_bytes}" ]]; then
    echo "Byte-count mismatch for ${target}: ${observed_bytes} != ${expected_bytes}" >&2
    exit 1
  fi
  if [[ "${observed_md5}" != "${expected_md5}" ]]; then
    echo "MD5 mismatch for ${target}: ${observed_md5} != ${expected_md5}" >&2
    exit 1
  fi
  echo "verified ${target} bytes=${observed_bytes} md5=${observed_md5}"
}

download_verified \
  "https://zenodo.org/records/3935737/files/ZymoBIOMICS.STD.refseq.v2.zip?download=1" \
  "${raw_dir}/ZymoBIOMICS.STD.refseq.v2.zip" \
  "21453226" \
  "be18fd9195379082096061b8249489b3"

download_verified \
  "https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR123/053/SRR12324253/SRR12324253_1.fastq.gz" \
  "${raw_dir}/SRR12324253_1.fastq.gz" \
  "3013365268" \
  "1485049c9792c7e43d267fe7bb84a2bd"

download_verified \
  "https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR123/053/SRR12324253/SRR12324253_2.fastq.gz" \
  "${raw_dir}/SRR12324253_2.fastq.gz" \
  "3010771717" \
  "fbd8494da79fc796d6725a4e242a9b9c"

echo "Article 58 source assets passed official byte-count and MD5 gates."
