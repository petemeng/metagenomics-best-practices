#!/usr/bin/env bash
set -euo pipefail

root="${1:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
out="${2:-${root}/data/raw/article55/fastq}"
mkdir -p "${out}"

download_fastq() {
  local run="$1"
  local mate="$2"
  local expected_md5="$3"
  local expected_bytes="$4"
  local rel="vol1/fastq/ERR103/${run: -3}/${run}/${run}_${mate}.fastq.gz"
  local url="https://ftp.sra.ebi.ac.uk/${rel}"
  local target="${out}/${run}_${mate}.fastq.gz"

  if [[ ! -f "${target}" ]]; then
    curl --fail --location \
      --retry 8 --retry-all-errors --retry-delay 5 \
      --connect-timeout 30 --speed-time 90 --speed-limit 1024 \
      --continue-at - --output "${target}.part" "${url}"
    mv "${target}.part" "${target}"
  fi

  local observed_md5 observed_bytes
  observed_md5="$(md5sum "${target}" | awk '{print $1}')"
  observed_bytes="$(stat -c %s "${target}")"
  if [[ "${observed_md5}" != "${expected_md5}" || "${observed_bytes}" != "${expected_bytes}" ]]; then
    echo "identity check failed for ${target}" >&2
    echo "expected md5=${expected_md5} bytes=${expected_bytes}" >&2
    echo "observed md5=${observed_md5} bytes=${observed_bytes}" >&2
    exit 1
  fi
  echo "verified ${run}_${mate}.fastq.gz md5=${observed_md5} bytes=${observed_bytes}"
}

# Cook et al. 2024, PRJEB56639; three unamplified MiSeq paired-end libraries.
# Total compressed size: 2,077,860,761 bytes (1.94 GiB).
download_fastq ERR10359653 1 205beb60c97acc1e44b74e3978b0d70b 316101153
download_fastq ERR10359653 2 fccc18690e9ef5db157df9083f24a3ec 358529570
download_fastq ERR10359656 1 65564ee5c7b50f821440e4d13b94a927 436349283
download_fastq ERR10359656 2 ab5c9db6afefcd4c8bcb0c0435518d08 499721263
download_fastq ERR10359658 1 c9c1d7fe6325f4764434e5178b100da6 221718556
download_fastq ERR10359658 2 74620028ebac16b18a3801c511a415c8 245440936

echo "Article 55 Illumina FASTQ download complete: ${out}"
