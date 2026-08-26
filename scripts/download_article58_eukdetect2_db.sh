#!/usr/bin/env bash
set -euo pipefail

project_root="${1:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
db_root="${2:-${project_root}/db/microbial-eukaryotes-cache/eukdetect2-v2026-03-16}"
mode="${3:-download}"
manifest="${project_root}/data/small/58-eukdetect2-database-manifest.tsv"

[[ "${mode}" == "download" || "${mode}" == "verify" ]] || {
  echo "usage: $0 PROJECT_ROOT DB_ROOT [download|verify]" >&2
  exit 2
}
[[ -r "${manifest}" ]] || {
  echo "database manifest is missing: ${manifest}" >&2
  exit 2
}
mkdir -p "${db_root}"

verified=0
total_bytes=0
while IFS=$'\t' read -r file release doi bytes expected_md5 url; do
  [[ "${file}" != "File" ]] || continue
  target="${db_root}/${file}"
  part="${target}.part"
  if [[ "${mode}" == "download" && ! -s "${target}" ]]; then
    echo "Downloading ${file} (${bytes} bytes)"
    if command -v aria2c >/dev/null 2>&1; then
      aria2c \
        --continue=true \
        --max-connection-per-server=8 \
        --split=8 \
        --min-split-size=16M \
        --piece-length=4M \
        --file-allocation=none \
        --auto-file-renaming=false \
        --allow-overwrite=true \
        --max-tries=0 \
        --retry-wait=5 \
        --connect-timeout=30 \
        --timeout=90 \
        --summary-interval=30 \
        --dir="$(dirname "${part}")" \
        --out="$(basename "${part}")" \
        "${url}"
    else
      curl --fail --location \
        --retry 50 --retry-all-errors --retry-delay 5 \
        --connect-timeout 30 --speed-time 180 --speed-limit 1024 \
        --continue-at - --output "${part}" "${url}"
    fi
    mv "${part}" "${target}"
  fi
  [[ -s "${target}" ]] || {
    echo "missing database file: ${target}" >&2
    exit 2
  }
  observed_bytes="$(stat -c '%s' "${target}")"
  observed_md5="$(md5sum "${target}" | awk '{print $1}')"
  [[ "${observed_bytes}" == "${bytes}" ]] || {
    echo "byte-count mismatch for ${file}: ${observed_bytes} != ${bytes}" >&2
    exit 2
  }
  [[ "${observed_md5}" == "${expected_md5}" ]] || {
    echo "MD5 mismatch for ${file}: ${observed_md5} != ${expected_md5}" >&2
    exit 2
  }
  echo "VERIFIED ${file} md5=${observed_md5}"
  verified=$((verified + 1))
  total_bytes=$((total_bytes + bytes))
done <"${manifest}"

[[ "${verified}" -eq 14 ]] || {
  echo "expected 14 database files, verified ${verified}" >&2
  exit 2
}
{
  printf 'Database\tRelease\tDOI\tFiles\tInstalledBytes\n'
  printf 'EukDetect2 marker database\t2026-03-16\t10.5281/zenodo.19056625\t%s\t%s\n' \
    "${verified}" "${total_bytes}"
} >"${db_root}/database-release.tsv.tmp"
mv "${db_root}/database-release.tsv.tmp" "${db_root}/database-release.tsv"
echo "PASS EukDetect2 database: ${verified} files, ${total_bytes} bytes"
