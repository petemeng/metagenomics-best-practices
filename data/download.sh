#!/usr/bin/env bash
set -euo pipefail

root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
raw_dir="${root}/data/raw"
metadata="${raw_dir}/meta_all.tsv"
url="https://zenodo.org/records/3517209/files/meta_all.tsv?download=1"
expected_md5="da18b10fdabae6308329e80b73991f84"

mkdir -p "${raw_dir}"
if [[ ! -f "${metadata}" ]]; then
  partial="${metadata}.part"
  curl --fail --location \
    --retry 5 --retry-all-errors --retry-delay 5 \
    --connect-timeout 30 --speed-time 60 --speed-limit 1024 \
    --continue-at - --output "${partial}" "${url}"
  partial_md5="$(md5sum "${partial}" | awk '{print $1}')"
  if [[ "${partial_md5}" != "${expected_md5}" ]]; then
    echo "MD5 mismatch for incomplete download ${partial}" >&2
    rm -f "${partial}"
    exit 1
  fi
  mv "${partial}" "${metadata}"
fi

observed_md5="$(md5sum "${metadata}" | awk '{print $1}')"
if [[ "${observed_md5}" != "${expected_md5}" ]]; then
  echo "MD5 mismatch for ${metadata}" >&2
  exit 1
fi

echo "verified ${metadata} md5=${observed_md5}"

if [[ "${DOWNLOAD_ARTICLE20_CMD:-0}" == "1" ]]; then
  article20_dir="${raw_dir}/article20-cmd"
  mkdir -p "${article20_dir}"

  download_sha256() {
    local url="$1"
    local target="$2"
    local expected_sha256="$3"
    if [[ ! -f "${target}" ]]; then
      curl --fail --location \
        --retry 5 --retry-all-errors --retry-delay 5 \
        --connect-timeout 30 --speed-time 60 --speed-limit 1024 \
        --output "${target}.part" "${url}"
      mv "${target}.part" "${target}"
    fi
    local observed_sha256
    observed_sha256="$(sha256sum "${target}" | awk '{print $1}')"
    if [[ "${observed_sha256}" != "${expected_sha256}" ]]; then
      echo "SHA-256 mismatch for ${target}" >&2
      exit 1
    fi
    echo "verified ${target} sha256=${observed_sha256}"
  }

  cmd_base="https://mghp.osn.xsede.org/bir190004-bucket01/ExperimentHub/curatedMetagenomicData/2021-10-14/AsnicarF_2017"
  download_sha256 \
    "${cmd_base}/2021-10-14.AsnicarF_2017.pathway_abundance.rda" \
    "${article20_dir}/AsnicarF_2017.pathway_abundance.rda" \
    "ead7c78c075fec92a7d641b731594e068b2ba2a47479151d081c338f615af121"
  download_sha256 \
    "${cmd_base}/2021-10-14.AsnicarF_2017.pathway_coverage.rda" \
    "${article20_dir}/AsnicarF_2017.pathway_coverage.rda" \
    "73a1b77b70f88e9028e8707ba3e99b93f0ff99cd91401a3966c4f7a31dbfc3a1"

  echo "Run scripts/prepare_article20_cmd_pathways.R to create the checksum-locked small exports."
fi

if [[ "${DOWNLOAD_ARTICLE30_ASSEMBLY:-0}" == "1" ]]; then
  bash "${root}/scripts/download_article30_assembly_reads.sh" \
    "${root}" \
    "${raw_dir}/article30"
fi

if [[ "${DOWNLOAD_ARTICLE31_LONG_READ:-0}" == "1" ]]; then
  bash "${root}/scripts/download_article31_long_reads.sh" \
    "${root}" \
    "${raw_dir}/article31"
fi

if [[ "${DOWNLOAD_ARTICLE32_HYBRID:-0}" == "1" ]]; then
  bash "${root}/scripts/download_article32_hybrid_sources.sh" "${root}"
fi

if [[ "${DOWNLOAD_ARTICLE55_ILLUMINA:-0}" == "1" ]]; then
  bash "${root}/scripts/download_article55_illumina_fastq.sh" \
    "${root}" \
    "${raw_dir}/article55/fastq"
fi

if [[ "${DOWNLOAD_ARTICLE57_PLASMIDS:-0}" == "1" ]]; then
  bash "${root}/scripts/download_article57_plasmid_sources.sh" \
    "${root}" \
    "data/raw/article57"
fi

if [[ "${DOWNLOAD_ARTICLE58_EUKARYOTES:-0}" == "1" ]]; then
  bash "${root}/scripts/download_article58_eukaryote_sources.sh" \
    "${root}" \
    "data/raw/article58"
fi
