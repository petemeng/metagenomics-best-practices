#!/usr/bin/env bash
set -euo pipefail

usage() {
  echo "Usage: download_article32_hybrid_sources.sh <project-root>" >&2
}
if [[ "$#" -ne 1 ]]; then usage; exit 2; fi

project_root="$(cd "$1" && pwd)"
cache_article30="${project_root}/data/raw/article30/full"
cache_article31="${project_root}/data/raw/article31/full"
article32="${project_root}/data/raw/article32"
sources="${article32}/sources"
mkdir -p "${sources}" "${article32}"

reuse_verified_cache() {
  local cache="$1" target="$2" expected_md5="$3" expected_bytes="$4"
  if [[ -e "${target}" || ! -f "${cache}" ]]; then
    return
  fi
  local bytes md5
  bytes="$(stat -c '%s' "${cache}")"
  md5="$(md5sum "${cache}" | cut -d ' ' -f 1)"
  if [[ "${bytes}" == "${expected_bytes}" && "${md5}" == "${expected_md5}" ]]; then
    ln "${cache}" "${target}"
    printf 'reused-verified-cache\t%s\t%s\n' "${cache}" "${target}"
  fi
}

download_one() {
  local url="$1" target="$2" expected_md5="$3" expected_bytes="$4"
  if [[ ! -f "${target}" ]]; then
    if command -v aria2c >/dev/null 2>&1; then
      aria2c --continue=true --allow-overwrite=true --auto-file-renaming=false \
        --file-allocation=none --max-connection-per-server=4 --split=4 \
        --min-split-size=20M --max-tries=8 --retry-wait=10 --timeout=60 \
        --dir "$(dirname "${target}")" --out "$(basename "${target}").part" \
        "${url}"
    else
      curl --fail --location --continue-at - --retry 8 --retry-all-errors \
        --retry-delay 5 --connect-timeout 30 --speed-time 120 --speed-limit 1024 \
        --output "${target}.part" "${url}"
    fi
    mv "${target}.part" "${target}"
  fi
  local bytes md5
  bytes="$(stat -c '%s' "${target}")"
  md5="$(md5sum "${target}" | cut -d ' ' -f 1)"
  if [[ "${bytes}" != "${expected_bytes}" || "${md5}" != "${expected_md5}" ]]; then
    echo "Source identity mismatch: ${target}" >&2
    echo "expected bytes=${expected_bytes} md5=${expected_md5}" >&2
    echo "observed bytes=${bytes} md5=${md5}" >&2
    exit 1
  fi
  printf 'verified\t%s\t%s\t%s\n' "${target}" "${bytes}" "${md5}"
}

reuse_verified_cache "${cache_article30}/ERR9765746_R1.fastq.gz" \
  "${sources}/ERR9765746_R1.fastq.gz" ed0e6e0ee846542531c742a45181cd6f 1740647656
reuse_verified_cache "${cache_article30}/ERR9765746_R2.fastq.gz" \
  "${sources}/ERR9765746_R2.fastq.gz" 5b60ac93cb69dff77ae38cfa501afd06 2104551765
reuse_verified_cache "${cache_article31}/ERR9765780.fastq.gz" \
  "${sources}/ERR9765780.fastq.gz" 33eb90ac7437b0039180f03e7a697269 3117261341
reuse_verified_cache "${cache_article31}/ERR9765783.fastq.gz" \
  "${sources}/ERR9765783.fastq.gz" 02ec4bc541b4e1ec5d0f58e4a519f2cb 3982506052

download_one \
  'https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR976/006/ERR9765746/ERR9765746_1.fastq.gz' \
  "${sources}/ERR9765746_R1.fastq.gz" ed0e6e0ee846542531c742a45181cd6f 1740647656
download_one \
  'https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR976/006/ERR9765746/ERR9765746_2.fastq.gz' \
  "${sources}/ERR9765746_R2.fastq.gz" 5b60ac93cb69dff77ae38cfa501afd06 2104551765
download_one \
  'https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR976/000/ERR9765780/ERR9765780.fastq.gz' \
  "${sources}/ERR9765780.fastq.gz" 33eb90ac7437b0039180f03e7a697269 3117261341
download_one \
  'https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR976/003/ERR9765783/ERR9765783.fastq.gz' \
  "${sources}/ERR9765783.fastq.gz" 02ec4bc541b4e1ec5d0f58e4a519f2cb 3982506052
download_one \
  'https://static-content.springer.com/esm/art%3A10.1038%2Fs41597-022-01762-z/MediaObjects/41597_2022_1762_MOESM2_ESM.xlsx' \
  "${article32}/Supplementary_Table_S2.xlsx" b279595390293d2f6b3e61f975576c14 14775

benchmark_repo="${article32}/benchmark_mock"
expected_commit="a429a3724d4593f35b8d7323b20252a6be90e1cd"
if [[ ! -d "${benchmark_repo}/.git" ]]; then
  git clone --filter=blob:none \
    https://forgemia.inra.fr/metagenopolis/benchmark_mock.git \
    "${benchmark_repo}"
fi
observed_commit="$(git -C "${benchmark_repo}" rev-parse HEAD)"
if [[ "${observed_commit}" != "${expected_commit}" ]]; then
  echo "Benchmark repository is not at the locked commit: ${observed_commit}" >&2
  echo "Expected ${expected_commit}; use a clean Article 32 raw directory." >&2
  exit 1
fi
for identity in \
  "script_r/Supplementary_Table_S1.xlsx:937653a56fea7fbfcbe35b3f35c721b4125072ba4ab04c44c9d454697240c6df" \
  "reference/MOCK_001.fasta.gz:ad2641ac155006387e722ae5ec8592fa077f9ab3cf411cc11f5757430b8e752f"; do
  relative="${identity%%:*}"
  expected_sha256="${identity##*:}"
  observed_sha256="$(sha256sum "${benchmark_repo}/${relative}" | cut -d ' ' -f 1)"
  if [[ "${observed_sha256}" != "${expected_sha256}" ]]; then
    echo "Benchmark payload SHA-256 mismatch: ${relative}" >&2
    exit 1
  fi
done
printf 'verified\tbenchmark_mock\t%s\n' "${observed_commit}"

echo "Article 32 official sources are checksum-identified and ready."
