#!/usr/bin/env bash
set -euo pipefail

usage() {
  printf '%s\n' \
    'Usage:' \
    '  bootstrap_article19_humann_databases.sh --project-root DIR [--cache-root DIR] download' \
    '  bootstrap_article19_humann_databases.sh --project-root DIR [--cache-root DIR] extract-metaphlan' \
    '  bootstrap_article19_humann_databases.sh --project-root DIR [--cache-root DIR] extract-chocophlan' \
    '  bootstrap_article19_humann_databases.sh --project-root DIR [--cache-root DIR] extract-uniref90' \
    '  bootstrap_article19_humann_databases.sh --project-root DIR [--cache-root DIR] extract' \
    '  bootstrap_article19_humann_databases.sh --project-root DIR [--cache-root DIR] inventory' \
    '  bootstrap_article19_humann_databases.sh --project-root DIR [--cache-root DIR] all' >&2
}

project_root=""
cache_root=""
while [[ "$#" -gt 0 ]]; do
  case "$1" in
    --project-root)
      project_root="$2"
      shift 2
      ;;
    --cache-root)
      cache_root="$2"
      shift 2
      ;;
    download|extract-metaphlan|extract-chocophlan|extract-uniref90|extract|inventory|all)
      command_name="$1"
      shift
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      printf 'Unknown argument: %s\n' "$1" >&2
      usage
      exit 2
      ;;
  esac
done

if [[ -z "${project_root}" || -z "${command_name:-}" ]]; then
  usage
  exit 2
fi

project_root="$(cd "${project_root}" && pwd)"
cache_root="${cache_root:-${project_root}/db/humann-cache}"
mkdir -p "${cache_root}"/{archives,headers,installed,logs,manifests}
cache_root="$(cd "${cache_root}" && pwd)"

export LC_ALL=C
export TZ=UTC

result_path="${cache_root}/manifests/bootstrap-results.tsv"

asset_ids=(
  metaphlan-vjun23-metadata
  metaphlan-vjun23-bowtie2
  humann-chocophlan-full
  humann-uniref90-full
)
asset_urls=(
  https://cmprod1.cibio.unitn.it/biobakery4/metaphlan_databases/mpa_vJun23_CHOCOPhlAnSGB_202403.tar
  https://cmprod1.cibio.unitn.it/biobakery4/metaphlan_databases/bowtie2_indexes/mpa_vJun23_CHOCOPhlAnSGB_202403_bt2.tar
  https://huttenhower.sph.harvard.edu/humann_data/chocophlan/full_chocophlan.v201901_v31.tar.gz
  https://huttenhower.sph.harvard.edu/humann_data/uniprot/uniref_annotated/uniref90_annotated_v201901b_full.tar.gz
)
asset_bytes=(3316336640 22989752320 16502062909 20579913329)
asset_algorithms=(md5 md5 none none)
asset_checksums=(
  d985de75a217cd319e721863f68e7d33
  8caae86b4d2931416cbdbb92f5985cef
  unpublished
  unpublished
)

write_header() {
  printf '%s\n' \
    $'database_id\tarchive_path\teffective_url\tremote_content_length\tbytes\tpublisher_checksum_algorithm\tpublisher_checksum\tobserved_publisher_checksum\tobserved_sha256\tretrieved_at_utc\tarchive_integrity' \
    > "${result_path}"
}

fetch_one() {
  local database_id="$1"
  local url="$2"
  local expected_bytes="$3"
  local algorithm="$4"
  local expected_checksum="$5"
  local archive_dir archive_name target part observed_bytes observed_checksum
  local observed_sha256 effective_url remote_content_length retrieved_at
  local retrieved_record integrity

  archive_dir="${cache_root}/archives/${database_id}"
  archive_name="$(basename "${url}")"
  target="${archive_dir}/${archive_name}"
  part="${target}.part"
  mkdir -p "${archive_dir}"

  curl --silent --show-error --location --head --max-time 120 \
    --dump-header "${cache_root}/headers/${database_id}.headers.txt" \
    --output /dev/null "${url}"
  effective_url="$(curl --silent --show-error --location --head --max-time 120 \
    --write-out '%{url_effective}' --output /dev/null "${url}")"
  remote_content_length="$(
    awk 'tolower($1) == "content-length:" {gsub("\\r", "", $2); value=$2} END {print value}' \
      "${cache_root}/headers/${database_id}.headers.txt"
  )"
  if [[ "${remote_content_length}" != "${expected_bytes}" ]]; then
    printf 'Remote Content-Length mismatch for %s: expected %s observed %s\n' \
      "${database_id}" "${expected_bytes}" "${remote_content_length:-missing}" >&2
    exit 1
  fi

  if [[ ! -s "${target}" ]]; then
    if [[ "${effective_url}" == *'.data.globus.org/'* ]]; then
      printf 'Downloading %s with explicit HTTP ranges to %s\n' \
        "${database_id}" "${target}"
      bash "${project_root}/scripts/download_ranged_archive.sh" \
        --url "${effective_url}" \
        --bytes "${expected_bytes}" \
        --output "${target}" \
        --segments 16
    else
      printf 'Downloading %s to %s\n' "${database_id}" "${part}"
      aria2c \
        --continue=true \
        --max-connection-per-server=16 \
        --split=16 \
        --min-split-size=16M \
        --piece-length=4M \
        --file-allocation=none \
        --auto-file-renaming=false \
        --allow-overwrite=true \
        --max-tries=0 \
        --retry-wait=5 \
        --connect-timeout=30 \
        --timeout=60 \
        --summary-interval=30 \
        --dir="${archive_dir}" \
        --out="${archive_name}.part" \
        "${effective_url}"
      mv "${part}" "${target}"
    fi
  fi

  observed_bytes="$(stat -c '%s' "${target}")"
  if [[ "${observed_bytes}" != "${expected_bytes}" ]]; then
    printf 'Byte-count mismatch for %s: expected %s observed %s\n' \
      "${database_id}" "${expected_bytes}" "${observed_bytes}" >&2
    exit 1
  fi

  case "${algorithm}" in
    md5)
      observed_checksum="$(md5sum "${target}" | awk '{print $1}')"
      if [[ "${observed_checksum}" != "${expected_checksum}" ]]; then
        printf 'Publisher MD5 mismatch for %s\n' "${database_id}" >&2
        exit 1
      fi
      ;;
    none)
      observed_checksum="unpublished"
      ;;
    *)
      printf 'Unsupported checksum algorithm: %s\n' "${algorithm}" >&2
      exit 1
      ;;
  esac

  tar -tf "${target}" > "${cache_root}/logs/${database_id}.tar-list.txt"
  integrity="tar-list-pass"
  observed_sha256="$(sha256sum "${target}" | awk '{print $1}')"
  retrieved_record="${target}.retrieved-at-utc"
  if [[ ! -s "${retrieved_record}" ]]; then
    date -u -r "${target}" +'%Y-%m-%dT%H:%M:%SZ' > "${retrieved_record}"
  fi
  retrieved_at="$(tr -d '[:space:]' < "${retrieved_record}")"
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "${database_id}" "${target}" "${effective_url}" "${remote_content_length}" \
    "${observed_bytes}" "${algorithm}" "${expected_checksum}" \
    "${observed_checksum}" "${observed_sha256}" \
    "${retrieved_at}" "${integrity}" >> "${result_path}"
  printf 'LOCKED %s bytes=%s sha256=%s\n' \
    "${database_id}" "${observed_bytes}" "${observed_sha256}"
}

download_all() {
  write_header
  local i
  for i in "${!asset_ids[@]}"; do
    fetch_one \
      "${asset_ids[$i]}" \
      "${asset_urls[$i]}" \
      "${asset_bytes[$i]}" \
      "${asset_algorithms[$i]}" \
      "${asset_checksums[$i]}"
  done
}

extract_one() {
  local database_id="$1"
  local archive_name="$2"
  local destination="$3"
  local allow_nonempty="${4:-no}"
  local target sentinel
  target="${cache_root}/archives/${database_id}/${archive_name}"
  # Keep extraction sentinels outside database directories: HUMAnN validates
  # every directory entry and rejects unrelated marker files.
  sentinel="${cache_root}/manifests/${database_id}.article19-extracted"
  if [[ ! -s "${target}" ]]; then
    printf 'Archive missing: %s\n' "${target}" >&2
    exit 1
  fi
  if [[ -e "${sentinel}" ]]; then
    printf 'Extraction already complete: %s\n' "${destination}"
    return
  fi
  if [[ -d "${destination}" ]] &&
     [[ -n "$(find "${destination}" -mindepth 1 -maxdepth 1 -print -quit)" ]] &&
     [[ "${allow_nonempty}" != "yes" ]]; then
    printf 'Refusing to extract into non-empty unsentinelled directory: %s\n' \
      "${destination}" >&2
    exit 1
  fi
  mkdir -p "${destination}"
  /usr/bin/time -v \
    -o "${cache_root}/logs/${database_id}.extract.resources.txt" \
    tar --no-same-owner --no-same-permissions -xf "${target}" -C "${destination}" \
    > "${cache_root}/logs/${database_id}.extract.log" 2>&1
  printf '%s\n' "${database_id}" > "${sentinel}"
}

extract_all() {
  extract_metaphlan
  extract_chocophlan
  extract_uniref90
}

extract_chocophlan() {
  extract_one \
    humann-chocophlan-full \
    full_chocophlan.v201901_v31.tar.gz \
    "${cache_root}/installed/humann-chocophlan-v201901-v31"
}

extract_uniref90() {
  extract_one \
    humann-uniref90-full \
    uniref90_annotated_v201901b_full.tar.gz \
    "${cache_root}/installed/humann-uniref90-v201901b"
}

extract_metaphlan() {
  extract_one \
    metaphlan-vjun23-metadata \
    mpa_vJun23_CHOCOPhlAnSGB_202403.tar \
    "${cache_root}/installed/metaphlan-vjun23"
  extract_one \
    metaphlan-vjun23-bowtie2 \
    mpa_vJun23_CHOCOPhlAnSGB_202403_bt2.tar \
    "${cache_root}/installed/metaphlan-vjun23" \
    yes
}

inventory_all() {
  local output="${cache_root}/manifests/installed-inventory.tsv"
  local checksum_output="${cache_root}/manifests/installed-files.sha256"
  printf '%s\n' $'database_id\tinstalled_path\tfiles\tinstalled_bytes' > "${output}"
  : > "${checksum_output}"
  local database_id destination files installed_bytes
  while IFS=$'\t' read -r database_id destination; do
    files="$(find "${destination}" -type f ! -name '.*.article19-extracted' | wc -l)"
    installed_bytes="$(find "${destination}" -type f ! -name '.*.article19-extracted' -printf '%s\n' | awk '{s += $1} END {print s + 0}')"
    printf '%s\t%s\t%s\t%s\n' \
      "${database_id}" "${destination}" "${files}" "${installed_bytes}" >> "${output}"
    find "${destination}" -type f ! -name '.*.article19-extracted' -print0 \
      | sort -z \
      | xargs -0 sha256sum >> "${checksum_output}"
  done <<EOF
metaphlan-vjun23	${cache_root}/installed/metaphlan-vjun23
humann-chocophlan-full	${cache_root}/installed/humann-chocophlan-v201901-v31
humann-uniref90-full	${cache_root}/installed/humann-uniref90-v201901b
EOF
  printf 'Inventory written to %s\n' "${output}"
}

case "${command_name}" in
  download)
    download_all
    ;;
  extract-metaphlan)
    extract_metaphlan
    ;;
  extract-chocophlan)
    extract_chocophlan
    ;;
  extract-uniref90)
    extract_uniref90
    ;;
  extract)
    extract_all
    ;;
  inventory)
    inventory_all
    ;;
  all)
    download_all
    extract_all
    inventory_all
    ;;
esac
