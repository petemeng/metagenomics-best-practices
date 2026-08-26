#!/usr/bin/env bash
set -euo pipefail

project_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
manifest="${project_root}/data/small/36-eggnog-database-manifest.tsv"

usage() {
  cat <<'EOF'
Usage:
  scripts/bootstrap_article36_eggnog_database.sh list
  scripts/bootstrap_article36_eggnog_database.sh audit --database-dir /absolute/path
  scripts/bootstrap_article36_eggnog_database.sh download --database-dir /absolute/path

The download is about 12.1 GB compressed and installs about 51.9 GB. Archives
are retained under DATABASE_DIR/archives; every archive and installed file is
checked against the Article 36 SHA-256 manifest. The official eggNOG 5.0.2
URLs use HTTP, so checksum verification is mandatory.
EOF
}

die() {
  echo "ERROR: $*" >&2
  exit 2
}

command="${1:-}"
shift || true
database_dir=""
while (($#)); do
  case "$1" in
    --database-dir)
      (($# >= 2)) || die "--database-dir requires a path"
      database_dir="$2"
      shift 2
      ;;
    *) die "unknown argument: $1" ;;
  esac
done

[[ -r "${manifest}" ]] || die "manifest is missing: ${manifest}"

case "${command}" in
  list)
    awk -F '\t' 'BEGIN{OFS="\t"} NR==1{print "Asset","Release","ArchiveFile","ArchiveBytes","InstalledFile","InstalledBytes"; next} {print $1,$2,$4,$5,$7,$8}' "${manifest}"
    exit 0
    ;;
  audit|download) ;;
  *) usage; exit 2 ;;
esac

[[ -n "${database_dir}" ]] || die "set --database-dir"
[[ "${database_dir}" = /* ]] || die "--database-dir must be absolute"
database_dir="$(realpath -m "${database_dir}")"
case "${database_dir}" in
  "${project_root}"|"${project_root}"/*) die "database directory must be outside the Git repository" ;;
esac
mkdir -p "${database_dir}" "${database_dir}/archives"

verify_file() {
  local path="$1" expected_bytes="$2" expected_sha="$3" label="$4"
  [[ -s "${path}" ]] || die "missing ${label}: ${path}"
  local observed_bytes observed_sha
  observed_bytes="$(stat -c '%s' "${path}")"
  [[ "${observed_bytes}" == "${expected_bytes}" ]] ||
    die "${label} byte mismatch: expected ${expected_bytes}, observed ${observed_bytes}"
  observed_sha="$(sha256sum "${path}" | awk '{print $1}')"
  [[ "${observed_sha}" == "${expected_sha}" ]] ||
    die "${label} SHA-256 mismatch: expected ${expected_sha}, observed ${observed_sha}"
  printf 'PASS\t%s\t%s\t%s\n' "${label}" "${observed_bytes}" "${observed_sha}"
}

download_archives() {
  while IFS=$'\t' read -r archive_file source_url archive_bytes archive_sha; do
    local_path="${database_dir}/archives/${archive_file}"
    if [[ -s "${local_path}" ]]; then
      verify_file "${local_path}" "${archive_bytes}" "${archive_sha}" "archive:${archive_file}"
      continue
    fi
    part="${local_path}.part"
    curl --fail --location --retry 8 --retry-delay 5 --continue-at - \
      --output "${part}" "${source_url}"
    verify_file "${part}" "${archive_bytes}" "${archive_sha}" "archive:${archive_file}"
    mv "${part}" "${local_path}"
  done < <(
    awk -F '\t' 'BEGIN{OFS="\t"} NR>1 && !seen[$4]++ {print $4,$3,$5,$6}' "${manifest}"
  )
}

install_archives() {
  gzip -cd "${database_dir}/archives/eggnog.db.gz" > "${database_dir}/eggnog.db.tmp"
  mv "${database_dir}/eggnog.db.tmp" "${database_dir}/eggnog.db"
  gzip -cd "${database_dir}/archives/eggnog_proteins.dmnd.gz" > "${database_dir}/eggnog_proteins.dmnd.tmp"
  mv "${database_dir}/eggnog_proteins.dmnd.tmp" "${database_dir}/eggnog_proteins.dmnd"
  tar -tzf "${database_dir}/archives/eggnog.taxa.tar.gz" | \
    awk '$0 != "eggnog.taxa.db" && $0 != "eggnog.taxa.db.traverse.pkl" {bad=1} END{exit bad}' ||
    die "taxonomy archive contains unexpected paths"
  tar -xzf "${database_dir}/archives/eggnog.taxa.tar.gz" \
    -C "${database_dir}" eggnog.taxa.db eggnog.taxa.db.traverse.pkl
}

audit_installed() {
  tail -n +2 "${manifest}" | while IFS=$'\t' read -r \
    asset release source_url archive_file archive_bytes archive_sha \
    installed_file installed_bytes installed_sha required notes; do
    [[ "${required}" == "yes" ]] || continue
    verify_file "${database_dir}/${installed_file}" "${installed_bytes}" "${installed_sha}" "installed:${asset}"
  done
  version="$(sqlite3 "${database_dir}/eggnog.db" 'SELECT version FROM version;')"
  [[ "${version}" == "5.0.2" ]] || die "SQLite release mismatch: ${version}"
  printf 'PASS\teggNOG-release\t%s\n' "${version}"
}

if [[ "${command}" == "download" ]]; then
  download_archives
  install_archives
fi
audit_installed
