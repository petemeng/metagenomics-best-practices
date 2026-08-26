#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
project_root="$(cd "${script_dir}/.." && pwd)"
manifest="${DB_MANIFEST:-${project_root}/data/small/11-database-manifest.tsv}"

usage() {
  cat <<'EOF'
Usage:
  db/download_db.sh list
  db/download_db.sh audit
  DB_ROOT=/absolute/path db/download_db.sh download DATABASE_ID
  DB_ROOT=/absolute/path db/download_db.sh verify DATABASE_ID
  DB_ROOT=/absolute/path db/download_db.sh antismash8-install
  DB_ROOT=/absolute/path db/download_db.sh antismash8-verify

Set DB_MANIFEST=data/small/15-database-manifest.tsv to use the two-archive
MetaPhlAn Article 15 registry (`marker_metadata` and `bowtie2_index`).
Set DB_MANIFEST=data/small/16-database-manifest.tsv to use the checksum-locked
Kraken2 Standard-8 Article 16 registry.
Set DB_MANIFEST=data/small/17-database-manifest.tsv to compare the same-date
Kraken2 Standard-8, Standard-16 and PlusPF-8 Article 17 registries.
Set DB_MANIFEST=data/small/18-database-manifest.tsv to retrieve the immutable
mOTUs marker-gene database 4.1 archive used by Article 18.
Set DB_MANIFEST=data/small/19-database-manifest.tsv to verify the Article 19
MetaPhlAn publisher checksums and the explicitly labelled HUMAnN retrieval
locks after their one-time bootstrap acquisition.
Set DB_MANIFEST=data/small/44-mag-qc-database-manifest.tsv to retrieve the
immutable CheckM2 v3, GUNC ProGenomes 2.1, and CheckM1 2015-01-16 archives.
Set DB_MANIFEST=data/small/46-gtdbtk-database-manifest.tsv to retrieve the
immutable 60,806,405,195-byte GTDB-Tk R232 full reference package for
GTDB-Tk 2.7.2. Extract it outside Git and set GTDBTK_DATA_PATH to that path.
The registry locks compressed byte counts and publisher checksums; extraction
is deliberately separate so database files remain outside Git.

For Article 40, create the antiSMASH 8.0.4 environment first, activate it,
and set DB_ROOT to the dedicated antiSMASH database directory. The install
command uses the official downloader; install and verify both fail unless the
PFAM 35.0, MIBiG 4.0, and MITE 1.3 files match the locked local manifest.

Routine QA runs only `list` and `audit`. `download` is an explicit,
multi-gigabyte operator action. A row is downloadable only when its gate is
`enabled` and it contains a release-specific URL plus an approved locked checksum.
For retrieval-locked assets, `enabled` is set only after the initial HTTPS,
byte-count, archive-integrity, and local SHA-256 evidence has been reviewed.
When aria2c is available, download uses 16 resumable ranged connections;
set DB_DOWNLOAD_ENGINE=curl to force the single-connection fallback.
EOF
}

die() {
  echo "ERROR: $*" >&2
  exit 2
}

[[ -r "${manifest}" ]] || die "database manifest is not readable: ${manifest}"

IFS= read -r manifest_header <"${manifest}"
case "${manifest_header%%$'\t'*}" in
  database_id)
    manifest_schema="database_registry"
    ;;
  AssetRole)
    manifest_schema="metaphlan_asset_registry"
    ;;
  *)
    die "unsupported database manifest schema: ${manifest}"
    ;;
esac

lookup_row() {
  local database_id="$1"
  awk -F '\t' -v id="${database_id}" '
    NR > 1 && $1 == id { print; found = 1; exit }
    END { if (!found) exit 3 }
  ' "${manifest}"
}

read_row() {
  local database_id="$1"
  local row
  row="$(lookup_row "${database_id}")" ||
    die "unknown database_id: ${database_id}"
  if [[ "${manifest_schema}" == "database_registry" ]]; then
    IFS=$'\t' read -r \
      DB_ID TOOL TOOL_VERSION RELEASE_ID RELEASE_DATE ARCHIVE_URL \
      CHECKSUM_ALGORITHM EXPECTED_CHECKSUM EXPECTED_COMPRESSED_BYTES \
      EXPECTED_INSTALLED_BYTES DOWNLOAD_GATE VALIDATION_STATUS NOTES <<<"${row}"
  else
    local asset_role index_name release asset_name official_last_modified
    local archive_sha256 extracted_location source_checked
    IFS=$'\t' read -r \
      asset_role index_name release TOOL_VERSION asset_name ARCHIVE_URL \
      EXPECTED_CHECKSUM EXPECTED_COMPRESSED_BYTES archive_sha256 \
      official_last_modified extracted_location VALIDATION_STATUS \
      source_checked <<<"${row}"
    DB_ID="${asset_role}"
    TOOL="MetaPhlAn"
    RELEASE_ID="${index_name}"
    RELEASE_DATE="${official_last_modified%%T*}"
    CHECKSUM_ALGORITHM="md5"
    EXPECTED_INSTALLED_BYTES="0"
    DOWNLOAD_GATE="enabled"
    NOTES="${release};${asset_name};sha256=${archive_sha256};${extracted_location};source_checked=${source_checked}"
  fi
}

require_db_root() {
  [[ -n "${DB_ROOT:-}" ]] || die "set DB_ROOT to an absolute database path"
  [[ "${DB_ROOT}" = /* ]] || die "DB_ROOT must be absolute"
  local resolved_root
  resolved_root="$(realpath -m "${DB_ROOT}")"
  case "${resolved_root}" in
    "${project_root}"|"${project_root}"/*)
      die "DB_ROOT must be outside the Git repository"
      ;;
  esac
  DB_ROOT="${resolved_root}"
}

checksum_command() {
  case "${CHECKSUM_ALGORITHM}" in
    md5) printf '%s\n' md5sum ;;
    sha256) printf '%s\n' sha256sum ;;
    *) die "${DB_ID} has no supported locked checksum" ;;
  esac
}

require_enabled_row() {
  [[ "${DOWNLOAD_GATE}" == "enabled" ]] ||
    die "${DB_ID} is fail-closed: ${DOWNLOAD_GATE}"
  [[ "${RELEASE_ID}" != "unknown" ]] ||
    die "${DB_ID} has no immutable release identifier"
  [[ "${ARCHIVE_URL}" != *"/latest/"* ]] ||
    die "${DB_ID} uses a mutable latest URL"
  [[ "${EXPECTED_CHECKSUM}" =~ ^[0-9a-fA-F]{32}$|^[0-9a-fA-F]{64}$ ]] ||
    die "${DB_ID} has no valid locked checksum"
}

archive_path() {
  local archive_name
  archive_name="$(basename "${ARCHIVE_URL%%\?*}")"
  printf '%s\n' \
    "${DB_ROOT}/archives/${DB_ID}/${archive_name}"
}

verify_archive() {
  local target="$1"
  [[ -s "${target}" ]] || die "archive is missing or empty: ${target}"
  local checker observed observed_bytes
  observed_bytes="$(stat -c '%s' "${target}")"
  if [[ "${EXPECTED_COMPRESSED_BYTES}" =~ ^[0-9]+$ ]] &&
     (( EXPECTED_COMPRESSED_BYTES > 0 )); then
    (( observed_bytes == EXPECTED_COMPRESSED_BYTES )) ||
      die "byte-count mismatch for ${DB_ID}: expected ${EXPECTED_COMPRESSED_BYTES}, observed ${observed_bytes}"
  fi
  checker="$(checksum_command)"
  observed="$("${checker}" "${target}" | awk '{print $1}')"
  [[ "${observed,,}" == "${EXPECTED_CHECKSUM,,}" ]] ||
    die "checksum mismatch for ${DB_ID}: expected ${EXPECTED_CHECKSUM}, observed ${observed}"

  local sentinel tmp
  sentinel="${target}.${CHECKSUM_ALGORITHM}.verified"
  tmp="${sentinel}.tmp.$$"
  {
    printf 'database_id\t%s\n' "${DB_ID}"
    printf 'release_id\t%s\n' "${RELEASE_ID}"
    printf 'algorithm\t%s\n' "${CHECKSUM_ALGORITHM}"
    printf 'checksum\t%s\n' "${observed,,}"
    printf 'bytes\t%s\n' "${observed_bytes}"
    printf 'archive\t%s\n' "$(basename "${target}")"
  } >"${tmp}"
  mv "${tmp}" "${sentinel}"
  echo "VERIFIED ${DB_ID} ${observed,,}"
}

verify_antismash8_database() {
  local locked_manifest="${project_root}/data/small/40-bgc-database-manifest.tsv"
  [[ -r "${locked_manifest}" ]] || die "Article 40 database manifest is missing"
  python3 - "${locked_manifest}" "${DB_ROOT}" <<'PY'
import csv
import hashlib
import sys
from pathlib import Path

manifest, database = Path(sys.argv[1]), Path(sys.argv[2])
rows = list(csv.DictReader(manifest.open(encoding="utf-8"), delimiter="\t"))
for row in rows:
    path = database / row["File"]
    if not path.is_file():
        raise SystemExit(f"MISSING {path}")
    observed_bytes = path.stat().st_size
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    observed_hash = digest.hexdigest()
    if observed_bytes != int(row["Bytes"]) or observed_hash != row["SHA256"]:
        raise SystemExit(f"MISMATCH {path}")
    print(f"PASS\t{row['Database']}\t{row['Release']}\t{row['File']}")
print(f"PASS antiSMASH 8 database audit: {len(rows)} locked files")
PY
}

command="${1:-}"
case "${command}" in
  list)
    if [[ "${manifest_schema}" == "database_registry" ]]; then
      awk -F '\t' '
        BEGIN {
          OFS = "\t"
          print "database_id", "tool", "tool_version", "release_id",
                "download_gate", "validation_status"
        }
        NR > 1 { print $1, $2, $3, $4, $11, $12 }
      ' "${manifest}"
    else
      awk -F '\t' '
        BEGIN {
          OFS = "\t"
          print "database_id", "tool", "tool_version", "release_id",
                "download_gate", "validation_status"
        }
        NR > 1 { print $1, "MetaPhlAn", $4, $2, "enabled", $12 }
      ' "${manifest}"
    fi
    ;;

  audit)
    python3 - "${manifest}" <<'PY'
import csv
import re
import sys
from pathlib import Path

path = Path(sys.argv[1])
rows = list(csv.DictReader(path.open(encoding="utf-8"), delimiter="\t"))
if not rows:
    raise SystemExit("manifest has no database rows")
schema = "database_registry" if "database_id" in rows[0] else "metaphlan_asset_registry"
id_field = "database_id" if schema == "database_registry" else "AssetRole"
ids = [row[id_field] for row in rows]
if len(ids) != len(set(ids)):
    raise SystemExit("database_id values are not unique")
for row in rows:
    if schema == "database_registry":
        database_id = row["database_id"]
        enabled = row["download_gate"] == "enabled"
        release = row["release_id"]
        url = row["archive_url"]
        algorithm = row["checksum_algorithm"]
        checksum = row["expected_checksum"]
        byte_count = row["expected_compressed_bytes"]
    else:
        database_id = row["AssetRole"]
        enabled = True
        release = row["IndexName"]
        url = row["AssetURL"]
        algorithm = "md5"
        checksum = row["OfficialMD5"]
        byte_count = row["ArchiveBytes"]
        if not re.fullmatch(r"[0-9a-fA-F]{64}", row["ArchiveSHA256"]):
            raise SystemExit(f"{database_id}: invalid local SHA-256")
    if enabled:
        if release == "unknown":
            raise SystemExit(f"{database_id}: release is unknown")
        if "/latest/" in url:
            raise SystemExit(f"{database_id}: mutable latest URL")
        if algorithm not in {"md5", "sha256"}:
            raise SystemExit(f"{database_id}: unsupported checksum")
        if not re.fullmatch(r"[0-9a-fA-F]{32}|[0-9a-fA-F]{64}", checksum):
            raise SystemExit(f"{database_id}: invalid checksum")
        if not byte_count.isdigit() or int(byte_count) <= 0:
            raise SystemExit(f"{database_id}: invalid archive byte count")
print(f"PASS database manifest audit: {len(rows)} rows")
PY
    ;;

  download)
    [[ "$#" -eq 2 ]] || die "download requires DATABASE_ID"
    read_row "$2"
    require_enabled_row
    require_db_root
    target="$(archive_path)"
    mkdir -p "$(dirname "${target}")"
    if [[ -s "${target}" ]]; then
      echo "Archive already exists; verifying instead of downloading."
      verify_archive "${target}"
      exit 0
    fi

    if [[ "${EXPECTED_COMPRESSED_BYTES}" =~ ^[0-9]+$ ]] &&
       (( EXPECTED_COMPRESSED_BYTES > 0 )); then
      available_bytes="$(df -PB1 "${DB_ROOT}" | awk 'NR == 2 {print $4}')"
      required_bytes=$(( EXPECTED_COMPRESSED_BYTES + EXPECTED_COMPRESSED_BYTES / 5 ))
      (( available_bytes >= required_bytes )) ||
        die "insufficient free space for archive plus 20% headroom"
    fi

    part="${target}.part"
    echo "Downloading ${DB_ID} ${RELEASE_ID} to ${part}"
    download_engine="${DB_DOWNLOAD_ENGINE:-auto}"
    if [[ "${download_engine}" != "curl" ]] && command -v aria2c >/dev/null 2>&1; then
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
        --dir="$(dirname "${part}")" \
        --out="$(basename "${part}")" \
        "${ARCHIVE_URL}"
    else
      curl --fail --location \
        --retry 50 --retry-all-errors --retry-delay 5 \
        --connect-timeout 30 --speed-time 180 --speed-limit 1024 \
        --continue-at - --output "${part}" "${ARCHIVE_URL}"
    fi
    mv "${part}" "${target}"
    verify_archive "${target}"
    ;;

  verify)
    [[ "$#" -eq 2 ]] || die "verify requires DATABASE_ID"
    read_row "$2"
    require_enabled_row
    require_db_root
    verify_archive "$(archive_path)"
    ;;

  antismash8-install)
    [[ "$#" -eq 1 ]] || die "antismash8-install takes no positional arguments"
    require_db_root
    downloader="${ANTISMASH_DOWNLOADER:-download-antismash-databases}"
    command -v "${downloader}" >/dev/null 2>&1 ||
      die "activate env/antismash8.yml (or set ANTISMASH_DOWNLOADER) first"
    mkdir -p "${DB_ROOT}"
    "${downloader}" --database-dir "${DB_ROOT}"
    verify_antismash8_database
    ;;

  antismash8-verify)
    [[ "$#" -eq 1 ]] || die "antismash8-verify takes no positional arguments"
    require_db_root
    verify_antismash8_database
    ;;

  -h|--help|help)
    usage
    ;;

  *)
    usage >&2
    exit 2
    ;;
esac
