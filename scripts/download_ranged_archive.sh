#!/usr/bin/env bash
set -euo pipefail

usage() {
  printf '%s\n' \
    'Usage: download_ranged_archive.sh --url URL --bytes N --output FILE [--segments N]' >&2
}

url=""
expected_bytes=""
output=""
segments=16
while [[ "$#" -gt 0 ]]; do
  case "$1" in
    --url)
      url="$2"
      shift 2
      ;;
    --bytes)
      expected_bytes="$2"
      shift 2
      ;;
    --output)
      output="$2"
      shift 2
      ;;
    --segments)
      segments="$2"
      shift 2
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

if [[ -z "${url}" || -z "${expected_bytes}" || -z "${output}" ]]; then
  usage
  exit 2
fi
if [[ ! "${expected_bytes}" =~ ^[1-9][0-9]*$ ]] ||
   [[ ! "${segments}" =~ ^[1-9][0-9]*$ ]]; then
  printf 'Bytes and segments must be positive integers\n' >&2
  exit 2
fi
if (( segments > 64 )); then
  printf 'Refusing more than 64 parallel range segments\n' >&2
  exit 2
fi

mkdir -p "$(dirname "${output}")"
if [[ -s "${output}" ]]; then
  observed="$(stat -c '%s' "${output}")"
  if [[ "${observed}" == "${expected_bytes}" ]]; then
    printf 'Range download already complete: %s bytes=%s\n' \
      "${output}" "${observed}"
    exit 0
  fi
  printf 'Refusing to overwrite wrong-sized output: %s bytes=%s expected=%s\n' \
    "${output}" "${observed}" "${expected_bytes}" >&2
  exit 1
fi

segment_root="${output}.range-parts"
mkdir -p "${segment_root}"
chunk_size=$(( (expected_bytes + segments - 1) / segments ))
pids=()

for (( index=0; index<segments; index++ )); do
  start=$(( index * chunk_size ))
  if (( start >= expected_bytes )); then
    break
  fi
  end=$(( start + chunk_size - 1 ))
  if (( end >= expected_bytes )); then
    end=$(( expected_bytes - 1 ))
  fi
  segment_bytes=$(( end - start + 1 ))
  printf -v label '%03d' "${index}"
  segment="${segment_root}/segment-${label}.bin"
  temporary="${segment}.downloading"
  (
    if [[ -s "${segment}" ]] &&
       [[ "$(stat -c '%s' "${segment}")" == "${segment_bytes}" ]]; then
      printf 'Range %s already complete: bytes=%s-%s\n' \
        "${label}" "${start}" "${end}"
      exit 0
    fi
    http_code="$(
      curl --fail --silent --show-error --location \
        --retry 100 --retry-all-errors --retry-delay 2 \
        --connect-timeout 30 --max-time 0 \
        --range "${start}-${end}" \
        --output "${temporary}" \
        --write-out '%{http_code}' \
        "${url}"
    )"
    if [[ "${http_code}" != "206" ]]; then
      printf 'Range %s returned HTTP %s instead of 206\n' \
        "${label}" "${http_code}" >&2
      exit 1
    fi
    observed="$(stat -c '%s' "${temporary}")"
    if [[ "${observed}" != "${segment_bytes}" ]]; then
      printf 'Range %s byte mismatch: expected %s observed %s\n' \
        "${label}" "${segment_bytes}" "${observed}" >&2
      exit 1
    fi
    mv "${temporary}" "${segment}"
    printf 'Range %s complete: bytes=%s-%s size=%s\n' \
      "${label}" "${start}" "${end}" "${observed}"
  ) &
  pids+=("$!")
done

while :; do
  active=0
  for pid in "${pids[@]}"; do
    if kill -0 "${pid}" 2>/dev/null; then
      active=$(( active + 1 ))
    fi
  done
  complete_bytes=0
  for path in "${segment_root}"/segment-*.bin \
              "${segment_root}"/segment-*.bin.downloading; do
    if [[ -f "${path}" ]]; then
      complete_bytes=$(( complete_bytes + $(stat -c '%s' "${path}") ))
    fi
  done
  printf 'Range download progress: %s/%s bytes; active=%s\n' \
    "${complete_bytes}" "${expected_bytes}" "${active}"
  if (( active == 0 )); then
    break
  fi
  sleep 30
done

failed=0
for pid in "${pids[@]}"; do
  if ! wait "${pid}"; then
    failed=$(( failed + 1 ))
  fi
done
if (( failed > 0 )); then
  printf 'Range download failed in %s worker(s); completed segments are reusable\n' \
    "${failed}" >&2
  exit 1
fi

assembling="${output}.assembling.$$"
: > "${assembling}"
for (( index=0; index<segments; index++ )); do
  start=$(( index * chunk_size ))
  if (( start >= expected_bytes )); then
    break
  fi
  printf -v label '%03d' "${index}"
  segment="${segment_root}/segment-${label}.bin"
  if [[ ! -s "${segment}" ]]; then
    printf 'Missing completed range segment: %s\n' "${segment}" >&2
    exit 1
  fi
  cat "${segment}" >> "${assembling}"
done

observed="$(stat -c '%s' "${assembling}")"
if [[ "${observed}" != "${expected_bytes}" ]]; then
  printf 'Assembled byte mismatch: expected %s observed %s\n' \
    "${expected_bytes}" "${observed}" >&2
  exit 1
fi
mv "${assembling}" "${output}"
printf 'Range download assembled: %s bytes=%s segments=%s\n' \
  "${output}" "${observed}" "${segments}"

