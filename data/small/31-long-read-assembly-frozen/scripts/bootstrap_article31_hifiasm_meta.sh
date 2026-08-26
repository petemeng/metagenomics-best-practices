#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat >&2 <<'USAGE'
Usage:
  bootstrap_article31_hifiasm_meta.sh \
    --project-root DIR --env-prefix DIR --raw-dir DIR
USAGE
}

project_root=""
env_prefix=""
raw_dir=""

while [[ "$#" -gt 0 ]]; do
  case "$1" in
    --project-root) project_root="$2"; shift 2 ;;
    --env-prefix) env_prefix="$2"; shift 2 ;;
    --raw-dir) raw_dir="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown argument: $1" >&2; usage; exit 2 ;;
  esac
done

if [[ -z "${project_root}" || -z "${env_prefix}" || -z "${raw_dir}" ]]; then
  usage
  exit 2
fi

project_root="$(cd "${project_root}" && pwd)"
env_prefix="$(cd "${env_prefix}" && pwd)"
mkdir -p "${raw_dir}/tools"
raw_dir="$(cd "${raw_dir}" && pwd)"

tag="hamtv0.3.5"
release="0.3.5-r81"
commit="e4e24f5158091ad901c1ff6f68278559bd41a6b5"
url="https://github.com/xfengnefx/hifiasm-meta/archive/refs/tags/${tag}.tar.gz"
expected_sha256="8c1c1f394e0d4d3be2c78cb76c4122dd0caf2d088a8b986c161d4d90c194f560"
archive="${raw_dir}/tools/hifiasm-meta-${tag}.tar.gz"
source_dir="${raw_dir}/tools/hifiasm-meta-${tag}"
binary_dir="${raw_dir}/tools/bin"
binary="${binary_dir}/hifiasm_meta"
record="${raw_dir}/tools/hifiasm-meta-source.tsv"

if [[ ! -f "${archive}" ]]; then
  if command -v aria2c >/dev/null 2>&1; then
    aria2c \
      --continue=true --allow-overwrite=true --auto-file-renaming=false \
      --file-allocation=none --max-connection-per-server=4 --split=4 \
      --max-tries=8 --retry-wait=10 --timeout=60 \
      --dir "$(dirname "${archive}")" \
      --out "$(basename "${archive}").part" \
      "${url}"
  else
    curl --fail --location --continue-at - \
      --retry 8 --retry-all-errors --retry-delay 5 \
      --connect-timeout 30 \
      --output "${archive}.part" "${url}"
  fi
  mv "${archive}.part" "${archive}"
fi

observed_sha256="$(sha256sum "${archive}" | cut -d ' ' -f 1)"
if [[ "${observed_sha256}" != "${expected_sha256}" ]]; then
  echo "hifiasm-meta source checksum mismatch: ${archive}" >&2
  exit 1
fi

if [[ ! -d "${source_dir}" ]]; then
  tar -xzf "${archive}" -C "${raw_dir}/tools"
fi
if ! grep -Fq '{ "tsne-seed", ko_required_argument, 425 }' \
    "${source_dir}/CommandLines.cpp"; then
  echo "Expected --tsne-seed parser entry is absent from locked source" >&2
  exit 1
fi
if ! grep -Fq 'tsne_randomseed = 42;' "${source_dir}/CommandLines.cpp"; then
  echo "Expected t-SNE default seed is absent from locked source" >&2
  exit 1
fi

cxx="${env_prefix}/bin/x86_64-conda-linux-gnu-c++"
cc="${env_prefix}/bin/x86_64-conda-linux-gnu-cc"
make_bin="${env_prefix}/bin/make"
if [[ ! -x "${cxx}" || ! -x "${cc}" || ! -x "${make_bin}" ]]; then
  echo "The locked compiler toolchain is incomplete under ${env_prefix}" >&2
  exit 1
fi

if [[ ! -x "${binary}" ]]; then
  "${make_bin}" -C "${source_dir}" -j 16 \
    CC="${cc}" \
    CXX="${cxx}" \
    INCLUDES="-I${env_prefix}/include" \
    LIBS="-L${env_prefix}/lib -Wl,-rpath,${env_prefix}/lib -lz -lpthread -lm"
  mkdir -p "${binary_dir}"
  install -m 0755 "${source_dir}/hifiasm_meta" "${binary}"
fi

version_output="$("${binary}" --version 2>&1 | paste -sd ';' -)"
compiler_output="$("${cxx}" --version | head -n 1)"
binary_sha256="$(sha256sum "${binary}" | cut -d ' ' -f 1)"

{
  printf 'Tool\tRelease\tTag\tCommit\tSourceURL\tSourceSHA256\tBinarySHA256\tCompiler\tVersionOutput\tTSNESeedOption\tTSNESeedDefault\n'
  printf 'hifiasm-meta\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "${release}" "${tag}" "${commit}" "${url}" "${expected_sha256}" \
    "${binary_sha256}" "${compiler_output}" "${version_output}" \
    '--tsne-seed' '42'
} > "${record}"

echo "verified hifiasm-meta ${release}: ${binary}"
