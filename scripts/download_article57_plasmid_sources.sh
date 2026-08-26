#!/usr/bin/env bash
set -euo pipefail

project_root="${1:-.}"
raw_dir="${2:-data/raw/article57}"
benchmark_url="https://forgemia.inra.fr/metagenopolis/benchmark_mock.git"
benchmark_commit="a429a3724d4593f35b8d7323b20252a6be90e1cd"
reference_sha256="5617cc377fc503141d7a27d7c52ce874e3393e3939a9fdcbbd43fe0268c6092c"

project_root="$(cd "$project_root" && pwd)"
mkdir -p "$project_root/$raw_dir"
repo="$project_root/$raw_dir/benchmark_mock"

if [[ ! -d "$repo/.git" ]]; then
  git clone "$benchmark_url" "$repo"
fi
git -C "$repo" fetch --tags origin
git -C "$repo" checkout --detach "$benchmark_commit"

observed_commit="$(git -C "$repo" rev-parse HEAD)"
[[ "$observed_commit" == "$benchmark_commit" ]] || {
  printf 'Unexpected benchmark commit: %s\n' "$observed_commit" >&2
  exit 1
}

reference="$repo/reference/MOCK_002.fasta.gz"
observed_sha256="$(sha256sum "$reference" | awk '{print $1}')"
[[ "$observed_sha256" == "$reference_sha256" ]] || {
  printf 'Reference SHA-256 mismatch: %s\n' "$observed_sha256" >&2
  exit 1
}

read -r records plasmids < <(
  gzip -cd "$reference" | awk '
    /^>/ {records += 1; if (tolower($0) ~ /plasmid/) plasmids += 1}
    END {print records + 0, plasmids + 0}'
)
[[ "$records" -eq 399 && "$plasmids" -eq 43 ]] || {
  printf 'Unexpected reference inventory: records=%s plasmids=%s\n' "$records" "$plasmids" >&2
  exit 1
}

printf 'verified\tcommit=%s\treplicons=%s\tplasmids=%s\tsha256=%s\n' \
  "$observed_commit" "$records" "$plasmids" "$observed_sha256"
