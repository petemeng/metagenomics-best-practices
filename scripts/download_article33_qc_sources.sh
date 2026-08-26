#!/usr/bin/env bash
set -euo pipefail

project_root="${1:-.}"
raw_dir="${2:-data/raw/article33}"
project_root="$(cd "$project_root" && pwd)"
mkdir -p "$raw_dir"
raw_dir="$(cd "$raw_dir" && pwd)"

repo="$raw_dir/benchmark_mock"
url="https://forgemia.inra.fr/metagenopolis/benchmark_mock.git"
commit="a429a3724d4593f35b8d7323b20252a6be90e1cd"
s1_sha="937653a56fea7fbfcbe35b3f35c721b4125072ba4ab04c44c9d454697240c6df"

if [[ ! -d "$repo/.git" ]]; then
  git clone "$url" "$repo"
fi
if ! git -C "$repo" cat-file -e "$commit^{commit}" 2>/dev/null; then
  git -C "$repo" fetch origin "$commit"
fi
git -C "$repo" checkout --detach "$commit"

observed_commit="$(git -C "$repo" rev-parse HEAD)"
[[ "$observed_commit" == "$commit" ]] || { printf 'Commit mismatch: %s\n' "$observed_commit" >&2; exit 1; }
printf '%s  %s\n' "$s1_sha" "$repo/script_r/Supplementary_Table_S1.xlsx" | sha256sum -c -

for required in \
  "$repo/profiling/MOCK_001.list" \
  "$repo/profiling/MOCK_002.list" \
  "$repo/reference/MOCK_001.fasta.gz" \
  "$repo/reference/MOCK_002.fasta.gz"; do
  [[ -s "$required" ]] || { printf 'Missing benchmark payload: %s\n' "$required" >&2; exit 1; }
done

mock1_count="$(wc -l < "$repo/profiling/MOCK_001.list")"
mock2_count="$(wc -l < "$repo/profiling/MOCK_002.list")"
genome_files="$(find "$repo/reference/all_genomes_listed" -maxdepth 1 -type f -name '*.fna.gz' | wc -l)"
[[ "$mock1_count" -eq 71 && "$mock2_count" -eq 87 && "$genome_files" -ge 87 ]] || {
  printf 'Unexpected truth inventory: MOCK1=%s MOCK2=%s genome_files=%s\n' "$mock1_count" "$mock2_count" "$genome_files" >&2
  exit 1
}

printf 'Article 33 reference sources verified at %s\n' "$repo"
printf 'Assemblies are read from checksum-locked data/small/30-32 frozen bundles; no FASTQ download is required.\n'
