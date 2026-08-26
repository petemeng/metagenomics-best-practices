#!/usr/bin/env bash
set -euo pipefail

project_root="${1:-.}"
raw_dir="${2:-data/raw/article34}"
benchmark_url="https://forgemia.inra.fr/metagenopolis/benchmark_mock.git"
benchmark_commit="a429a3724d4593f35b8d7323b20252a6be90e1cd"
paper_url="https://www.ebi.ac.uk/europepmc/webservices/rest/PMC9074274/fullTextXML"

project_root="$(cd "$project_root" && pwd)"
mkdir -p "$project_root/$raw_dir"
repo="$project_root/$raw_dir/benchmark_mock"
paper_xml="$project_root/$raw_dir/PMC9074274.fullTextXML"

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

curl -fL --retry 5 --retry-all-errors --connect-timeout 20 \
  --max-time 180 "$paper_url" -o "$paper_xml"

mock1_count="$(wc -l < "$repo/profiling/MOCK_001.list")"
mock2_count="$(wc -l < "$repo/profiling/MOCK_002.list")"
genome_files="$(find "$repo/reference/all_genomes_listed" -maxdepth 1 -type f -name '*.fna.gz' | wc -l)"

[[ "$mock1_count" -eq 71 && "$mock2_count" -eq 87 && "$genome_files" -eq 91 ]] || {
  printf 'Unexpected truth inventory: MOCK1=%s MOCK2=%s genome_files=%s\n' \
    "$mock1_count" "$mock2_count" "$genome_files" >&2
  exit 1
}

grep -q '10.1186/s40168-022-01259-2' "$paper_xml"
grep -q 'Prodigal' "$paper_xml"
grep -q 'MMseqs2' "$paper_xml"

sha256sum \
  "$repo/script_r/Supplementary_Table_S1.xlsx" \
  "$repo/reference/MOCK_002.fasta.gz" \
  "$paper_xml"
printf 'Locked benchmark commit: %s\n' "$observed_commit"
