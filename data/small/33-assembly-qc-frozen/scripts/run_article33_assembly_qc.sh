#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage: scripts/run_article33_assembly_qc.sh \
  --project-root DIR --qc-prefix DIR --benchmark-repo DIR --work-dir DIR [--threads N] [--force-prepare]
EOF
}

project_root=""
qc_prefix=""
benchmark_repo=""
work_dir=""
threads=32
force_prepare=0
while [[ $# -gt 0 ]]; do
  case "$1" in
    --project-root) project_root="$2"; shift 2 ;;
    --qc-prefix) qc_prefix="$2"; shift 2 ;;
    --benchmark-repo) benchmark_repo="$2"; shift 2 ;;
    --work-dir) work_dir="$2"; shift 2 ;;
    --threads) threads="$2"; shift 2 ;;
    --force-prepare) force_prepare=1; shift ;;
    -h|--help) usage; exit 0 ;;
    *) printf 'Unknown argument: %s\n' "$1" >&2; usage >&2; exit 2 ;;
  esac
done

for value in project_root qc_prefix benchmark_repo work_dir; do
  if [[ -z "${!value}" ]]; then
    printf 'Missing required argument: %s\n' "$value" >&2
    usage >&2
    exit 2
  fi
done

project_root="$(cd "$project_root" && pwd)"
qc_prefix="$(cd "$qc_prefix" && pwd)"
benchmark_repo="$(cd "$benchmark_repo" && pwd)"
mkdir -p "$work_dir"
work_dir="$(cd "$work_dir" && pwd)"

python="$qc_prefix/bin/python"
quast="$qc_prefix/bin/quast.py"
metaquast="$qc_prefix/bin/metaquast.py"
for executable in "$python" "$quast" "$metaquast"; do
  [[ -x "$executable" ]] || { printf 'Missing executable: %s\n' "$executable" >&2; exit 1; }
done
[[ "$threads" =~ ^[1-9][0-9]*$ ]] || { printf 'Threads must be a positive integer\n' >&2; exit 2; }

mkdir -p "$work_dir/logs" "$work_dir/resources" "$work_dir/summary"

if [[ "$force_prepare" -eq 1 || ! -s "$work_dir/summary/prepare-summary.json" ]]; then
  /usr/bin/time -v -o "$work_dir/resources/prepare-inputs.txt" \
    "$python" "$project_root/scripts/prepare_article33_qc_inputs.py" \
      --project-root "$project_root" \
      --benchmark-repo "$benchmark_repo" \
      --work-dir "$work_dir" \
      > "$work_dir/logs/prepare-inputs.log" 2>&1
else
  printf 'Reusing checksum-audited prepared inputs: %s\n' "$work_dir/summary/prepare-summary.json"
fi

lineage="$work_dir/summary/input-lineage.tsv"
[[ -s "$lineage" ]] || { printf 'Missing input lineage: %s\n' "$lineage" >&2; exit 1; }

mapfile -t all_branches < <("$python" -c 'import csv,sys; print("\n".join(r["Branch"] for r in csv.DictReader(open(sys.argv[1]), delimiter="\t")))' "$lineage")
assemblies=()
for branch in "${all_branches[@]}"; do
  fasta="$work_dir/assemblies/$branch.fasta"
  [[ -s "$fasta" ]] || { printf 'Missing prepared assembly: %s\n' "$fasta" >&2; exit 1; }
  assemblies+=("$fasta")
done
labels="$(IFS=,; printf '%s' "${all_branches[*]}")"

if [[ ! -f "$work_dir/quast/.article33-complete" ]]; then
  rm -rf "$work_dir/quast.partial"
  /usr/bin/time -v -o "$work_dir/resources/quast-reference-free.txt" \
    "$quast" \
      --threads "$threads" --min-contig 1000 --no-icarus --no-plots \
      --space-efficient --report-all-metrics \
      --labels "$labels" \
      "${assemblies[@]}" \
      -o "$work_dir/quast.partial" \
      > "$work_dir/logs/quast-reference-free.log" 2>&1
  rm -rf "$work_dir/quast"
  mv "$work_dir/quast.partial" "$work_dir/quast"
  touch "$work_dir/quast/.article33-complete"
else
  printf 'Reusing completed reference-free QUAST output\n'
fi

for evaluation_set in MOCK1 MOCK2 'MOCK1+MOCK2'; do
  safe_set="${evaluation_set//+/_}"
  output="$work_dir/metaquast/$safe_set"
  if [[ -f "$output/.article33-complete" ]]; then
    printf 'Reusing completed MetaQUAST output: %s\n' "$evaluation_set"
    continue
  fi
  mapfile -t group_branches < <("$python" -c 'import csv,sys; print("\n".join(r["Branch"] for r in csv.DictReader(open(sys.argv[1]), delimiter="\t") if r["EvaluationSet"] == sys.argv[2]))' "$lineage" "$evaluation_set")
  group_assemblies=()
  for branch in "${group_branches[@]}"; do
    group_assemblies+=("$work_dir/assemblies/$branch.fasta")
  done
  group_labels="$(IFS=,; printf '%s' "${group_branches[*]}")"
  reference_csv="$(find "$work_dir/truth/references/$evaluation_set" -maxdepth 1 -type l -name '*.fna.gz' -print | sort | paste -sd, -)"
  [[ -n "$reference_csv" ]] || { printf 'No references for %s\n' "$evaluation_set" >&2; exit 1; }
  mkdir -p "$work_dir/metaquast"
  rm -rf "$output.partial"
  /usr/bin/time -v -o "$work_dir/resources/metaquast-$safe_set.txt" \
    "$metaquast" \
      --min-alignment 500 --fragmented --min-identity 97 \
      --split-scaffolds --threads "$threads" --min-contig 1000 \
      --no-icarus --no-plots --space-efficient --report-all-metrics \
      -r "$reference_csv" \
      "${group_assemblies[@]}" \
      --labels "$group_labels" \
      -o "$output.partial" \
      > "$work_dir/logs/metaquast-$safe_set.log" 2>&1
  mv "$output.partial" "$output"
  touch "$output/.article33-complete"
done

{
  printf 'Tool\tVersion\tExecutable\n'
  printf 'QUAST/MetaQUAST\t%s\t%s\n' "$($quast --version 2>&1 | sed 's/^QUAST v//')" '${QC_ENV_PREFIX}/bin/quast.py'
  printf 'Python\t%s\t%s\n' "$($python --version 2>&1 | sed 's/^Python //')" '${QC_ENV_PREFIX}/bin/python'
} > "$work_dir/tool-versions.tsv"

"$python" "$project_root/scripts/summarize_article33_assembly_qc.py" \
  --work-dir "$work_dir" \
  --output-dir "$work_dir/summary"

printf 'Article 33 compute completed: %s\n' "$work_dir/summary/run-summary.json"
