#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat >&2 <<'USAGE'
Usage:
  run_article32_hybrid_assembly.sh \
    --project-root DIR \
    --hybrid-prefix DIR \
    --raw-dir DIR

The script is restart-aware. SPAdes branches use --continue after an interrupted
run; atomic temporary outputs protect non-resumable alignment and MetaQUAST
steps. Raw reads, references, SAM files, and work trees remain Git-ignored.
USAGE
}

project_root=""
hybrid_prefix=""
raw_dir=""
while [[ "$#" -gt 0 ]]; do
  case "$1" in
    --project-root) project_root="$2"; shift 2 ;;
    --hybrid-prefix) hybrid_prefix="$2"; shift 2 ;;
    --raw-dir) raw_dir="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown argument: $1" >&2; usage; exit 2 ;;
  esac
done
if [[ -z "${project_root}" || -z "${hybrid_prefix}" || -z "${raw_dir}" ]]; then
  usage
  exit 2
fi

project_root="$(cd "${project_root}" && pwd)"
hybrid_prefix="$(cd "${hybrid_prefix}" && pwd)"
raw_dir="$(cd "${raw_dir}" && pwd)"
work_dir="${raw_dir}/work"
selected_dir="${raw_dir}/selected"
clean_dir="${raw_dir}/clean"
truth_dir="${work_dir}/truth"
assemblies_dir="${work_dir}/assemblies"
normalized_base="${work_dir}/normalized/base"
normalized_final="${work_dir}/normalized/final"
polish_dir="${work_dir}/polish"
logs_dir="${work_dir}/logs"
resources_dir="${work_dir}/resources"
tmp_dir="${work_dir}/tmp"
mkdir -p "${selected_dir}" "${clean_dir}" "${truth_dir}" "${assemblies_dir}" \
  "${normalized_base}" "${normalized_final}" "${polish_dir}" "${logs_dir}" \
  "${resources_dir}" "${tmp_dir}"

export LC_ALL=C
export TZ=UTC
export PYTHONHASHSEED=0
export PYTHONNOUSERSITE=1
export PATH="${hybrid_prefix}/bin:/usr/bin:/bin"
export TMPDIR="${tmp_dir}"
export XDG_CACHE_HOME="${work_dir}/.cache"
export MPLCONFIGDIR="${work_dir}/.matplotlib"
export MIMALLOC_VERBOSE=0
mkdir -p "${XDG_CACHE_HOME}" "${MPLCONFIGDIR}"

python="${hybrid_prefix}/bin/python"
fastp="${hybrid_prefix}/bin/fastp"
spades="${hybrid_prefix}/bin/spades.py"
bwa="${hybrid_prefix}/bin/bwa"
polypolish="${hybrid_prefix}/bin/polypolish"
metaquast="${hybrid_prefix}/bin/metaquast.py"
for executable in "${python}" "${fastp}" "${spades}" "${bwa}" "${polypolish}" "${metaquast}"; do
  if [[ ! -x "${executable}" ]]; then
    echo "Missing executable: ${executable}" >&2
    exit 1
  fi
done

full_r1="${raw_dir}/sources/ERR9765746_R1.fastq.gz"
full_r2="${raw_dir}/sources/ERR9765746_R2.fastq.gz"
ont_reads="${raw_dir}/sources/ERR9765780.fastq.gz"
hifi_reads="${raw_dir}/sources/ERR9765783.fastq.gz"
ont_flye="${project_root}/data/small/31-long-read-assembly-frozen/assemblies/flye-ont-r9.ge1000.fna.gz"
hifi_flye="${project_root}/data/small/31-long-read-assembly-frozen/assemblies/flye-hifi.ge1000.fna.gz"
benchmark_repo="${raw_dir}/benchmark_mock"
supplement_s2="${raw_dir}/Supplementary_Table_S2.xlsx"
for path in "${full_r1}" "${full_r2}" "${ont_reads}" "${hifi_reads}" \
  "${ont_flye}" "${hifi_flye}" "${supplement_s2}"; do
  if [[ ! -s "${path}" ]]; then
    echo "Missing required input: ${path}" >&2
    exit 1
  fi
done

selected_r1="${selected_dir}/ERR9765746_selected10m_R1.fastq.gz"
selected_r2="${selected_dir}/ERR9765746_selected10m_R2.fastq.gz"
selection_summary="${selected_dir}/ERR9765746_selection-summary.json"
if [[ -s "${selected_r1}" && -s "${selected_r2}" && -s "${selection_summary}" ]]; then
  echo "skip exact 10M-pair selection: complete outputs exist"
elif [[ -e "${selected_r1}" || -e "${selected_r2}" || -e "${selection_summary}" ]]; then
  echo "Refusing partial read-selection outputs" >&2
  exit 1
else
  /usr/bin/time -v -o "${resources_dir}/select-illumina-10m.txt" \
    "${python}" "${project_root}/scripts/select_article32_read_pairs.py" \
      --r1 "${full_r1}" --r2 "${full_r2}" \
      --output-r1 "${selected_r1}" --output-r2 "${selected_r2}" \
      --summary "${selection_summary}" --total-pairs 20597525 \
      --target-pairs 10000000 --seed 20260732 \
    > "${logs_dir}/select-illumina-10m.log" 2>&1
fi

clean_r1="${clean_dir}/ERR9765746_clean10m_R1.fastq.gz"
clean_r2="${clean_dir}/ERR9765746_clean10m_R2.fastq.gz"
fastp_json="${work_dir}/fastp.json"
if [[ -s "${clean_r1}" && -s "${clean_r2}" && -s "${fastp_json}" ]]; then
  echo "skip fastp: complete outputs exist"
elif [[ -e "${clean_r1}" || -e "${clean_r2}" || -e "${fastp_json}" ]]; then
  echo "Refusing partial fastp outputs" >&2
  exit 1
else
  /usr/bin/time -v -o "${resources_dir}/fastp-illumina-10m.txt" \
    "${fastp}" --in1 "${selected_r1}" --in2 "${selected_r2}" \
      --out1 "${clean_r1}" --out2 "${clean_r2}" \
      --json "${fastp_json}" --html "${work_dir}/fastp.html" \
      --report_title "ERR9765746 deterministic ten-million-pair subset" \
      --thread 32 --compression 6 --detect_adapter_for_pe \
      --qualified_quality_phred 20 --unqualified_percent_limit 40 \
      --n_base_limit 5 --length_required 50 --disable_trim_poly_g \
      --overrepresentation_analysis --overrepresentation_sampling 20 \
      --dont_overwrite \
    > "${logs_dir}/fastp-illumina-10m.log" 2>&1
fi

"${python}" "${project_root}/scripts/prepare_article32_truth.py" \
  --benchmark-repo "${benchmark_repo}" --supplement-s2 "${supplement_s2}" \
  --output-dir "${truth_dir}" \
  > "${logs_dir}/prepare-truth.log" 2>&1

next_attempt() {
  local stem="$1"
  local count
  count="$(find "${resources_dir}" -maxdepth 1 -type f -name "${stem}.attempt*.txt" | wc -l)"
  printf '%03d' "$((count + 1))"
}

run_spades() {
  local branch="$1"
  local long_kind="$2"
  local long_path="$3"
  local output="${assemblies_dir}/${branch}"
  if [[ -s "${output}/scaffolds.fasta" && -f "${output}/.article32-complete" ]]; then
    echo "skip SPAdes ${branch}: complete output exists"
    return
  fi
  local attempt
  attempt="$(next_attempt "assemble-${branch}")"
  if [[ -d "${output}" ]]; then
    /usr/bin/time -v -o "${resources_dir}/assemble-${branch}.attempt${attempt}.txt" \
      "${spades}" --continue -o "${output}" \
      > "${logs_dir}/assemble-${branch}.attempt${attempt}.log" 2>&1
  else
    local long_args=()
    if [[ "${long_kind}" == "ont" ]]; then
      long_args=(--nanopore "${long_path}")
    elif [[ "${long_kind}" == "hifi" ]]; then
      long_args=(--pacbio "${long_path}")
    elif [[ "${long_kind}" != "none" ]]; then
      echo "Unknown long-read kind: ${long_kind}" >&2
      exit 1
    fi
    /usr/bin/time -v -o "${resources_dir}/assemble-${branch}.attempt${attempt}.txt" \
      "${spades}" --meta --only-assembler \
        -1 "${clean_r1}" -2 "${clean_r2}" "${long_args[@]}" \
        -k 21,33,55,77 -t 32 -m 256 -o "${output}" \
      > "${logs_dir}/assemble-${branch}.attempt${attempt}.log" 2>&1
  fi
  if [[ ! -s "${output}/scaffolds.fasta" ]]; then
    echo "SPAdes did not produce scaffolds.fasta: ${branch}" >&2
    exit 1
  fi
  touch "${output}/.article32-complete"
}

run_spades spades-short-only none ""
run_spades spades-illumina-ont ont "${ont_reads}"
run_spades spades-illumina-hifi hifi "${hifi_reads}"

"${python}" "${project_root}/scripts/prepare_article32_assemblies.py" \
  --assembly "spades-short-only=${assemblies_dir}/spades-short-only/scaffolds.fasta" \
  --assembly "spades-illumina-ont=${assemblies_dir}/spades-illumina-ont/scaffolds.fasta" \
  --assembly "spades-illumina-hifi=${assemblies_dir}/spades-illumina-hifi/scaffolds.fasta" \
  --assembly "flye-ont=${ont_flye}" --assembly "flye-hifi=${hifi_flye}" \
  --output-dir "${normalized_base}" --min-length 1000

polish_draft="${normalized_base}/flye-ont.ge1000.fasta"
index_marker="${polish_dir}/bwa-index.complete"
if [[ ! -f "${index_marker}" ]]; then
  /usr/bin/time -v -o "${resources_dir}/bwa-index-flye-ont.txt" \
    "${bwa}" index "${polish_draft}" \
    > "${logs_dir}/bwa-index-flye-ont.log" 2>&1
  touch "${index_marker}"
fi

align_mate() {
  local mate="$1"
  local reads="$2"
  local final="${polish_dir}/alignments_R${mate}.sam"
  local marker="${polish_dir}/alignments_R${mate}.complete"
  if [[ -s "${final}" && -f "${marker}" ]]; then
    echo "skip BWA R${mate}: complete SAM exists"
    return
  fi
  local temporary="${final}.partial"
  rm -f "${temporary}"
  /usr/bin/time -v -o "${resources_dir}/bwa-align-R${mate}.txt" \
    "${bwa}" mem -t 32 -a "${polish_draft}" "${reads}" \
    > "${temporary}" 2> "${logs_dir}/bwa-align-R${mate}.log"
  mv "${temporary}" "${final}"
  touch "${marker}"
}
align_mate 1 "${clean_r1}"
align_mate 2 "${clean_r2}"

filtered_r1="${polish_dir}/filtered_R1.sam"
filtered_r2="${polish_dir}/filtered_R2.sam"
if [[ ! -s "${filtered_r1}" || ! -s "${filtered_r2}" || ! -f "${polish_dir}/filter.complete" ]]; then
  rm -f "${filtered_r1}.partial" "${filtered_r2}.partial"
  /usr/bin/time -v -o "${resources_dir}/polypolish-filter.txt" \
    "${polypolish}" filter \
      --in1 "${polish_dir}/alignments_R1.sam" \
      --in2 "${polish_dir}/alignments_R2.sam" \
      --out1 "${filtered_r1}.partial" --out2 "${filtered_r2}.partial" \
    > "${logs_dir}/polypolish-filter.log" 2>&1
  mv "${filtered_r1}.partial" "${filtered_r1}"
  mv "${filtered_r2}.partial" "${filtered_r2}"
  touch "${polish_dir}/filter.complete"
fi

run_polypolish() {
  local mode="$1"
  local careful=()
  if [[ "${mode}" == "careful" ]]; then careful=(--careful); fi
  local final="${polish_dir}/flye-ont-polypolish-${mode}.fasta"
  if [[ -s "${final}" && -f "${polish_dir}/polypolish-${mode}.complete" ]]; then
    echo "skip Polypolish ${mode}: complete output exists"
    return
  fi
  local temporary="${final}.partial"
  rm -f "${temporary}"
  /usr/bin/time -v -o "${resources_dir}/polypolish-${mode}.txt" \
    "${polypolish}" polish "${careful[@]}" "${polish_draft}" \
      "${filtered_r1}" "${filtered_r2}" \
    > "${temporary}" 2> "${logs_dir}/polypolish-${mode}.log"
  mv "${temporary}" "${final}"
  touch "${polish_dir}/polypolish-${mode}.complete"
}
run_polypolish default
run_polypolish careful

"${python}" "${project_root}/scripts/prepare_article32_assemblies.py" \
  --assembly "spades-short-only=${assemblies_dir}/spades-short-only/scaffolds.fasta" \
  --assembly "spades-illumina-ont=${assemblies_dir}/spades-illumina-ont/scaffolds.fasta" \
  --assembly "spades-illumina-hifi=${assemblies_dir}/spades-illumina-hifi/scaffolds.fasta" \
  --assembly "flye-ont=${ont_flye}" \
  --assembly "flye-ont-polypolish-default=${polish_dir}/flye-ont-polypolish-default.fasta" \
  --assembly "flye-ont-polypolish-careful=${polish_dir}/flye-ont-polypolish-careful.fasta" \
  --assembly "flye-hifi=${hifi_flye}" \
  --output-dir "${normalized_final}" --min-length 1000

metaquast_output="${work_dir}/metaquast"
if [[ -s "${metaquast_output}/combined_reference/transposed_report.tsv" && -f "${metaquast_output}/.article32-complete" ]]; then
  echo "skip MetaQUAST: complete output exists"
else
  if [[ -e "${metaquast_output}" ]]; then
    echo "Refusing incomplete MetaQUAST output: ${metaquast_output}" >&2
    echo "Inspect it, then remove or archive it before restarting." >&2
    exit 1
  fi
  reference_csv="$(find "${truth_dir}/references/MOCK1" -maxdepth 1 -type l -name '*.fna.gz' -print | sort | paste -sd, -)"
  if [[ -z "${reference_csv}" ]]; then
    echo "No MOCK1 reference files" >&2
    exit 1
  fi
  partial="${metaquast_output}.partial"
  rm -rf "${partial}"
  /usr/bin/time -v -o "${resources_dir}/metaquast-seven-branches.txt" \
    "${metaquast}" --min-alignment 500 --fragmented --min-identity 97 \
      --split-scaffolds --threads 32 --min-contig 1000 --no-icarus \
      -r "${reference_csv}" \
      "${normalized_final}/spades-short-only.ge1000.fasta" \
      "${normalized_final}/spades-illumina-ont.ge1000.fasta" \
      "${normalized_final}/spades-illumina-hifi.ge1000.fasta" \
      "${normalized_final}/flye-ont.ge1000.fasta" \
      "${normalized_final}/flye-ont-polypolish-default.ge1000.fasta" \
      "${normalized_final}/flye-ont-polypolish-careful.ge1000.fasta" \
      "${normalized_final}/flye-hifi.ge1000.fasta" \
      --labels spades-short-only,spades-illumina-ont,spades-illumina-hifi,flye-ont,flye-ont-polypolish-default,flye-ont-polypolish-careful,flye-hifi \
      -o "${partial}" \
    > "${logs_dir}/metaquast-seven-branches.log" 2>&1
  mv "${partial}" "${metaquast_output}"
  touch "${metaquast_output}/.article32-complete"
fi

"${python}" "${project_root}/scripts/summarize_article32_metaquast.py" \
  --metaquast-dir "${metaquast_output}" \
  --normalized-dir "${normalized_final}" \
  --truth-manifest "${truth_dir}/truth-manifest.tsv" \
  --resource-dir "${resources_dir}" --log-dir "${logs_dir}" \
  --output-dir "${work_dir}/summary" \
  > "${logs_dir}/summarize-metaquast.log" 2>&1

{
  printf 'Tool\tVersion\tExecutable\n'
  printf 'fastp\t%s\t%s\n' "$("${fastp}" --version 2>&1 | sed 's/^fastp //')" '${HYBRID_ENV_PREFIX}/bin/fastp'
  printf 'SPAdes\t%s\t%s\n' "$("${spades}" --version 2>&1 | sed -n 's/^SPAdes genome assembler v\([^ ]*\).*/\1/p')" '${HYBRID_ENV_PREFIX}/bin/spades.py'
  printf 'BWA-MEM\t%s\t%s\n' "$("${bwa}" 2>&1 | sed -n 's/^Version: //p')" '${HYBRID_ENV_PREFIX}/bin/bwa'
  printf 'Polypolish\t%s\t%s\n' "$("${polypolish}" --version 2>&1 | awk '{print $2}')" '${HYBRID_ENV_PREFIX}/bin/polypolish'
  printf 'MetaQUAST\t%s\t%s\n' "$("${metaquast}" --version 2>&1 | sed 's/^QUAST v//; s/ (MetaQUAST mode)//')" '${HYBRID_ENV_PREFIX}/bin/metaquast.py'
  printf 'Python\t%s\t%s\n' "$("${python}" --version 2>&1 | sed 's/^Python //')" '${HYBRID_ENV_PREFIX}/bin/python'
} > "${work_dir}/tool-versions.tsv"

echo "Article 32 hybrid assembly and polishing computation completed."
