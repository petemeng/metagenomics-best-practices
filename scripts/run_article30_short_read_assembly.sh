#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat >&2 <<'USAGE'
Usage:
  run_article30_short_read_assembly.sh \
    --project-root DIR \
    --assembly-prefix DIR \
    --read-qc-prefix DIR \
    --raw-dir DIR \
    --frozen-dir DIR

The four full ENA archives, selected/clean FASTQs, assembler work directories,
Bowtie2 indices, and transient SAM streams remain under the Git-ignored raw
directory. Compact contigs, metrics, normalized logs, and checksums are frozen.
USAGE
}

project_root=""
assembly_prefix=""
read_qc_prefix=""
raw_dir=""
frozen_dir=""

while [[ "$#" -gt 0 ]]; do
  case "$1" in
    --project-root) project_root="$2"; shift 2 ;;
    --assembly-prefix) assembly_prefix="$2"; shift 2 ;;
    --read-qc-prefix) read_qc_prefix="$2"; shift 2 ;;
    --raw-dir) raw_dir="$2"; shift 2 ;;
    --frozen-dir) frozen_dir="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown argument: $1" >&2; usage; exit 2 ;;
  esac
done

if [[ -z "${project_root}" || -z "${assembly_prefix}" || -z "${read_qc_prefix}" ||
      -z "${raw_dir}" || -z "${frozen_dir}" ]]; then
  usage
  exit 2
fi

project_root="$(cd "${project_root}" && pwd)"
assembly_prefix="$(cd "${assembly_prefix}" && pwd)"
read_qc_prefix="$(cd "${read_qc_prefix}" && pwd)"
raw_dir="$(cd "${raw_dir}" && pwd)"
if [[ -e "${frozen_dir}" ]]; then
  if [[ -n "$(find "${frozen_dir}" -mindepth 1 -maxdepth 1 -print -quit)" ]]; then
    echo "Refusing to overwrite non-empty frozen directory: ${frozen_dir}" >&2
    exit 1
  fi
fi

work_dir="${raw_dir}/work"
mkdir -p \
  "${raw_dir}/clean" \
  "${work_dir}/fastp" \
  "${work_dir}/assemblies" \
  "${work_dir}/indexes" \
  "${work_dir}/mapping" \
  "${work_dir}/logs" \
  "${work_dir}/resources" \
  "${work_dir}/tmp"

export LC_ALL=C
export TZ=UTC
export PYTHONHASHSEED=0
export PYTHONNOUSERSITE=1
export PATH="${assembly_prefix}/bin:${read_qc_prefix}/bin:/usr/bin:/bin"
export TMPDIR="${work_dir}/tmp"
export XDG_CACHE_HOME="${work_dir}/.cache"
export MPLCONFIGDIR="${work_dir}/.matplotlib"
mkdir -p "${XDG_CACHE_HOME}" "${MPLCONFIGDIR}"

fastp="${read_qc_prefix}/bin/fastp"
megahit="${assembly_prefix}/bin/megahit"
metaspades="${assembly_prefix}/bin/metaspades.py"
bowtie2="${assembly_prefix}/bin/bowtie2"
bowtie2_build="${assembly_prefix}/bin/bowtie2-build"
python="${assembly_prefix}/bin/python"
sam_summary="${project_root}/scripts/summarize_article30_sam.py"
freezer="${project_root}/scripts/freeze_article30_short_read_assembly.py"

for executable in "${fastp}" "${megahit}" "${metaspades}" "${bowtie2}" \
  "${bowtie2_build}" "${python}"; do
  if [[ ! -x "${executable}" ]]; then
    echo "Missing executable: ${executable}" >&2
    exit 1
  fi
done

selected_r1_mock1="${raw_dir}/selected/ERR9765746_selected2m_R1.fastq.gz"
selected_r2_mock1="${raw_dir}/selected/ERR9765746_selected2m_R2.fastq.gz"
selected_r1_mock2="${raw_dir}/selected/ERR9765747_selected2m_R1.fastq.gz"
selected_r2_mock2="${raw_dir}/selected/ERR9765747_selected2m_R2.fastq.gz"
for path in "${selected_r1_mock1}" "${selected_r2_mock1}" \
  "${selected_r1_mock2}" "${selected_r2_mock2}"; do
  if [[ ! -s "${path}" ]]; then
    echo "Missing selected FASTQ: ${path}" >&2
    exit 1
  fi
done

clean_r1_mock1="${raw_dir}/clean/ERR9765746_clean_R1.fastq.gz"
clean_r2_mock1="${raw_dir}/clean/ERR9765746_clean_R2.fastq.gz"
clean_r1_mock2="${raw_dir}/clean/ERR9765747_clean_R1.fastq.gz"
clean_r2_mock2="${raw_dir}/clean/ERR9765747_clean_R2.fastq.gz"

run_fastp() {
  local mock="$1" run="$2" input1="$3" input2="$4" output1="$5" output2="$6"
  local json="${work_dir}/fastp/${mock}.json"
  if [[ -s "${output1}" && -s "${output2}" && -s "${json}" ]]; then
    echo "skip fastp ${mock}: complete outputs exist"
    return
  fi
  if [[ -e "${output1}" || -e "${output2}" || -e "${json}" ]]; then
    echo "Refusing partial fastp outputs for ${mock}" >&2
    exit 1
  fi
  /usr/bin/time -v \
    -o "${work_dir}/resources/fastp-${mock}.txt" \
    "${fastp}" \
      --in1 "${input1}" --in2 "${input2}" \
      --out1 "${output1}" --out2 "${output2}" \
      --json "${json}" \
      --html "${work_dir}/fastp/${mock}.html" \
      --report_title "${run} deterministic two-million-pair subset" \
      --thread 16 --compression 6 --detect_adapter_for_pe \
      --qualified_quality_phred 20 --unqualified_percent_limit 40 \
      --n_base_limit 5 --length_required 50 --disable_trim_poly_g \
      --overrepresentation_analysis --overrepresentation_sampling 20 \
      --dont_overwrite \
      > "${work_dir}/logs/fastp-${mock}.log" 2>&1
}

run_fastp MOCK1 ERR9765746 \
  "${selected_r1_mock1}" "${selected_r2_mock1}" \
  "${clean_r1_mock1}" "${clean_r2_mock1}"
run_fastp MOCK2 ERR9765747 \
  "${selected_r1_mock2}" "${selected_r2_mock2}" \
  "${clean_r1_mock2}" "${clean_r2_mock2}"

run_assembly() {
  local branch="$1" final_name="$2"
  shift 2
  local output_dir="${work_dir}/assemblies/${branch}"
  if [[ -s "${output_dir}/${final_name}" && -s "${work_dir}/resources/assemble-${branch}.txt" ]]; then
    echo "skip assembly ${branch}: complete output exists"
    return
  fi
  if [[ -e "${output_dir}" ]]; then
    echo "Refusing partial assembly directory: ${output_dir}" >&2
    exit 1
  fi
  /usr/bin/time -v \
    -o "${work_dir}/resources/assemble-${branch}.txt" \
    "$@" \
    > "${work_dir}/logs/assemble-${branch}.log" 2>&1
  if [[ ! -s "${output_dir}/${final_name}" ]]; then
    echo "Assembly did not produce ${final_name}: ${branch}" >&2
    exit 1
  fi
}

run_assembly megahit-single-MOCK1 final.contigs.fa \
  "${megahit}" -1 "${clean_r1_mock1}" -2 "${clean_r2_mock1}" \
    --presets meta-sensitive --min-contig-len 500 --num-cpu-threads 16 \
    --memory 68719476736 -o "${work_dir}/assemblies/megahit-single-MOCK1"
run_assembly megahit-single-MOCK2 final.contigs.fa \
  "${megahit}" -1 "${clean_r1_mock2}" -2 "${clean_r2_mock2}" \
    --presets meta-sensitive --min-contig-len 500 --num-cpu-threads 16 \
    --memory 68719476736 -o "${work_dir}/assemblies/megahit-single-MOCK2"
run_assembly megahit-coassembly final.contigs.fa \
  "${megahit}" -1 "${clean_r1_mock1},${clean_r1_mock2}" \
    -2 "${clean_r2_mock1},${clean_r2_mock2}" \
    --presets meta-sensitive --min-contig-len 500 --num-cpu-threads 16 \
    --memory 68719476736 -o "${work_dir}/assemblies/megahit-coassembly"

run_assembly metaspades-single-MOCK1 contigs.fasta \
  "${metaspades}" --only-assembler \
    -1 "${clean_r1_mock1}" -2 "${clean_r2_mock1}" \
    -t 16 -m 64 -o "${work_dir}/assemblies/metaspades-single-MOCK1"
run_assembly metaspades-single-MOCK2 contigs.fasta \
  "${metaspades}" --only-assembler \
    -1 "${clean_r1_mock2}" -2 "${clean_r2_mock2}" \
    -t 16 -m 64 -o "${work_dir}/assemblies/metaspades-single-MOCK2"
run_assembly metaspades-coassembly contigs.fasta \
  "${metaspades}" --only-assembler \
    -1 "${clean_r1_mock1}" -1 "${clean_r1_mock2}" \
    -2 "${clean_r2_mock1}" -2 "${clean_r2_mock2}" \
    -t 16 -m 64 -o "${work_dir}/assemblies/metaspades-coassembly"

build_index() {
  local branch="$1" fasta="$2"
  local prefix="${work_dir}/indexes/${branch}"
  local count
  count="$(find "${work_dir}/indexes" -maxdepth 1 -type f -name "${branch}*.bt2*" -size +0c | wc -l)"
  if [[ "${count}" == "6" ]]; then
    echo "skip Bowtie2 index ${branch}: six files exist"
    return
  fi
  if [[ "${count}" != "0" ]]; then
    echo "Refusing partial Bowtie2 index: ${branch} (${count}/6)" >&2
    exit 1
  fi
  /usr/bin/time -v \
    -o "${work_dir}/resources/index-${branch}.txt" \
    "${bowtie2_build}" --threads 16 "${fasta}" "${prefix}" \
    > "${work_dir}/logs/index-${branch}.log" 2>&1
}

build_index megahit-single-MOCK1 "${work_dir}/assemblies/megahit-single-MOCK1/final.contigs.fa"
build_index megahit-single-MOCK2 "${work_dir}/assemblies/megahit-single-MOCK2/final.contigs.fa"
build_index megahit-coassembly "${work_dir}/assemblies/megahit-coassembly/final.contigs.fa"
build_index metaspades-single-MOCK1 "${work_dir}/assemblies/metaspades-single-MOCK1/contigs.fasta"
build_index metaspades-single-MOCK2 "${work_dir}/assemblies/metaspades-single-MOCK2/contigs.fasta"
build_index metaspades-coassembly "${work_dir}/assemblies/metaspades-coassembly/contigs.fasta"

clean_pairs_mock1="$(${python} -c 'import json,sys; print(json.load(open(sys.argv[1]))["summary"]["after_filtering"]["total_reads"]//2)' "${work_dir}/fastp/MOCK1.json")"
clean_pairs_mock2="$(${python} -c 'import json,sys; print(json.load(open(sys.argv[1]))["summary"]["after_filtering"]["total_reads"]//2)' "${work_dir}/fastp/MOCK2.json")"

map_reads() {
  local sample="$1" branch="$2" reads1="$3" reads2="$4" expected_pairs="$5"
  local stem="${sample}__${branch}"
  local output="${work_dir}/mapping/${stem}.json"
  if [[ -s "${output}" && -s "${work_dir}/resources/map-${stem}.txt" ]]; then
    echo "skip mapping ${stem}: complete summary exists"
    return
  fi
  if [[ -e "${output}" ]]; then
    echo "Refusing partial mapping summary: ${output}" >&2
    exit 1
  fi
  /usr/bin/time -v \
    -o "${work_dir}/resources/map-${stem}.txt" \
    "${bowtie2}" --very-sensitive -p 16 \
      --seed 20260730 \
      -x "${work_dir}/indexes/${branch}" -1 "${reads1}" -2 "${reads2}" \
      2> "${work_dir}/logs/map-${stem}.log" \
    | "${python}" "${sam_summary}" \
        --output "${output}" --sample "${sample}" \
        --assembly-branch "${branch}" --expected-pairs "${expected_pairs}"
}

map_reads MOCK1 megahit-single-MOCK1 "${clean_r1_mock1}" "${clean_r2_mock1}" "${clean_pairs_mock1}"
map_reads MOCK1 megahit-coassembly "${clean_r1_mock1}" "${clean_r2_mock1}" "${clean_pairs_mock1}"
map_reads MOCK1 metaspades-single-MOCK1 "${clean_r1_mock1}" "${clean_r2_mock1}" "${clean_pairs_mock1}"
map_reads MOCK1 metaspades-coassembly "${clean_r1_mock1}" "${clean_r2_mock1}" "${clean_pairs_mock1}"
map_reads MOCK2 megahit-single-MOCK2 "${clean_r1_mock2}" "${clean_r2_mock2}" "${clean_pairs_mock2}"
map_reads MOCK2 megahit-coassembly "${clean_r1_mock2}" "${clean_r2_mock2}" "${clean_pairs_mock2}"
map_reads MOCK2 metaspades-single-MOCK2 "${clean_r1_mock2}" "${clean_r2_mock2}" "${clean_pairs_mock2}"
map_reads MOCK2 metaspades-coassembly "${clean_r1_mock2}" "${clean_r2_mock2}" "${clean_pairs_mock2}"

{
  printf 'Tool\tVersion\tExecutable\n'
  printf 'fastp\t%s\t%s\n' "$("${fastp}" --version 2>&1 | sed 's/^fastp //')" '${READ_QC_ENV_PREFIX}/bin/fastp'
  printf 'MEGAHIT\t%s\t%s\n' "$("${megahit}" --version 2>&1 | sed 's/^MEGAHIT v//')" '${ASSEMBLY_ENV_PREFIX}/bin/megahit'
  printf 'metaSPAdes\t%s\t%s\n' "$("${metaspades}" --version 2>&1 | sed -n 's/^SPAdes genome assembler v\([^ ]*\).*/\1/p')" '${ASSEMBLY_ENV_PREFIX}/bin/metaspades.py'
  printf 'Bowtie2\t%s\t%s\n' "$("${bowtie2}" --version 2>&1 | sed -n 's/.* version \([0-9.]*\).*/\1/p' | head -n 1)" '${ASSEMBLY_ENV_PREFIX}/bin/bowtie2'
  printf 'Python\t%s\t%s\n' "$("${python}" --version 2>&1 | sed 's/^Python //')" '${ASSEMBLY_ENV_PREFIX}/bin/python'
} > "${work_dir}/tool-versions.tsv"

"${python}" "${freezer}" \
  --project-root "${project_root}" \
  --raw-dir "${raw_dir}" \
  --work-dir "${work_dir}" \
  --frozen-dir "${frozen_dir}"

echo "Article 30 short-read assembly run and freeze completed."
