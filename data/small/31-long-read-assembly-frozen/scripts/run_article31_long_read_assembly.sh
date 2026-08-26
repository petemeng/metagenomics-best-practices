#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat >&2 <<'USAGE'
Usage:
  run_article31_long_read_assembly.sh \
    --project-root DIR --env-prefix DIR --raw-dir DIR --frozen-dir DIR

The complete ENA archives, assembler work directories, PAF files, and source-
built hifiasm-meta binary remain under the Git-ignored raw directory. Compact
assemblies, normalized evidence tables, logs, resource records, and checksums
are frozen only after every branch has completed.
USAGE
}

project_root=""
env_prefix=""
raw_dir=""
frozen_dir=""

while [[ "$#" -gt 0 ]]; do
  case "$1" in
    --project-root) project_root="$2"; shift 2 ;;
    --env-prefix) env_prefix="$2"; shift 2 ;;
    --raw-dir) raw_dir="$2"; shift 2 ;;
    --frozen-dir) frozen_dir="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown argument: $1" >&2; usage; exit 2 ;;
  esac
done

if [[ -z "${project_root}" || -z "${env_prefix}" || -z "${raw_dir}" || -z "${frozen_dir}" ]]; then
  usage
  exit 2
fi

project_root="$(cd "${project_root}" && pwd)"
env_prefix="$(cd "${env_prefix}" && pwd)"
raw_dir="$(cd "${raw_dir}" && pwd)"
if [[ -e "${frozen_dir}" ]] && [[ -n "$(find "${frozen_dir}" -mindepth 1 -maxdepth 1 -print -quit)" ]]; then
  echo "Refusing to overwrite non-empty frozen directory: ${frozen_dir}" >&2
  exit 1
fi

work_dir="${raw_dir}/work"
mkdir -p \
  "${work_dir}/assemblies" \
  "${work_dir}/logs" \
  "${work_dir}/resources" \
  "${work_dir}/mapping" \
  "${work_dir}/junction-mapping" \
  "${work_dir}/normalized" \
  "${work_dir}/tmp"

export LC_ALL=C
export TZ=UTC
export PYTHONHASHSEED=0
export PYTHONNOUSERSITE=1
export OMP_NUM_THREADS=32
export TMPDIR="${work_dir}/tmp"
export XDG_CACHE_HOME="${work_dir}/.cache"
export MPLCONFIGDIR="${work_dir}/.matplotlib"
export PATH="${env_prefix}/bin:/usr/bin:/bin"
mkdir -p "${XDG_CACHE_HOME}" "${MPLCONFIGDIR}"

flye="${env_prefix}/bin/flye"
metamdbg="${env_prefix}/bin/metaMDBG"
minimap2="${env_prefix}/bin/minimap2"
samtools="${env_prefix}/bin/samtools"
seqkit="${env_prefix}/bin/seqkit"
python="${env_prefix}/bin/python"
hifiasm="${raw_dir}/tools/bin/hifiasm_meta"
prepare="${project_root}/scripts/prepare_article31_assemblies.py"
paf_summary="${project_root}/scripts/summarize_article31_paf.py"
junction_summary="${project_root}/scripts/summarize_article31_junctions.py"
freezer="${project_root}/scripts/freeze_article31_long_read_assembly.py"

for executable in "${flye}" "${metamdbg}" "${minimap2}" "${samtools}" \
  "${seqkit}" "${python}" "${hifiasm}"; do
  if [[ ! -x "${executable}" ]]; then
    echo "Missing executable: ${executable}" >&2
    exit 1
  fi
done
for script in "${prepare}" "${paf_summary}" "${junction_summary}" "${freezer}"; do
  if [[ ! -f "${script}" ]]; then
    echo "Missing Article 31 script: ${script}" >&2
    exit 1
  fi
done

ont_reads="${raw_dir}/full/ERR9765780.fastq.gz"
hifi_reads="${raw_dir}/full/ERR9765783.fastq.gz"
for fastq in "${ont_reads}" "${hifi_reads}"; do
  if [[ ! -s "${fastq}" ]]; then
    echo "Missing complete FASTQ archive: ${fastq}" >&2
    exit 1
  fi
done

next_attempt() {
  local label="$1"
  local count
  count="$(find "${work_dir}/resources" -maxdepth 1 -type f -name "${label}.attempt*.txt" | wc -l)"
  printf '%03d' "$((count + 1))"
}

run_timed() {
  local label="$1"
  shift
  local attempt resource log
  attempt="$(next_attempt "${label}")"
  resource="${work_dir}/resources/${label}.attempt${attempt}.txt"
  log="${work_dir}/logs/${label}.attempt${attempt}.log"
  echo "start ${label} attempt ${attempt}"
  /usr/bin/time -v -o "${resource}" "$@" > "${log}" 2>&1
  echo "finish ${label} attempt ${attempt}"
}

run_flye() {
  local branch="$1"
  local read_flag="$2"
  local reads="$3"
  local iterations="$4"
  local output_dir="${work_dir}/assemblies/${branch}"
  local complete="${output_dir}/.article31-complete"
  if [[ -s "${output_dir}/assembly.fasta" && -s "${output_dir}/assembly_info.txt" && -f "${complete}" ]]; then
    echo "skip assembly ${branch}: complete sentinel and outputs exist"
    return
  fi
  local resume=()
  if [[ -s "${output_dir}/flye.log" ]]; then
    resume=(--resume)
  fi
  run_timed "assemble-${branch}" \
    "${flye}" "${read_flag}" "${reads}" \
      --out-dir "${output_dir}" \
      --threads 32 \
      --meta \
      --min-overlap 2000 \
      --iterations "${iterations}" \
      --deterministic \
      "${resume[@]}"
  test -s "${output_dir}/assembly.fasta"
  test -s "${output_dir}/assembly_info.txt"
  touch "${complete}"
}

run_flye flye-ont-r9 --nano-raw "${ont_reads}" 2
run_flye flye-hifi --pacbio-hifi "${hifi_reads}" 1

hifiasm_dir="${work_dir}/assemblies/hifiasm-meta-hifi"
hifiasm_prefix="${hifiasm_dir}/asm"
if [[ -s "${hifiasm_prefix}.p_ctg.gfa" && -f "${hifiasm_dir}/.article31-complete" ]]; then
  echo "skip assembly hifiasm-meta-hifi: complete sentinel and primary GFA exist"
else
  mkdir -p "${hifiasm_dir}"
  run_timed "assemble-hifiasm-meta-hifi" \
    "${hifiasm}" \
      -t 32 \
      --force-rs \
      --tsne-seed 20260731 \
      -o "${hifiasm_prefix}" \
      "${hifi_reads}"
  test -s "${hifiasm_prefix}.p_ctg.gfa"
  touch "${hifiasm_dir}/.article31-complete"
fi

metamdbg_dir="${work_dir}/assemblies/metamdbg-hifi"
if [[ -s "${metamdbg_dir}/contigs.fasta.gz" && -f "${metamdbg_dir}/.article31-complete" ]]; then
  echo "skip assembly metamdbg-hifi: complete sentinel and contigs exist"
else
  run_timed "assemble-metamdbg-hifi" \
    "${metamdbg}" asm \
      --out-dir "${metamdbg_dir}" \
      --in-hifi "${hifi_reads}" \
      --threads 32 \
      --min-contig-length 50 \
      --min-contig-coverage 1
  test -s "${metamdbg_dir}/contigs.fasta.gz"
  touch "${metamdbg_dir}/.article31-complete"
fi

"${python}" "${prepare}" --work-dir "${work_dir}"

map_reads() {
  local branch="$1"
  local assembler="$2"
  local platform="$3"
  local preset="$4"
  local assembly="$5"
  local reads="$6"
  local expected_reads="$7"
  local expected_bases="$8"
  local paf="${work_dir}/mapping/${branch}.paf"
  local summary="${work_dir}/mapping/${branch}.json"
  if [[ -s "${summary}" && -f "${work_dir}/mapping/${branch}.complete" ]]; then
    local reference_bytes reference_sha256
    reference_bytes="$(stat -c '%s' "${assembly}")"
    reference_sha256="$(sha256sum "${assembly}" | cut -d ' ' -f 1)"
    if "${python}" -c \
      'import json,sys; x=json.load(open(sys.argv[1])); ok=(int(x.get("ReferenceThresholdBp",-1))==1000 and int(x.get("ReferenceBytes",-1))==int(sys.argv[2]) and x.get("ReferenceSHA256")==sys.argv[3] and x.get("IdentityComputation")=="minimap2-c-paf-matches-over-block" and int(x.get("BaseLevelCigarPAFRecords",-1))==int(x.get("PrimaryPAFRecords",-2))); raise SystemExit(0 if ok else 1)' \
      "${summary}" "${reference_bytes}" "${reference_sha256}"; then
      echo "skip mapping ${branch}: complete summary and reference identity match"
      return
    fi
    echo "rerun mapping ${branch}: stale summary/reference identity" >&2
  fi
  run_timed "map-${branch}" \
    "${minimap2}" \
      -c -x "${preset}" -t 32 -K 1g -I 8G \
      --secondary=no \
      -o "${paf}" \
      "${assembly}" "${reads}"
  "${python}" "${paf_summary}" \
    --paf "${paf}" \
    --output "${summary}" \
    --branch "${branch}" \
    --assembler "${assembler}" \
    --platform "${platform}" \
    --expected-reads "${expected_reads}" \
    --expected-bases "${expected_bases}" \
    --reference-threshold-bp 1000 \
    --reference-fasta "${assembly}"
  touch "${work_dir}/mapping/${branch}.complete"
}

map_reads flye-ont-r9 Flye "ONT R9" map-ont \
  "${work_dir}/normalized/flye-ont-r9.ge1000.fasta" \
  "${ont_reads}" 696944 3125920499
map_reads flye-hifi Flye "PacBio HiFi" map-hifi \
  "${work_dir}/normalized/flye-hifi.ge1000.fasta" \
  "${hifi_reads}" 524805 5400038744
map_reads hifiasm-meta-hifi hifiasm-meta "PacBio HiFi" map-hifi \
  "${work_dir}/normalized/hifiasm-meta-hifi.ge1000.fasta" \
  "${hifi_reads}" 524805 5400038744
map_reads metamdbg-hifi metaMDBG "PacBio HiFi" map-hifi \
  "${work_dir}/normalized/metamdbg-hifi.ge1000.fasta" \
  "${hifi_reads}" 524805 5400038744

map_junctions() {
  local branch="$1"
  local preset="$2"
  local reads="$3"
  local junctions="${work_dir}/normalized/${branch}.junctions.fasta"
  local paf="${work_dir}/junction-mapping/${branch}.paf"
  if [[ -f "${work_dir}/junction-mapping/${branch}.complete" ]]; then
    echo "skip junction mapping ${branch}: complete sentinel exists"
    return
  fi
  if [[ -s "${junctions}" ]]; then
    run_timed "junction-${branch}" \
      "${minimap2}" \
        -c -x "${preset}" -t 32 -K 1g -I 8G \
        --secondary=no \
        -o "${paf}" \
        "${junctions}" "${reads}"
  else
    : > "${paf}"
  fi
  touch "${work_dir}/junction-mapping/${branch}.complete"
}

map_junctions flye-ont-r9 map-ont "${ont_reads}"
map_junctions flye-hifi map-hifi "${hifi_reads}"
map_junctions hifiasm-meta-hifi map-hifi "${hifi_reads}"
map_junctions metamdbg-hifi map-hifi "${hifi_reads}"
"${python}" "${junction_summary}" --work-dir "${work_dir}"

{
  printf 'Tool\tVersion\tExecutable\n'
  printf 'Flye\t%s\t%s\n' "$("${flye}" --version 2>&1)" '${LONG_READ_ENV_PREFIX}/bin/flye'
  printf 'metaMDBG\t1.4\t%s\n' '${LONG_READ_ENV_PREFIX}/bin/metaMDBG'
  printf 'minimap2\t%s\t%s\n' "$("${minimap2}" --version 2>&1)" '${LONG_READ_ENV_PREFIX}/bin/minimap2'
  printf 'samtools\t%s\t%s\n' "$("${samtools}" --version | sed -n '1s/^samtools //p')" '${LONG_READ_ENV_PREFIX}/bin/samtools'
  printf 'SeqKit\t%s\t%s\n' "$("${seqkit}" version | sed 's/^seqkit v//')" '${LONG_READ_ENV_PREFIX}/bin/seqkit'
  printf 'Python\t%s\t%s\n' "$("${python}" --version 2>&1 | sed 's/^Python //')" '${LONG_READ_ENV_PREFIX}/bin/python'
  printf 'hifiasm-meta release\t0.3.5-r81\t%s\n' '${ARTICLE31_RAW_DIR}/tools/bin/hifiasm_meta'
  printf 'hifiasm-meta embedded\t%s\t%s\n' "$("${hifiasm}" --version 2>&1 | paste -sd ';' -)" '${ARTICLE31_RAW_DIR}/tools/bin/hifiasm_meta'
} > "${work_dir}/tool-versions.tsv"

"${python}" "${freezer}" \
  --project-root "${project_root}" \
  --env-prefix "${env_prefix}" \
  --raw-dir "${raw_dir}" \
  --work-dir "${work_dir}" \
  --frozen-dir "${frozen_dir}"

echo "Article 31 long-read assembly run and freeze completed."
