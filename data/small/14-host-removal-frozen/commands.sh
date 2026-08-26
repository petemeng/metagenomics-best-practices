#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat >&2 <<'USAGE'
Usage:
  run_article14_host_removal.sh \
    --project-root DIR \
    --environment-prefix DIR \
    --raw-dir DIR \
    --index-archive FILE \
    --index-dir DIR \
    --frozen-dir DIR

The input and derived FASTQ files remain in the Git-ignored raw directory.
The large Hostile archive and extracted Bowtie2 index remain in the Git-ignored
index directory. Only aggregate reports, normalized logs, and checksums are
written to the frozen directory.
USAGE
}

project_root=""
environment_prefix=""
raw_dir=""
index_archive=""
index_dir=""
frozen_dir=""

while [[ "$#" -gt 0 ]]; do
  case "$1" in
    --project-root)
      project_root="$2"
      shift 2
      ;;
    --environment-prefix)
      environment_prefix="$2"
      shift 2
      ;;
    --raw-dir)
      raw_dir="$2"
      shift 2
      ;;
    --index-archive)
      index_archive="$2"
      shift 2
      ;;
    --index-dir)
      index_dir="$2"
      shift 2
      ;;
    --frozen-dir)
      frozen_dir="$2"
      shift 2
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "Unknown argument: $1" >&2
      usage
      exit 2
      ;;
  esac
done

if [[ -z "${project_root}" ||
      -z "${environment_prefix}" ||
      -z "${raw_dir}" ||
      -z "${index_archive}" ||
      -z "${index_dir}" ||
      -z "${frozen_dir}" ]]; then
  usage
  exit 2
fi

project_root="$(cd "${project_root}" && pwd)"
environment_prefix="$(cd "${environment_prefix}" && pwd)"
mkdir -p "${raw_dir}" "${index_dir}"
raw_dir="$(cd "${raw_dir}" && pwd)"
index_dir="$(cd "${index_dir}" && pwd)"
index_archive="$(realpath "${index_archive}")"

if [[ -e "${frozen_dir}" ]] &&
   [[ -n "$(find "${frozen_dir}" -mindepth 1 -maxdepth 1 -print -quit)" ]]; then
  echo "Refusing to overwrite non-empty frozen directory: ${frozen_dir}" >&2
  exit 1
fi
mkdir -p "${frozen_dir}"
frozen_dir="$(cd "${frozen_dir}" && pwd)"

work_dir="${raw_dir}/work"
if [[ -e "${work_dir}" ]] &&
   [[ -n "$(find "${work_dir}" -mindepth 1 -maxdepth 1 -print -quit)" ]]; then
  echo "Refusing to overwrite non-empty Article 14 work directory: ${work_dir}" >&2
  exit 1
fi
mkdir -p "${work_dir}"

hostile="${environment_prefix}/bin/hostile"
bowtie2="${environment_prefix}/bin/bowtie2"
samtools="${environment_prefix}/bin/samtools"
fastp="${environment_prefix}/bin/fastp"
seqkit="${environment_prefix}/bin/seqkit"
python="${environment_prefix}/bin/python"

for executable in \
  "${hostile}" "${bowtie2}" "${samtools}" "${fastp}" "${seqkit}" "${python}"; do
  if [[ ! -x "${executable}" ]]; then
    echo "Missing executable: ${executable}" >&2
    exit 1
  fi
done

export PATH="${environment_prefix}/bin:${PATH}"
export LC_ALL=C
export TZ=UTC
export MPLCONFIGDIR="${work_dir}/.matplotlib"
export XDG_CACHE_HOME="${work_dir}/.cache"
export HOSTILE_CACHE_DIR="${index_dir}"
mkdir -p "${MPLCONFIGDIR}" "${XDG_CACHE_HOME}"

expected_archive_sha256="5b584f5c28abeec5dba78bd37b53fa476dd42af57051d2fb7d2f2098e3a2df13"
expected_archive_bytes="3934284979"

if [[ ! -s "${index_archive}" ]]; then
  echo "Hostile index archive is missing or empty: ${index_archive}" >&2
  exit 1
fi
observed_archive_bytes="$(stat -c '%s' "${index_archive}")"
observed_archive_sha256="$(sha256sum "${index_archive}" | awk '{print $1}')"
if [[ "${observed_archive_bytes}" != "${expected_archive_bytes}" ]]; then
  echo "Hostile index size mismatch: ${observed_archive_bytes}" >&2
  exit 1
fi
if [[ "${observed_archive_sha256}" != "${expected_archive_sha256}" ]]; then
  echo "Hostile index SHA-256 mismatch: ${observed_archive_sha256}" >&2
  exit 1
fi

if [[ ! -s "${index_dir}/human-t2t-hla.1.bt2" ]]; then
  tar --extract \
    --file "${index_archive}" \
    --directory "${index_dir}" \
    --no-same-owner \
    --no-same-permissions
fi

for suffix in 1 2 3 4 rev.1 rev.2; do
  if [[ ! -s "${index_dir}/human-t2t-hla.${suffix}.bt2" ]]; then
    echo "Extracted Bowtie2 index is incomplete: ${suffix}" >&2
    exit 1
  fi
done

source_manifest="${project_root}/data/small/14-source-manifest.tsv"
controls_summary="${raw_dir}/controls-summary.json"
human_r1="${raw_dir}/ERR194147_prefix20k_R1.fastq.gz"
human_r2="${raw_dir}/ERR194147_prefix20k_R2.fastq.gz"
mock_r1="${raw_dir}/ERR9765746_prefix20k_R1.fastq.gz"
mock_r2="${raw_dir}/ERR9765746_prefix20k_R2.fastq.gz"

if [[ ! -s "${human_r1}" ||
      ! -s "${human_r2}" ||
      ! -s "${mock_r1}" ||
      ! -s "${mock_r2}" ||
      ! -s "${controls_summary}" ]]; then
  "${python}" "${project_root}/scripts/build_article14_fastq_controls.py" \
    --manifest "${source_manifest}" \
    --output-dir "${raw_dir}"
fi

mkdir -p \
  "${frozen_dir}/fastp" \
  "${frozen_dir}/hostile" \
  "${frozen_dir}/logs" \
  "${work_dir}/hostile/human" \
  "${work_dir}/hostile/mock" \
  "${work_dir}/fastp"

cp "${controls_summary}" "${frozen_dir}/controls-summary.json"
cp "$0" "${frozen_dir}/commands.sh"

{
  printf 'Tool\tVersion\tExecutable\n'
  printf 'Hostile\t%s\t%s\n' \
    "$("${hostile}" --version | sed 's/^hostile //')" \
    '${HOST_REMOVAL_ENV_PREFIX}/bin/hostile'
  printf 'Bowtie2\t%s\t%s\n' \
    "$("${bowtie2}" --version | sed -n '1s/.*version //p')" \
    '${HOST_REMOVAL_ENV_PREFIX}/bin/bowtie2'
  printf 'Samtools\t%s\t%s\n' \
    "$("${samtools}" --version | sed -n '1s/^samtools //p')" \
    '${HOST_REMOVAL_ENV_PREFIX}/bin/samtools'
  printf 'fastp\t%s\t%s\n' \
    "$("${fastp}" --version | sed 's/^fastp //')" \
    '${HOST_REMOVAL_ENV_PREFIX}/bin/fastp'
  printf 'SeqKit\t%s\t%s\n' \
    "$("${seqkit}" version | sed 's/^seqkit v//')" \
    '${HOST_REMOVAL_ENV_PREFIX}/bin/seqkit'
  printf 'Python\t%s\t%s\n' \
    "$("${python}" --version | sed 's/^Python //')" \
    '${HOST_REMOVAL_ENV_PREFIX}/bin/python'
} > "${frozen_dir}/tool-versions.tsv"

{
  printf 'IndexName\tAssetName\tExpectedBytes\tObservedBytes\tExpectedSHA256\tObservedSHA256\tStatus\n'
  printf 'human-t2t-hla\thuman-t2t-hla.tar\t%s\t%s\t%s\t%s\tVERIFIED\n' \
    "${expected_archive_bytes}" \
    "${observed_archive_bytes}" \
    "${expected_archive_sha256}" \
    "${observed_archive_sha256}"
} > "${frozen_dir}/index-archive-audit.tsv"

run_hostile() {
  local control="$1"
  local read1="$2"
  local read2="$3"
  local output_dir="${work_dir}/hostile/${control}"

  /usr/bin/time -v \
    -o "${frozen_dir}/logs/hostile-${control}.resources.txt" \
    "${hostile}" clean \
    --fastq1 "${read1}" \
    --fastq2 "${read2}" \
    --aligner bowtie2 \
    --index human-t2t-hla \
    --rename \
    --reorder \
    --threads 8 \
    --airplane \
    --output "${output_dir}" \
    --force \
    > "${frozen_dir}/hostile/${control}-hostile.json" \
    2> "${frozen_dir}/logs/hostile-${control}.log"
}

run_fastp_branch() {
  local control="$1"
  local read1="$2"
  local read2="$3"
  local branch="$4"
  local extra_flag="${5:-}"
  local extra_value="${6:-}"
  local output1="${work_dir}/fastp/${control}-${branch}_R1.fastq.gz"
  local output2="${work_dir}/fastp/${control}-${branch}_R2.fastq.gz"
  local json="${frozen_dir}/fastp/${control}-${branch}.json"
  local html="${work_dir}/fastp/${control}-${branch}.html"
  local resource="${frozen_dir}/logs/fastp-${control}-${branch}.resources.txt"
  local log="${frozen_dir}/logs/fastp-${control}-${branch}.log"
  local optional_args=()

  if [[ -n "${extra_flag}" ]]; then
    optional_args+=("${extra_flag}")
  fi
  if [[ -n "${extra_value}" ]]; then
    optional_args+=("${extra_value}")
  fi

  /usr/bin/time -v \
    -o "${resource}" \
    "${fastp}" \
    --in1 "${read1}" \
    --in2 "${read2}" \
    --out1 "${output1}" \
    --out2 "${output2}" \
    --json "${json}" \
    --html "${html}" \
    --report_title "${control} ${branch} sensitivity branch" \
    --thread 4 \
    --compression 6 \
    --disable_adapter_trimming \
    --disable_quality_filtering \
    --disable_length_filtering \
    --disable_trim_poly_g \
    --dont_overwrite \
    "${optional_args[@]}" \
    > "${log}" 2>&1
}

run_hostile human "${human_r1}" "${human_r2}"
run_hostile mock "${mock_r1}" "${mock_r2}"

for control in human mock; do
  if [[ "${control}" == "human" ]]; then
    control_r1="${human_r1}"
    control_r2="${human_r2}"
  else
    control_r1="${mock_r1}"
    control_r2="${mock_r2}"
  fi
  run_fastp_branch "${control}" "${control_r1}" "${control_r2}" baseline
  for threshold in 20 30 40; do
    run_fastp_branch \
      "${control}" \
      "${control_r1}" \
      "${control_r2}" \
      "complexity-${threshold}" \
      --low_complexity_filter \
      "--complexity_threshold=${threshold}"
  done
  run_fastp_branch \
    "${control}" \
    "${control_r1}" \
    "${control_r2}" \
    dedup \
    --dedup \
    "--dup_calc_accuracy=1"
done

"${python}" "${project_root}/scripts/validate_article14_host_removal.py" \
  --project-root "${project_root}" \
  --environment-prefix "${environment_prefix}" \
  --frozen-dir "${frozen_dir}" \
  --raw-dir "${raw_dir}" \
  --work-dir "${work_dir}" \
  --index-archive "${index_archive}" \
  --index-dir "${index_dir}" \
  --initialize-frozen

echo "Article 14 one-time run completed: ${frozen_dir}"
