#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat >&2 <<'USAGE'
Usage:
  run_article13_read_qc.sh \
    --project-root DIR \
    --environment-prefix DIR \
    --raw-dir DIR \
    --frozen-dir DIR

The raw directory is excluded from Git. The frozen directory receives only
reports, logs, manifests, and checksums; no raw or cleaned FASTQ is copied.
USAGE
}

project_root=""
environment_prefix=""
raw_dir=""
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

if [[ -z "${project_root}" || -z "${environment_prefix}" || -z "${raw_dir}" || -z "${frozen_dir}" ]]; then
  usage
  exit 2
fi

project_root="$(cd "${project_root}" && pwd)"
environment_prefix="$(cd "${environment_prefix}" && pwd)"
mkdir -p "${raw_dir}"
raw_dir="$(cd "${raw_dir}" && pwd)"

if [[ -e "${frozen_dir}" ]] && [[ -n "$(find "${frozen_dir}" -mindepth 1 -maxdepth 1 -print -quit)" ]]; then
  echo "Refusing to overwrite non-empty frozen directory: ${frozen_dir}" >&2
  exit 1
fi
mkdir -p "${frozen_dir}"
frozen_dir="$(cd "${frozen_dir}" && pwd)"

fastqc="${environment_prefix}/bin/fastqc"
fastp="${environment_prefix}/bin/fastp"
multiqc="${environment_prefix}/bin/multiqc"
python="${environment_prefix}/bin/python"

for executable in "${fastqc}" "${fastp}" "${multiqc}" "${python}"; do
  if [[ ! -x "${executable}" ]]; then
    echo "Missing executable: ${executable}" >&2
    exit 1
  fi
done

export LC_ALL=C
export TZ=UTC
export MPLCONFIGDIR="${raw_dir}/.matplotlib"
export XDG_CACHE_HOME="${raw_dir}/.cache"
mkdir -p "${MPLCONFIGDIR}" "${XDG_CACHE_HOME}"

source_manifest="${project_root}/data/small/13-source-manifest.tsv"
subset_summary="${raw_dir}/subset-summary.json"
raw_r1="${raw_dir}/ERR9765746_prefix100k_R1.fastq.gz"
raw_r2="${raw_dir}/ERR9765746_prefix100k_R2.fastq.gz"
clean_r1="${raw_dir}/ERR9765746_clean_R1.fastq.gz"
clean_r2="${raw_dir}/ERR9765746_clean_R2.fastq.gz"

if [[ ! -s "${raw_r1}" || ! -s "${raw_r2}" || ! -s "${subset_summary}" ]]; then
  "${python}" "${project_root}/scripts/build_article13_fastq_subset.py" \
    --manifest "${source_manifest}" \
    --output-dir "${raw_dir}"
fi

mkdir -p \
  "${frozen_dir}/raw_fastqc" \
  "${frozen_dir}/fastp" \
  "${frozen_dir}/clean_fastqc" \
  "${frozen_dir}/multiqc" \
  "${frozen_dir}/logs"

cp "${subset_summary}" "${frozen_dir}/subset-summary.json"
cp "$0" "${frozen_dir}/commands.sh"

{
  printf 'Tool\tVersion\tExecutable\n'
  printf 'FastQC\t%s\t%s\n' \
    "$("${fastqc}" --version | sed 's/^FastQC v//')" \
    '${READ_QC_ENV_PREFIX}/bin/fastqc'
  printf 'fastp\t%s\t%s\n' \
    "$("${fastp}" --version | sed 's/^fastp //')" \
    '${READ_QC_ENV_PREFIX}/bin/fastp'
  printf 'MultiQC\t%s\t%s\n' \
    "$("${multiqc}" --version | sed 's/^multiqc, version //')" \
    '${READ_QC_ENV_PREFIX}/bin/multiqc'
  printf 'Python\t%s\t%s\n' \
    "$("${python}" --version | sed 's/^Python //')" \
    '${READ_QC_ENV_PREFIX}/bin/python'
} > "${frozen_dir}/tool-versions.tsv"

/usr/bin/time -v \
  -o "${frozen_dir}/logs/fastqc-before.resources.txt" \
  "${fastqc}" \
  --noextract \
  --threads 2 \
  --memory 512 \
  --svg \
  --outdir "${frozen_dir}/raw_fastqc" \
  "${raw_r1}" "${raw_r2}" \
  > "${frozen_dir}/logs/fastqc-before.log" 2>&1

/usr/bin/time -v \
  -o "${frozen_dir}/logs/fastp.resources.txt" \
  "${fastp}" \
  --in1 "${raw_r1}" \
  --in2 "${raw_r2}" \
  --out1 "${clean_r1}" \
  --out2 "${clean_r2}" \
  --json "${frozen_dir}/fastp/ERR9765746_fastp.json" \
  --html "${frozen_dir}/fastp/ERR9765746_fastp.html" \
  --report_title "ERR9765746 first 100,000 pairs" \
  --thread 4 \
  --compression 6 \
  --detect_adapter_for_pe \
  --qualified_quality_phred 20 \
  --unqualified_percent_limit 40 \
  --n_base_limit 5 \
  --length_required 50 \
  --disable_trim_poly_g \
  --overrepresentation_analysis \
  --overrepresentation_sampling 20 \
  --dont_overwrite \
  > "${frozen_dir}/logs/fastp.log" 2>&1

/usr/bin/time -v \
  -o "${frozen_dir}/logs/fastqc-after.resources.txt" \
  "${fastqc}" \
  --noextract \
  --threads 2 \
  --memory 512 \
  --svg \
  --outdir "${frozen_dir}/clean_fastqc" \
  "${clean_r1}" "${clean_r2}" \
  > "${frozen_dir}/logs/fastqc-after.log" 2>&1

/usr/bin/time -v \
  -o "${frozen_dir}/logs/multiqc.resources.txt" \
  "${multiqc}" \
  --force \
  --no-version-check \
  --module fastqc \
  --module fastp \
  --require-logs \
  --dirs \
  --dirs-depth 1 \
  --fullnames \
  --data-dir \
  --data-format json \
  --title "ERR9765746 read QC: raw versus clean" \
  --filename "13-multiqc-report.html" \
  --outdir "${frozen_dir}/multiqc" \
  "${frozen_dir}/raw_fastqc" \
  "${frozen_dir}/clean_fastqc" \
  "${frozen_dir}/fastp" \
  > "${frozen_dir}/logs/multiqc.log" 2>&1

"${python}" "${project_root}/scripts/validate_article13_read_qc.py" \
  --project-root "${project_root}" \
  --environment-prefix "${environment_prefix}" \
  --frozen-dir "${frozen_dir}" \
  --raw-dir "${raw_dir}" \
  --initialize-frozen

echo "Article 13 one-time QC completed: ${frozen_dir}"
