#!/usr/bin/env bash
#SBATCH --job-name=mg-compute-smoke
#SBATCH --array=1-3%3
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=00:10:00
#SBATCH --output=logs/10-smoke-%A_%a.out
#SBATCH --error=logs/10-smoke-%A_%a.err
#SBATCH --signal=B:USR1@120

set -euo pipefail

: "${PROJECT_ROOT:?Set PROJECT_ROOT to the repository root}"
RUN_MODE="${RUN_MODE:-container}"
TASK_TABLE="${TASK_TABLE:-${PROJECT_ROOT}/data/small/10-job-array.tsv}"
METRICS="${METRICS:-${PROJECT_ROOT}/data/small/08-read-prefix-metrics.tsv}"
OUT_DIR="${OUT_DIR:-${PROJECT_ROOT}/results/10-computing-hpc-cloud/array-output}"
TASK_ID="${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is required}"
SEED=20260719

if [[ ! -r "${TASK_TABLE}" || ! -r "${METRICS}" ]]; then
  echo "Required task table or metrics input is not readable" >&2
  exit 2
fi

task_line="$(
  awk -F '\t' -v task_id="${TASK_ID}" \
    'NR > 1 && $1 == task_id {print; found = 1} END {if (!found) exit 3}' \
    "${TASK_TABLE}"
)"
IFS=$'\t' read -r task_id platform platform_label run_accession mate \
  expected_rows source_bytes input_scope <<< "${task_line}"

safe_platform="$(printf '%s' "${platform}" | tr '[:upper:]' '[:lower:]')"
final_json="${OUT_DIR}/$(printf '%02d' "${task_id}")-${safe_platform}.json"
done_file="${final_json}.done"
mkdir -p "${OUT_DIR}"

if [[ -s "${final_json}" && -s "${done_file}" ]]; then
  expected_sha="$(cut -d ' ' -f 1 "${done_file}")"
  observed_sha="$(sha256sum "${final_json}" | cut -d ' ' -f 1)"
  if [[ "${expected_sha}" == "${observed_sha}" ]]; then
    echo "ACTION=SKIPPED TASK_ID=${task_id} PLATFORM=${platform} CHECKSUM=${observed_sha}"
    exit 0
  fi
  echo "Existing output checksum does not match sentinel: ${final_json}" >&2
  exit 4
fi

if [[ -n "${SLURM_TMPDIR:-}" ]]; then
  scratch_root="${SLURM_TMPDIR}"
  mkdir -p "${scratch_root}"
  task_scratch="${scratch_root}/article10-${SLURM_JOB_ID:-local}-${task_id}"
  mkdir -p "${task_scratch}"
else
  task_scratch="$(mktemp -d "${TMPDIR:-/tmp}/article10-${task_id}.XXXXXX")"
fi

cleanup() {
  rm -rf "${task_scratch}"
}
on_signal() {
  echo "ACTION=INTERRUPTED TASK_ID=${task_id}; no final sentinel written" >&2
  exit 99
}
trap cleanup EXIT
trap on_signal USR1 TERM INT

cp "${METRICS}" "${task_scratch}/08-read-prefix-metrics.tsv"
scratch_json="${task_scratch}/task-summary.json"

case "${RUN_MODE}" in
  container)
    : "${APPTAINER_IMAGE:?Set APPTAINER_IMAGE for container mode}"
    : "${APPTAINER_IMAGE_SHA256:?Set APPTAINER_IMAGE_SHA256 for container mode}"
    printf '%s  %s\n' "${APPTAINER_IMAGE_SHA256}" "${APPTAINER_IMAGE}" |
      sha256sum -c -
    command -v apptainer >/dev/null 2>&1 || {
      echo "Apptainer is required in container mode" >&2
      exit 5
    }
    apptainer exec \
      --cleanenv \
      --bind "${PROJECT_ROOT}:/work:ro,${task_scratch}:/scratch:rw" \
      --pwd /work \
      "${APPTAINER_IMAGE}" \
      python /work/scripts/article10_task.py \
        --metrics /scratch/08-read-prefix-metrics.tsv \
        --platform "${platform}" \
        --run-accession "${run_accession}" \
        --expected-rows "${expected_rows}" \
        --seed "${SEED}" \
        --output /scratch/task-summary.json
    ;;
  native-smoke)
    python3 "${PROJECT_ROOT}/scripts/article10_task.py" \
      --metrics "${task_scratch}/08-read-prefix-metrics.tsv" \
      --platform "${platform}" \
      --run-accession "${run_accession}" \
      --expected-rows "${expected_rows}" \
      --seed "${SEED}" \
      --output "${scratch_json}"
    ;;
  *)
    echo "RUN_MODE must be container or native-smoke" >&2
    exit 6
    ;;
esac

python3 - "${scratch_json}" "${platform}" "${run_accession}" "${expected_rows}" <<'PY'
import json
import sys
from pathlib import Path

path = Path(sys.argv[1])
platform = sys.argv[2]
run_accession = sys.argv[3]
expected_rows = int(sys.argv[4])
payload = json.loads(path.read_text(encoding="utf-8"))
assert payload["status"] == "passed"
assert payload["platform"] == platform
assert payload["run_accession"] == run_accession
assert payload["rows"] == expected_rows
PY

tmp_final="${final_json}.tmp.${SLURM_JOB_ID:-local}.${task_id}"
cp "${scratch_json}" "${tmp_final}"
mv "${tmp_final}" "${final_json}"
final_sha="$(sha256sum "${final_json}" | cut -d ' ' -f 1)"
printf '%s  %s\n' "${final_sha}" "$(basename "${final_json}")" > "${done_file}"

echo "ACTION=COMPLETED TASK_ID=${task_id} PLATFORM=${platform} ROWS=${expected_rows} CHECKSUM=${final_sha}"
