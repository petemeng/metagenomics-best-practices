#!/usr/bin/env bash
set -euo pipefail

usage() {
  printf '%s\n' \
    'Usage: download_article35_gene_abundance_sources.sh PROJECT_ROOT [RAW_DIR] [READ_QC_ENV_PREFIX]' \
    '' \
    'Downloads the two checksum-locked ENA runs through the Article 30 downloader,' \
    'creates the exact paired fastp inputs when absent, and verifies the Article 35' \
    'clean-read identities. The 36.3-GB UniRef90 database is managed separately by' \
    'scripts/bootstrap_article19_humann_databases.sh.' >&2
}

if [[ "$#" -lt 1 || "$#" -gt 3 ]]; then
  usage
  exit 2
fi

project_root="$(cd "$1" && pwd)"
raw_dir="${2:-${project_root}/data/raw/article30}"
read_qc_env="${3:-${HOME}/miniconda3/envs/metagenome-read-qc-2026.07}"

bash "${project_root}/scripts/download_article30_assembly_reads.sh" \
  "${project_root}" "${raw_dir}"

fastp="${read_qc_env}/bin/fastp"
if [[ ! -x "${fastp}" ]]; then
  printf 'Missing executable: %s\n' "${fastp}" >&2
  printf 'Create the locked environment from env/read-qc-linux-64.lock first.\n' >&2
  exit 1
fi

mkdir -p "${raw_dir}/clean" "${raw_dir}/work/fastp"

clean_one() {
  local sample="$1"
  local run="$2"
  local out1="${raw_dir}/clean/${run}_clean_R1.fastq.gz"
  local out2="${raw_dir}/clean/${run}_clean_R2.fastq.gz"
  if [[ -s "${out1}" && -s "${out2}" ]]; then
    printf 'clean reads already exist: %s\n' "${run}"
    return
  fi
  if [[ -e "${out1}" || -e "${out2}" ]]; then
    printf 'Refusing incomplete clean-read pair for %s\n' "${run}" >&2
    exit 1
  fi
  "${fastp}" \
    --in1 "${raw_dir}/selected/${run}_selected2m_R1.fastq.gz" \
    --in2 "${raw_dir}/selected/${run}_selected2m_R2.fastq.gz" \
    --out1 "${out1}" --out2 "${out2}" \
    --json "${raw_dir}/work/fastp/${sample}.json" \
    --html "${raw_dir}/work/fastp/${sample}.html" \
    --thread 16 --compression 6 --detect_adapter_for_pe \
    --qualified_quality_phred 20 --unqualified_percent_limit 40 \
    --n_base_limit 5 --length_required 50 --disable_trim_poly_g \
    --overrepresentation_analysis --overrepresentation_sampling 20 \
    --dont_overwrite
}

clean_one MOCK1 ERR9765746
clean_one MOCK2 ERR9765747

python3 - "${project_root}" "${raw_dir}" <<'PY'
import hashlib
import json
import sys
from pathlib import Path

root = Path(sys.argv[1])
raw = Path(sys.argv[2])
summary = json.loads(
    (root / "data/small/30-short-read-assembly-frozen/run-summary.json").read_text()
)
specs = (("MOCK1", "ERR9765746"), ("MOCK2", "ERR9765747"))
for sample, run in specs:
    for mate in ("R1", "R2"):
        path = raw / "clean" / f"{run}_clean_{mate}.fastq.gz"
        hasher = hashlib.sha256()
        with path.open("rb") as handle:
            for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
                hasher.update(block)
        digest = hasher.hexdigest()
        expected = summary["clean_fastq_audit"][sample][mate]
        if path.stat().st_size != expected["CompressedBytes"]:
            raise SystemExit(f"byte-count mismatch: {path}")
        if digest != expected["CompressedSHA256"]:
            raise SystemExit(f"SHA-256 mismatch: {path}")
        print(f"verified\t{sample}\t{mate}\t{path.stat().st_size}\t{digest}")
PY

catalog="${project_root}/data/small/34-nonredundant-gene-catalog-frozen"
for relative in \
  catalog/megahit-mix-primary.fna.gz \
  catalog/megahit-mix-primary.faa.gz \
  primary-catalog-representatives.tsv.gz; do
  test -s "${catalog}/${relative}"
done

printf '%s\n' \
  'Article 35 read and catalog sources are ready.' \
  'For de novo UniRef90 annotation, download/extract HUMAnN UniRef90 v201901b with:' \
  "  bash ${project_root}/scripts/bootstrap_article19_humann_databases.sh --project-root ${project_root} --cache-root ${project_root}/db/humann-cache download" \
  "  bash ${project_root}/scripts/bootstrap_article19_humann_databases.sh --project-root ${project_root} --cache-root ${project_root}/db/humann-cache extract-uniref90"
