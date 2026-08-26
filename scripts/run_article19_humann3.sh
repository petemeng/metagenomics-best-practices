#!/usr/bin/env bash
set -euo pipefail

usage() {
  printf '%s\n' \
    'Usage:' \
    '  run_article19_humann3.sh --project-root DIR --environment-prefix DIR' \
    '    --cache-root DIR --raw-dir DIR --frozen-dir DIR' \
    '    [--profile-only] [--resume]' >&2
}

project_root=""
environment_prefix=""
cache_root=""
raw_dir=""
frozen_dir=""
resume="no"
profile_only="no"
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
    --cache-root)
      cache_root="$2"
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
    --resume)
      resume="yes"
      shift
      ;;
    --profile-only)
      profile_only="yes"
      shift
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      printf 'Unknown argument: %s\n' "$1" >&2
      usage
      exit 2
      ;;
  esac
done

if [[ -z "${project_root}" || -z "${environment_prefix}" ||
      -z "${cache_root}" || -z "${raw_dir}" || -z "${frozen_dir}" ]]; then
  usage
  exit 2
fi

project_root="$(cd "${project_root}" && pwd)"
environment_prefix="$(cd "${environment_prefix}" && pwd)"
cache_root="$(cd "${cache_root}" && pwd)"
mkdir -p "${raw_dir}"
raw_dir="$(cd "${raw_dir}" && pwd)"

if [[ -e "${frozen_dir}" ]] &&
   [[ -n "$(find "${frozen_dir}" -mindepth 1 -maxdepth 1 -print -quit)" ]]; then
  printf 'Refusing to overwrite non-empty frozen directory: %s\n' "${frozen_dir}" >&2
  exit 1
fi

work_dir="${raw_dir}/work"
if [[ "${resume}" == "no" && -e "${work_dir}" ]] &&
   [[ -n "$(find "${work_dir}" -mindepth 1 -maxdepth 1 -print -quit)" ]]; then
  printf 'Refusing to overwrite non-empty Article 19 work directory: %s\n' "${work_dir}" >&2
  exit 1
fi
mkdir -p "${work_dir}"/{metaphlan,humann,logs,tmp}

export LC_ALL=C
export TZ=UTC
export PYTHONHASHSEED=0
export PYTHONNOUSERSITE=1
export PATH="${environment_prefix}/bin:/usr/bin:/bin"
export MPLCONFIGDIR="${work_dir}/.matplotlib"
export XDG_CACHE_HOME="${work_dir}/.cache"
mkdir -p "${MPLCONFIGDIR}" "${XDG_CACHE_HOME}"

humann="${environment_prefix}/bin/humann"
metaphlan="${environment_prefix}/bin/metaphlan"
bowtie2="${environment_prefix}/bin/bowtie2"
diamond="${environment_prefix}/bin/diamond"
renorm="${environment_prefix}/bin/humann_renorm_table"
regroup="${environment_prefix}/bin/humann_regroup_table"
python="${environment_prefix}/bin/python"
metaphlan_compat_dir="${project_root}/scripts/humann39-compat-bin"
for executable in \
  "${humann}" "${metaphlan}" "${bowtie2}" "${diamond}" \
  "${renorm}" "${regroup}" "${python}"; do
  if [[ ! -x "${executable}" ]]; then
    printf 'Missing executable: %s\n' "${executable}" >&2
    exit 1
  fi
done
if [[ ! -x "${metaphlan_compat_dir}/metaphlan" ]]; then
  printf 'Missing Article 19 MetaPhlAn compatibility wrapper\n' >&2
  exit 1
fi
export ARTICLE19_METAPHLAN_REAL="${metaphlan}"

clean_r1="${project_root}/data/raw/article13/ERR9765746_clean_R1.fastq.gz"
clean_r2="${project_root}/data/raw/article13/ERR9765746_clean_R2.fastq.gz"
for input_file in "${clean_r1}" "${clean_r2}"; do
  if [[ ! -s "${input_file}" ]]; then
    printf 'Missing Article 13 clean FASTQ: %s\n' "${input_file}" >&2
    exit 1
  fi
done

manifest="${project_root}/data/small/19-database-manifest.tsv"
manifest_row() {
  local database_id="$1"
  awk -F '\t' -v id="${database_id}" 'NR > 1 && $1 == id {print; found=1; exit} END {if (!found) exit 3}' "${manifest}"
}

verify_archive() {
  local database_id="$1"
  local archive="$2"
  local row tool tool_version release release_date url algorithm checksum bytes
  local installed gate status notes observed_bytes observed_checksum
  row="$(manifest_row "${database_id}")"
  IFS=$'\t' read -r \
    _ tool tool_version release release_date url algorithm checksum bytes \
    installed gate status notes <<< "${row}"
  if [[ "${gate}" != "enabled" ]]; then
    printf '%s remains fail-closed: %s\n' "${database_id}" "${gate}" >&2
    exit 1
  fi
  if [[ ! -s "${archive}" ]]; then
    printf 'Missing database archive: %s\n' "${archive}" >&2
    exit 1
  fi
  observed_bytes="$(stat -c '%s' "${archive}")"
  if [[ "${observed_bytes}" != "${bytes}" ]]; then
    printf 'Archive byte mismatch for %s\n' "${database_id}" >&2
    exit 1
  fi
  case "${algorithm}" in
    md5)
      observed_checksum="$(md5sum "${archive}" | awk '{print $1}')"
      ;;
    sha256)
      observed_checksum="$(sha256sum "${archive}" | awk '{print $1}')"
      ;;
    *)
      printf 'Unsupported locked checksum for %s: %s\n' "${database_id}" "${algorithm}" >&2
      exit 1
      ;;
  esac
  if [[ "${observed_checksum}" != "${checksum}" ]]; then
    printf 'Archive checksum mismatch for %s\n' "${database_id}" >&2
    exit 1
  fi
}

mpa_metadata_archive="${cache_root}/archives/metaphlan-vjun23-metadata/mpa_vJun23_CHOCOPhlAnSGB_202403.tar"
mpa_bowtie_archive="${cache_root}/archives/metaphlan-vjun23-bowtie2/mpa_vJun23_CHOCOPhlAnSGB_202403_bt2.tar"
choco_archive="${cache_root}/archives/humann-chocophlan-full/full_chocophlan.v201901_v31.tar.gz"
uniref_archive="${cache_root}/archives/humann-uniref90-full/uniref90_annotated_v201901b_full.tar.gz"
verify_archive metaphlan-vjun23-metadata "${mpa_metadata_archive}"
verify_archive metaphlan-vjun23-bowtie2 "${mpa_bowtie_archive}"
if [[ "${profile_only}" == "no" ]]; then
  verify_archive humann-chocophlan-full "${choco_archive}"
  verify_archive humann-uniref90-full "${uniref_archive}"
fi

mpa_dir="${cache_root}/installed/metaphlan-vjun23"
choco_root="${cache_root}/installed/humann-chocophlan-v201901-v31"
uniref_root="${cache_root}/installed/humann-uniref90-v201901b"
index_name="mpa_vJun23_CHOCOPhlAnSGB_202403"

if [[ ! -s "${mpa_dir}/${index_name}.pkl" ]]; then
  printf 'Missing MetaPhlAn vJun23 metadata pickle\n' >&2
  exit 1
fi
mpa_index_count="$(find "${mpa_dir}" -maxdepth 1 -type f -name "${index_name}*.bt2*" -size +0c | wc -l)"
if [[ "${mpa_index_count}" != "6" ]]; then
  printf 'Expected six MetaPhlAn vJun23 Bowtie2 index files; observed %s\n' "${mpa_index_count}" >&2
  exit 1
fi

combined_fastq="${work_dir}/ERR9765746_clean_all_reads.fastq.gz"
if [[ ! -s "${combined_fastq}" ]]; then
  gzip -cd "${clean_r1}" "${clean_r2}" | gzip -n > "${combined_fastq}"
fi

"${python}" - "${combined_fastq}" <<'PY'
import gzip
import sys
from pathlib import Path

path = Path(sys.argv[1])
records = 0
ids = set()
with gzip.open(path, "rt", encoding="utf-8") as handle:
    while True:
        header = handle.readline()
        if not header:
            break
        sequence = handle.readline()
        plus = handle.readline()
        quality = handle.readline()
        if not sequence or not plus or not quality or not header.startswith("@") or not plus.startswith("+"):
            raise SystemExit("Malformed concatenated FASTQ")
        records += 1
        ids.add(header[1:].strip().replace(" ", ""))
if records != 199982 or len(ids) != records:
    raise SystemExit(f"Concatenated FASTQ invariant failed: records={records}, unique_ids={len(ids)}")
print(f"PASS concatenated FASTQ: records={records}, unique_ids={len(ids)}")
PY

vjun_profile="${work_dir}/metaphlan/ERR9765746-vJun23-profile.tsv"
vjun_mapout="${work_dir}/metaphlan/ERR9765746-vJun23.mapout.bz2"
if [[ ! -s "${vjun_profile}" ]]; then
  /usr/bin/time -v \
    -o "${work_dir}/logs/metaphlan-vjun23.resources.txt" \
    "${metaphlan}" \
      "${clean_r1},${clean_r2}" \
      --input_type fastq \
      --db_dir "${mpa_dir}" \
      --index "${index_name}" \
      --offline \
      --bowtie2_exe "${bowtie2}" \
      --mapout "${vjun_mapout}" \
      --tmp_dir "${work_dir}/tmp" \
      --nproc 8 \
      --read_min_len 70 \
      --perc_nonzero 0.33 \
      --stat tavg_g \
      --stat_q 0.2 \
      -t rel_ab \
      --sample_id ERR9765746_MOCK1_vJun23 \
      --output_file "${vjun_profile}" \
      > "${work_dir}/logs/metaphlan-vjun23.log" 2>&1
fi
if ! grep -q 'vJun23' "${vjun_profile}"; then
  printf 'Compatible profile does not declare vJun23\n' >&2
  exit 1
fi
if [[ "${profile_only}" == "yes" ]]; then
  printf 'Article 19 compatible MetaPhlAn profile complete: %s\n' "${vjun_profile}"
  exit 0
fi

mapfile -t choco_dirs < <(
  find "${choco_root}" -type f -name 'g__*.v201901_v31.ffn.gz' -printf '%h\n' | sort -u
)
mapfile -t uniref_dirs < <(
  find "${uniref_root}" -type f -name '*uniref90*201901b*.dmnd' -printf '%h\n' | sort -u
)
if [[ "${#choco_dirs[@]}" != "1" || "${#uniref_dirs[@]}" != "1" ]]; then
  printf 'Could not resolve exactly one ChocoPhlAn and one UniRef90 database directory\n' >&2
  printf 'Choco directories: %s; UniRef directories: %s\n' \
    "${#choco_dirs[@]}" "${#uniref_dirs[@]}" >&2
  exit 1
fi
choco_dir="${choco_dirs[0]}"
uniref_dir="${uniref_dirs[0]}"

incompatible_log="${work_dir}/logs/humann-vJan26-expected-rejection.log"
incompatible_status="${work_dir}/logs/humann-vJan26-expected-rejection.status"
if [[ ! -s "${incompatible_status}" ]] ||
   [[ "$(tr -d '[:space:]' < "${incompatible_status}")" == "0" ]] ||
   ! grep -q 'v3 or vJun23' "${incompatible_log}"; then
  incompatible_dir="${work_dir}/incompatible-vJan26-attempt-$$"
  mkdir -p "${incompatible_dir}"
  set +e
  "${humann}" \
    --input "${combined_fastq}" \
    --output "${incompatible_dir}" \
    --threads 1 \
    --taxonomic-profile "${project_root}/data/small/15-metaphlan-frozen/profile-all.tsv" \
    --metaphlan "${metaphlan_compat_dir}" \
    --nucleotide-database "${choco_dir}" \
    --protein-database "${uniref_dir}" \
    --bypass-translated-search \
    --output-basename incompatible_vJan26 \
    > "${incompatible_log}" 2>&1
  status="$?"
  set -e
  printf '%s\n' "${status}" > "${incompatible_status}"
fi
if [[ "$(tr -d '[:space:]' < "${incompatible_status}")" == "0" ]] ||
   ! grep -q 'v3 or vJun23' "${incompatible_log}"; then
  printf 'HUMAnN vJan26 compatibility rejection did not match the expected contract\n' >&2
  exit 1
fi

humann_output="${work_dir}/humann"
humann_log="${work_dir}/logs/ERR9765746-humann3.log"
humann_args=(
  --input "${combined_fastq}"
  --output "${humann_output}"
  --threads 8
  --taxonomic-profile "${vjun_profile}"
  --metaphlan "${metaphlan_compat_dir}"
  --nucleotide-database "${choco_dir}"
  --protein-database "${uniref_dir}"
  --search-mode uniref90
  --prescreen-threshold 0.01
  --bowtie-options '--very-sensitive --seed 20260722'
  --memory-use minimum
  --pathways metacyc
  --minpath on
  --gap-fill on
  --log-level DEBUG
  --o-log "${humann_log}"
  --output-basename ERR9765746_MOCK1
  --output-format tsv
)
if [[ "${resume}" == "yes" ]]; then
  humann_args+=(--resume)
fi

if [[ ! -s "${humann_output}/ERR9765746_MOCK1_pathcoverage.tsv" ]]; then
  /usr/bin/time -v \
    -o "${work_dir}/logs/humann3.resources.txt" \
    "${humann}" "${humann_args[@]}" \
    > "${work_dir}/logs/humann3.stdout-stderr.log" 2>&1
fi

gene_rpk="${humann_output}/ERR9765746_MOCK1_genefamilies.tsv"
path_rpk="${humann_output}/ERR9765746_MOCK1_pathabundance.tsv"
path_coverage="${humann_output}/ERR9765746_MOCK1_pathcoverage.tsv"
for output_file in "${gene_rpk}" "${path_rpk}" "${path_coverage}"; do
  if [[ ! -s "${output_file}" ]]; then
    printf 'Missing HUMAnN output: %s\n' "${output_file}" >&2
    exit 1
  fi
done

if ! grep -q '^bypass prescreen = False$' "${humann_log}"; then
  printf 'HUMAnN run did not preserve taxonomic prescreening\n' >&2
  exit 1
fi
if ! grep -q 'Total species selected from prescreen: 89$' "${humann_log}"; then
  printf 'HUMAnN prescreen species count does not match the locked run\n' >&2
  exit 1
fi
custom_choco_files="$(
  awk '/Adding file to database:/ {count += 1} END {print count + 0}' "${humann_log}"
)"
if [[ "${custom_choco_files}" != "37" ]]; then
  printf 'Expected 37 sample-specific ChocoPhlAn files; observed %s\n' \
    "${custom_choco_files}" >&2
  exit 1
fi

gene_cpm="${work_dir}/humann/ERR9765746_MOCK1_genefamilies-cpm.tsv"
gene_relab="${work_dir}/humann/ERR9765746_MOCK1_genefamilies-relab.tsv"
path_cpm="${work_dir}/humann/ERR9765746_MOCK1_pathabundance-cpm.tsv"
path_relab="${work_dir}/humann/ERR9765746_MOCK1_pathabundance-relab.tsv"
reaction_rpk="${work_dir}/humann/ERR9765746_MOCK1_reactions-rpk.tsv"
reaction_cpm="${work_dir}/humann/ERR9765746_MOCK1_reactions-cpm.tsv"
reaction_relab="${work_dir}/humann/ERR9765746_MOCK1_reactions-relab.tsv"
"${renorm}" --input "${gene_rpk}" --units cpm --mode community --special y --update-snames --output "${gene_cpm}"
"${renorm}" --input "${gene_rpk}" --units relab --mode community --special y --update-snames --output "${gene_relab}"
"${renorm}" --input "${path_rpk}" --units cpm --mode community --special y --update-snames --output "${path_cpm}"
"${renorm}" --input "${path_rpk}" --units relab --mode community --special y --update-snames --output "${path_relab}"
"${regroup}" \
  --input "${gene_rpk}" \
  --groups uniref90_rxn \
  --function sum \
  --ungrouped Y \
  --protected Y \
  --output "${reaction_rpk}" \
  > "${work_dir}/logs/humann-regroup-uniref90-rxn.log" 2>&1
"${renorm}" --input "${reaction_rpk}" --units cpm --mode community --special y --update-snames --output "${reaction_cpm}"
"${renorm}" --input "${reaction_rpk}" --units relab --mode community --special y --update-snames --output "${reaction_relab}"

stage="${frozen_dir}.staging.$$"
if [[ -e "${stage}" ]]; then
  printf 'Frozen staging directory already exists: %s\n' "${stage}" >&2
  exit 1
fi
mkdir -p "${stage}"
cp "$0" "${stage}/commands.sh"
cp "${metaphlan_compat_dir}/metaphlan" "${stage}/metaphlan-version-compat"

"${python}" "${project_root}/scripts/validate_article19_humann3.py" \
  --project-root "${project_root}" \
  --environment-prefix "${environment_prefix}" \
  --frozen-dir "${stage}" \
  --initialize-frozen \
  --cache-root "${cache_root}" \
  --work-dir "${work_dir}" \
  --combined-fastq "${combined_fastq}" \
  --vjun-profile "${vjun_profile}" \
  --incompatible-log "${incompatible_log}" \
  --incompatible-status "${incompatible_status}" \
  --gene-rpk "${gene_rpk}" \
  --gene-cpm "${gene_cpm}" \
  --gene-relab "${gene_relab}" \
  --reaction-rpk "${reaction_rpk}" \
  --reaction-cpm "${reaction_cpm}" \
  --reaction-relab "${reaction_relab}" \
  --regroup-log "${work_dir}/logs/humann-regroup-uniref90-rxn.log" \
  --path-rpk "${path_rpk}" \
  --path-cpm "${path_cpm}" \
  --path-relab "${path_relab}" \
  --path-coverage "${path_coverage}"

mkdir -p "$(dirname "${frozen_dir}")"
mv "${stage}" "${frozen_dir}"
printf 'Article 19 frozen evidence initialized: %s\n' "${frozen_dir}"
