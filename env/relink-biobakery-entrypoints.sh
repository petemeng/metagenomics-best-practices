#!/usr/bin/env bash
set -euo pipefail

environment_name="${1:-metagenome-biobakery-2026.07}"

if ! command -v mamba >/dev/null 2>&1; then
  echo "ERROR: mamba is required to relink the exact executable packages." >&2
  exit 1
fi
if ! command -v conda >/dev/null 2>&1; then
  echo "ERROR: conda is required to audit the repaired environment." >&2
  exit 1
fi

# HUMAnN 3.9's source distribution installs bundled files named bowtie2 and
# diamond. Bioconda also installs current standalone packages, so the HUMAnN
# link step can leave older bundled executables at those paths even though
# conda-meta reports the requested modern versions. Relinking the two exact
# packages after HUMAnN restores the executable/package identity.
mamba install -y \
  -n "${environment_name}" \
  -c conda-forge \
  -c bioconda \
  --override-channels \
  --strict-channel-priority \
  --force-reinstall \
  "bowtie2=2.5.5" \
  "diamond=2.2.4"

bowtie2_report="$(
  conda run -n "${environment_name}" \
    env PYTHONNOUSERSITE=1 PYTHONPATH= \
    bowtie2 --version 2>&1
)"
diamond_report="$(
  conda run -n "${environment_name}" \
    env PYTHONNOUSERSITE=1 PYTHONPATH= \
    diamond version 2>&1
)"

if ! grep -Fq "version 2.5.5" <<<"${bowtie2_report}"; then
  echo "ERROR: Bowtie2 executable is not version 2.5.5 after relink." >&2
  exit 1
fi
if ! grep -Fq "diamond version 2.2.4" <<<"${diamond_report}"; then
  echo "ERROR: DIAMOND executable is not version 2.2.4 after relink." >&2
  exit 1
fi

echo "PASS Bowtie2 executable: 2.5.5"
echo "PASS DIAMOND executable: 2.2.4"
