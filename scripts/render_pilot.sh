#!/usr/bin/env bash
set -euo pipefail

if [[ "$#" -ne 2 ]]; then
  echo "Usage: render_pilot.sh <project-root> <run-root>" >&2
  exit 2
fi

project_root="$(cd "$1" && pwd)"
run_root="$(cd "$2" && pwd)"
stage="${run_root}/pilot"

quarto_bin="${QUARTO_BIN:-}"
if [[ -z "${quarto_bin}" ]]; then
  quarto_bin="$(command -v quarto || true)"
fi
if [[ -z "${quarto_bin}" ]]; then
  local_quarto="${project_root}/../16s/.tools/quarto-1.9.38/bin/quarto"
  if [[ -x "${local_quarto}" ]]; then
    quarto_bin="${local_quarto}"
  fi
fi
if [[ -z "${quarto_bin}" ]]; then
  echo "Quarto executable not found; set QUARTO_BIN." >&2
  exit 1
fi

mkdir -p \
  "${stage}/chapters" \
  "${stage}/data/small" \
  "${stage}/data/small/07-decontam" \
  "${stage}/env" \
  "${stage}/figures" \
  "${stage}/results" \
  "${run_root}/.cache/quarto"
cp "${project_root}/qa/pilot/_quarto.yml" "${stage}/_quarto.yml"
cp "${project_root}/index.qmd" "${stage}/index.qmd"
cp \
  "${project_root}/chapters/02-three-analysis-layers.qmd" \
  "${project_root}/chapters/03-study-design-power.qmd" \
  "${project_root}/chapters/04-sequencing-depth.qmd" \
  "${project_root}/chapters/05-library-prep-absolute-quantification.qmd" \
  "${project_root}/chapters/06-host-depletion-low-biomass.qmd" \
  "${project_root}/chapters/07-contamination-controls.qmd" \
  "${project_root}/chapters/08-fastq-short-long-reads.qmd" \
  "${project_root}/chapters/09-wsl2-conda.qmd" \
  "${project_root}/chapters/10-computing-hpc-cloud.qmd" \
  "${project_root}/chapters/11-install-biobakery-assembly.qmd" \
  "${project_root}/chapters/12-install-r-cmd.qmd" \
  "${project_root}/chapters/13-read-qc-fastp.qmd" \
  "${project_root}/chapters/14-host-removal-complexity-duplicates.qmd" \
  "${project_root}/chapters/15-metaphlan4.qmd" \
  "${project_root}/chapters/16-kraken2-bracken.qmd" \
  "${project_root}/chapters/17-kraken2-database-confidence.qmd" \
  "${project_root}/chapters/18-profiler-benchmark.qmd" \
  "${project_root}/chapters/19-humann3.qmd" \
  "${project_root}/chapters/20-functional-profile-normalization.qmd" \
  "${project_root}/chapters/21-metagenomic-table-semantics.qmd" \
  "${project_root}/chapters/22-alpha-beta-diversity.qmd" \
  "${project_root}/chapters/23-pcoa-cap-permanova.qmd" \
  "${project_root}/chapters/24-differential-abundance.qmd" \
  "${project_root}/chapters/25-composition-core-microbiome.qmd" \
  "${project_root}/chapters/26-curated-metagenomic-data-lineage.qmd" \
  "${project_root}/chapters/27-machine-learning-roc.qmd" \
  "${project_root}/chapters/28-cross-cohort-validation.qmd" \
  "${project_root}/chapters/29-cooccurrence-network.qmd" \
  "${project_root}/chapters/30-short-read-assembly.qmd" \
  "${project_root}/chapters/31-long-read-assembly.qmd" \
  "${project_root}/chapters/32-hybrid-assembly-polishing.qmd" \
  "${project_root}/chapters/33-assembly-qc.qmd" \
  "${project_root}/chapters/34-nonredundant-gene-catalog.qmd" \
  "${project_root}/chapters/35-gene-abundance.qmd" \
  "${project_root}/chapters/36-eggnog-functional-annotation.qmd" \
  "${project_root}/chapters/37-cazymes-dbcan.qmd" \
  "${project_root}/chapters/38-resistome.qmd" \
  "${project_root}/chapters/39-virulome.qmd" \
  "${project_root}/chapters/40-bgc-natural-products.qmd" \
  "${project_root}/chapters/41-read-mapping-depth.qmd" \
  "${project_root}/chapters/42-binning.qmd" \
  "${project_root}/chapters/43-bin-refinement.qmd" \
  "${project_root}/chapters/44-mag-qc-checkm2-gunc-mimag.qmd" \
  "${stage}/chapters/"
cp "${project_root}/styles.scss" "${stage}/styles.scss"
cp "${project_root}/references.bib" "${stage}/references.bib"
cp \
  "${project_root}/data/small/01-assay-evidence.tsv" \
  "${stage}/data/small/"
cp \
  "${run_root}/data/small/01-crc-cohort-summary.tsv" \
  "${stage}/data/small/"
cp \
  "${project_root}/data/small/02-layer-capabilities.tsv" \
  "${project_root}/data/small/02-published-anchors.tsv" \
  "${stage}/data/small/"
cp \
  "${run_root}/data/small/03-crc-design-audit.tsv" \
  "${stage}/data/small/"
cp \
  "${project_root}/data/small/04-lake-lanier.npo" \
  "${project_root}/data/small/04-depth-evidence.tsv" \
  "${stage}/data/small/"
cp \
  "${run_root}/data/small/04-crc-library-size.tsv" \
  "${run_root}/data/small/04-lake-lanier-coverage.tsv" \
  "${stage}/data/small/"
cp \
  "${run_root}/data/small/05-costea-mock-profiles.tsv" \
  "${run_root}/data/small/05-costea-bias-summary.tsv" \
  "${stage}/data/small/"
cp \
  "${project_root}/data/small/05-syndna-mock-benchmark.tsv" \
  "${project_root}/data/small/05-control-placement.tsv" \
  "${stage}/data/small/"
cp \
  "${project_root}/data/small/06-marotz-host-depletion-source.tsv" \
  "${project_root}/data/small/06-longhi-saponin-source.tsv" \
  "${project_root}/data/small/06-host-filter-contract.tsv" \
  "${project_root}/data/small/06-data-NOTICE.txt" \
  "${stage}/data/small/"
cp \
  "${run_root}/data/small/06-host-depletion-read-budget.tsv" \
  "${run_root}/data/small/06-saponin-tradeoff.tsv" \
  "${stage}/data/small/"
cp \
  "${run_root}/figures/06-host-depletion-efficiency.pdf" \
  "${run_root}/figures/06-host-depletion-efficiency.png" \
  "${run_root}/figures/06-host-depletion-efficiency.tiff" \
  "${run_root}/figures/06-saponin-tradeoff.pdf" \
  "${run_root}/figures/06-saponin-tradeoff.png" \
  "${run_root}/figures/06-saponin-tradeoff.tiff" \
  "${run_root}/figures/06-host-depletion-decision.pdf" \
  "${run_root}/figures/06-host-depletion-decision.png" \
  "${run_root}/figures/06-host-depletion-decision.tiff" \
  "${stage}/figures/"
cp \
  "${project_root}/data/small/07-decontam/otutab.tsv" \
  "${project_root}/data/small/07-decontam/taxonomy.tsv" \
  "${project_root}/data/small/07-decontam/metadata.tsv" \
  "${project_root}/data/small/07-decontam/source_summary.json" \
  "${stage}/data/small/07-decontam/"
cp \
  "${project_root}/data/small/07-salter-shotgun-evidence.tsv" \
  "${project_root}/data/small/07-index-hopping-evidence.tsv" \
  "${project_root}/data/small/07-data-NOTICE.txt" \
  "${stage}/data/small/"
cp \
  "${run_root}/figures/07-control-library-size.pdf" \
  "${run_root}/figures/07-control-library-size.png" \
  "${run_root}/figures/07-control-library-size.tiff" \
  "${run_root}/figures/07-contaminant-prevalence.pdf" \
  "${run_root}/figures/07-contaminant-prevalence.png" \
  "${run_root}/figures/07-contaminant-prevalence.tiff" \
  "${run_root}/figures/07-contaminant-burden.pdf" \
  "${run_root}/figures/07-contaminant-burden.png" \
  "${run_root}/figures/07-contaminant-burden.tiff" \
  "${run_root}/figures/07-index-hopping-evidence.pdf" \
  "${run_root}/figures/07-index-hopping-evidence.png" \
  "${run_root}/figures/07-index-hopping-evidence.tiff" \
  "${stage}/figures/"
cp \
  "${project_root}/data/small/08-ena-fastq-sources.tsv" \
  "${project_root}/data/small/08-platform-benchmark.tsv" \
  "${project_root}/data/small/08-native-format-contract.tsv" \
  "${project_root}/data/small/08-read-prefix-metrics.tsv" \
  "${project_root}/data/small/08-fastq-anatomy.tsv" \
  "${project_root}/data/small/08-prefix-source-summary.json" \
  "${project_root}/data/small/08-data-NOTICE.txt" \
  "${stage}/data/small/"
cp \
  "${project_root}/data/small/09-environment-contract.tsv" \
  "${project_root}/data/small/09-data-NOTICE.txt" \
  "${stage}/data/small/"
cp \
  "${project_root}/data/small/10-job-array.tsv" \
  "${project_root}/data/small/10-runtime-contract.tsv" \
  "${project_root}/data/small/10-container-smoke.log" \
  "${project_root}/data/small/10-data-NOTICE.txt" \
  "${stage}/data/small/"
cp \
  "${project_root}/env/platform-smoke.yml" \
  "${stage}/env/"
cp \
  "${project_root}/env/assembly.yml" \
  "${project_root}/env/assembly-linux-64.lock" \
  "${project_root}/env/biobakery.yml" \
  "${project_root}/env/biobakery-linux-64.lock" \
  "${project_root}/env/relink-biobakery-entrypoints.sh" \
  "${stage}/env/"
cp \
  "${run_root}/figures/09-wsl2-layer-map.pdf" \
  "${run_root}/figures/09-wsl2-layer-map.png" \
  "${run_root}/figures/09-wsl2-layer-map.tiff" \
  "${run_root}/figures/09-environment-validation.pdf" \
  "${run_root}/figures/09-environment-validation.png" \
  "${run_root}/figures/09-environment-validation.tiff" \
  "${stage}/figures/"
cp \
  "${run_root}/figures/10-input-resource-budget.pdf" \
  "${run_root}/figures/10-input-resource-budget.png" \
  "${run_root}/figures/10-input-resource-budget.tiff" \
  "${run_root}/figures/10-resource-control-loop.pdf" \
  "${run_root}/figures/10-resource-control-loop.png" \
  "${run_root}/figures/10-resource-control-loop.tiff" \
  "${run_root}/figures/10-restart-safe-array.pdf" \
  "${run_root}/figures/10-restart-safe-array.png" \
  "${run_root}/figures/10-restart-safe-array.tiff" \
  "${stage}/figures/"
cp \
  "${project_root}/data/small/11-environment-evidence.tsv" \
  "${project_root}/data/small/11-install-self-tests.log" \
  "${project_root}/data/small/11-solver-audit.tsv" \
  "${project_root}/data/small/11-database-manifest.tsv" \
  "${project_root}/data/small/11-data-NOTICE.txt" \
  "${stage}/data/small/"
cp \
  "${run_root}/figures/11-environment-boundaries.pdf" \
  "${run_root}/figures/11-environment-boundaries.png" \
  "${run_root}/figures/11-environment-boundaries.tiff" \
  "${run_root}/figures/11-toolchain-entrypoints.pdf" \
  "${run_root}/figures/11-toolchain-entrypoints.png" \
  "${run_root}/figures/11-toolchain-entrypoints.tiff" \
  "${run_root}/figures/11-database-storage-contract.pdf" \
  "${run_root}/figures/11-database-storage-contract.png" \
  "${run_root}/figures/11-database-storage-contract.tiff" \
  "${stage}/figures/"
cp \
  "${project_root}/env/renv.lock" \
  "${project_root}/env/machine-learning-renv.lock" \
  "${project_root}/env/cross-cohort-renv.lock" \
  "${project_root}/env/network-renv.lock" \
  "${stage}/env/"
cp \
  "${project_root}/data/small/12-package-contract.tsv" \
  "${project_root}/data/small/12-cmd-resource-manifest.tsv" \
  "${project_root}/data/small/12-cmd-asnicarf-2017-relative-abundance.rds" \
  "${project_root}/data/small/12-resource-retrieval.log" \
  "${project_root}/data/small/12-data-NOTICE.txt" \
  "${stage}/data/small/"
cp \
  "${run_root}/figures/12-r-data-access-boundaries.pdf" \
  "${run_root}/figures/12-r-data-access-boundaries.png" \
  "${run_root}/figures/12-r-data-access-boundaries.tiff" \
  "${run_root}/figures/12-package-role-contract.pdf" \
  "${run_root}/figures/12-package-role-contract.png" \
  "${run_root}/figures/12-package-role-contract.tiff" \
  "${run_root}/figures/12-cmd-object-contract.pdf" \
  "${run_root}/figures/12-cmd-object-contract.png" \
  "${run_root}/figures/12-cmd-object-contract.tiff" \
  "${stage}/figures/"
cp \
  "${project_root}/env/read-qc.yml" \
  "${project_root}/env/read-qc-linux-64.lock" \
  "${stage}/env/"
cp \
  "${project_root}/data/small/13-source-manifest.tsv" \
  "${project_root}/data/small/13-data-NOTICE.txt" \
  "${stage}/data/small/"
cp -R \
  "${project_root}/data/small/13-qc-frozen" \
  "${stage}/data/small/"
cp \
  "${run_root}/figures/13-per-cycle-quality.pdf" \
  "${run_root}/figures/13-per-cycle-quality.png" \
  "${run_root}/figures/13-per-cycle-quality.tiff" \
  "${run_root}/figures/13-read-pair-fate.pdf" \
  "${run_root}/figures/13-read-pair-fate.png" \
  "${run_root}/figures/13-read-pair-fate.tiff" \
  "${run_root}/figures/13-fastqc-module-states.pdf" \
  "${run_root}/figures/13-fastqc-module-states.png" \
  "${run_root}/figures/13-fastqc-module-states.tiff" \
  "${stage}/figures/"
cp \
  "${project_root}/env/host-removal.yml" \
  "${project_root}/env/host-removal-linux-64.lock" \
  "${stage}/env/"
cp \
  "${project_root}/data/small/14-source-manifest.tsv" \
  "${project_root}/data/small/14-index-manifest.tsv" \
  "${project_root}/data/small/14-data-NOTICE.txt" \
  "${stage}/data/small/"
cp -R \
  "${project_root}/data/small/14-host-removal-frozen" \
  "${stage}/data/small/"
cp \
  "${run_root}/figures/14-host-removal-recall-retention.pdf" \
  "${run_root}/figures/14-host-removal-recall-retention.png" \
  "${run_root}/figures/14-host-removal-recall-retention.tiff" \
  "${run_root}/figures/14-complexity-sensitivity.pdf" \
  "${run_root}/figures/14-complexity-sensitivity.png" \
  "${run_root}/figures/14-complexity-sensitivity.tiff" \
  "${run_root}/figures/14-duplicate-decision.pdf" \
  "${run_root}/figures/14-duplicate-decision.png" \
  "${run_root}/figures/14-duplicate-decision.tiff" \
  "${stage}/figures/"
cp \
  "${project_root}/data/small/15-source-manifest.tsv" \
  "${project_root}/data/small/15-database-manifest.tsv" \
  "${project_root}/data/small/15-data-NOTICE.txt" \
  "${stage}/data/small/"
cp -R \
  "${project_root}/data/small/15-metaphlan-frozen" \
  "${stage}/data/small/"
cp \
  "${run_root}/figures/15-metaphlan-composition.pdf" \
  "${run_root}/figures/15-metaphlan-composition.png" \
  "${run_root}/figures/15-metaphlan-composition.tiff" \
  "${run_root}/figures/15-sgb-marker-support.pdf" \
  "${run_root}/figures/15-sgb-marker-support.png" \
  "${run_root}/figures/15-sgb-marker-support.tiff" \
  "${run_root}/figures/15-detection-quantification-sensitivity.pdf" \
  "${run_root}/figures/15-detection-quantification-sensitivity.png" \
  "${run_root}/figures/15-detection-quantification-sensitivity.tiff" \
  "${stage}/figures/"
cp \
  "${project_root}/env/kraken.yml" \
  "${project_root}/env/kraken-linux-64.lock" \
  "${stage}/env/"
cp \
  "${project_root}/data/small/16-source-manifest.tsv" \
  "${project_root}/data/small/16-database-manifest.tsv" \
  "${project_root}/data/small/16-standard8-files.md5" \
  "${project_root}/data/small/16-data-NOTICE.txt" \
  "${stage}/data/small/"
cp -R \
  "${project_root}/data/small/16-kraken-bracken-frozen" \
  "${stage}/data/small/"
cp \
  "${run_root}/figures/16-kraken-classification-ledger.pdf" \
  "${run_root}/figures/16-kraken-classification-ledger.png" \
  "${run_root}/figures/16-kraken-classification-ledger.tiff" \
  "${run_root}/figures/16-bracken-redistribution.pdf" \
  "${run_root}/figures/16-bracken-redistribution.png" \
  "${run_root}/figures/16-bracken-redistribution.tiff" \
  "${run_root}/figures/16-bracken-parameter-sensitivity.pdf" \
  "${run_root}/figures/16-bracken-parameter-sensitivity.png" \
  "${run_root}/figures/16-bracken-parameter-sensitivity.tiff" \
  "${stage}/figures/"
cp \
  "${project_root}/data/small/17-source-manifest.tsv" \
  "${project_root}/data/small/17-database-manifest.tsv" \
  "${project_root}/data/small/17-standard16-files.md5" \
  "${project_root}/data/small/17-pluspf8-files.md5" \
  "${project_root}/data/small/17-mock1-truth.tsv" \
  "${project_root}/data/small/17-ncbi-genome-snapshot.jsonl" \
  "${project_root}/data/small/17-ncbi-sequence-snapshot.jsonl" \
  "${project_root}/data/small/17-truth-provenance.tsv" \
  "${project_root}/data/small/17-data-NOTICE.txt" \
  "${stage}/data/small/"
cp -R \
  "${project_root}/data/small/17-kraken-database-confidence-frozen" \
  "${stage}/data/small/"
cp \
  "${run_root}/figures/17-reference-coverage.pdf" \
  "${run_root}/figures/17-reference-coverage.png" \
  "${run_root}/figures/17-reference-coverage.tiff" \
  "${run_root}/figures/17-confidence-tradeoff.pdf" \
  "${run_root}/figures/17-confidence-tradeoff.png" \
  "${run_root}/figures/17-confidence-tradeoff.tiff" \
  "${run_root}/figures/17-database-stability.pdf" \
  "${run_root}/figures/17-database-stability.png" \
  "${run_root}/figures/17-database-stability.tiff" \
  "${run_root}/figures/17-hit-group-control.pdf" \
  "${run_root}/figures/17-hit-group-control.png" \
  "${run_root}/figures/17-hit-group-control.tiff" \
  "${stage}/figures/"
cp \
  "${project_root}/env/motus.yml" \
  "${project_root}/env/motus-linux-64.lock" \
  "${stage}/env/"
cp \
  "${project_root}/data/small/18-source-manifest.tsv" \
  "${project_root}/data/small/18-database-manifest.tsv" \
  "${project_root}/data/small/18-data-NOTICE.txt" \
  "${stage}/data/small/"
cp -R \
  "${project_root}/data/small/18-profiler-benchmark-frozen" \
  "${stage}/data/small/"
cp \
  "${run_root}/figures/18-feature-space-crosswalk.pdf" \
  "${run_root}/figures/18-feature-space-crosswalk.png" \
  "${run_root}/figures/18-feature-space-crosswalk.tiff" \
  "${run_root}/figures/18-recovery-nontruth.pdf" \
  "${run_root}/figures/18-recovery-nontruth.png" \
  "${run_root}/figures/18-recovery-nontruth.tiff" \
  "${run_root}/figures/18-composition-agreement.pdf" \
  "${run_root}/figures/18-composition-agreement.png" \
  "${run_root}/figures/18-composition-agreement.tiff" \
  "${run_root}/figures/18-resource-footprint.pdf" \
  "${run_root}/figures/18-resource-footprint.png" \
  "${run_root}/figures/18-resource-footprint.tiff" \
  "${stage}/figures/"
cp \
  "${project_root}/data/small/19-source-manifest.tsv" \
  "${project_root}/data/small/19-database-manifest.tsv" \
  "${project_root}/data/small/19-data-NOTICE.txt" \
  "${stage}/data/small/"
cp -R \
  "${project_root}/data/small/19-humann3-frozen" \
  "${stage}/data/small/"
cp \
  "${run_root}/figures/19-read-flow.pdf" \
  "${run_root}/figures/19-read-flow.png" \
  "${run_root}/figures/19-read-flow.tiff" \
  "${run_root}/figures/19-gene-family-stratification.pdf" \
  "${run_root}/figures/19-gene-family-stratification.png" \
  "${run_root}/figures/19-gene-family-stratification.tiff" \
  "${run_root}/figures/19-pathway-contributions.pdf" \
  "${run_root}/figures/19-pathway-contributions.png" \
  "${run_root}/figures/19-pathway-contributions.tiff" \
  "${run_root}/figures/19-abundance-coverage.pdf" \
  "${run_root}/figures/19-abundance-coverage.png" \
  "${run_root}/figures/19-abundance-coverage.tiff" \
  "${stage}/figures/"
cp -R \
  "${project_root}/data/small/20-cmd-pathway" \
  "${project_root}/data/small/20-functional-normalization-frozen" \
  "${stage}/data/small/"
cp \
  "${run_root}/figures/20-normalization-denominators.pdf" \
  "${run_root}/figures/20-normalization-denominators.png" \
  "${run_root}/figures/20-normalization-denominators.tiff" \
  "${run_root}/figures/20-special-feature-budget.pdf" \
  "${run_root}/figures/20-special-feature-budget.png" \
  "${run_root}/figures/20-special-feature-budget.tiff" \
  "${run_root}/figures/20-pathway-contributions.pdf" \
  "${run_root}/figures/20-pathway-contributions.png" \
  "${run_root}/figures/20-pathway-contributions.tiff" \
  "${run_root}/figures/20-prevalence-zero-sensitivity.pdf" \
  "${run_root}/figures/20-prevalence-zero-sensitivity.png" \
  "${run_root}/figures/20-prevalence-zero-sensitivity.tiff" \
  "${stage}/figures/"
cp -R \
  "${project_root}/data/small/21-table-semantics-frozen" \
  "${stage}/data/small/"
cp \
  "${project_root}/data/small/21-data-NOTICE.txt" \
  "${stage}/data/small/"
cp \
  "${run_root}/figures/21-table-unit-map.pdf" \
  "${run_root}/figures/21-table-unit-map.png" \
  "${run_root}/figures/21-table-unit-map.tiff" \
  "${run_root}/figures/21-denominator-closure.pdf" \
  "${run_root}/figures/21-denominator-closure.png" \
  "${run_root}/figures/21-denominator-closure.tiff" \
  "${run_root}/figures/21-genome-equivalent-calibration.pdf" \
  "${run_root}/figures/21-genome-equivalent-calibration.png" \
  "${run_root}/figures/21-genome-equivalent-calibration.tiff" \
  "${run_root}/figures/21-zero-strata-semantics.pdf" \
  "${run_root}/figures/21-zero-strata-semantics.png" \
  "${run_root}/figures/21-zero-strata-semantics.tiff" \
  "${stage}/figures/"
cp -R \
  "${project_root}/data/small/22-diversity-inputs" \
  "${stage}/data/small/"
cp \
  "${project_root}/data/small/22-data-NOTICE.txt" \
  "${stage}/data/small/"
cp \
  "${project_root}/figures/22-korchagina-fig5-original.png" \
  "${run_root}/figures/22-resolution-boundaries.pdf" \
  "${run_root}/figures/22-resolution-boundaries.png" \
  "${run_root}/figures/22-resolution-boundaries.tiff" \
  "${run_root}/figures/22-alpha-hill-numbers.pdf" \
  "${run_root}/figures/22-alpha-hill-numbers.png" \
  "${run_root}/figures/22-alpha-hill-numbers.tiff" \
  "${run_root}/figures/22-beta-ordination.pdf" \
  "${run_root}/figures/22-beta-ordination.png" \
  "${run_root}/figures/22-beta-ordination.tiff" \
  "${run_root}/figures/22-sensitivity-recruitment.pdf" \
  "${run_root}/figures/22-sensitivity-recruitment.png" \
  "${run_root}/figures/22-sensitivity-recruitment.tiff" \
  "${stage}/figures/"
cp -R \
  "${project_root}/data/small/23-ordination-permanova" \
  "${stage}/data/small/"
cp \
  "${project_root}/data/small/23-data-NOTICE.txt" \
  "${stage}/data/small/"
cp \
  "${project_root}/figures/23-korchagina-fig1-original.png" \
  "${run_root}/figures/23-design-permutation-space.pdf" \
  "${run_root}/figures/23-design-permutation-space.png" \
  "${run_root}/figures/23-design-permutation-space.tiff" \
  "${run_root}/figures/23-pcoa-cap.pdf" \
  "${run_root}/figures/23-pcoa-cap.png" \
  "${run_root}/figures/23-pcoa-cap.tiff" \
  "${run_root}/figures/23-permanova-dispersion.pdf" \
  "${run_root}/figures/23-permanova-dispersion.png" \
  "${run_root}/figures/23-permanova-dispersion.tiff" \
  "${run_root}/figures/23-pairwise-sensitivity.pdf" \
  "${run_root}/figures/23-pairwise-sensitivity.png" \
  "${run_root}/figures/23-pairwise-sensitivity.tiff" \
  "${stage}/figures/"
cp -R \
  "${project_root}/data/small/24-differential-abundance" \
  "${stage}/data/small/"
cp -R \
  "${run_root}/results/24-differential-abundance" \
  "${stage}/results/"
cp \
  "${project_root}/figures/24-zeller-fig1-original.png" \
  "${run_root}/figures/24-zeller-marker-redraw.pdf" \
  "${run_root}/figures/24-zeller-marker-redraw.png" \
  "${run_root}/figures/24-zeller-marker-redraw.tiff" \
  "${run_root}/figures/24-maaslin3-two-part.pdf" \
  "${run_root}/figures/24-maaslin3-two-part.png" \
  "${run_root}/figures/24-maaslin3-two-part.tiff" \
  "${run_root}/figures/24-functional-associations.pdf" \
  "${run_root}/figures/24-functional-associations.png" \
  "${run_root}/figures/24-functional-associations.tiff" \
  "${run_root}/figures/24-method-concordance.pdf" \
  "${run_root}/figures/24-method-concordance.png" \
  "${run_root}/figures/24-method-concordance.tiff" \
  "${stage}/figures/"
cp -R \
  "${run_root}/data/small/25-composition-core" \
  "${stage}/data/small/"
cp \
  "${run_root}/data/small/25-data-NOTICE.txt" \
  "${stage}/data/small/"
cp -R \
  "${run_root}/results/25-composition-core" \
  "${stage}/results/"
cp \
  "${project_root}/figures/25-hmp-fig3-original.jpg" \
  "${run_root}/figures/25-multirank-mean-composition.pdf" \
  "${run_root}/figures/25-multirank-mean-composition.png" \
  "${run_root}/figures/25-multirank-mean-composition.tiff" \
  "${run_root}/figures/25-individual-stool-composition.pdf" \
  "${run_root}/figures/25-individual-stool-composition.png" \
  "${run_root}/figures/25-individual-stool-composition.tiff" \
  "${run_root}/figures/25-prevalence-abundance.pdf" \
  "${run_root}/figures/25-prevalence-abundance.png" \
  "${run_root}/figures/25-prevalence-abundance.tiff" \
  "${run_root}/figures/25-core-membership-sensitivity.pdf" \
  "${run_root}/figures/25-core-membership-sensitivity.png" \
  "${run_root}/figures/25-core-membership-sensitivity.tiff" \
  "${stage}/figures/"
cp -R \
  "${run_root}/data/small/26-cmd-lineage" \
  "${stage}/data/small/"
cp \
  "${run_root}/data/small/26-data-NOTICE.txt" \
  "${stage}/data/small/"
cp -R \
  "${run_root}/results/26-cmd-lineage" \
  "${stage}/results/"
cp \
  "${project_root}/figures/26-cmd3-fig1-original.png" \
  "${run_root}/figures/26-resource-release-audit.pdf" \
  "${run_root}/figures/26-resource-release-audit.png" \
  "${run_root}/figures/26-resource-release-audit.tiff" \
  "${run_root}/figures/26-metadata-completeness.pdf" \
  "${run_root}/figures/26-metadata-completeness.png" \
  "${run_root}/figures/26-metadata-completeness.tiff" \
  "${run_root}/figures/26-query-attrition.pdf" \
  "${run_root}/figures/26-query-attrition.png" \
  "${run_root}/figures/26-query-attrition.tiff" \
  "${run_root}/figures/26-lineage-compatibility.pdf" \
  "${run_root}/figures/26-lineage-compatibility.png" \
  "${run_root}/figures/26-lineage-compatibility.tiff" \
  "${stage}/figures/"
cp -R \
  "${run_root}/data/small/27-machine-learning" \
  "${stage}/data/small/"
cp \
  "${run_root}/data/small/27-data-NOTICE.txt" \
  "${stage}/data/small/"
cp -R \
  "${run_root}/results/27-machine-learning" \
  "${stage}/results/"
cp \
  "${project_root}/figures/27-zeller-fig1-original.png" \
  "${run_root}/figures/27-nested-roc.pdf" \
  "${run_root}/figures/27-nested-roc.png" \
  "${run_root}/figures/27-nested-roc.tiff" \
  "${run_root}/figures/27-model-performance.pdf" \
  "${run_root}/figures/27-model-performance.png" \
  "${run_root}/figures/27-model-performance.tiff" \
  "${run_root}/figures/27-permutation-importance.pdf" \
  "${run_root}/figures/27-permutation-importance.png" \
  "${run_root}/figures/27-permutation-importance.tiff" \
  "${run_root}/figures/27-leakage-permutation-audit.pdf" \
  "${run_root}/figures/27-leakage-permutation-audit.png" \
  "${run_root}/figures/27-leakage-permutation-audit.tiff" \
  "${stage}/figures/"
cp -R \
  "${run_root}/data/small/28-cross-cohort" \
  "${stage}/data/small/"
cp \
  "${run_root}/data/small/28-data-NOTICE.txt" \
  "${stage}/data/small/"
cp -R \
  "${run_root}/results/28-cross-cohort" \
  "${stage}/results/"
cp \
  "${project_root}/figures/28-wirbel-fig3-original.png" \
  "${run_root}/figures/28-meta-signature.pdf" \
  "${run_root}/figures/28-meta-signature.png" \
  "${run_root}/figures/28-meta-signature.tiff" \
  "${run_root}/figures/28-single-study-transfer.pdf" \
  "${run_root}/figures/28-single-study-transfer.png" \
  "${run_root}/figures/28-single-study-transfer.tiff" \
  "${run_root}/figures/28-lodo-forest.pdf" \
  "${run_root}/figures/28-lodo-forest.png" \
  "${run_root}/figures/28-lodo-forest.tiff" \
  "${run_root}/figures/28-lodo-roc.pdf" \
  "${run_root}/figures/28-lodo-roc.png" \
  "${run_root}/figures/28-lodo-roc.tiff" \
  "${run_root}/figures/28-validation-stability.pdf" \
  "${run_root}/figures/28-validation-stability.png" \
  "${run_root}/figures/28-validation-stability.tiff" \
  "${stage}/figures/"
cp -R \
  "${run_root}/data/small/29-network" \
  "${stage}/data/small/"
cp \
  "${run_root}/data/small/29-data-NOTICE.txt" \
  "${stage}/data/small/"
cp -R \
  "${run_root}/results/29-network" \
  "${stage}/results/"
cp \
  "${project_root}/figures/29-kurtz-fig2-original.png" \
  "${run_root}/figures/29-feature-contract.pdf" \
  "${run_root}/figures/29-feature-contract.png" \
  "${run_root}/figures/29-feature-contract.tiff" \
  "${run_root}/figures/29-conditional-network.pdf" \
  "${run_root}/figures/29-conditional-network.png" \
  "${run_root}/figures/29-conditional-network.tiff" \
  "${run_root}/figures/29-zippi-roles.pdf" \
  "${run_root}/figures/29-zippi-roles.png" \
  "${run_root}/figures/29-zippi-roles.tiff" \
  "${run_root}/figures/29-edge-stability-sensitivity.pdf" \
  "${run_root}/figures/29-edge-stability-sensitivity.png" \
  "${run_root}/figures/29-edge-stability-sensitivity.tiff" \
  "${run_root}/figures/29-robustness-null.pdf" \
  "${run_root}/figures/29-robustness-null.png" \
  "${run_root}/figures/29-robustness-null.tiff" \
  "${stage}/figures/"
cp \
  "${project_root}/data/small/30-source-manifest.tsv" \
  "${project_root}/data/small/30-data-NOTICE.txt" \
  "${stage}/data/small/"
cp -R \
  "${project_root}/data/small/30-short-read-assembly-frozen" \
  "${stage}/data/small/"
cp -R \
  "${run_root}/results/30-short-read-assembly" \
  "${stage}/results/"
cp \
  "${run_root}/figures/30-assembly-branch-design.pdf" \
  "${run_root}/figures/30-assembly-branch-design.png" \
  "${run_root}/figures/30-assembly-branch-design.tiff" \
  "${run_root}/figures/30-contiguity-output.pdf" \
  "${run_root}/figures/30-contiguity-output.png" \
  "${run_root}/figures/30-contiguity-output.tiff" \
  "${run_root}/figures/30-recruitment-tradeoff.pdf" \
  "${run_root}/figures/30-recruitment-tradeoff.png" \
  "${run_root}/figures/30-recruitment-tradeoff.tiff" \
  "${run_root}/figures/30-resource-footprint.pdf" \
  "${run_root}/figures/30-resource-footprint.png" \
  "${run_root}/figures/30-resource-footprint.tiff" \
  "${stage}/figures/"
cp \
  "${project_root}/data/small/31-source-manifest.tsv" \
  "${project_root}/data/small/31-data-NOTICE.txt" \
  "${project_root}/data/small/31-software-releases.tsv" \
  "${stage}/data/small/"
cp -R \
  "${project_root}/data/small/31-long-read-assembly-frozen" \
  "${stage}/data/small/"
cp -R \
  "${run_root}/results/31-long-read-assembly" \
  "${stage}/results/"
cp \
  "${run_root}/figures/31-long-read-branch-design.pdf" \
  "${run_root}/figures/31-long-read-branch-design.png" \
  "${run_root}/figures/31-long-read-branch-design.tiff" \
  "${run_root}/figures/31-long-read-contiguity.pdf" \
  "${run_root}/figures/31-long-read-contiguity.png" \
  "${run_root}/figures/31-long-read-contiguity.tiff" \
  "${run_root}/figures/31-long-read-readback.pdf" \
  "${run_root}/figures/31-long-read-readback.png" \
  "${run_root}/figures/31-long-read-readback.tiff" \
  "${run_root}/figures/31-circular-resource-audit.pdf" \
  "${run_root}/figures/31-circular-resource-audit.png" \
  "${run_root}/figures/31-circular-resource-audit.tiff" \
  "${stage}/figures/"
cp \
  "${project_root}/data/small/32-source-manifest.tsv" \
  "${project_root}/data/small/32-reference-manifest.tsv" \
  "${project_root}/data/small/32-software-releases.tsv" \
  "${project_root}/data/small/32-branch-contract.tsv" \
  "${project_root}/data/small/32-data-NOTICE.txt" \
  "${stage}/data/small/"
cp -R \
  "${project_root}/data/small/32-hybrid-assembly-polishing-frozen" \
  "${stage}/data/small/"
cp -R \
  "${run_root}/results/32-hybrid-assembly-polishing" \
  "${stage}/results/"
cp \
  "${run_root}/figures/32-hybrid-branch-design.pdf" \
  "${run_root}/figures/32-hybrid-branch-design.png" \
  "${run_root}/figures/32-hybrid-branch-design.tiff" \
  "${run_root}/figures/32-recovery-contiguity.pdf" \
  "${run_root}/figures/32-recovery-contiguity.png" \
  "${run_root}/figures/32-recovery-contiguity.tiff" \
  "${run_root}/figures/32-consensus-error.pdf" \
  "${run_root}/figures/32-consensus-error.png" \
  "${run_root}/figures/32-consensus-error.tiff" \
  "${run_root}/figures/32-abundance-resource-audit.pdf" \
  "${run_root}/figures/32-abundance-resource-audit.png" \
  "${run_root}/figures/32-abundance-resource-audit.tiff" \
  "${stage}/figures/"
cp \
  "${project_root}/data/small/33-source-manifest.tsv" \
  "${project_root}/data/small/33-software-releases.tsv" \
  "${project_root}/data/small/33-data-NOTICE.txt" \
  "${stage}/data/small/"
cp -R \
  "${project_root}/data/small/33-assembly-qc-frozen" \
  "${stage}/data/small/"
cp -R \
  "${run_root}/results/33-assembly-qc" \
  "${stage}/results/"
cp \
  "${run_root}/figures/33-qc-evidence-ladder.pdf" \
  "${run_root}/figures/33-qc-evidence-ladder.png" \
  "${run_root}/figures/33-qc-evidence-ladder.tiff" \
  "${run_root}/figures/33-n50-na50-correctness.pdf" \
  "${run_root}/figures/33-n50-na50-correctness.png" \
  "${run_root}/figures/33-n50-na50-correctness.tiff" \
  "${run_root}/figures/33-recovery-error-tradeoff.pdf" \
  "${run_root}/figures/33-recovery-error-tradeoff.png" \
  "${run_root}/figures/33-recovery-error-tradeoff.tiff" \
  "${run_root}/figures/33-diagnostic-task-gates.pdf" \
  "${run_root}/figures/33-diagnostic-task-gates.png" \
  "${run_root}/figures/33-diagnostic-task-gates.tiff" \
  "${stage}/figures/"
cp -R \
  "${project_root}/data/small/34-nonredundant-gene-catalog-frozen" \
  "${stage}/data/small/"
cp -R \
  "${run_root}/results/34-nonredundant-gene-catalog" \
  "${stage}/results/"
cp \
  "${run_root}/figures/34-gene-catalog-workflow.pdf" \
  "${run_root}/figures/34-gene-catalog-workflow.png" \
  "${run_root}/figures/34-gene-catalog-workflow.tiff" \
  "${run_root}/figures/34-gene-length-distributions.pdf" \
  "${run_root}/figures/34-gene-length-distributions.png" \
  "${run_root}/figures/34-gene-length-distributions.tiff" \
  "${run_root}/figures/34-strategy-truth-audit.pdf" \
  "${run_root}/figures/34-strategy-truth-audit.png" \
  "${run_root}/figures/34-strategy-truth-audit.tiff" \
  "${run_root}/figures/34-threshold-method-sensitivity.pdf" \
  "${run_root}/figures/34-threshold-method-sensitivity.png" \
  "${run_root}/figures/34-threshold-method-sensitivity.tiff" \
  "${stage}/figures/"
cp -R \
  "${project_root}/data/small/35-gene-abundance-frozen" \
  "${stage}/data/small/"
cp -R \
  "${run_root}/results/35-gene-abundance" \
  "${stage}/results/"
cp \
  "${run_root}/figures/35-read-to-gene-ledger.pdf" \
  "${run_root}/figures/35-read-to-gene-ledger.png" \
  "${run_root}/figures/35-read-to-gene-ledger.tiff" \
  "${run_root}/figures/35-mapping-policy-sensitivity.pdf" \
  "${run_root}/figures/35-mapping-policy-sensitivity.png" \
  "${run_root}/figures/35-mapping-policy-sensitivity.tiff" \
  "${run_root}/figures/35-unit-normalization-audit.pdf" \
  "${run_root}/figures/35-unit-normalization-audit.png" \
  "${run_root}/figures/35-unit-normalization-audit.tiff" \
  "${run_root}/figures/35-functional-aggregation.pdf" \
  "${run_root}/figures/35-functional-aggregation.png" \
  "${run_root}/figures/35-functional-aggregation.tiff" \
  "${stage}/figures/"
cp -R \
  "${project_root}/data/small/36-eggnog-functional-annotation-frozen" \
  "${stage}/data/small/"
cp -R \
  "${run_root}/results/36-eggnog-functional-annotation" \
  "${stage}/results/"
cp \
  "${run_root}/figures/36-field-coverage.pdf" \
  "${run_root}/figures/36-field-coverage.png" \
  "${run_root}/figures/36-field-coverage.tiff" \
  "${run_root}/figures/36-functional-dark-matter.pdf" \
  "${run_root}/figures/36-functional-dark-matter.png" \
  "${run_root}/figures/36-functional-dark-matter.tiff" \
  "${run_root}/figures/36-cog-category-profile.pdf" \
  "${run_root}/figures/36-cog-category-profile.png" \
  "${run_root}/figures/36-cog-category-profile.tiff" \
  "${run_root}/figures/36-go-evidence-policy.pdf" \
  "${run_root}/figures/36-go-evidence-policy.png" \
  "${run_root}/figures/36-go-evidence-policy.tiff" \
  "${stage}/figures/"
for article_dir in \
  37-cazymes-dbcan-frozen \
  38-resistome-card-rgi-frozen \
  39-virulome-vfdb-abricate-frozen \
  40-bgc-natural-products-frozen \
  41-read-mapping-depth-frozen \
  42-binning-comparison-frozen \
  43-bin-refinement-frozen \
  44-mag-qc-mimag-graph-frozen; do
  cp -R "${project_root}/data/small/${article_dir}" "${stage}/data/small/"
done
for result_dir in \
  37-cazymes-dbcan \
  38-resistome-card-rgi \
  39-virulome-vfdb-abricate \
  40-bgc-natural-products \
  41-read-mapping-depth \
  42-binning-comparison \
  43-bin-refinement \
  44-mag-qc-mimag-graph; do
  cp -R "${run_root}/results/${result_dir}" "${stage}/results/"
done
for stem in \
  37-tool-consensus \
  37-cazyme-class-profile \
  37-family-abundance \
  37-cgc-substrate \
  38-evidence-tiers \
  38-primary-hit-quality \
  38-drug-class-profile \
  38-positive-controls \
  39-threshold-database-sensitivity \
  39-hit-quality \
  39-vfc-category-profile \
  39-context-controls \
  40-tool-bgc-yield \
  40-fragmentation-sensitivity \
  40-bgc-type-profile \
  40-mibig-similarity-abundance \
  41-mapping-fate \
  41-depth-concordance \
  41-depth-breadth \
  42-binner-quality-yield \
  42-recovery-purity \
  42-single-vs-multisample \
  42-taxonomy-coverage \
  43-refinement-yield \
  43-quality-landscape \
  43-refinement-provenance \
  43-method-selection \
  44-quality-landscape \
  44-mimag-requirements \
  44-assembly-graph-audit \
  44-checkm-audit; do
  cp \
    "${run_root}/figures/${stem}.pdf" \
    "${run_root}/figures/${stem}.png" \
    "${run_root}/figures/${stem}.tiff" \
    "${stage}/figures/"
done

(
  cd "${stage}"
  export R_LIBS_USER="${project_root}/.r-lib:${HOME}/R/library${R_LIBS_USER:+:${R_LIBS_USER}}"
  XDG_CACHE_HOME="${run_root}/.cache/quarto" "${quarto_bin}" render
)

test -f "${stage}/_site/index.html"
test -f "${stage}/_site/chapters/02-three-analysis-layers.html"
test -f "${stage}/_site/chapters/03-study-design-power.html"
test -f "${stage}/_site/chapters/04-sequencing-depth.html"
test -f "${stage}/_site/chapters/05-library-prep-absolute-quantification.html"
test -f "${stage}/_site/chapters/06-host-depletion-low-biomass.html"
test -f "${stage}/_site/chapters/07-contamination-controls.html"
test -f "${stage}/_site/chapters/08-fastq-short-long-reads.html"
test -f "${stage}/_site/chapters/09-wsl2-conda.html"
test -f "${stage}/_site/chapters/10-computing-hpc-cloud.html"
test -f "${stage}/_site/chapters/11-install-biobakery-assembly.html"
test -f "${stage}/_site/chapters/12-install-r-cmd.html"
test -f "${stage}/_site/chapters/13-read-qc-fastp.html"
test -f "${stage}/_site/chapters/14-host-removal-complexity-duplicates.html"
test -f "${stage}/_site/chapters/15-metaphlan4.html"
test -f "${stage}/_site/chapters/16-kraken2-bracken.html"
test -f "${stage}/_site/chapters/17-kraken2-database-confidence.html"
test -f "${stage}/_site/chapters/18-profiler-benchmark.html"
test -f "${stage}/_site/chapters/19-humann3.html"
test -f "${stage}/_site/chapters/20-functional-profile-normalization.html"
test -f "${stage}/_site/chapters/21-metagenomic-table-semantics.html"
test -f "${stage}/_site/chapters/22-alpha-beta-diversity.html"
test -f "${stage}/figures/22-korchagina-fig5-original.png"
test -f "${stage}/_site/chapters/23-pcoa-cap-permanova.html"
test -f "${stage}/figures/23-korchagina-fig1-original.png"
test -f "${stage}/_site/chapters/24-differential-abundance.html"
test -f "${stage}/figures/24-zeller-fig1-original.png"
test -f "${stage}/_site/chapters/25-composition-core-microbiome.html"
test -f "${stage}/figures/25-hmp-fig3-original.jpg"
test -f "${stage}/_site/chapters/26-curated-metagenomic-data-lineage.html"
test -f "${stage}/figures/26-cmd3-fig1-original.png"
test -f "${stage}/_site/chapters/27-machine-learning-roc.html"
test -f "${stage}/figures/27-zeller-fig1-original.png"
test -f "${stage}/_site/chapters/28-cross-cohort-validation.html"
test -f "${stage}/figures/28-wirbel-fig3-original.png"
test -f "${stage}/_site/chapters/29-cooccurrence-network.html"
test -f "${stage}/figures/29-kurtz-fig2-original.png"
test -f "${stage}/_site/chapters/30-short-read-assembly.html"
test -f "${stage}/_site/chapters/31-long-read-assembly.html"
test -f "${stage}/_site/chapters/32-hybrid-assembly-polishing.html"
test -f "${stage}/_site/chapters/33-assembly-qc.html"
test -f "${stage}/_site/chapters/34-nonredundant-gene-catalog.html"
test -f "${stage}/_site/chapters/35-gene-abundance.html"
test -f "${stage}/_site/chapters/36-eggnog-functional-annotation.html"
test -f "${stage}/_site/chapters/37-cazymes-dbcan.html"
test -f "${stage}/_site/chapters/38-resistome.html"
test -f "${stage}/_site/chapters/39-virulome.html"
test -f "${stage}/_site/chapters/40-bgc-natural-products.html"
test -f "${stage}/_site/chapters/41-read-mapping-depth.html"
test -f "${stage}/_site/chapters/42-binning.html"
test -f "${stage}/_site/chapters/43-bin-refinement.html"
test -f "${stage}/_site/chapters/44-mag-qc-checkm2-gunc-mimag.html"
for stem in \
  01-assay-boundaries \
  01-crc-cohorts \
  02-layer-decision \
  02-workflow-anchors \
  03-cohort-balance \
  03-covariate-completeness \
  03-power-sensitivity \
  04-crc-library-depth \
  04-nonpareil-saturation \
  04-endpoint-depth-anchors \
  05-extraction-bias \
  05-protocol-bias-range \
  05-syndna-quantification \
  06-host-depletion-efficiency \
  06-saponin-tradeoff \
  06-host-depletion-decision \
  07-control-library-size \
  07-contaminant-prevalence \
  07-contaminant-burden \
  07-index-hopping-evidence \
  08-read-geometry \
  08-span-survival \
  08-assembly-impact \
  09-wsl2-layer-map \
  09-environment-validation \
  10-input-resource-budget \
  10-resource-control-loop \
  10-restart-safe-array \
  11-environment-boundaries \
  11-toolchain-entrypoints \
  11-database-storage-contract \
  12-r-data-access-boundaries \
  12-package-role-contract \
  12-cmd-object-contract \
  13-per-cycle-quality \
  13-read-pair-fate \
  13-fastqc-module-states \
  14-host-removal-recall-retention \
  14-complexity-sensitivity \
  14-duplicate-decision \
  15-metaphlan-composition \
  15-sgb-marker-support \
  15-detection-quantification-sensitivity \
  16-kraken-classification-ledger \
  16-bracken-redistribution \
  16-bracken-parameter-sensitivity \
  17-reference-coverage \
  17-confidence-tradeoff \
  17-database-stability \
  17-hit-group-control \
  18-feature-space-crosswalk \
  18-recovery-nontruth \
  18-composition-agreement \
  18-resource-footprint \
  19-read-flow \
  19-gene-family-stratification \
  19-pathway-contributions \
  19-abundance-coverage \
  20-normalization-denominators \
  20-special-feature-budget \
  20-pathway-contributions \
  20-prevalence-zero-sensitivity \
  21-table-unit-map \
  21-denominator-closure \
  21-genome-equivalent-calibration \
  21-zero-strata-semantics \
  22-resolution-boundaries \
  22-alpha-hill-numbers \
  22-beta-ordination \
  22-sensitivity-recruitment \
  23-design-permutation-space \
  23-pcoa-cap \
  23-permanova-dispersion \
  23-pairwise-sensitivity \
  24-zeller-marker-redraw \
  24-maaslin3-two-part \
  24-functional-associations \
  24-method-concordance \
  25-multirank-mean-composition \
  25-individual-stool-composition \
  25-prevalence-abundance \
  25-core-membership-sensitivity \
  26-resource-release-audit \
  26-metadata-completeness \
  26-query-attrition \
  26-lineage-compatibility \
  27-nested-roc \
  27-model-performance \
  27-permutation-importance \
  27-leakage-permutation-audit \
  28-meta-signature \
  28-single-study-transfer \
  28-lodo-forest \
  28-lodo-roc \
  28-validation-stability \
  29-feature-contract \
  29-conditional-network \
  29-zippi-roles \
  29-edge-stability-sensitivity \
  29-robustness-null \
  30-assembly-branch-design \
  30-contiguity-output \
  30-recruitment-tradeoff \
  30-resource-footprint \
  31-long-read-branch-design \
  31-long-read-contiguity \
  31-long-read-readback \
  31-circular-resource-audit \
  32-hybrid-branch-design \
  32-recovery-contiguity \
  32-consensus-error \
  32-abundance-resource-audit \
  33-qc-evidence-ladder \
  33-n50-na50-correctness \
  33-recovery-error-tradeoff \
  33-diagnostic-task-gates \
  34-gene-catalog-workflow \
  34-gene-length-distributions \
  34-strategy-truth-audit \
  34-threshold-method-sensitivity \
  35-read-to-gene-ledger \
  35-mapping-policy-sensitivity \
  35-unit-normalization-audit \
  35-functional-aggregation \
  36-field-coverage \
  36-functional-dark-matter \
  36-cog-category-profile \
  36-go-evidence-policy \
  37-tool-consensus \
  37-cazyme-class-profile \
  37-family-abundance \
  37-cgc-substrate \
  38-evidence-tiers \
  38-primary-hit-quality \
  38-drug-class-profile \
  38-positive-controls \
  39-threshold-database-sensitivity \
  39-hit-quality \
  39-vfc-category-profile \
  39-context-controls \
  40-tool-bgc-yield \
  40-fragmentation-sensitivity \
  40-bgc-type-profile \
  40-mibig-similarity-abundance \
  41-mapping-fate \
  41-depth-concordance \
  41-depth-breadth \
  42-binner-quality-yield \
  42-recovery-purity \
  42-single-vs-multisample \
  42-taxonomy-coverage \
  43-refinement-yield \
  43-quality-landscape \
  43-refinement-provenance \
  43-method-selection \
  44-quality-landscape \
  44-mimag-requirements \
  44-assembly-graph-audit \
  44-checkm-audit; do
  test -f "${stage}/figures/${stem}.pdf"
  test -f "${stage}/figures/${stem}.png"
  test -f "${stage}/figures/${stem}.tiff"
done
