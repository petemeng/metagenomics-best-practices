#!/usr/bin/env python3
"""Validate the 77-article contract and all currently executable articles."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any

import yaml


SHARED_ROOT = Path(__file__).resolve().parents[2]
if str(SHARED_ROOT) not in sys.path:
    sys.path.insert(0, str(SHARED_ROOT))

from tutorial_automation.manifest import load_manifest


EXPECTED_COUNT = 77
EXPECTED_WECHAT_TITLE_PREFIX = "宏基因组最佳实践"
EXPECTED_WECHAT_TITLE_MAX_CHARS = 64
EXECUTABLE_TOKEN_NUMBERS = set(range(1, 45))
UPSTREAM_EVAL_FALSE = {
    6, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20,
    30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44,
}
REQUIRED_SECTIONS = (
    "这一步对应论文里的哪张图",
    "理论",
    "准备工作",
    "可复制代码",
    "审计与升级",
    "出版级美化",
    "常见坑",
    "这段 Methods 怎么写",
    "换成你自己的数据怎么做",
    "参考",
)
REQUIRED_ARTICLE_TOKENS = {
    1: (
        "data/small/01-assay-evidence.tsv",
        "data/small/01-crc-cohort-summary.tsv",
        "set.seed(",
        "01-assay-boundaries",
        "01-crc-cohorts",
        "save_pub(",
    ),
    2: (
        "data/small/02-layer-capabilities.tsv",
        "data/small/02-published-anchors.tsv",
        "set.seed(",
        "02-layer-decision",
        "02-workflow-anchors",
        "save_pub(",
    ),
    3: (
        "data/small/03-crc-design-audit.tsv",
        "set.seed(",
        "03-cohort-balance",
        "03-covariate-completeness",
        "03-power-sensitivity",
        "save_pub(",
    ),
    4: (
        "data/small/04-crc-library-size.tsv",
        "data/small/04-lake-lanier-coverage.tsv",
        "data/small/04-depth-evidence.tsv",
        "set.seed(",
        "04-crc-library-depth",
        "04-nonpareil-saturation",
        "04-endpoint-depth-anchors",
        "save_pub(",
    ),
    5: (
        "data/small/05-costea-mock-profiles.tsv",
        "data/small/05-costea-bias-summary.tsv",
        "data/small/05-syndna-mock-benchmark.tsv",
        "data/small/05-control-placement.tsv",
        "set.seed(",
        "05-extraction-bias",
        "05-protocol-bias-range",
        "05-syndna-quantification",
        "save_pub(",
    ),
    6: (
        "data/small/06-marotz-host-depletion-source.tsv",
        "data/small/06-longhi-saponin-source.tsv",
        "data/small/06-host-filter-contract.tsv",
        "data/small/06-host-depletion-read-budget.tsv",
        "data/small/06-saponin-tradeoff.tsv",
        "set.seed(",
        "Hostile 2.0.2",
        "human-t2t-hla",
        "06-host-depletion-efficiency",
        "06-saponin-tradeoff",
        "06-host-depletion-decision",
        "save_pub(",
    ),
    7: (
        'data_dir <- "data/small/07-decontam"',
        '"otutab.tsv"',
        '"taxonomy.tsv"',
        '"metadata.tsv"',
        "data/small/07-salter-shotgun-evidence.tsv",
        "data/small/07-index-hopping-evidence.tsv",
        "data/small/07-data-NOTICE.txt",
        "set.seed(",
        "decontam 1.24.0",
        "07-control-library-size",
        "07-contaminant-prevalence",
        "07-contaminant-burden",
        "07-index-hopping-evidence",
        "save_pub(",
    ),
    8: (
        "data/small/08-ena-fastq-sources.tsv",
        "data/small/08-platform-benchmark.tsv",
        "data/small/08-native-format-contract.tsv",
        "data/small/08-read-prefix-metrics.tsv",
        "data/small/08-fastq-anatomy.tsv",
        "data/small/08-prefix-source-summary.json",
        "data/small/08-data-NOTICE.txt",
        "set.seed(",
        "ERR9765746",
        "ERR9765780",
        "ERR9765783",
        "08-read-geometry",
        "08-span-survival",
        "08-assembly-impact",
        "save_pub(",
    ),
    9: (
        "env/platform-smoke.yml",
        "data/small/09-environment-contract.tsv",
        "data/small/09-data-NOTICE.txt",
        "data/small/08-ena-fastq-sources.tsv",
        "data/small/08-prefix-source-summary.json",
        "data/small/08-read-prefix-metrics.tsv",
        "scripts/validate_article09_environment.py",
        "metagenome-platform-smoke-2026.07",
        "Miniforge 26.3.2-2",
        "8f758b6ffdcf1561ece7d187ff34bc3f5a174fd8c6da66a101b206fcc869d20c",
        "wsl.exe --list --verbose",
        "Python 3.12.13",
        "seqkit v2.10.0",
        "pigz 2.8",
        "Matplotlib 3.10.5",
        "random.seed(20260719)",
        "09-wsl2-layer-map",
        "09-environment-validation",
        "save_pub(",
    ),
    10: (
        "data/small/10-job-array.tsv",
        "data/small/10-runtime-contract.tsv",
        "data/small/10-container-smoke.log",
        "data/small/10-data-NOTICE.txt",
        "scripts/article10_task.py",
        "scripts/10_slurm_array_smoke.sh",
        "scripts/validate_article10_compute.py",
        "10,944,966,814",
        "Apptainer 1.5.2",
        "sha256:eafc1edb577d2e9b458664a15f23ea1c370214193226069eb22921169fc7e43f",
        "random.seed(20260719)",
        "10-input-resource-budget",
        "10-resource-control-loop",
        "10-restart-safe-array",
        "save_pub(",
    ),
    11: (
        "env/biobakery.yml",
        "env/biobakery-linux-64.lock",
        "env/relink-biobakery-entrypoints.sh",
        "env/assembly.yml",
        "env/assembly-linux-64.lock",
        "data/small/11-environment-evidence.tsv",
        "data/small/11-install-self-tests.log",
        "data/small/11-solver-audit.tsv",
        "data/small/11-database-manifest.tsv",
        "db/download_db.sh",
        "scripts/validate_article11_installation.py",
        "mpa_vJan26_CHOCOPhlAnSGB_202605",
        "7162b0c3493663dce9abef08ccc06aea",
        "R232",
        "25a59e0352b1fd150c589f56559767d4",
        "random.seed(20260719)",
        "11-environment-boundaries",
        "11-toolchain-entrypoints",
        "11-database-storage-contract",
        "save_pub <- function",
    ),
    12: (
        "env/renv.lock",
        "data/small/12-package-contract.tsv",
        "data/small/12-cmd-resource-manifest.tsv",
        "data/small/12-cmd-asnicarf-2017-relative-abundance.rds",
        "data/small/12-resource-retrieval.log",
        "data/small/12-data-NOTICE.txt",
        "scripts/validate_article12_r_cmd.R",
        "R 4.4.1",
        "Bioconductor 3.19",
        "curatedMetagenomicData 3.12.0",
        "AsnicarF_2017.relative_abundance",
        "EH7091",
        "MetaPhlAn 3",
        "HUMAnN 3",
        "2952774730bff2af9e13c9c40058320aed524dc2dc7408b44dc6e697c06564b2",
        "set.seed(20260720)",
        "12-r-data-access-boundaries",
        "12-package-role-contract",
        "12-cmd-object-contract",
        "save_pub <- function",
    ),
    13: (
        "env/read-qc.yml",
        "env/read-qc-linux-64.lock",
        "data/small/13-source-manifest.tsv",
        "data/small/13-qc-frozen/run-summary.json",
        "data/small/13-qc-frozen/file-checksums.sha256",
        "data/small/13-data-NOTICE.txt",
        "scripts/build_article13_fastq_subset.py",
        "scripts/run_article13_read_qc.sh",
        "scripts/validate_article13_read_qc.py",
        "FastQC 0.12.1",
        "fastp 1.3.6",
        "MultiQC 1.35",
        "ERR9765746",
        "99.991%",
        "set.seed(20260720)",
        "13-per-cycle-quality",
        "13-read-pair-fate",
        "13-fastqc-module-states",
        "save_pub <- function",
    ),
    14: (
        "env/host-removal.yml",
        "env/host-removal-linux-64.lock",
        "data/small/14-source-manifest.tsv",
        "data/small/14-index-manifest.tsv",
        "data/small/14-host-removal-frozen/run-summary.json",
        "data/small/14-host-removal-frozen/file-checksums.sha256",
        "data/small/14-data-NOTICE.txt",
        "scripts/build_article14_fastq_controls.py",
        "scripts/run_article14_host_removal.sh",
        "scripts/validate_article14_host_removal.py",
        "Hostile 2.0.2",
        "Bowtie2 2.5.4",
        "ERR194147",
        "ERR9765746",
        "99.750%",
        "100.000%",
        "set.seed(20260720)",
        "14-host-removal-recall-retention",
        "14-complexity-sensitivity",
        "14-duplicate-decision",
        "save_pub <- function",
    ),
    15: (
        "env/biobakery.yml",
        "env/biobakery-linux-64.lock",
        "data/small/15-source-manifest.tsv",
        "data/small/15-database-manifest.tsv",
        "data/small/15-metaphlan-frozen/run-summary.json",
        "data/small/15-metaphlan-frozen/file-checksums.sha256",
        "data/small/15-data-NOTICE.txt",
        "db/download_db.sh",
        "scripts/run_article15_metaphlan4.sh",
        "scripts/validate_article15_metaphlan4.py",
        "MetaPhlAn 4.2.5",
        "Bowtie2 2.5.5",
        "mpa_vJan26_CHOCOPhlAnSGB_202605",
        "ViralMarkerNamesInMapout",
        "MetadataExcludedSGBMarkerNamesInMapout",
        "199,929",
        "set.seed(20260721)",
        "15-metaphlan-composition",
        "15-sgb-marker-support",
        "15-detection-quantification-sensitivity",
        "save_pub <- function",
    ),
    16: (
        "env/kraken.yml",
        "env/kraken-linux-64.lock",
        "data/small/16-source-manifest.tsv",
        "data/small/16-database-manifest.tsv",
        "data/small/16-kraken-bracken-frozen/run-summary.json",
        "data/small/16-kraken-bracken-frozen/file-checksums.sha256",
        "data/small/16-data-NOTICE.txt",
        "db/download_db.sh",
        "scripts/run_article16_kraken2_bracken.sh",
        "scripts/validate_article16_kraken2_bracken.py",
        "Kraken2 2.17.1",
        "Bracken package 3.1p1",
        "kraken2-standard8-20260626",
        "86,373",
        "set.seed(20260721)",
        "16-kraken-classification-ledger",
        "16-bracken-redistribution",
        "16-bracken-parameter-sensitivity",
        "save_pub <- function",
    ),
    17: (
        "env/kraken.yml",
        "env/kraken-linux-64.lock",
        "data/small/17-source-manifest.tsv",
        "data/small/17-database-manifest.tsv",
        "data/small/17-mock1-truth.tsv",
        "data/small/17-kraken-database-confidence-frozen/run-summary.json",
        "data/small/17-kraken-database-confidence-frozen/file-checksums.sha256",
        "data/small/17-data-NOTICE.txt",
        "db/download_db.sh",
        "scripts/run_article17_kraken2_database_confidence.sh",
        "scripts/validate_article17_kraken2_database_confidence.py",
        "Kraken2 2.17.1",
        "Bracken package 3.1p1",
        "Standard-8",
        "Standard-16",
        "PlusPF-8",
        "55/63",
        "15,585 MiB",
        "set.seed(20260721)",
        "17-reference-coverage",
        "17-confidence-tradeoff",
        "17-database-stability",
        "17-hit-group-control",
        "save_pub <- function",
    ),
    18: (
        "env/motus.yml",
        "env/motus-linux-64.lock",
        "data/small/18-source-manifest.tsv",
        "data/small/18-database-manifest.tsv",
        "data/small/18-profiler-benchmark-frozen/run-summary.json",
        "data/small/18-profiler-benchmark-frozen/truth-feature-crosswalk.tsv",
        "data/small/18-profiler-benchmark-frozen/file-checksums.sha256",
        "data/small/18-data-NOTICE.txt",
        "db/download_db.sh",
        "scripts/run_article18_profiler_benchmark.sh",
        "scripts/validate_article18_profiler_benchmark.py",
        "MetaPhlAn 4.2.5",
        "Kraken2 2.17.1",
        "mOTUs 4.1.0",
        "GTDB R226",
        "52 个物种",
        "70.738%",
        "18-feature-space-crosswalk",
        "18-recovery-nontruth",
        "18-composition-agreement",
        "18-resource-footprint",
        "save_pub <- function",
    ),
    19: (
        "env/biobakery.yml",
        "env/biobakery-linux-64.lock",
        "data/small/19-source-manifest.tsv",
        "data/small/19-database-manifest.tsv",
        "data/small/19-humann3-frozen/run-summary.json",
        "data/small/19-humann3-frozen/database-audit.tsv",
        "data/small/19-humann3-frozen/regroup-audit.tsv",
        "data/small/19-humann3-frozen/reactions-rpk.tsv",
        "data/small/19-humann3-frozen/file-checksums.sha256",
        "data/small/19-data-NOTICE.txt",
        "db/download_db.sh",
        "scripts/bootstrap_article19_humann_databases.sh",
        "scripts/download_ranged_archive.sh",
        "scripts/run_article19_humann3.sh",
        "scripts/validate_article19_humann3.py",
        "HUMAnN 3.9",
        "MetaPhlAn 4.2.5",
        "mpa_vJun23_CHOCOPhlAnSGB_202403",
        "v201901_v31",
        "v201901b",
        "humann_regroup_table",
        "uniref90_rxn",
        "1,402",
        "1.935",
        "199,982",
        "--seed 20260722",
        "19-read-flow",
        "19-gene-family-stratification",
        "19-pathway-contributions",
        "19-abundance-coverage",
        "save_pub <- function",
    ),
    20: (
        "env/biobakery-linux-64.lock",
        "env/renv.lock",
        "data/small/19-humann3-frozen/genefamilies-rpk.tsv",
        "data/small/19-humann3-frozen/reactions-rpk.tsv",
        "data/small/19-humann3-frozen/pathabundance-rpk.tsv",
        "data/small/19-humann3-frozen/pathcoverage.tsv",
        "data/small/20-cmd-pathway/pathway-abundance.tsv.gz",
        "data/small/20-cmd-pathway/pathway-coverage.tsv.gz",
        "data/small/20-cmd-pathway/sample-metadata.tsv",
        "data/small/20-functional-normalization-frozen/file-checksums.sha256",
        "scripts/prepare_article20_cmd_pathways.R",
        "scripts/run_article20_humann_normalization.sh",
        "scripts/validate_article20_functional_normalization.py",
        "HUMAnN 3.9",
        "curatedMetagenomicData 3.12.0",
        "community",
        "levelwise",
        "UNMAPPED",
        "UNINTEGRATED",
        "UNGROUPED",
        "set.seed(20260722)",
        "20-normalization-denominators",
        "20-special-feature-budget",
        "20-pathway-contributions",
        "20-prevalence-zero-sensitivity",
        "save_pub <- function",
    ),
    21: (
        "env/biobakery-linux-64.lock",
        "env/renv.lock",
        "data/small/13-qc-frozen",
        "data/small/15-metaphlan-frozen/profile-all.tsv",
        "data/small/16-kraken-bracken-frozen/bracken-species-r150-t10.tsv",
        "data/small/19-humann3-frozen/genefamilies-rpk.tsv",
        "data/small/19-humann3-frozen/pathcoverage.tsv",
        "data/small/20-cmd-pathway/pathway-abundance.tsv.gz",
        "data/small/20-cmd-pathway/pathway-coverage.tsv.gz",
        "data/small/21-table-semantics-frozen/microbecensus-read-length-sensitivity.tsv",
        "scripts/run_article21_microbecensus.sh",
        "dfc42d356bfd7943633cde6c0fbfc0b116f29ae2",
        "199,982",
        "29,809,773",
        "11.5711948566",
        "set.seed(20260722)",
        "21-table-unit-map",
        "21-denominator-closure",
        "21-genome-equivalent-calibration",
        "21-zero-strata-semantics",
        "save_pub <- function",
    ),
    22: (
        "env/renv.lock",
        "data/small/22-diversity-inputs",
        "species-relative-abundance.tsv.gz",
        "gene-family-prevalence10.tsv.gz",
        "mag-relative-abundance.tsv.gz",
        "hot-spring-sample-metadata.tsv",
        "mag-recruitment.tsv",
        "scripts/prepare_article22_diversity_data.R",
        "EH7091",
        "EH7086",
        "30284068 v2",
        "415,581",
        "178,928",
        "12.98%",
        "set.seed(20260722)",
        "22-resolution-boundaries",
        "22-alpha-hill-numbers",
        "22-beta-ordination",
        "22-sensitivity-recruitment",
        "save_pub <- function",
    ),
    23: (
        "env/renv.lock",
        "data/small/23-ordination-permanova",
        "spring-mag-relative-abundance.tsv.gz",
        "spring-mag-recruitment-weighted.tsv.gz",
        "spring-metadata.tsv",
        "analysis-contract.tsv",
        "scripts/prepare_article23_ordination_data.R",
        "scripts/validate_article23_ordination_permanova.R",
        "30284068 v2",
        "set.seed(20260723)",
        "adonis2",
        "betadisper",
        "23-design-permutation-space",
        "23-pcoa-cap",
        "23-permanova-dispersion",
        "23-pairwise-sensitivity",
        "save_pub <- function",
    ),
    24: (
        "env/differential-abundance-linux-64.lock",
        "env/install-differential-abundance-r.R",
        "data/small/24-differential-abundance",
        "species-relative-abundance.tsv.gz",
        "species-pseudocounts.tsv.gz",
        "pathway-relative-abundance.tsv.gz",
        "scripts/prepare_article24_differential_data.R",
        "scripts/validate_article24_differential_abundance.R",
        "curatedMetagenomicData 3.12.0",
        "MaAsLin3 1.5.3",
        "ANCOM-BC2 2.12.0",
        "ALDEx2 1.42.0",
        "set.seed(20260724)",
        "24-zeller-marker-redraw",
        "24-maaslin3-two-part",
        "24-functional-associations",
        "24-method-concordance",
        "save_pub <- function",
    ),
    25: (
        "data/small/25-composition-core",
        "phylum-relative-abundance.tsv.gz",
        "genus-relative-abundance.tsv.gz",
        "species-relative-abundance.tsv.gz",
        "scripts/prepare_article25_composition_core.R",
        "scripts/validate_article25_composition_core.R",
        "curatedMetagenomicData 3.12.0",
        "set.seed(20260725)",
        "25-multirank-mean-composition",
        "25-individual-stool-composition",
        "25-prevalence-abundance",
        "25-core-membership-sensitivity",
        "save_pub <- function",
    ),
    26: (
        "data/small/26-cmd-lineage",
        "resource-catalog.tsv",
        "sample-catalog.tsv.gz",
        "sample-id-collisions.tsv",
        "merged-species-relative-abundance.tsv.gz",
        "lineage-cards.tsv",
        "scripts/prepare_article26_cmd_lineage.R",
        "scripts/validate_article26_cmd_lineage.R",
        "curatedMetagenomicData 3.12.0",
        "MetaPhlAn 3",
        "HUMAnN 3",
        "CHOCOPhlAn 201901",
        "22,588",
        "22,710",
        "set.seed(20260725)",
        "26-resource-release-audit",
        "26-metadata-completeness",
        "26-query-attrition",
        "26-lineage-compatibility",
        "save_pub <- function",
    ),
    27: (
        "env/machine-learning-renv.lock",
        "data/small/27-machine-learning",
        "species-relative-abundance.tsv.gz",
        "analysis-contract.tsv",
        "scripts/prepare_article27_machine_learning.R",
        "scripts/validate_article27_machine_learning.R",
        "curatedMetagenomicData 3.12.0",
        "MetaPhlAn 3",
        "CHOCOPhlAn 201901",
        "nested cross-validation",
        "Random forest",
        "XGBoost",
        "AUROC",
        "AUPRC",
        "Brier",
        "permutation importance",
        "data leakage",
        "set.seed(20260727)",
        "27-nested-roc",
        "27-model-performance",
        "27-permutation-importance",
        "27-leakage-permutation-audit",
        "save_pub <- function",
    ),
    28: (
        "env/cross-cohort-renv.lock",
        "data/small/28-cross-cohort",
        "analysis-contract.tsv",
        "scripts/prepare_article28_cross_cohort.R",
        "scripts/validate_article28_cross_cohort.R",
        "curatedMetagenomicData 3.12.0",
        "MetaPhlAn 3",
        "CHOCOPhlAn 201901",
        "leave-one-training-cohort-out",
        "Random-effects",
        "Hedges",
        "hierarchical bootstrap",
        "ComBat",
        "set.seed(20260728)",
        "28-meta-signature",
        "28-single-study-transfer",
        "28-lodo-forest",
        "28-lodo-roc",
        "28-validation-stability",
        "save_pub <- function",
    ),
    29: (
        "env/network-renv.lock",
        "data/small/29-network",
        "analysis-contract.tsv",
        "scripts/prepare_article29_network_data.R",
        "scripts/validate_article29_network.R",
        "Figshare 30284068 v2",
        "huge 1.5.1",
        "graphical lasso",
        "StARS",
        "BroadRegion",
        "Zi",
        "Pi",
        "set.seed(20260729)",
        "29-feature-contract",
        "29-conditional-network",
        "29-zippi-roles",
        "29-edge-stability-sensitivity",
        "29-robustness-null",
        "save_pub <- function",
    ),
    30: (
        "data/small/30-short-read-assembly-frozen",
        "scripts/run_article30_short_read_assembly.sh",
        "scripts/validate_article30_short_read_assembly.py",
        "PRJEB52977",
        "ERR9765746",
        "ERR9765747",
        "MEGAHIT 1.2.9",
        "metaSPAdes 4.3.0",
        "Bowtie2 2.5.5",
        "20260730",
        "不是算法单因素",
        "不代表普适赢家",
        "30-assembly-branch-design",
        "30-contiguity-output",
        "30-recruitment-tradeoff",
        "30-resource-footprint",
        "save_pub <- function",
    ),
    31: (
        "data/small/31-long-read-assembly-frozen",
        "scripts/run_article31_long_read_assembly.sh",
        "scripts/validate_article31_long_read_assembly.py",
        "PRJEB52977",
        "ERR9765780",
        "ERR9765783",
        "Flye 2.9.6",
        "hifiasm-meta 0.3.5",
        "metaMDBG 1.4",
        "minimap2 2.31",
        "20260731",
        "不是平台单因素",
        "不代表普适赢家",
        "31-long-read-branch-design",
        "31-long-read-contiguity",
        "31-long-read-readback",
        "31-circular-resource-audit",
        "save_pub <- function",
    ),
    32: (
        "data/small/32-hybrid-assembly-polishing-frozen",
        "data/small/32-branch-contract.tsv",
        "scripts/run_article32_hybrid_assembly.sh",
        "scripts/validate_article32_hybrid_assembly.py",
        "PRJEB52977",
        "ERR9765746",
        "ERR9765780",
        "ERR9765783",
        "SPAdes 4.3.0",
        "Polypolish 0.6.1",
        "BWA-MEM v0.7.19",
        "MetaQUAST 5.3.0",
        "minimap2 2.28-r1209",
        "20260732",
        "不是平台单因素实验",
        "不能写“Polypolish 保证提升准确性”",
        "32-hybrid-branch-design",
        "32-recovery-contiguity",
        "32-consensus-error",
        "32-abundance-resource-audit",
        "save_pub <- function",
    ),
    33: (
        "data/small/33-assembly-qc-frozen",
        "data/small/33-source-manifest.tsv",
        "scripts/prepare_article33_qc_inputs.py",
        "scripts/run_article33_assembly_qc.sh",
        "scripts/summarize_article33_assembly_qc.py",
        "scripts/validate_article33_assembly_qc.py",
        "PRJEB52977",
        "QUAST 5.3.0",
        "MetaQUAST 5.3.0",
        "minimap2 2.28-r1209",
        "20260733",
        "1,271",
        "不构成平台或算法的单因素因果比较",
        "no universal cutoff",
        "33-qc-evidence-ladder",
        "33-n50-na50-correctness",
        "33-recovery-error-tradeoff",
        "33-diagnostic-task-gates",
        "save_pub <- function",
    ),
    34: (
        "data/small/34-nonredundant-gene-catalog-frozen",
        "env/gene-catalog-linux-64.lock",
        "scripts/prepare_article34_gene_catalog_inputs.py",
        "scripts/run_article34_gene_catalog.py",
        "scripts/summarize_article34_gene_catalog.py",
        "scripts/validate_article34_gene_catalog.py",
        "PRJEB52977",
        "Prodigal 2.6.3",
        "MMseqs2 9.d36de",
        "CD-HIT 4.8.1",
        "20260734",
        "441,407",
        "270,679",
        "260,868",
        "93,782",
        "216,191",
        "不能据此宣布某种 assembly strategy 在自然群落中普遍最优",
        "不是多数表决序列",
        "34-gene-catalog-workflow",
        "34-gene-length-distributions",
        "34-strategy-truth-audit",
        "34-threshold-method-sensitivity",
        "save_pub <- function",
    ),
    35: (
        "data/small/35-gene-abundance-frozen",
        "env/gene-abundance-linux-64.lock",
        "scripts/download_article35_gene_abundance_sources.sh",
        "scripts/prepare_article35_gene_abundance_inputs.py",
        "scripts/parse_article35_sam.py",
        "scripts/run_article35_gene_abundance.py",
        "scripts/summarize_article35_gene_abundance.py",
        "scripts/freeze_article35_gene_abundance.py",
        "scripts/validate_article35_gene_abundance.py",
        "PRJEB52977",
        "Bowtie2 2.5.5",
        "SAMtools 1.23.1",
        "HTSeq 2.1.2",
        "DIAMOND 2.2.4",
        "UniRef90 v201901b",
        "20260735",
        "2,784,234",
        "2,777,443",
        "NO_UNIREF90_HIT",
        "UNIREF90_NO_REACTION",
        "不是生物学重复",
        "不是 absolute abundance",
        "35-read-to-gene-ledger",
        "35-mapping-policy-sensitivity",
        "35-unit-normalization-audit",
        "35-functional-aggregation",
        "save_pub <- function",
    ),
    36: (
        "data/small/36-eggnog-functional-annotation-frozen",
        "env/eggnog-annotation-linux-64.lock",
        "scripts/bootstrap_article36_eggnog_database.sh",
        "scripts/prepare_article36_eggnog_inputs.py",
        "scripts/run_article36_eggnog_annotation.py",
        "scripts/summarize_article36_eggnog_annotation.py",
        "scripts/freeze_article36_eggnog_annotation.py",
        "scripts/validate_article36_eggnog_annotation.py",
        "eggNOG-mapper 2.1.15",
        "eggNOG 5.0.2",
        "DIAMOND 2.0.15",
        "20260736",
        "non-electronic",
        "all-evidence",
        "No seed ortholog",
        "fractional allocation",
        "36-field-coverage",
        "36-functional-dark-matter",
        "36-cog-category-profile",
        "36-go-evidence-policy",
        "save_pub <- function",
    ),
    37: (
        "data/small/37-cazymes-dbcan-frozen",
        "env/cazyme-linux-64.lock",
        "scripts/prepare_article37_cazymes.py",
        "scripts/run_article37_cazymes.py",
        "scripts/summarize_article37_cazymes.py",
        "scripts/freeze_article37_cazymes.py",
        "scripts/validate_article37_cazymes.py",
        "dbCAN 5.2.9",
        "db_v5-2-9_5-5-2026",
        "20260737",
        "93,782",
        "2,050",
        "37-tool-consensus",
        "37-cazyme-class-profile",
        "37-family-abundance",
        "37-cgc-substrate",
        "save_pub <- function",
    ),
    38: (
        "data/small/38-resistome-card-rgi-frozen",
        "env/resistome-linux-64.lock",
        "scripts/prepare_article38_resistome.py",
        "scripts/run_article38_resistome.py",
        "scripts/summarize_article38_resistome.py",
        "scripts/freeze_article38_resistome.py",
        "scripts/validate_article38_resistome.py",
        "RGI 6.0.8",
        "CARD 4.0.1",
        "20260738",
        "Perfect",
        "Strict",
        "Loose",
        "38-evidence-tiers",
        "38-primary-hit-quality",
        "38-drug-class-profile",
        "38-positive-controls",
        "save_pub <- function",
    ),
    39: (
        "data/small/39-virulome-vfdb-abricate-frozen",
        "env/virulome-linux-64.lock",
        "scripts/prepare_article39_virulome.py",
        "scripts/run_article39_virulome.py",
        "scripts/summarize_article39_virulome.py",
        "scripts/freeze_article39_virulome.py",
        "scripts/validate_article39_virulome.py",
        "ABRicate 1.4.0",
        "VFDB 2026-07-24",
        "20260739",
        "core set A",
        "full set B",
        "39-threshold-database-sensitivity",
        "39-hit-quality",
        "39-vfc-category-profile",
        "39-context-controls",
        "save_pub <- function",
    ),
    40: (
        "data/small/40-bgc-natural-products-frozen",
        "env/bgc-gecco-linux-64.lock",
        "env/antismash8-linux-64.lock",
        "scripts/prepare_article40_bgc.py",
        "scripts/run_article40_bgc.py",
        "scripts/summarize_article40_bgc.py",
        "scripts/freeze_article40_bgc.py",
        "scripts/validate_article40_bgc.py",
        "antiSMASH v8.0.4",
        "GECCO v0.10.3",
        "MIBiG v4.0",
        "20260740",
        "25% reciprocal",
        "40-tool-bgc-yield",
        "40-fragmentation-sensitivity",
        "40-bgc-type-profile",
        "40-mibig-similarity-abundance",
        "save_pub <- function",
    ),
    41: (
        "data/small/41-read-mapping-depth-frozen",
        "env/assembly-linux-64.lock",
        "scripts/prepare_article41_mapping_depth.py",
        "scripts/run_article41_mapping_depth.py",
        "scripts/summarize_article41_mapping_depth.py",
        "scripts/freeze_article41_mapping_depth.py",
        "scripts/validate_article41_mapping_depth.py",
        "Bowtie2 v2.5.5",
        "SAMtools v1.23.1",
        "MetaBAT2 v2.18",
        "20260741",
        "41-mapping-fate",
        "41-depth-concordance",
        "41-depth-breadth",
        "save_pub <- function",
    ),
    42: (
        "data/small/42-binning-comparison-frozen",
        "env/binning-linux-64.lock",
        "scripts/prepare_article42_binning.py",
        "scripts/run_article42_binning.py",
        "scripts/summarize_article42_binning.py",
        "scripts/freeze_article42_binning.py",
        "scripts/validate_article42_binning.py",
        "MetaBAT2 v2.18",
        "SemiBin2 v2.3.0",
        "VAMB v5.0.4",
        "20260742",
        "42-binner-quality-yield",
        "42-recovery-purity",
        "42-single-vs-multisample",
        "42-taxonomy-coverage",
        "save_pub <- function",
    ),
    43: (
        "data/small/43-bin-refinement-frozen",
        "env/mag-qc-linux-64.lock",
        "scripts/prepare_article43_refinement.py",
        "scripts/run_article43_refinement.py",
        "scripts/summarize_article43_refinement.py",
        "scripts/freeze_article43_refinement.py",
        "scripts/validate_article43_refinement.py",
        "DAS Tool 1.1.7",
        "Binette 1.2.1",
        "CheckM2 1.1.0",
        "20260743",
        "43-refinement-yield",
        "43-quality-landscape",
        "43-refinement-provenance",
        "43-method-selection",
        "save_pub <- function",
    ),
    44: (
        "data/small/44-mag-qc-mimag-graph-frozen",
        "env/checkm1-linux-64.lock",
        "scripts/prepare_article44_mag_qc.py",
        "scripts/run_article44_mag_qc.py",
        "scripts/summarize_article44_mag_qc.py",
        "scripts/freeze_article44_mag_qc.py",
        "scripts/validate_article44_mag_qc.py",
        "CheckM2 1.1.0",
        "GUNC 1.1.0",
        "CheckM1 1.2.5",
        "barrnap 1.10.5",
        "tRNAscan-SE 2.0.13",
        "20260744",
        "44-quality-landscape",
        "44-mimag-requirements",
        "44-assembly-graph-audit",
        "44-checkm-audit",
        "save_pub <- function",
    ),
}
PROHIBITED_PUBLIC_PATTERNS = (
    "作者代码通常长这样",
    "本篇可独立跑通",
    "这体现全系列",
    "即本文",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def frontmatter(path: Path) -> dict[str, Any]:
    text = path.read_text(encoding="utf-8")
    if not text.startswith("---\n"):
        return {}
    end = text.find("\n---\n", 4)
    if end < 0:
        return {}
    payload = yaml.safe_load(text[4:end])
    return payload if isinstance(payload, dict) else {}


def flatten_book_chapters(items: list[Any]) -> list[str]:
    paths: list[str] = []
    for item in items:
        if isinstance(item, str):
            paths.append(item)
        elif isinstance(item, dict):
            nested = item.get("chapters", [])
            if isinstance(nested, list):
                paths.extend(flatten_book_chapters(nested))
    return paths


def main() -> int:
    args = parse_args()
    root = args.project_root.resolve()
    manifest_path = args.manifest.resolve()
    normalized_manifest = load_manifest(manifest_path)
    manifest = yaml.safe_load(manifest_path.read_text(encoding="utf-8"))
    chapters = manifest.get("series", {}).get("chapters", [])
    errors: list[str] = []
    warnings: list[str] = []

    numbers = [item.get("number") for item in chapters]
    expected_numbers = list(range(1, EXPECTED_COUNT + 1))
    if numbers != expected_numbers:
        errors.append("series.chapters must contain consecutive numbers 1..77")

    validated_articles = manifest.get("series", {}).get("validated_articles", [])
    if validated_articles != expected_numbers:
        errors.append("series.validated_articles must contain consecutive numbers 1..77")

    if manifest.get("series", {}).get("total_articles") != EXPECTED_COUNT:
        errors.append("series.total_articles must equal 77")

    wechat_contract = manifest.get("publication", {}).get("wechat", {})
    if not isinstance(wechat_contract, dict):
        errors.append("publication.wechat must be a mapping")
        wechat_contract = {}
    title_prefix = str(wechat_contract.get("title_prefix", "")).strip()
    if title_prefix != EXPECTED_WECHAT_TITLE_PREFIX:
        errors.append(
            "publication.wechat.title_prefix must equal "
            f"{EXPECTED_WECHAT_TITLE_PREFIX}"
        )
    try:
        title_max_chars = int(wechat_contract.get("title_max_chars", 0))
    except (TypeError, ValueError):
        title_max_chars = 0
    if title_max_chars != EXPECTED_WECHAT_TITLE_MAX_CHARS:
        errors.append(
            "publication.wechat.title_max_chars must equal "
            f"{EXPECTED_WECHAT_TITLE_MAX_CHARS}"
        )

    wechat_titles: list[str] = []
    for item in chapters:
        number = int(item["number"])
        topic = str(item.get("wechat_title", item.get("title", ""))).strip()
        if not topic:
            errors.append(f"article {number:02d} has an empty WeChat topic title")
            continue
        if topic.startswith(EXPECTED_WECHAT_TITLE_PREFIX):
            errors.append(
                f"article {number:02d} WeChat topic repeats the series prefix"
            )
        public_title = f"{EXPECTED_WECHAT_TITLE_PREFIX}｜{number}. {topic}"
        wechat_titles.append(public_title)
        if len(public_title) > EXPECTED_WECHAT_TITLE_MAX_CHARS:
            errors.append(
                f"article {number:02d} WeChat title has {len(public_title)} "
                f"characters; maximum is {EXPECTED_WECHAT_TITLE_MAX_CHARS}"
            )
        if f"｜0{number}. " in public_title or f"｜第 {number} 篇" in public_title:
            errors.append(f"article {number:02d} uses a legacy WeChat title form")
    if len(wechat_titles) != EXPECTED_COUNT:
        errors.append("WeChat title sequence must contain exactly 77 titles")
    if len(set(wechat_titles)) != len(wechat_titles):
        errors.append("WeChat public titles must be unique")

    files = [item.get("file") for item in chapters]
    if len(files) != len(set(files)):
        errors.append("series.chapters contains duplicate file paths")

    quarto_path = root / "_quarto.yml"
    if not quarto_path.exists():
        errors.append("_quarto.yml is missing")
        quarto_chapters: list[str] = []
    else:
        quarto = yaml.safe_load(quarto_path.read_text(encoding="utf-8"))
        quarto_chapters = flatten_book_chapters(
            quarto.get("book", {}).get("chapters", [])
        )
        if quarto_chapters != files:
            errors.append("_quarto.yml chapter order does not match tutorial.yaml")

    for item in chapters:
        number = int(item["number"])
        chapter_path = root / item["file"]
        if not chapter_path.exists():
            errors.append(f"missing chapter file: {item['file']}")
            continue
        metadata = frontmatter(chapter_path)
        title = str(metadata.get("title", ""))
        if item["title"] not in title:
            errors.append(
                f"title mismatch for article {number:02d}: {item['file']}"
            )
        text = chapter_path.read_text(encoding="utf-8")
        if metadata.get("draft") is True:
            errors.append(f"article {number:02d} is still marked draft")
        execute = metadata.get("execute", {})
        if not isinstance(execute, dict) or not isinstance(execute.get("eval"), bool):
            errors.append(f"article {number:02d} must declare boolean execute.eval")
        elif number in EXECUTABLE_TOKEN_NUMBERS:
            expected_eval = number not in UPSTREAM_EVAL_FALSE
            if execute.get("eval") is not expected_eval:
                errors.append(
                    f"article {number:02d} must use eval: "
                    f"{str(expected_eval).lower()}"
                )
        if not isinstance(execute, dict) or execute.get("freeze") != "auto":
            errors.append(f"article {number:02d} must use freeze: auto")
        for section in REQUIRED_SECTIONS:
            if section not in text:
                errors.append(
                    f"article {number:02d} is missing section: {section}"
                )
        if number in EXECUTABLE_TOKEN_NUMBERS:
            for token in REQUIRED_ARTICLE_TOKENS[number]:
                if token not in text:
                    errors.append(
                        f"article {number:02d} is missing executable token: {token}"
                    )
        for pattern in PROHIBITED_PUBLIC_PATTERNS:
            if pattern in text:
                errors.append(
                    f"article {number:02d} contains prohibited public text: {pattern}"
                )

    for old_path in (
        "样板-28-跨队列验证-curatedMetagenomicData.qmd",
        "样板-44-MAG质控-CheckM2-GUNC-MIMAG.qmd",
    ):
        if (root / old_path).exists():
            errors.append(f"obsolete sample path still exists: {old_path}")

    status = "passed" if not errors else "failed"
    payload = {
        "status": status,
        "manifest_hash": normalized_manifest["manifest_hash"],
        "chapter_count": len(chapters),
        "expected_chapter_count": EXPECTED_COUNT,
        "validated_article_count": len(validated_articles),
        "structurally_validated_articles": expected_numbers,
        "executable_token_articles": sorted(EXECUTABLE_TOKEN_NUMBERS),
        "quarto_chapter_count": len(quarto_chapters),
        "wechat_title_count": len(wechat_titles),
        "wechat_title_max_chars_observed": max(
            (len(title) for title in wechat_titles), default=0
        ),
        "errors": errors,
        "warnings": warnings,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(payload, ensure_ascii=False, indent=2))
    return 0 if status == "passed" else 1


if __name__ == "__main__":
    raise SystemExit(main())
