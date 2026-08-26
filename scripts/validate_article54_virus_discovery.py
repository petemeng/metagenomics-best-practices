#!/usr/bin/env python3
"""Fail-closed validation for Article 54 virus-discovery evidence."""

from __future__ import annotations

import argparse
import json
import math
import re
from collections import Counter
from pathlib import Path

from article42_44_validation_utils import (
    Audit,
    as_bool,
    audit_checksums,
    audit_figures,
    finish,
    read_tsv,
    sha256,
)


FIGURES = (
    "54-discovery-consensus",
    "54-checkv-quality",
    "54-votu-threshold",
    "54-virus-evidence-ladder",
    "54-library-design-bias",
)
EXPECTED_RAW = {
    "PMC6871006.xml": (
        183_966,
        "289c221ece18f17a0a8930a354cfc3e01b712751abfdd52ff7b054c6f31c36fa",
    ),
    "checkv-test-sequences.fna": (
        910_376,
        "4347d1f5f52d4e2b6986845afb8e401fd146bb5a9628768e0ba125ec380144f5",
    ),
    "checkv-upstream-ground-truth-quality-summary.tsv": (
        5_518,
        "f0d692ab02446bca007722b0f97a928b4cc893c94bcd311e151196bb9afa5a76",
    ),
}
EXPECTED_DATABASE_IDS = {
    "genomad-db-v1.9",
    "checkv-db-v1.5-mirror",
    "virsorter2-db-core",
    "virsorter2-combined-hmm-00",
    "virsorter2-combined-hmm-01",
    "virsorter2-combined-hmm-02",
    "virsorter2-pfam-archaea",
    "virsorter2-pfam-bacteria",
    "virsorter2-pfam-eukaryota",
    "virsorter2-pfam-mixed",
    "virsorter2-pfam-viruses",
    "virsorter2-pfam-accessions",
}


def near(value: object, expected: float, tolerance: float = 1e-8) -> bool:
    try:
        return math.isclose(float(value), expected, abs_tol=tolerance)
    except (TypeError, ValueError):
        return False


def fasta_ids(path: Path) -> list[str]:
    return [
        line[1:].split()[0]
        for line in path.read_text(encoding="utf-8").splitlines()
        if line.startswith(">")
    ]


def valid_tilde_fences(text: str) -> bool:
    open_fence = False
    for line in text.splitlines():
        if not line.startswith("~~~"):
            continue
        if open_fence:
            if line.strip() != "~~~":
                return False
            open_fence = False
        else:
            if line.strip() not in {"~~~bash", "~~~r", "~~~text"}:
                return False
            open_fence = True
    return not open_fence


def audit_chapter(chapter: Path, frozen: Path, audit: Audit) -> None:
    text = chapter.read_text(encoding="utf-8")
    plot = (frozen / "scripts/plot_article54_virus_discovery.R").read_text(
        encoding="utf-8"
    )
    checks = {
        "draft-false": "draft: false" in text,
        "eval-false": re.search(r"execute:\s*\n\s+eval:\s+false", text) is not None,
        "freeze-auto": "freeze: auto" in text,
        "bibliography": "bibliography: ../references.bib" in text,
        "seed": "20260754" in text,
        "inline-theme": all(
            token in text
            for token in (
                "pal_pub <-",
                "scale_color_pub <-",
                "scale_fill_pub",
                "theme_pub <-",
                "save_pub <-",
            )
        ),
        "resource-contract": all(
            token in text for token in ("CPU", "RAM", "磁盘", "秒")
        ),
        "real-fixture": all(
            token in text
            for token in (
                "46 条真实回归序列",
                "PMC6871006",
                "6a118f20e895105ce0e4f10257955494c60f1293",
            )
        ),
        "frozen-input": "data/small/54-virus-discovery-quality-frozen" in text,
        "methods-results": all(
            token in text
            for token in (
                "Discovery、quality 与 vOTU 的英文模板",
                "Results template",
                "Read depth from amplified virome libraries was treated as nonquantitative",
            )
        ),
        "versions": all(
            token in text
            for token in (
                "geNomad v1.12.0",
                "VirSorter2 v2.2.4",
                "CheckV v1.1.1",
                "BLASTN v2.17.0",
            )
        ),
        "database-identity": all(
            token in text
            for token in (
                "v1.9 / ICTV MSL39",
                "v1.5 / 2023-01-10",
                "Zenodo 4269607 v0.4",
                "UNOFFICIAL MIRROR",
            )
        ),
        "hard-gates": all(
            token in text
            for token in (
                "--min-score 0.5",
                "--min_ani 95",
                "--min_tcov 85",
                "较短序列 ≥85%",
            )
        ),
        "result-counts": all(
            token in text
            for token in (
                "geNomad 检出 **41**",
                "VirSorter2 默认 `score ≥ 0.5` 检出 **43**",
                "交集 **40**",
                "并集 **44**",
                "只保留 **38**",
                "95.16% ANI、84.76% AF",
            )
        ),
        "interpretation-boundaries": all(
            token in text
            for token in (
                "不能替代病毒发现工具",
                "不是自动“完整基因组”",
                "Not-determined",
                "nonquantitative",
                "不能直接并表比较",
            )
        ),
        "database-build-complete": all(
            token in text
            for token in (
                "conda_envs/328974fc",
                "virsorter2-pfam-archaea:Pfam-A-Archaea.hmm",
                "a67245759c1529aee6485825ae5c3912",
            )
        ),
        "citations": all(
            token in text
            for token in (
                "@camargo2024genomad",
                "@guo2021virsorter2",
                "@nayfach2021checkv",
                "@roux2019miuvig",
                "@parrasmolto2018viromebias",
                "@conceicaoneto2015netovir",
            )
        ),
        "figures": all(f"../figures/{stem}.png" in text for stem in FIGURES),
        "vector-raster-export": all(
            token in text
            for token in (".pdf", ".png", ".tiff", 'compression = "lzw"')
        ),
        "plot-labels-english": re.search(r"[\u3400-\u9fff]", plot) is None,
        "code-fences": valid_tilde_fences(text),
        "no-duplicate-project-root": re.search(
            r'--project-root "\$ROOT"\s*\\\s*\n\s*--project-root', text
        )
        is None,
        "no-source-theme": 'source("R/theme_pub.R")' not in text
        and "source('R/theme_pub.R')" not in text,
        "no-placeholders": re.search(
            r"__[A-Z][A-Z0-9_]*__|\bTODO\b|\bTBD\b|\bNNN\b", text
        )
        is None,
        "no-meta-prose": not any(
            token in text
            for token in (
                "本篇可独立",
                "本文可独立",
                "全系列约定",
                "接口只学一次",
                "作者代码通常长这样",
                "（即本文）",
                "无头服务器",
            )
        ),
    }
    sections = (
        "这一步对应论文里的哪张图",
        "理论：",
        "准备工作",
        "可复制代码",
        "审计与升级",
        "出版级美化",
        "常见坑",
        "这段 Methods 怎么写",
        "换成你自己的数据怎么做",
        "参考",
    )
    for section in sections:
        checks[f"section-{section}"] = (
            re.search(rf"(?m)^##\s+{re.escape(section)}", text) is not None
        )
    for check, status in checks.items():
        audit.add("Chapter", check, status, check)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--frozen-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--figure-dir", type=Path, required=True)
    parser.add_argument("--chapter", type=Path, required=True)
    args = parser.parse_args()
    frozen = args.frozen_dir.resolve()
    audit = Audit()
    audit_checksums(frozen, audit)

    summary = json.loads((frozen / "summary.json").read_text(encoding="utf-8"))
    exact = {
        "article": 54,
        "seed": 20260754,
        "input_sequences": 46,
        "input_total_bp": 894_215,
        "input_n50_bp": 27_099,
        "genomad_detected": 41,
        "virsorter2_default_detected": 43,
        "virsorter2_high_confidence": 38,
        "detected_by_both": 40,
        "genomad_only": 1,
        "virsorter2_only": 3,
        "detected_union": 44,
        "detected_by_neither": 2,
        "checkv_proviruses": 3,
        "checkv_low_confidence_dtr": 5,
        "genomad_dtr": 3,
        "votu_clusters": 46,
        "votu_multi_member_clusters": 0,
        "votu_nonself_pairs": 10,
        "votu_pairs_passing_both": 0,
    }
    for key, expected in exact.items():
        audit.add(
            "Summary",
            key,
            summary.get(key) == expected,
            {"expected": expected, "observed": summary.get(key)},
        )
    audit.add(
        "Summary",
        "quality-counts",
        summary.get("checkv_quality_counts")
        == {
            "High-quality": 1,
            "Medium-quality": 4,
            "Low-quality": 33,
            "Not-determined": 8,
        },
        summary.get("checkv_quality_counts"),
    )
    audit.add(
        "Summary",
        "checkv-fixture-exact",
        summary.get("checkv_fixture_exact_match") is True,
        summary.get("checkv_fixture_exact_match"),
    )
    audit.add(
        "Summary",
        "boundary-pair",
        set(summary.get("boundary_pair", []))
        == {"UHGV-0001702", "UHGV-0001715"}
        and near(summary.get("boundary_pair_ani_pct"), 95.16)
        and near(summary.get("boundary_pair_shorter_af_pct"), 84.76),
        summary.get("boundary_pair"),
    )
    audit.add(
        "Summary",
        "no-random-output",
        summary.get("random_output_requested") is False,
        summary.get("random_output_requested"),
    )

    assets = read_tsv(frozen / "asset-check-audit.tsv")
    audit.add(
        "Input",
        "three-checksum-gated-assets",
        len(assets) == 3 and all(as_bool(row["ChecksumPass"]) for row in assets),
        len(assets),
    )
    raw_results = []
    for name, (expected_bytes, expected_hash) in EXPECTED_RAW.items():
        path = frozen / "raw" / name
        raw_results.append(
            path.is_file()
            and path.stat().st_size == expected_bytes
            and sha256(path) == expected_hash
        )
    audit.add("Input", "raw-source-snapshots", all(raw_results), raw_results)
    ids = fasta_ids(frozen / "raw/checkv-test-sequences.fna")
    audit.add(
        "Input",
        "46-unique-fasta-records",
        len(ids) == 46 and len(set(ids)) == 46,
        len(ids),
    )
    fixture = frozen / "raw/checkv-upstream-ground-truth-quality-summary.tsv"
    observed = frozen / "upstream/checkv-quality-summary.tsv"
    audit.add(
        "Input",
        "byte-identical-checkv-regression",
        fixture.read_bytes() == observed.read_bytes(),
        sha256(observed),
    )
    manifest = read_tsv(frozen / "asset-manifest.tsv")
    audit.add(
        "Input",
        "source-provenance",
        len(manifest) == 3
        and {row["DOI"] for row in manifest}
        == {"10.1038/nbt.4306", "10.1038/s41587-020-00774-7"}
        and {row["License"] for row in manifest} == {"CC BY 4.0", "BSD-3-Clause"}
        and all(row["SourceIdentity"] for row in manifest),
        Counter(row["DOI"] for row in manifest),
    )

    database = read_tsv(frozen / "54-virus-database-manifest.tsv")
    database_by_id = {row["database_id"]: row for row in database}
    audit.add(
        "Database",
        "twelve-locked-assets",
        len(database) == 12
        and set(database_by_id) == EXPECTED_DATABASE_IDS
        and len(database_by_id) == len(database),
        len(database),
    )
    audit.add(
        "Database",
        "immutable-enabled-checksums",
        all(
            row["download_gate"] == "enabled"
            and row["validation_status"].startswith("VERIFIED_")
            and re.fullmatch(r"[0-9a-f]{32}", row["expected_checksum"])
            and int(row["expected_compressed_bytes"]) > 0
            and "latest" not in row["archive_url"]
            for row in database
        ),
        "12 publisher MD5 and byte-count gates",
    )
    audit.add(
        "Database",
        "release-identities",
        database_by_id["genomad-db-v1.9"]["release_id"] == "v1.9-ICTV-MSL39"
        and database_by_id["checkv-db-v1.5-mirror"]["release_id"] == "v1.5"
        and database_by_id["virsorter2-db-core"]["release_id"]
        == "Zenodo-4269607-v0.4"
        and "UNOFFICIAL MIRROR"
        in database_by_id["checkv-db-v1.5-mirror"]["notes"],
        "geNomad v1.9; CheckV v1.5 mirror; VirSorter2 v0.4",
    )

    evidence = read_tsv(frozen / "virus-evidence-matrix.tsv")
    pattern_counts = Counter(row["DiscoveryPattern"] for row in evidence)
    audit.add(
        "Discovery",
        "evidence-coordinate",
        len(evidence) == 46
        and len({row["ContigID"] for row in evidence}) == 46
        and {row["ContigID"] for row in evidence} == set(ids),
        len(evidence),
    )
    audit.add(
        "Discovery",
        "overlap-counts",
        pattern_counts
        == Counter({"Both": 40, "VirSorter2 only": 3, "Neither": 2, "geNomad only": 1}),
        pattern_counts,
    )
    discordant = {
        row["ContigID"]: row["DiscoveryPattern"]
        for row in evidence
        if row["DiscoveryPattern"] != "Both"
    }
    audit.add(
        "Discovery",
        "discordant-identities",
        discordant
        == {
            "2013338001_____MeugFOFF_C1475": "geNomad only",
            "2013843002_____DCKB1_C2382": "VirSorter2 only",
            "2013843002_____DCKB1_C3067": "VirSorter2 only",
            "2014031003_____YNP3_C2743": "Neither",
            "2014031003_____YNP3_C3059": "VirSorter2 only",
            "2014031003_____YNP3_C3265": "Neither",
        },
        discordant,
    )
    audit.add(
        "Discovery",
        "high-confidence-rule",
        sum(as_bool(row["VirSorter2HighConfidence"]) for row in evidence) == 38
        and all(
            not as_bool(row["VirSorter2HighConfidence"])
            or (
                float(row["VirSorter2FinalScore"]) >= 0.9
                or (
                    float(row["VirSorter2FinalScore"]) >= 0.7
                    and int(row["VirSorter2Hallmarks"]) >= 1
                )
            )
            for row in evidence
        ),
        "38 calls",
    )
    genomad = read_tsv(frozen / "upstream/genomad-virus-summary.tsv")
    genomad_all = read_tsv(
        frozen / "upstream/genomad-aggregated-classification.tsv"
    )
    vs2_final = read_tsv(frozen / "upstream/virsorter2-final-viral-score.tsv")
    vs2_all = read_tsv(frozen / "upstream/virsorter2-all-fullseq-proba.tsv")
    audit.add(
        "Discovery",
        "upstream-output-counts",
        [len(genomad), len(genomad_all), len(vs2_final), len(vs2_all)]
        == [41, 46, 43, 46],
        [len(genomad), len(genomad_all), len(vs2_final), len(vs2_all)],
    )

    quality = Counter(row["CheckVQuality"] for row in evidence)
    audit.add(
        "CheckV",
        "quality-tiers",
        quality
        == Counter(
            {
                "High-quality": 1,
                "Medium-quality": 4,
                "Low-quality": 33,
                "Not-determined": 8,
            }
        ),
        quality,
    )
    proviruses = [row for row in evidence if row["CheckVProvirus"] == "Yes"]
    special = next(row for row in evidence if row["ContigID"] == "UHGV-0000346")
    audit.add(
        "CheckV",
        "provirus-boundary",
        len(proviruses) == 3
        and special["CheckVProviralLength"] == "4899"
        and special["CheckVQuality"] == "High-quality"
        and near(special["CompletenessPct"], 100.0)
        and near(special["ContaminationPct"], 30.29),
        [row["ContigID"] for row in proviruses],
    )
    terminal = read_tsv(frozen / "terminal-repeat-audit.tsv")
    audit.add(
        "CheckV",
        "low-confidence-terminal-repeats",
        len(terminal) == 5
        and sum(row["geNomadTopology"] == "DTR" for row in terminal) == 3
        and all(row["CheckVPrediction"] == "DTR" for row in terminal)
        and all(row["CheckVConfidence"] == "low" for row in terminal)
        and all("ambiguous bases" in row["CheckVReason"] for row in terminal)
        and all(not as_bool(row["MIUViGFinished"]) for row in terminal),
        len(terminal),
    )

    pairs = read_tsv(frozen / "votu-pairwise-threshold-audit.tsv")
    boundary = [row for row in pairs if as_bool(row["BoundaryPair"])]
    audit.add(
        "vOTU",
        "ten-pairs-no-merge",
        len(pairs) == 10
        and all(not as_bool(row["SameVOTU"]) for row in pairs)
        and sum(as_bool(row["PassANI95"]) for row in pairs) == 5
        and sum(as_bool(row["PassAF85"]) for row in pairs) == 0,
        len(pairs),
    )
    audit.add(
        "vOTU",
        "boundary-two-gate-decision",
        len(boundary) == 1
        and near(boundary[0]["ANIpct"], 95.16)
        and near(boundary[0]["ShorterAlignmentFractionPct"], 84.76)
        and as_bool(boundary[0]["PassANI95"])
        and not as_bool(boundary[0]["PassAF85"]),
        boundary,
    )
    clusters = read_tsv(frozen / "votu-cluster-summary.tsv")
    raw_clusters = [
        line.split("\t")
        for line in (frozen / "upstream/votu-clusters.tsv")
        .read_text(encoding="utf-8")
        .splitlines()
        if line
    ]
    audit.add(
        "vOTU",
        "46-singletons",
        len(clusters) == 46
        and all(int(row["Members"]) == 1 for row in clusters)
        and len(raw_clusters) == 46
        and all(len(row) == 2 and row[0] == row[1] for row in raw_clusters),
        len(clusters),
    )

    assertions = read_tsv(frozen / "miuvig-source-assertions.tsv")
    audit.add(
        "MIUViG",
        "five-source-assertions",
        len(assertions) == 5
        and all(as_bool(row["Pass"]) for row in assertions)
        and all(row["Source"] == "PMC6871006" for row in assertions),
        len(assertions),
    )
    library = read_tsv(frozen / "library-design-bias-audit.tsv")
    audit.add(
        "MIUViG",
        "library-design-boundary",
        len(library) == 4
        and {row["Library"] for row in library}
        == {"Total metagenome", "Virus-enriched virome", "Amplified virome"}
        and any(row["QuantitativeUse"] == "Nonquantitative read depth" for row in library)
        and any("0.22 um" in row["KnownBlindSpot"] for row in library),
        Counter(row["Library"] for row in library),
    )
    ladder = read_tsv(frozen / "virus-evidence-ladder.tsv")
    audit.add(
        "MIUViG",
        "five-stage-evidence-ladder",
        [int(row["Order"]) for row in ladder] == [1, 2, 3, 4, 5]
        and ladder[-1]["Decision"] == "MIUViG reporting",
        [row["Decision"] for row in ladder],
    )

    resources = read_tsv(frozen / "resource-summary.tsv")
    resource_by_tool = {row["Tool"]: row for row in resources}
    audit.add(
        "Execution",
        "measured-resources",
        set(resource_by_tool) == {"geNomad", "CheckV", "VirSorter2"}
        and all(
            int(row["Threads"]) == 16
            and int(row["ExitStatus"]) == 0
            and float(row["WallSeconds"]) > 0
            and float(row["PeakRAMGiB"]) > 0
            and int(row["OutputBytes"]) > 0
            for row in resources
        )
        and near(resource_by_tool["geNomad"]["WallSeconds"], 131.88)
        and near(resource_by_tool["CheckV"]["PeakRAMGiB"], 1.945)
        and near(resource_by_tool["VirSorter2"]["WallSeconds"], 881.5),
        resources,
    )
    determinism = read_tsv(frozen / "determinism-audit.tsv")
    audit.add(
        "Execution",
        "determinism-contract",
        len(determinism) == 5
        and all(not as_bool(row["RandomProcess"]) for row in determinism)
        and all(row["Status"] == "PASS" for row in determinism),
        len(determinism),
    )
    versions = {row["Tool"]: row["Version"] for row in read_tsv(frozen / "tool-versions.tsv")}
    audit.add(
        "Version",
        "analysis-stack",
        versions
        == {
            "geNomad": "1.12.0",
            "VirSorter2": "2.2.4",
            "CheckV": "1.1.1",
            "BLASTN / anicalc / aniclust": "BLAST 2.17.0 / CheckV 1.1.1",
        },
        versions,
    )
    plot_versions = {
        row["Package"]: row["Version"]
        for row in read_tsv(frozen / "plot-software-versions.tsv")
    }
    audit.add(
        "Version",
        "plot-stack",
        plot_versions
        == {
            "R": "4.4.1",
            "ggplot2": "3.5.2",
            "dplyr": "1.1.4",
            "readr": "2.1.5",
            "tidyr": "1.3.1",
            "patchwork": "1.3.2",
        },
        plot_versions,
    )
    lock_lengths = {}
    for name, minimum in (
        ("virus-discovery-linux-64.lock", 180),
        ("virsorter2-linux-64.lock", 140),
        ("virsorter2-runtime-linux-64.lock", 110),
    ):
        lines = (frozen / "env" / name).read_text(encoding="utf-8").splitlines()
        lock_lengths[name] = len(lines)
        audit.add(
            "Version",
            f"exact-lock-{name}",
            "@EXPLICIT" in lines and len(lines) >= minimum,
            len(lines),
        )

    audit_chapter(args.chapter.resolve(), frozen, audit)
    audit_figures(args.figure_dir.resolve(), audit, FIGURES)
    return finish(
        article=54,
        audit=audit,
        output=args.output_dir.resolve(),
        payload={
            "input_sequences": summary["input_sequences"],
            "detected_union": summary["detected_union"],
            "checkv_quality_counts": summary["checkv_quality_counts"],
            "votu_clusters": summary["votu_clusters"],
        },
    )


if __name__ == "__main__":
    raise SystemExit(main())
