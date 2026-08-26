#!/usr/bin/env python3
"""Fail-closed validation for Article 55 virus taxonomy and abundance evidence."""

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
    "55-breadth-heatmap",
    "55-threshold-prevalence",
    "55-abundance-multimapping",
    "55-votu-taxonomy",
    "55-lifecycle-evidence",
)
EXPECTED_RAW = {
    "PMC10926689.xml": (
        140_775,
        "74d28ef0c7d21a14eaf80a8f517f5f9eadf5385f586b8f20cf280535491d4c4e",
    ),
    "PMC11521287.xml": (
        157_363,
        "a28b8f2478ace607d497eb61f23e5c097107e30810d2ee077473cfd36b6f44eb",
    ),
    "mgen-10-1198-s002.xlsx": (
        376_869,
        "4f7b44c6c276a6b2b60142543c5f2515917370ddac340a40bd85b37f13475642",
    ),
    "ena-prjeb56639-filereport.tsv": (
        3_158,
        "c56c6b024f43dccbb56b22b8aeccdd7b3345a00081b7b19afce713b099eeb096",
    ),
}


def near(value: object, expected: float, tolerance: float = 1e-6) -> bool:
    try:
        return math.isclose(float(value), expected, abs_tol=tolerance)
    except (TypeError, ValueError):
        return False


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
    plot = (frozen / "scripts/plot_article55_virus_abundance.R").read_text(
        encoding="utf-8"
    )
    checks = {
        "draft-false": "draft: false" in text,
        "eval-false": re.search(r"execute:\s*\n\s+eval:\s+false", text) is not None,
        "freeze-auto": "freeze: auto" in text,
        "bibliography": "bibliography: ../references.bib" in text,
        "seed": "20260755" in text,
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
        "resource-contract": all(token in text for token in ("CPU", "RAM", "磁盘", "秒")),
        "real-data": all(
            token in text
            for token in (
                "PRJEB56639",
                "ERR10359653",
                "ERR10359656",
                "ERR10359658",
                "Supplementary Table S4",
                "15-phage mock virome",
            )
        ),
        "frozen-input": "data/small/55-virus-taxonomy-abundance-frozen" in text,
        "threshold-contract": all(
            token in text
            for token in (
                "minid=0.90",
                "ambiguous=all",
                "ambiguous=toss",
                "mean depth ≥1×",
                "breadth ≥75%",
                "--min_ani 95",
                "--min_tcov 85",
            )
        ),
        "result-counts": all(
            token in text
            for token in (
                "**9、11、10**",
                "合并 reads 后检出 **11**",
                "**12 个 vOTU**",
                "12/12",
                "9/12",
                "少 127 bp",
            )
        ),
        "interpretation-boundaries": all(
            token in text
            for token in (
                "technical occupancy",
                "不能当三名受试者",
                "逐成员 abundance 不可辨识",
                "family 未分类不等于“非病毒”",
                "Physical state",
                "不能一一对应",
            )
        ),
        "versions": all(
            token in text
            for token in (
                "geNomad v1.12.0",
                "geNomad DB v1.9 / ICTV MSL39",
                "BBMap v38.69",
                "SAMtools 1.24",
                "BLASTN v2.17.0",
            )
        ),
        "methods-results": all(
            token in text
            for token in (
                "Abundance、taxonomy 与 lifestyle 英文模板",
                "Results template",
                "Technical-library occupancy was reported separately from biological prevalence",
            )
        ),
        "citations": all(
            token in text
            for token in (
                "@cook2024viromicsbenchmark",
                "@roux2019miuvig",
                "@camargo2024genomad",
                "@bushnell2014bbmap",
                "@li2009samtools",
                "@zhang2024deeppl",
                "@ictv2024msl39",
            )
        ),
        "figures": all(f"../figures/{stem}.png" in text for stem in FIGURES),
        "vector-raster-export": all(
            token in text for token in (".pdf", ".png", ".tiff", 'compression = "lzw"')
        ),
        "plot-labels-english": re.search(r"[\u3400-\u9fff]", plot) is None,
        "code-fences": valid_tilde_fences(text),
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
        "article": 55,
        "seed": 20260755,
        "mock_phages": 15,
        "illumina_libraries": 3,
        "illumina_read_pairs": 2_062_012,
        "primary_presence_min_depth_x": 1.0,
        "primary_presence_min_breadth_pct": 75.0,
        "genomad_detected": 15,
        "votu_clusters": 12,
        "votu_multi_member_clusters": 1,
        "largest_votu_members": 4,
        "j1_j2_ani_pct": 100.0,
        "j1_j2_shorter_af_pct": 100.0,
        "deeppl_confirmed_correct": 12,
        "phatyp_confirmed_correct": 9,
        "slur29_reference_length_delta_bp": -127,
        "bbmap_smoke_j1_depth_x": 0.0,
        "bbmap_smoke_j2_depth_x": 2.0,
        "random_output_requested": False,
    }
    for key, expected in exact.items():
        observed = summary.get(key)
        status = near(observed, expected) if isinstance(expected, float) else observed == expected
        audit.add("Summary", key, status, {"expected": expected, "observed": observed})
    audit.add(
        "Summary",
        "detected-at-75pct",
        summary.get("detected_at_75pct")
        == {
            "lib1_illumina": 9,
            "lib2_illumina": 11,
            "lib3_illumina": 10,
            "pooled_illumina": 11,
        },
        summary.get("detected_at_75pct"),
    )
    audit.add(
        "Summary",
        "replicate-prevalence-counts",
        summary.get("replicate_prevalence_counts")
        == {"0": 4, "1": 1, "2": 1, "3": 9},
        summary.get("replicate_prevalence_counts"),
    )
    audit.add(
        "Summary",
        "universal-nondetections",
        set(summary.get("universal_nondetections", []))
        == {"CDMH1", "HP1", "phix174", "vB_Eco_mar005P1"},
        summary.get("universal_nondetections"),
    )
    audit.add(
        "Summary",
        "largest-votu-members",
        set(summary.get("largest_votu_member_ids", []))
        == {
            "vB_EcoS_swan01",
            "vB_Eco_mar002J2",
            "vB_Eco_mar001J1",
            "vB_Eco_SLUR29",
        },
        summary.get("largest_votu_member_ids"),
    )
    audit.add(
        "Summary",
        "lifecycle-counts",
        summary.get("confirmed_lifecycle_counts")
        == {"Not reported": 1, "Temperate": 5, "Unknown": 2, "Virulent": 7},
        summary.get("confirmed_lifecycle_counts"),
    )
    family_expected = {
        "Autographiviridae": 2,
        "Demerecviridae": 1,
        "Drexlerviridae": 4,
        "Kyanoviridae": 1,
        "Microviridae": 1,
        "Naomviridae": 1,
        "Straboviridae": 1,
        "Unclassified at family": 4,
    }
    audit.add(
        "Summary",
        "family-counts",
        summary.get("genomad_family_counts") == family_expected,
        summary.get("genomad_family_counts"),
    )

    for filename, (expected_bytes, expected_hash) in EXPECTED_RAW.items():
        path = frozen / "raw" / filename
        audit.add("Raw", f"{filename}-exists", path.is_file(), str(path))
        if path.is_file():
            audit.add("Raw", f"{filename}-bytes", path.stat().st_size == expected_bytes, path.stat().st_size)
            audit.add("Raw", f"{filename}-sha256", sha256(path) == expected_hash, sha256(path))

    abundance = read_tsv(frozen / "illumina-abundance-evidence.tsv")
    audit.add("Abundance", "rows", len(abundance) == 60, len(abundance))
    audit.add("Abundance", "phages", len({row["PhageID"] for row in abundance}) == 15, "15 expected")
    audit.add(
        "Abundance",
        "presence-gates",
        all(near(row["PresenceDepthGateX"], 1) and near(row["PresenceBreadthGatePct"], 75) for row in abundance),
        "1x/75%",
    )

    threshold = read_tsv(frozen / "breadth-threshold-sensitivity.tsv")
    expected_threshold = {
        ("lib1_illumina", "50"): 11,
        ("lib1_illumina", "70"): 10,
        ("lib1_illumina", "75"): 9,
        ("lib1_illumina", "90"): 9,
        ("lib2_illumina", "50"): 11,
        ("lib2_illumina", "70"): 11,
        ("lib2_illumina", "75"): 11,
        ("lib2_illumina", "90"): 11,
        ("lib3_illumina", "50"): 11,
        ("lib3_illumina", "70"): 10,
        ("lib3_illumina", "75"): 10,
        ("lib3_illumina", "90"): 9,
        ("pooled_illumina", "50"): 11,
        ("pooled_illumina", "70"): 11,
        ("pooled_illumina", "75"): 11,
        ("pooled_illumina", "90"): 11,
    }
    observed_threshold = {
        (row["Dataset"], row["BreadthThresholdPct"]): int(row["PhagesDetected"])
        for row in threshold
    }
    audit.add(
        "Abundance",
        "threshold-table",
        observed_threshold == expected_threshold,
        [
            {"Dataset": dataset, "BreadthThresholdPct": threshold_value, "PhagesDetected": count}
            for (dataset, threshold_value), count in sorted(observed_threshold.items())
        ],
    )

    prevalence = read_tsv(frozen / "replicate-prevalence.tsv")
    prevalence_count = Counter(int(row["ReplicatesPresent"]) for row in prevalence)
    audit.add("Abundance", "prevalence-rows", len(prevalence) == 15, len(prevalence))
    audit.add("Abundance", "prevalence-distribution", prevalence_count == Counter({3: 9, 0: 4, 1: 1, 2: 1}), dict(prevalence_count))

    clusters = read_tsv(frozen / "votu-cluster-summary.tsv")
    large = [row for row in clusters if int(row["MemberCount"]) > 1]
    audit.add("vOTU", "cluster-count", len(clusters) == 12, len(clusters))
    audit.add("vOTU", "single-multimember", len(large) == 1 and int(large[0]["MemberCount"]) == 4, large)
    if large:
        audit.add(
            "vOTU",
            "four-member-identity",
            set(large[0]["Members"].split(","))
            == {"vB_EcoS_swan01", "vB_Eco_mar002J2", "vB_Eco_mar001J1", "vB_Eco_SLUR29"},
            large[0]["Members"],
        )

    pairs = read_tsv(frozen / "votu-pairwise-threshold-audit.tsv")
    j12 = next(
        row
        for row in pairs
        if {row["PhageA"], row["PhageB"]}
        == {"vB_Eco_mar001J1", "vB_Eco_mar002J2"}
    )
    audit.add("vOTU", "alignable-pairs", len(pairs) == 7, len(pairs))
    audit.add(
        "vOTU",
        "j1-j2-exact",
        near(j12["ANIPct"], 100) and near(j12["ShorterAlignmentFractionPct"], 100) and as_bool(j12["SameVOTU"]),
        j12,
    )

    taxonomy = read_tsv(frozen / "taxonomy-votu-ledger.tsv")
    audit.add("Taxonomy", "rows", len(taxonomy) == 15, len(taxonomy))
    audit.add(
        "Taxonomy",
        "family-unclassified",
        sum(row["geNomadFamily"] == "Unclassified at family" for row in taxonomy) == 4,
        Counter(row["geNomadFamily"] for row in taxonomy),
    )
    audit.add(
        "Taxonomy",
        "database-release",
        all(row["TaxonomyDatabase"] == "geNomad DB v1.9; ICTV MSL39" for row in taxonomy),
        "v1.9/MSL39",
    )

    lifecycle = read_tsv(frozen / "lifecycle-evidence-ledger.tsv")
    confirmed = [row for row in lifecycle if row["ConfirmedLifecycle"] in {"Temperate", "Virulent"}]
    audit.add("Lifecycle", "rows", len(lifecycle) == 15, len(lifecycle))
    audit.add("Lifecycle", "confirmed-count", len(confirmed) == 12, len(confirmed))
    audit.add("Lifecycle", "deeppl-matches", sum(as_bool(row["DeepPLMatchesConfirmed"]) for row in lifecycle) == 12, "12")
    audit.add("Lifecycle", "phatyp-matches", sum(as_bool(row["PhaTYPMatchesConfirmed"]) for row in lifecycle) == 9, "9")
    audit.add(
        "Lifecycle",
        "unknowns-not-promoted",
        all(row["EvidenceStatus"] == "Prediction only or not reported" for row in lifecycle if row["ConfirmedLifecycle"] in {"Unknown", "Not reported"}),
        "prediction boundary",
    )

    smoke = read_tsv(frozen / "bbmap-duplicate-reference-smoke.tsv")
    smoke_by_id = {row["PhageID"]: row for row in smoke}
    audit.add("Smoke", "rows", len(smoke) == 15, len(smoke))
    audit.add("Smoke", "j1-zero", near(smoke_by_id["vB_Eco_mar001J1"]["MeanDepthX"], 0), smoke_by_id["vB_Eco_mar001J1"])
    audit.add("Smoke", "j2-double", near(smoke_by_id["vB_Eco_mar002J2"]["MeanDepthX"], 2), smoke_by_id["vB_Eco_mar002J2"])

    reference = read_tsv(frozen / "source/reference-audit.tsv")
    audit.add("Reference", "rows", len(reference) == 15, len(reference))
    audit.add("Reference", "all-identity-pass", all(as_bool(row["ChecksumAndLengthPass"]) for row in reference), "15/15")
    slur = next(row for row in reference if row["PhageID"] == "vB_Eco_SLUR29")
    audit.add("Reference", "slur29-delta", int(slur["LengthDeltaBp"]) == -127, slur)

    for assertion_file in ("author-source-assertions.tsv", "lifecycle-source-assertions.tsv"):
        rows = read_tsv(frozen / assertion_file)
        audit.add("Source", assertion_file, bool(rows) and all(as_bool(row["Pass"]) for row in rows), rows)

    resources = read_tsv(frozen / "resource-summary.tsv")
    audit.add("Resource", "rows", len(resources) == 6, len(resources))
    audit.add("Resource", "all-exit-zero", all(int(row["ExitStatus"]) == 0 for row in resources), resources)
    genomad = next(row for row in resources if row["Tool"] == "geNomad")
    audit.add("Resource", "genomad-time", near(genomad["WallSeconds"], 145.36), genomad)
    audit.add("Resource", "genomad-ram", near(genomad["PeakRAMGiB"], 2.654), genomad)
    bbmap = next(row for row in resources if row["Tool"] == "BBMap smoke test")
    audit.add("Resource", "bbmap-time", near(bbmap["WallSeconds"], 29.29), bbmap)
    audit.add("Resource", "bbmap-ram", near(bbmap["PeakRAMGiB"], 1.395), bbmap)

    versions = {row["Tool"]: row["Version"] for row in read_tsv(frozen / "tool-versions.tsv")}
    audit.add("Environment", "versions", versions == {
        "geNomad": "1.12.0",
        "BLASTN": "2.17.0",
        "anicalc / aniclust": "CheckV 1.1.1 distribution",
        "BBMap": "38.69 (study-reported and smoke-tested)",
        "SAMtools": "1.24",
        "SeqKit": "2.13.0",
    }, versions)
    lock = (frozen / "env/virus-abundance-linux-64.lock").read_text(encoding="utf-8")
    audit.add("Environment", "lock-explicit", "@EXPLICIT" in lock, "@EXPLICIT")
    audit.add("Environment", "lock-bbmap", "bbmap-38.69-h516909a_0" in lock, "BBMap 38.69")
    audit.add("Environment", "lock-samtools", "samtools-1.24-h9dcdb79_1" in lock, "SAMtools 1.24")
    audit.add("Environment", "lock-seqkit", "seqkit-2.13.0-he881be0_0" in lock, "SeqKit 2.13.0")

    determinism = read_tsv(frozen / "determinism-audit.tsv")
    audit.add("Determinism", "all-pass", bool(determinism) and all(row["Status"] == "PASS" for row in determinism), determinism)
    audit.add("Determinism", "plot-seed", any(row["Component"] == "plots" and row["Seed"] == "20260755" for row in determinism), determinism)

    audit_figures(args.figure_dir.resolve(), audit, FIGURES)
    audit_chapter(args.chapter.resolve(), frozen, audit)
    return finish(
        article=55,
        audit=audit,
        output=args.output_dir.resolve(),
        payload={
            "frozen_dir": str(frozen),
            "chapter": str(args.chapter.resolve()),
            "figures": list(FIGURES),
        },
    )


if __name__ == "__main__":
    raise SystemExit(main())
