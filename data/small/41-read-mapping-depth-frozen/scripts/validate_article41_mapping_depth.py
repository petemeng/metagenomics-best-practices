#!/usr/bin/env python3
"""Fail-closed validation for Article 41 mapping and contig-depth evidence."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import re
from dataclasses import dataclass, field
from pathlib import Path

from PIL import Image


FIGURES = ("41-mapping-fate", "41-depth-concordance", "41-depth-breadth")
REQUIRED_SECTIONS = (
    "这一步对应论文里的哪张图", "理论", "准备工作", "可复制代码",
    "审计与升级", "出版级美化", "常见坑", "这段 Methods 怎么写",
    "换成你自己的数据怎么做", "参考",
)


def read_tsv(path: Path) -> list[dict[str, str]]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def close(left: float, right: float, tolerance: float = 1e-9) -> bool:
    return math.isclose(left, right, rel_tol=tolerance, abs_tol=tolerance)


@dataclass
class Audit:
    rows: list[dict[str, object]] = field(default_factory=list)

    def add(self, category: str, check: str, status: bool, detail: object) -> None:
        self.rows.append(
            {
                "Category": category,
                "CheckID": check,
                "Status": "PASS" if status else "FAIL",
                "Detail": detail if isinstance(detail, str) else json.dumps(detail, ensure_ascii=False, sort_keys=True),
            }
        )

    @property
    def passed(self) -> int:
        return sum(row["Status"] == "PASS" for row in self.rows)

    @property
    def failed(self) -> int:
        return sum(row["Status"] == "FAIL" for row in self.rows)


def audit_checksums(frozen: Path, audit: Audit) -> None:
    manifest = frozen / "file-checksums.sha256"
    audit.add("Checksum", "manifest-exists", manifest.is_file(), str(manifest))
    expected: dict[str, str] = {}
    for line in manifest.read_text(encoding="utf-8").splitlines():
        digest, relative = line.split("  ", 1)
        expected[relative] = digest
        path = frozen / relative
        observed = sha256(path) if path.is_file() else "MISSING"
        audit.add("Checksum", relative, observed == digest, observed)
    actual = {
        path.relative_to(frozen).as_posix()
        for path in frozen.rglob("*")
        if path.is_file() and path.name != manifest.name
    }
    audit.add("Checksum", "manifest-complete", set(expected) == actual, {"listed": len(expected), "actual": len(actual)})


def audit_chapter(chapter: Path, audit: Audit) -> None:
    text = chapter.read_text(encoding="utf-8")
    expectations = {
        "draft-false": "draft: false" in text,
        "eval-false": re.search(r"execute:\s*\n\s+eval:\s+false", text) is not None,
        "freeze-auto": "freeze: auto" in text,
        "bibliography": "bibliography: ../references.bib" in text,
        "seed": "20260741" in text,
        "inline-theme": all(token in text for token in ("pal_pub <-", "theme_pub <-", "save_pub <-")),
        "resource-contract": all(token in text for token in ("RAM", "CPU", "磁盘", "耗时")),
        "methods-template": "Methods template" in text,
        "results-template": "Results template" in text,
        "real-study": all(token in text for token in ("PRJEB52977", "ERR9765746", "ERR9765747")),
        "versions": all(token in text for token in ("Bowtie2 2.5.5", "SAMtools 1.23.1", "MetaBAT2")),
        "frozen-input": "data/small/41-read-mapping-depth-frozen" in text,
        "three-figures": all(f"../figures/{stem}.png" in text for stem in FIGURES),
        "vector-raster-export": all(token in text for token in (".pdf", ".png", ".tiff", "compression = \"lzw\"")),
        "no-source-theme": 'source("R/theme_pub.R")' not in text and "source('R/theme_pub.R')" not in text,
        "no-placeholders": re.search(r"__[A-Z0-9_]+__|TODO|TBD", text) is None,
        "no-meta-prose": not any(
            token in text for token in
            ("本篇可独立", "本文可独立", "全系列约定", "接口只学一次", "作者代码通常长这样", "（即本文）", "无头服务器")
        ),
    }
    for section in REQUIRED_SECTIONS:
        expectations[f"section-{section}"] = re.search(rf"(?m)^##\s+{re.escape(section)}", text) is not None
    for check, status in expectations.items():
        audit.add("Chapter", check, status, check)


def audit_figures(figure_dir: Path, audit: Audit) -> None:
    for stem in FIGURES:
        for extension in ("pdf", "png", "tiff"):
            path = figure_dir / f"{stem}.{extension}"
            exists = path.is_file() and path.stat().st_size > 0
            audit.add("Figure", f"{stem}-{extension}-exists", exists, str(path))
            if not exists or extension == "pdf":
                continue
            with Image.open(path) as image:
                dpi = tuple(float(value) for value in image.info.get("dpi", (0, 0)))
                width, height = image.size
                audit.add("Figure", f"{stem}-{extension}-dimensions", width >= 1400 and height >= 900, f"{width}x{height}")
                audit.add("Figure", f"{stem}-{extension}-dpi", min(dpi) >= 340, dpi)
                audit.add("Figure", f"{stem}-{extension}-mode", image.mode in {"RGB", "RGBA"}, image.mode)
                if extension == "tiff":
                    compression = image.info.get("compression", "")
                    audit.add("Figure", f"{stem}-tiff-lzw", str(compression).lower() in {"tiff_lzw", "lzw"}, compression)


def audit_science(frozen: Path, audit: Audit) -> dict[str, object]:
    summary = json.loads((frozen / "run-summary.json").read_text(encoding="utf-8"))
    contract = json.loads((frozen / "run-contract.json").read_text(encoding="utf-8"))
    assembly = read_tsv(frozen / "assembly-summary.tsv")
    mapping = read_tsv(frozen / "mapping-summary.tsv")
    fate = read_tsv(frozen / "mapping-fate-long.tsv")
    depths = read_tsv(frozen / "depth-summary.tsv")
    correlation = read_tsv(frozen / "depth-correlation.tsv")
    wide = read_tsv(frozen / "contig-depth-wide.tsv.gz")
    long = read_tsv(frozen / "contig-depth-long.tsv.gz")
    tools = read_tsv(frozen / "tool-versions.tsv")
    inputs = read_tsv(frozen / "input-audit.tsv")
    commands = read_tsv(frozen / "command-log.tsv")
    resources = read_tsv(frozen / "resource-summary.tsv")
    raw = read_tsv(frozen / "raw/jgi-depth.tsv")
    pairs = read_tsv(frozen / "raw/paired-contigs.tsv")

    audit.add("Identity", "article", summary.get("article") == 41 and contract.get("article") == 41, {"summary": summary.get("article"), "contract": contract.get("article")})
    audit.add("Identity", "seed", contract.get("seed") == 20260741, contract.get("seed"))
    audit.add("Identity", "coordinate-system", "contigs >=1000 bp" in contract.get("coordinate_system", ""), contract.get("coordinate_system"))
    audit.add("Identity", "alignment-preset", contract["alignment"]["preset"] == "--very-sensitive", contract["alignment"])
    audit.add("Identity", "alignment-flags", contract["alignment"]["bam_exclude_flag"] == 3588, contract["alignment"])
    audit.add("Identity", "jgi-identity", contract["jgi_depth"]["minimum_end_to_end_identity_percent"] == 97, contract["jgi_depth"])
    audit.add("Identity", "jgi-mapq", contract["jgi_depth"]["minimum_mapq"] == 0, contract["jgi_depth"])
    audit.add("Identity", "jgi-edge-rule", contract["jgi_depth"]["edge_bases_excluded"] is True, contract["jgi_depth"])

    observed_tools = {row["Tool"]: row["Version"] for row in tools}
    audit.add("Version", "bowtie2", observed_tools.get("Bowtie2") == "2.5.5", observed_tools)
    audit.add("Version", "bowtie2-build", observed_tools.get("Bowtie2-build") == "2.5.5", observed_tools)
    audit.add("Version", "samtools", observed_tools.get("SAMtools") == "1.23.1", observed_tools)
    audit.add("Version", "jgi-depth", observed_tools.get("MetaBAT2 depth") == "2.18", observed_tools)
    audit.add("Input", "all-inputs-pass", len(inputs) == 6 and all(row["Status"] == "PASS" for row in inputs), len(inputs))
    audit.add("Input", "sha256-shaped", all(re.fullmatch(r"[0-9a-f]{64}", row["SHA256"]) for row in inputs), [row["SHA256"] for row in inputs])
    audit.add("Input", "source-runs", {"ERR9765746", "ERR9765747"}.issubset({row["Path"].split("/")[-1].split("_")[0] for row in inputs}), [row["Path"] for row in inputs])

    row = assembly[0]
    assembly_expected = {"Contigs": 18354, "TotalBp": 84811518, "MinimumBp": 1000, "MaximumBp": 1064594, "N50Bp": 15928}
    for key, expected in assembly_expected.items():
        audit.add("Assembly", key, int(float(row[key])) == expected, row[key])
    audit.add("Assembly", "gc-percent", close(float(row["GCPct"]), 48.86623064570074), row["GCPct"])

    expected_samples = {
        "MOCK1": {"pairs": 1_999_853, "zero": 730_631, "once": 1_250_638, "multiple": 18_584, "overall": 87.48, "breadth": 95.29653743492717, "jgi": 5.812685370068148},
        "MOCK2": {"pairs": 1_999_888, "zero": 650_532, "once": 1_329_470, "multiple": 19_886, "overall": 87.25, "breadth": 95.09637830088126, "jgi": 5.786746469696133},
    }
    mapping_by_sample = {row["Sample"]: row for row in mapping}
    depth_by_sample = {row["Sample"]: row for row in depths}
    for sample, expected in expected_samples.items():
        observed = mapping_by_sample[sample]
        counts = [int(observed[key]) for key in ("ConcordantZero", "ConcordantOnce", "ConcordantMultiple")]
        audit.add("Mapping", f"{sample}-pairs", int(observed["ReadPairs"]) == expected["pairs"], observed["ReadPairs"])
        audit.add("Mapping", f"{sample}-pair-conservation", sum(counts) == int(observed["ReadPairs"]), counts)
        audit.add("Mapping", f"{sample}-zero", counts[0] == expected["zero"], counts[0])
        audit.add("Mapping", f"{sample}-once", counts[1] == expected["once"], counts[1])
        audit.add("Mapping", f"{sample}-multiple", counts[2] == expected["multiple"], counts[2])
        audit.add("Mapping", f"{sample}-overall", close(float(observed["OverallAlignmentPct"]), expected["overall"]), observed["OverallAlignmentPct"])
        audit.add("Mapping", f"{sample}-bam-mapped", int(observed["BAMRecords"]) == int(observed["MappedRecords"]) > 0, {"bam": observed["BAMRecords"], "mapped": observed["MappedRecords"]})
        audit.add("Mapping", f"{sample}-proper-positive", 0 < int(observed["ProperlyPairedRecords"]) <= int(observed["BAMRecords"]), observed["ProperlyPairedRecords"])
        audit.add("Mapping", f"{sample}-breadth", close(float(observed["AssemblyBreadthPct"]), expected["breadth"]), observed["AssemblyBreadthPct"])
        depth_row = depth_by_sample[sample]
        audit.add("Depth", f"{sample}-contigs", int(depth_row["Contigs"]) == 18354, depth_row["Contigs"])
        audit.add("Depth", f"{sample}-positive", int(depth_row["PositiveDepthContigs"]) == (18353 if sample == "MOCK1" else 18354), depth_row["PositiveDepthContigs"])
        audit.add("Depth", f"{sample}-zero", int(depth_row["ZeroDepthContigs"]) == (1 if sample == "MOCK1" else 0), depth_row["ZeroDepthContigs"])
        audit.add("Depth", f"{sample}-jgi-weighted", close(float(depth_row["LengthWeightedJgiMeanDepth"]), expected["jgi"]), depth_row["LengthWeightedJgiMeanDepth"])
        audit.add("Depth", f"{sample}-quantile-order", 0 < float(depth_row["MedianPositiveDepth"]) <= float(depth_row["Q90PositiveDepth"]) <= float(depth_row["Q99PositiveDepth"]), depth_row)

    fate_totals: dict[str, int] = {}
    fate_percent: dict[str, float] = {}
    for row in fate:
        fate_totals[row["Sample"]] = fate_totals.get(row["Sample"], 0) + int(row["ReadPairs"])
        fate_percent[row["Sample"]] = fate_percent.get(row["Sample"], 0.0) + float(row["Percent"])
    for sample in expected_samples:
        audit.add("Ledger", f"{sample}-fate-rows", sum(row["Sample"] == sample for row in fate) == 3, sample)
        audit.add("Ledger", f"{sample}-fate-count-close", fate_totals[sample] == expected_samples[sample]["pairs"], fate_totals[sample])
        audit.add("Ledger", f"{sample}-fate-percent-close", close(fate_percent[sample], 100.0, 1e-7), fate_percent[sample])

    names = [row["Contig"] for row in wide]
    audit.add("Matrix", "wide-row-count", len(wide) == 18354, len(wide))
    audit.add("Matrix", "wide-unique-contigs", len(set(names)) == len(names), len(set(names)))
    audit.add("Matrix", "wide-total-bp", sum(int(row["LengthBp"]) for row in wide) == 84811518, sum(int(row["LengthBp"]) for row in wide))
    audit.add("Matrix", "wide-minimum-bp", min(int(row["LengthBp"]) for row in wide) == 1000, min(int(row["LengthBp"]) for row in wide))
    audit.add("Matrix", "long-row-count", len(long) == 2 * len(wide), len(long))
    audit.add("Matrix", "long-samples", {row["Sample"] for row in long} == {"MOCK1", "MOCK2"}, sorted({row["Sample"] for row in long}))
    audit.add("Matrix", "long-coordinate-product", len({(row["Contig"], row["Sample"]) for row in long}) == len(long), len(long))
    audit.add("Matrix", "depth-nonnegative", all(float(row["JgiMeanDepth"]) >= 0 for row in long), min(float(row["JgiMeanDepth"]) for row in long))
    audit.add("Matrix", "breadth-bounded", all(0 <= float(row["BreadthPct"]) <= 100 for row in long), {"min": min(float(row["BreadthPct"]) for row in long), "max": max(float(row["BreadthPct"]) for row in long)})
    audit.add("Matrix", "detected-both", sum(int(row["DetectedInBoth"]) for row in wide) == 18353, sum(int(row["DetectedInBoth"]) for row in wide))

    corr = correlation[0]
    audit.add("Correlation", "coordinate-set", corr["CoordinateSet"] == "All contigs >=1000 bp", corr["CoordinateSet"])
    audit.add("Correlation", "contigs", int(corr["Contigs"]) == 18354, corr["Contigs"])
    audit.add("Correlation", "detected-both", int(corr["DetectedInBoth"]) == 18353, corr["DetectedInBoth"])
    audit.add("Correlation", "pearson", close(float(corr["PearsonLog1p"]), 0.7916776134046771), corr["PearsonLog1p"])
    audit.add("Correlation", "spearman", close(float(corr["SpearmanLog1p"]), 0.48337113278398963), corr["SpearmanLog1p"])
    audit.add("Correlation", "coefficients-bounded", all(-1 <= float(corr[key]) <= 1 for key in ("PearsonLog1p", "SpearmanLog1p")), corr)

    audit.add("Raw", "jgi-row-count", len(raw) == 18354, len(raw))
    audit.add("Raw", "jgi-coordinate-identity", {row["contigName"] for row in raw} == set(names), len(raw))
    audit.add("Raw", "jgi-length-identity", all(int(row["contigLen"]) == int(wide_row["LengthBp"]) for row, wide_row in zip(raw, wide)), "rowwise")
    audit.add("Raw", "paired-linkage-nonempty", len(pairs) > 18354, len(pairs))
    audit.add("Raw", "paired-linkage-coverage-positive", all(float(row["AvgCoverage"]) >= 0 for row in pairs), len(pairs))

    audit.add("Execution", "commands", len(commands) == 15, len(commands))
    audit.add("Execution", "command-return-codes", all(int(row["ReturnCode"]) == 0 for row in commands), [row["ReturnCode"] for row in commands])
    audit.add("Execution", "map-seed", all("--seed 20260741" in row["Command"] for row in commands if row["Label"].startswith("map-")), "map commands")
    audit.add("Execution", "jgi-policy", any("--percentIdentity 97 --minMapQual 0" in row["Command"] for row in commands), "jgi command")
    audit.add("Execution", "resource-count", len(resources) == len(commands), {"resources": len(resources), "commands": len(commands)})
    audit.add("Execution", "resource-exit-status", all(int(row["ExitStatus"]) == 0 for row in resources), [row["ExitStatus"] for row in resources])
    audit.add("Execution", "resource-wall-positive", all(float(row["WallSeconds"]) >= 0 for row in resources), min(float(row["WallSeconds"]) for row in resources))
    audit.add("Execution", "resource-ram-positive", all(float(row["PeakRAMGiB"]) > 0 for row in resources), min(float(row["PeakRAMGiB"]) for row in resources))

    ledgers = summary.get("ledger_checks", {})
    audit.add("Boundary", "pair-conservation-declared", ledgers.get("bowtie_pair_conservation") is True, ledgers)
    audit.add("Boundary", "coverage-coordinate-declared", ledgers.get("coverage_coordinate_identity") is True, ledgers)
    audit.add("Boundary", "jgi-coordinate-declared", ledgers.get("jgi_coordinate_identity") is True, ledgers)
    audit.add("Boundary", "resource-status-declared", ledgers.get("resource_exit_status") is True, ledgers)
    return summary


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--frozen-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--figure-dir", type=Path, required=True)
    parser.add_argument("--chapter", type=Path, required=True)
    args = parser.parse_args()
    frozen = args.frozen_dir.resolve()
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    audit = Audit()
    audit_checksums(frozen, audit)
    summary = audit_science(frozen, audit)
    audit_chapter(args.chapter.resolve(), audit)
    audit_figures(args.figure_dir.resolve(), audit)
    write_tsv(output / "validation-checks.tsv", audit.rows)
    result = {
        "article": 41,
        "status": "passed" if audit.failed == 0 else "failed",
        "checks_passed": audit.passed,
        "checks_failed": audit.failed,
        "assembly": summary["assembly"],
        "samples": summary["samples"],
        "correlation": summary["correlation"],
    }
    (output / "validation-summary.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    (output / "validation.log").write_text(
        f"Article 41 validation: {result['status'].upper()}\nPASS={audit.passed}\nFAIL={audit.failed}\n"
        + "\n".join(f"{row['Status']}\t{row['Category']}\t{row['CheckID']}\t{row['Detail']}" for row in audit.rows)
        + "\n",
        encoding="utf-8",
    )
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0 if audit.failed == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
