#!/usr/bin/env python3
"""Create compact mapping/depth evidence and conservation audits for Article 41."""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
import statistics
from pathlib import Path

from article41_44_utils import (
    dump_json,
    fasta_summary,
    parse_time,
    read_tsv,
    write_tsv,
)


SAMPLES = ("MOCK1", "MOCK2")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--work-dir", type=Path, required=True)
    return parser.parse_args()


def parse_bowtie(path: Path) -> dict[str, object]:
    text = path.read_text(encoding="utf-8", errors="replace")
    patterns = {
        "ReadPairs": r"(?m)^\s*(\d+) reads; of these:",
        "ConcordantZero": r"(?m)^\s*(\d+) \([^\n]+ aligned concordantly 0 times$",
        "ConcordantOnce": r"(?m)^\s*(\d+) \([^\n]+ aligned concordantly exactly 1 time$",
        "ConcordantMultiple": r"(?m)^\s*(\d+) \([^\n]+ aligned concordantly >1 times$",
        "OverallAlignmentPct": r"(?m)^([0-9.]+)% overall alignment rate$",
    }
    values: dict[str, object] = {}
    for key, pattern in patterns.items():
        matches = re.findall(pattern, text)
        if not matches:
            raise ValueError(f"Cannot parse {key} from {path}")
        value = matches[-1]
        values[key] = float(value) if key == "OverallAlignmentPct" else int(value)
    if values["ConcordantZero"] + values["ConcordantOnce"] + values["ConcordantMultiple"] != values["ReadPairs"]:
        raise ValueError(f"Bowtie2 paired-read ledger does not conserve in {path}")
    return values


def parse_flagstat(path: Path) -> dict[str, int]:
    text = path.read_text(encoding="utf-8", errors="replace")
    values = {}
    for label, pattern in {
        "BAMRecords": r"(?m)^(\d+) \+ \d+ in total",
        "MappedRecords": r"(?m)^(\d+) \+ \d+ mapped \(",
        "ProperlyPairedRecords": r"(?m)^(\d+) \+ \d+ properly paired \(",
        "SingletonRecords": r"(?m)^(\d+) \+ \d+ singletons \(",
    }.items():
        match = re.search(pattern, text)
        if not match:
            raise ValueError(f"Cannot parse {label} from {path}")
        values[label] = int(match.group(1))
    return values


def quantile(values: list[float], probability: float) -> float:
    if not values:
        return math.nan
    ordered = sorted(values)
    position = (len(ordered) - 1) * probability
    lower = math.floor(position)
    upper = math.ceil(position)
    if lower == upper:
        return ordered[lower]
    weight = position - lower
    return ordered[lower] * (1 - weight) + ordered[upper] * weight


def pearson(x: list[float], y: list[float]) -> float:
    if len(x) != len(y) or len(x) < 2:
        return math.nan
    mx, my = statistics.mean(x), statistics.mean(y)
    numerator = sum((a - mx) * (b - my) for a, b in zip(x, y))
    denominator = math.sqrt(sum((a - mx) ** 2 for a in x) * sum((b - my) ** 2 for b in y))
    return numerator / denominator if denominator else math.nan


def ranks(values: list[float]) -> list[float]:
    order = sorted(range(len(values)), key=lambda index: values[index])
    output = [0.0] * len(values)
    cursor = 0
    while cursor < len(order):
        end = cursor + 1
        while end < len(order) and values[order[end]] == values[order[cursor]]:
            end += 1
        rank = (cursor + 1 + end) / 2
        for index in order[cursor:end]:
            output[index] = rank
        cursor = end
    return output


def read_coverage(path: Path) -> dict[str, dict[str, float]]:
    rows: dict[str, dict[str, float]] = {}
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError(path)
        first = reader.fieldnames[0]
        for raw in reader:
            contig = raw[first]
            rows[contig] = {
                "NumReads": int(raw["numreads"]),
                "CoveredBases": int(raw["covbases"]),
                "BreadthPct": float(raw["coverage"]),
                "SamtoolsMeanDepth": float(raw["meandepth"]),
                "MeanMapQ": float(raw["meanmapq"]),
            }
    return rows


def depth_columns(header: list[str]) -> tuple[dict[str, str], dict[str, str]]:
    means: dict[str, str] = {}
    variances: dict[str, str] = {}
    for column in header:
        base = Path(column.removesuffix("-var")).name
        for sample in SAMPLES:
            if base.startswith(sample) and column.endswith("-var"):
                variances[sample] = column
            elif base.startswith(sample) and not column.endswith("-var"):
                means[sample] = column
    if set(means) != set(SAMPLES) or set(variances) != set(SAMPLES):
        raise ValueError(f"Unexpected JGI depth header: {header}")
    return means, variances


def main() -> int:
    args = parse_args()
    work = args.work_dir.resolve()
    if not (work / ".article41-run-complete").is_file():
        raise FileNotFoundError("Run run_article41_mapping_depth.py first")
    summary_dir = work / "summary"
    summary_dir.mkdir(exist_ok=True)
    assembly_summary, contigs = fasta_summary(work / "inputs/megahit-coassembly.ge1000.fna")
    write_tsv(summary_dir / "assembly-summary.tsv", [assembly_summary])

    sample_contract = {row["Sample"]: row for row in read_tsv(work / "inputs/samples.tsv")}
    mapping_rows = []
    fate_rows = []
    coverage_by_sample: dict[str, dict[str, dict[str, float]]] = {}
    for sample in SAMPLES:
        bowtie = parse_bowtie(work / f"logs/map-{sample}.bowtie2.log")
        flags = parse_flagstat(work / f"depth/{sample}.flagstat.tsv")
        expected = int(sample_contract[sample]["ReadPairs"])
        if bowtie["ReadPairs"] != expected:
            raise ValueError(f"{sample} input-pair mismatch: {bowtie['ReadPairs']} != {expected}")
        coverage = read_coverage(work / f"depth/{sample}.coverage.tsv")
        coverage_by_sample[sample] = coverage
        if set(coverage) != set(contigs):
            raise ValueError(f"{sample} coverage coordinate mismatch")
        total_bp = sum(int(contigs[name]["LengthBp"]) for name in contigs)
        covered_bp = sum(int(row["CoveredBases"]) for row in coverage.values())
        depth_mass = sum(float(row["SamtoolsMeanDepth"]) * int(contigs[name]["LengthBp"]) for name, row in coverage.items())
        mapping_rows.append(
            {
                "Sample": sample,
                **bowtie,
                **flags,
                "AssemblyCoveredBp": covered_bp,
                "AssemblyBreadthPct": 100 * covered_bp / total_bp,
                "LengthWeightedSamtoolsMeanDepth": depth_mass / total_bp,
            }
        )
        for key, label in (
            ("ConcordantZero", "No concordant alignment"),
            ("ConcordantOnce", "Unique concordant alignment"),
            ("ConcordantMultiple", "Multiple concordant alignments"),
        ):
            count = int(bowtie[key])
            fate_rows.append(
                {
                    "Sample": sample,
                    "ReadFate": label,
                    "ReadPairs": count,
                    "Percent": 100 * count / expected,
                }
            )
    write_tsv(summary_dir / "mapping-summary.tsv", mapping_rows)
    write_tsv(summary_dir / "mapping-fate-long.tsv", fate_rows)

    jgi_path = work / "depth/jgi-depth.tsv"
    with jgi_path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError(jgi_path)
        means, variances = depth_columns(reader.fieldnames)
        jgi_rows = list(reader)
    if len(jgi_rows) != int(assembly_summary["Contigs"]):
        raise ValueError("JGI output lost assembly contigs")

    wide_rows = []
    long_rows = []
    for raw in jgi_rows:
        contig = raw["contigName"]
        if contig not in contigs:
            raise ValueError(f"JGI unknown contig: {contig}")
        length = int(raw["contigLen"])
        if length != int(contigs[contig]["LengthBp"]):
            raise ValueError(f"Length mismatch for {contig}")
        row: dict[str, object] = {
            "Contig": contig,
            "LengthBp": length,
            "GCPct": contigs[contig]["GCPct"],
            "TotalAverageDepth": float(raw["totalAvgDepth"]),
        }
        for sample in SAMPLES:
            depth = float(raw[means[sample]])
            variance = float(raw[variances[sample]])
            coverage = coverage_by_sample[sample][contig]
            row[f"{sample}MeanDepth"] = depth
            row[f"{sample}DepthVariance"] = variance
            row[f"{sample}BreadthPct"] = coverage["BreadthPct"]
            row[f"{sample}SamtoolsMeanDepth"] = coverage["SamtoolsMeanDepth"]
            row[f"{sample}MeanMapQ"] = coverage["MeanMapQ"]
            long_rows.append(
                {
                    "Contig": contig,
                    "Sample": sample,
                    "LengthBp": length,
                    "GCPct": contigs[contig]["GCPct"],
                    "JgiMeanDepth": depth,
                    "JgiDepthVariance": variance,
                    "BreadthPct": coverage["BreadthPct"],
                    "SamtoolsMeanDepth": coverage["SamtoolsMeanDepth"],
                    "MeanMapQ": coverage["MeanMapQ"],
                }
            )
        row["DetectedInBoth"] = int(row["MOCK1MeanDepth"] > 0 and row["MOCK2MeanDepth"] > 0)
        row["Log2DepthRatioMOCK2vsMOCK1"] = math.log2((float(row["MOCK2MeanDepth"]) + 0.01) / (float(row["MOCK1MeanDepth"]) + 0.01))
        wide_rows.append(row)
    write_tsv(summary_dir / "contig-depth-wide.tsv.gz", wide_rows)
    write_tsv(summary_dir / "contig-depth-long.tsv.gz", long_rows)

    depth_summary = []
    for sample in SAMPLES:
        values = [float(row[f"{sample}MeanDepth"]) for row in wide_rows]
        positive = [value for value in values if value > 0]
        weighted = sum(value * int(row["LengthBp"]) for value, row in zip(values, wide_rows)) / int(assembly_summary["TotalBp"])
        depth_summary.append(
            {
                "Sample": sample,
                "Contigs": len(values),
                "PositiveDepthContigs": len(positive),
                "ZeroDepthContigs": len(values) - len(positive),
                "PositiveDepthPct": 100 * len(positive) / len(values),
                "MedianPositiveDepth": statistics.median(positive),
                "Q90PositiveDepth": quantile(positive, 0.90),
                "Q99PositiveDepth": quantile(positive, 0.99),
                "LengthWeightedJgiMeanDepth": weighted,
            }
        )
    write_tsv(summary_dir / "depth-summary.tsv", depth_summary)
    x = [math.log1p(float(row["MOCK1MeanDepth"])) for row in wide_rows]
    y = [math.log1p(float(row["MOCK2MeanDepth"])) for row in wide_rows]
    both = [row for row in wide_rows if float(row["MOCK1MeanDepth"]) > 0 and float(row["MOCK2MeanDepth"]) > 0]
    write_tsv(
        summary_dir / "depth-correlation.tsv",
        [
            {
                "CoordinateSet": "All contigs >=1000 bp",
                "Contigs": len(wide_rows),
                "PearsonLog1p": pearson(x, y),
                "SpearmanLog1p": pearson(ranks(x), ranks(y)),
                "DetectedInBoth": len(both),
            }
        ],
    )
    resource_rows = [parse_time(path) for path in sorted((work / "logs").glob("*.time.txt"))]
    if any(int(row["ExitStatus"]) != 0 for row in resource_rows):
        raise RuntimeError("Non-zero command in resource logs")
    write_tsv(summary_dir / "resource-summary.tsv", resource_rows)
    dump_json(
        summary_dir / "run-summary.json",
        {
            "article": 41,
            "assembly": assembly_summary,
            "samples": mapping_rows,
            "depth": depth_summary,
            "correlation": read_tsv(summary_dir / "depth-correlation.tsv")[0],
            "ledger_checks": {
                "bowtie_pair_conservation": True,
                "coverage_coordinate_identity": True,
                "jgi_coordinate_identity": True,
                "resource_exit_status": True,
            },
        },
    )
    (work / ".article41-summary-complete").write_text("PASS\n", encoding="utf-8")
    print(json.dumps(json.loads((summary_dir / "run-summary.json").read_text()), indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
