#!/usr/bin/env python3
"""Stream deterministic FASTQ prefixes from ENA without storing raw FASTQ.

This one-time provenance builder writes only per-read metrics, hashes, and
short display prefixes. It intentionally does not create a FASTQ file under
data/small, because raw FASTQ stays outside Git in this project.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
import math
import statistics
import urllib.request
from pathlib import Path
from typing import Iterator


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--sources", type=Path, required=True)
    parser.add_argument("--metrics", type=Path, required=True)
    parser.add_argument("--anatomy", type=Path, required=True)
    parser.add_argument("--summary", type=Path, required=True)
    parser.add_argument("--timeout", type=int, default=120)
    return parser.parse_args()


def sha256_text(value: str) -> str:
    return hashlib.sha256(value.encode("ascii")).hexdigest()


def normalized_read_id(header: str) -> str:
    token = header[1:].split()[0]
    if token.endswith("/1") or token.endswith("/2"):
        return token[:-2]
    return token


def read_fastq_prefix(
    url: str,
    records: int,
    timeout: int,
) -> Iterator[tuple[int, str, str, str, str]]:
    request = urllib.request.Request(
        url,
        headers={"User-Agent": "metagenomics-best-practices/0.1"},
    )
    with urllib.request.urlopen(request, timeout=timeout) as response:
        with gzip.GzipFile(fileobj=response) as compressed:
            with io.TextIOWrapper(
                compressed,
                encoding="ascii",
                newline="",
            ) as text:
                for index in range(1, records + 1):
                    lines = [text.readline().rstrip("\r\n") for _ in range(4)]
                    if any(line == "" for line in lines):
                        raise ValueError(
                            f"FASTQ ended before record {index}: {url}"
                        )
                    header, sequence, plus, quality = lines
                    if not header.startswith("@"):
                        raise ValueError(
                            f"Record {index} header does not start with @"
                        )
                    if not plus.startswith("+"):
                        raise ValueError(
                            f"Record {index} separator does not start with +"
                        )
                    if len(sequence) != len(quality):
                        raise ValueError(
                            f"Record {index} sequence/quality lengths differ"
                        )
                    if not set(sequence.upper()) <= set("ACGTN"):
                        raise ValueError(
                            f"Record {index} contains unexpected base symbols"
                        )
                    yield index, header, sequence, plus, quality


def metric_row(
    source: dict[str, str],
    index: int,
    header: str,
    sequence: str,
    quality: str,
) -> dict[str, object]:
    q_values = [ord(char) - 33 for char in quality]
    if min(q_values) < 0 or max(q_values) > 93:
        raise ValueError(
            f"Non-Phred+33-compatible character in {source['RunAccession']}"
        )
    error_probabilities = [10 ** (-q / 10) for q in q_values]
    mean_error = statistics.fmean(error_probabilities)
    error_probability_mean_q = (
        -10 * math.log10(mean_error) if mean_error > 0 else math.inf
    )
    sequence_upper = sequence.upper()
    called_bases = len(sequence_upper) - sequence_upper.count("N")
    gc = sequence_upper.count("G") + sequence_upper.count("C")
    return {
        "PlatformKey": source["PlatformKey"],
        "PlatformLabel": source["PlatformLabel"],
        "RunAccession": source["RunAccession"],
        "SampleAccession": source["SampleAccession"],
        "Mate": source["Mate"],
        "ReadIndex": index,
        "ReadID": normalized_read_id(header),
        "ReadLength": len(sequence),
        "GCPercent": 100 * gc / called_bases if called_bases else "",
        "NCount": sequence_upper.count("N"),
        "ArithmeticMeanBaseQ": statistics.fmean(q_values),
        "ErrorProbabilityMeanQ": error_probability_mean_q,
        "ExpectedErrors": sum(error_probabilities),
        "ExpectedAccuracyPercent": 100 * (1 - mean_error),
        "MinimumBaseQ": min(q_values),
        "MaximumBaseQ": max(q_values),
        "SequenceSHA256": sha256_text(sequence),
        "QualitySHA256": sha256_text(quality),
    }


def write_tsv(
    path: Path,
    rows: list[dict[str, object]],
    fieldnames: list[str],
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> int:
    args = parse_args()
    with args.sources.open(encoding="utf-8", newline="") as handle:
        sources = list(csv.DictReader(handle, delimiter="\t"))

    expected_rows = {
        ("Illumina", "R1"),
        ("Illumina", "R2"),
        ("ONT", "SE"),
        ("PacBio", "SE"),
    }
    observed_rows = {
        (row["PlatformKey"], row["Mate"]) for row in sources
    }
    if observed_rows != expected_rows:
        raise ValueError(
            f"Unexpected source rows: {sorted(observed_rows)}"
        )
    if {row["SampleAccession"] for row in sources} != {"SAMEA14435832"}:
        raise ValueError("The platform runs are not linked to one sample")

    metric_rows: list[dict[str, object]] = []
    anatomy_rows: list[dict[str, object]] = []
    read_ids: dict[tuple[str, str], list[str]] = {}
    stream_audit: list[dict[str, object]] = []

    for source in sources:
        record_count = int(source["PrefixRecords"])
        source_ids: list[str] = []
        first_header = ""
        first_sequence = ""
        first_plus = ""
        first_quality = ""
        observed = 0
        for index, header, sequence, plus, quality in read_fastq_prefix(
            source["HTTPSURL"],
            record_count,
            args.timeout,
        ):
            observed += 1
            source_ids.append(normalized_read_id(header))
            if index == 1:
                first_header = header
                first_sequence = sequence
                first_plus = plus
                first_quality = quality
            if source["IncludeInMetrics"].upper() == "TRUE":
                metric_rows.append(
                    metric_row(
                        source,
                        index,
                        header,
                        sequence,
                        quality,
                    )
                )
        if observed != record_count:
            raise ValueError(
                f"Expected {record_count} records, observed {observed}"
            )
        if len(source_ids) != len(set(source_ids)):
            raise ValueError(
                f"Duplicate IDs in prefix for {source['RunAccession']} "
                f"{source['Mate']}"
            )
        read_ids[(source["PlatformKey"], source["Mate"])] = source_ids
        first_q = [ord(char) - 33 for char in first_quality[:12]]
        anatomy_rows.append(
            {
                "PlatformKey": source["PlatformKey"],
                "PlatformLabel": source["PlatformLabel"],
                "RunAccession": source["RunAccession"],
                "Mate": source["Mate"],
                "Header": first_header,
                "SequencePrefix60": first_sequence[:60],
                "Separator": first_plus,
                "QualityPrefix60": first_quality[:60],
                "DecodedPhredFirst12": ",".join(map(str, first_q)),
                "FullReadLength": len(first_sequence),
                "FullSequenceSHA256": sha256_text(first_sequence),
                "FullQualitySHA256": sha256_text(first_quality),
            }
        )
        stream_audit.append(
            {
                "platform": source["PlatformKey"],
                "mate": source["Mate"],
                "run_accession": source["RunAccession"],
                "records": observed,
                "first_read_id": source_ids[0],
                "last_read_id": source_ids[-1],
            }
        )

    illumina_r1 = read_ids[("Illumina", "R1")]
    illumina_r2 = read_ids[("Illumina", "R2")]
    paired_ids_identical = illumina_r1 == illumina_r2
    if not paired_ids_identical:
        mismatch = next(
            index
            for index, pair in enumerate(zip(illumina_r1, illumina_r2), 1)
            if pair[0] != pair[1]
        )
        raise ValueError(f"Illumina mates diverge at prefix record {mismatch}")

    metric_fields = [
        "PlatformKey",
        "PlatformLabel",
        "RunAccession",
        "SampleAccession",
        "Mate",
        "ReadIndex",
        "ReadID",
        "ReadLength",
        "GCPercent",
        "NCount",
        "ArithmeticMeanBaseQ",
        "ErrorProbabilityMeanQ",
        "ExpectedErrors",
        "ExpectedAccuracyPercent",
        "MinimumBaseQ",
        "MaximumBaseQ",
        "SequenceSHA256",
        "QualitySHA256",
    ]
    anatomy_fields = [
        "PlatformKey",
        "PlatformLabel",
        "RunAccession",
        "Mate",
        "Header",
        "SequencePrefix60",
        "Separator",
        "QualityPrefix60",
        "DecodedPhredFirst12",
        "FullReadLength",
        "FullSequenceSHA256",
        "FullQualitySHA256",
    ]
    write_tsv(args.metrics, metric_rows, metric_fields)
    write_tsv(args.anatomy, anatomy_rows, anatomy_fields)

    summary = {
        "source_project": "PRJEB52977",
        "sample_accession": "SAMEA14435832",
        "selection": (
            "first 5000 complete FASTQ records from each "
            "checksum-identified ENA file"
        ),
        "raw_fastq_stored": False,
        "metrics_rows": len(metric_rows),
        "anatomy_rows": len(anatomy_rows),
        "illumina_prefix_mates_synchronized": paired_ids_identical,
        "streams": stream_audit,
        "metrics_sha256": file_sha256(args.metrics),
        "anatomy_sha256": file_sha256(args.anatomy),
    }
    args.summary.parent.mkdir(parents=True, exist_ok=True)
    args.summary.write_text(
        json.dumps(summary, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
