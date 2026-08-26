#!/usr/bin/env python3
"""Stream deterministic paired FASTQ controls for metagenomics Article 14.

Complete ENA archives remain remote. The script writes only the requested
complete-record prefixes into Git-ignored scratch using deterministic gzip
metadata. Human read identifiers are never written to the frozen evidence:
only SHA-256 digests and aggregate sequence metrics are recorded.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
import statistics
import urllib.request
from collections import defaultdict
from pathlib import Path
from typing import BinaryIO, Iterator, TextIO


USER_AGENT = "metagenomics-best-practices/article14"
EXPECTED_CONTROLS = ("human_positive", "mock_retention")
EXPECTED_MATES = ("R1", "R2")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--timeout", type=int, default=300)
    return parser.parse_args()


def sha256_bytes(value: str) -> str:
    return hashlib.sha256(value.encode("ascii")).hexdigest()


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def normalized_read_id(header: str) -> str:
    token = header[1:].split()[0]
    if token.endswith("/1") or token.endswith("/2"):
        token = token[:-2]
    return token


def fastq_records(
    handle: TextIO,
    source: str,
) -> Iterator[tuple[str, str, str, str]]:
    record_index = 0
    while True:
        header = handle.readline()
        if header == "":
            return
        sequence = handle.readline()
        plus = handle.readline()
        quality = handle.readline()
        record_index += 1
        if "" in (sequence, plus, quality):
            raise ValueError(
                f"Truncated FASTQ record {record_index} in {source}"
            )
        values = [
            value.rstrip("\r\n")
            for value in (header, sequence, plus, quality)
        ]
        header_text, sequence_text, plus_text, quality_text = values
        if not header_text.startswith("@"):
            raise ValueError(
                f"Record {record_index} header does not start with @ in {source}"
            )
        if not plus_text.startswith("+"):
            raise ValueError(
                f"Record {record_index} separator does not start with + in {source}"
            )
        if len(sequence_text) != len(quality_text):
            raise ValueError(
                f"Record {record_index} sequence/quality lengths differ in {source}"
            )
        if not set(sequence_text.upper()) <= set("ACGTN"):
            raise ValueError(
                f"Record {record_index} has unexpected base symbols in {source}"
            )
        q_values = [ord(character) - 33 for character in quality_text]
        if min(q_values) < 0 or max(q_values) > 93:
            raise ValueError(
                f"Record {record_index} is not Phred+33 compatible in {source}"
            )
        yield header_text, sequence_text, plus_text, quality_text


def open_remote_fastq(
    url: str,
    timeout: int,
) -> tuple[BinaryIO, gzip.GzipFile, TextIO, str | None]:
    request = urllib.request.Request(
        url,
        headers={"User-Agent": USER_AGENT},
    )
    response = urllib.request.urlopen(request, timeout=timeout)
    content_length = response.headers.get("Content-Length")
    compressed = gzip.GzipFile(fileobj=response)
    text = io.TextIOWrapper(compressed, encoding="ascii", newline="")
    return response, compressed, text, content_length


def write_prefix(
    row: dict[str, str],
    output_path: Path,
    timeout: int,
) -> tuple[dict[str, object], list[str]]:
    requested_records = int(row["PrefixRecords"])
    output_path.parent.mkdir(parents=True, exist_ok=True)
    read_ids: list[str] = []
    uncompressed_digest = hashlib.sha256()
    sequence_digest = hashlib.sha256()
    length_values: list[int] = []
    total_bases = 0
    q30_bases = 0
    n_bases = 0
    first_id = ""
    last_id = ""
    observed_records = 0

    response, compressed, text, content_length = open_remote_fastq(
        row["HTTPSURL"],
        timeout,
    )
    try:
        with output_path.open("wb") as raw_output:
            with gzip.GzipFile(
                filename="",
                mode="wb",
                fileobj=raw_output,
                compresslevel=6,
                mtime=0,
            ) as gzip_output:
                for record in fastq_records(text, row["HTTPSURL"]):
                    header, sequence, plus, quality = record
                    observed_records += 1
                    read_id = normalized_read_id(header)
                    read_ids.append(read_id)
                    if observed_records == 1:
                        first_id = read_id
                    last_id = read_id
                    length_values.append(len(sequence))
                    total_bases += len(sequence)
                    q30_bases += sum(
                        ord(character) - 33 >= 30 for character in quality
                    )
                    n_bases += sequence.upper().count("N")
                    sequence_digest.update(sequence.upper().encode("ascii"))
                    sequence_digest.update(b"\n")
                    serialized = (
                        f"{header}\n{sequence}\n{plus}\n{quality}\n"
                    ).encode("ascii")
                    uncompressed_digest.update(serialized)
                    gzip_output.write(serialized)
                    if observed_records == requested_records:
                        break
    finally:
        text.close()
        compressed.close()
        response.close()

    if observed_records != requested_records:
        raise ValueError(
            f"{row['ControlID']} {row['Mate']} expected "
            f"{requested_records} records, observed {observed_records}"
        )
    if len(read_ids) != len(set(read_ids)):
        raise ValueError(
            f"Duplicate normalized IDs in {row['ControlID']} {row['Mate']}"
        )

    remote_size_matches = (
        content_length is not None
        and int(content_length) == int(row["ENABytes"])
    )
    summary: dict[str, object] = {
        "control_id": row["ControlID"],
        "expected_class": row["ExpectedClass"],
        "project_accession": row["ProjectAccession"],
        "sample_accession": row["SampleAccession"],
        "run_accession": row["RunAccession"],
        "mate": row["Mate"],
        "records": observed_records,
        "first_read_id_sha256": sha256_bytes(first_id),
        "last_read_id_sha256": sha256_bytes(last_id),
        "total_bases": total_bases,
        "minimum_read_length": min(length_values),
        "median_read_length": statistics.median(length_values),
        "maximum_read_length": max(length_values),
        "q30_base_fraction": q30_bases / total_bases,
        "n_base_fraction": n_bases / total_bases,
        "sequence_only_sha256": sequence_digest.hexdigest(),
        "uncompressed_fastq_sha256": uncompressed_digest.hexdigest(),
        "compressed_fastq_sha256": file_sha256(output_path),
        "compressed_subset_bytes": output_path.stat().st_size,
        "archive_reported_md5": row["ENAReportedMD5"],
        "archive_reported_bytes": int(row["ENABytes"]),
        "http_content_length": (
            int(content_length) if content_length is not None else None
        ),
        "http_content_length_matches_archive": remote_size_matches,
        "complete_archive_md5_verified": False,
        "output_file": output_path.name,
    }
    return summary, read_ids


def validate_manifest(
    rows: list[dict[str, str]],
) -> dict[str, list[dict[str, str]]]:
    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        grouped[row["ControlID"]].append(row)
    if tuple(grouped) != EXPECTED_CONTROLS:
        raise ValueError(
            f"Expected ordered controls {EXPECTED_CONTROLS}, observed {tuple(grouped)}"
        )
    for control_id, control_rows in grouped.items():
        if [row["Mate"] for row in control_rows] != list(EXPECTED_MATES):
            raise ValueError(
                f"{control_id} must contain ordered R1 and R2 rows"
            )
        for field in (
            "ExpectedClass",
            "ProjectAccession",
            "SampleAccession",
            "RunAccession",
            "Layout",
            "PrefixRecords",
        ):
            if len({row[field] for row in control_rows}) != 1:
                raise ValueError(f"{control_id} mates disagree for {field}")
        if control_rows[0]["Layout"] != "PAIRED":
            raise ValueError(f"{control_id} is not PAIRED")
    return dict(grouped)


def main() -> int:
    args = parse_args()
    with args.manifest.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    grouped = validate_manifest(rows)
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    all_summaries: list[dict[str, object]] = []
    control_summaries: dict[str, dict[str, object]] = {}
    for control_id, control_rows in grouped.items():
        ids_by_mate: dict[str, list[str]] = {}
        mate_summaries: list[dict[str, object]] = []
        for row in control_rows:
            output_path = output_dir / (
                f"{row['RunAccession']}_prefix"
                f"{int(row['PrefixRecords']) // 1000}k_"
                f"{row['Mate']}.fastq.gz"
            )
            summary, read_ids = write_prefix(
                row,
                output_path,
                args.timeout,
            )
            mate_summaries.append(summary)
            all_summaries.append(summary)
            ids_by_mate[row["Mate"]] = read_ids

        mates_synchronized = ids_by_mate["R1"] == ids_by_mate["R2"]
        if not mates_synchronized:
            mismatch = next(
                (
                    index
                    for index, (read1_id, read2_id) in enumerate(
                        zip(
                            ids_by_mate["R1"],
                            ids_by_mate["R2"],
                        ),
                        start=1,
                    )
                    if read1_id != read2_id
                ),
                None,
            )
            raise ValueError(
                f"{control_id} mates diverge at record {mismatch}"
            )
        pair_id_digest = hashlib.sha256()
        for read_id in ids_by_mate["R1"]:
            pair_id_digest.update(read_id.encode("ascii"))
            pair_id_digest.update(b"\n")
        control_summaries[control_id] = {
            "expected_class": control_rows[0]["ExpectedClass"],
            "project_accession": control_rows[0]["ProjectAccession"],
            "sample_accession": control_rows[0]["SampleAccession"],
            "run_accession": control_rows[0]["RunAccession"],
            "input_pairs": int(control_rows[0]["PrefixRecords"]),
            "selection_rule": control_rows[0]["SelectionRule"],
            "mates_synchronized": True,
            "normalized_pair_id_sha256": pair_id_digest.hexdigest(),
            "mate_summaries": mate_summaries,
        }

    payload = {
        "status": "passed",
        "privacy_boundary": "aggregate_only_no_human_sequences_or_read_ids_frozen",
        "complete_archive_md5_verified": False,
        "controls": control_summaries,
        "total_control_pairs": sum(
            int(value["input_pairs"]) for value in control_summaries.values()
        ),
        "mate_files": all_summaries,
    }
    summary_path = output_dir / "controls-summary.json"
    summary_path.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    print(
        json.dumps(
            {
                "status": "passed",
                "controls": list(control_summaries),
                "pairs_per_control": {
                    key: value["input_pairs"]
                    for key, value in control_summaries.items()
                },
                "summary": str(summary_path),
            },
            ensure_ascii=False,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
