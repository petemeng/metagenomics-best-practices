#!/usr/bin/env python3
"""Build a deterministic paired FASTQ prefix from checksum-identified ENA files.

The complete ENA files remain remote and their archive MD5 values are metadata
identifiers, not locally reverified checksums. Only the requested complete
record prefix is streamed. Output gzip files use mtime=0 so the same prefix has
stable compressed and uncompressed SHA-256 values.
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
from pathlib import Path
from typing import BinaryIO, Iterator, TextIO


USER_AGENT = "metagenomics-best-practices/article13"
EXPECTED_MATES = ("R1", "R2")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--timeout", type=int, default=180)
    return parser.parse_args()


def normalized_read_id(header: str) -> str:
    token = header[1:].split()[0]
    if token.endswith("/1") or token.endswith("/2"):
        token = token[:-2]
    return token


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def fastq_records(handle: TextIO, source: str) -> Iterator[tuple[str, str, str, str]]:
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
        lines = [
            value.rstrip("\r\n")
            for value in (header, sequence, plus, quality)
        ]
        header_text, sequence_text, plus_text, quality_text = lines
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
    ids: list[str] = []
    uncompressed_digest = hashlib.sha256()
    length_values: list[int] = []
    total_bases = 0
    total_q30_bases = 0
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
                    ids.append(read_id)
                    if observed_records == 1:
                        first_id = read_id
                    last_id = read_id
                    length_values.append(len(sequence))
                    total_bases += len(sequence)
                    total_q30_bases += sum(
                        ord(character) - 33 >= 30 for character in quality
                    )
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
            f"{row['Mate']} expected {requested_records} records, "
            f"observed {observed_records}"
        )
    if len(ids) != len(set(ids)):
        raise ValueError(f"Duplicate normalized IDs in {row['Mate']} prefix")

    remote_size_matches = (
        content_length is not None
        and int(content_length) == int(row["ENABytes"])
    )
    summary = {
        "mate": row["Mate"],
        "records": observed_records,
        "first_read_id": first_id,
        "last_read_id": last_id,
        "total_bases": total_bases,
        "minimum_read_length": min(length_values),
        "median_read_length": statistics.median(length_values),
        "maximum_read_length": max(length_values),
        "q30_base_fraction": total_q30_bases / total_bases,
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
    return summary, ids


def main() -> int:
    args = parse_args()
    with args.manifest.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if [row["Mate"] for row in rows] != list(EXPECTED_MATES):
        raise ValueError("Source manifest must contain ordered R1 and R2 rows")
    invariant_fields = (
        "ProjectAccession",
        "SampleAccession",
        "RunAccession",
        "Layout",
        "PrefixRecords",
    )
    for field in invariant_fields:
        if len({row[field] for row in rows}) != 1:
            raise ValueError(f"R1/R2 disagree for {field}")
    if rows[0]["Layout"] != "PAIRED":
        raise ValueError("Article 13 requires a PAIRED run")

    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    summaries: list[dict[str, object]] = []
    ids_by_mate: dict[str, list[str]] = {}
    for row in rows:
        output_path = output_dir / (
            f"{row['RunAccession']}_prefix"
            f"{int(row['PrefixRecords']) // 1000}k_{row['Mate']}.fastq.gz"
        )
        summary, read_ids = write_prefix(row, output_path, args.timeout)
        summaries.append(summary)
        ids_by_mate[row["Mate"]] = read_ids

    mates_synchronized = ids_by_mate["R1"] == ids_by_mate["R2"]
    if not mates_synchronized:
        mismatch = next(
            (
                index
                for index, (read1_id, read2_id) in enumerate(
                    zip(ids_by_mate["R1"], ids_by_mate["R2"]),
                    start=1,
                )
                if read1_id != read2_id
            ),
            None,
        )
        raise ValueError(f"R1/R2 normalized IDs diverge at record {mismatch}")

    id_digest = hashlib.sha256()
    for read_id in ids_by_mate["R1"]:
        id_digest.update(read_id.encode("ascii"))
        id_digest.update(b"\n")
    payload = {
        "status": "passed",
        "source_project": rows[0]["ProjectAccession"],
        "source_sample": rows[0]["SampleAccession"],
        "run_accession": rows[0]["RunAccession"],
        "layout": rows[0]["Layout"],
        "selection": rows[0]["SelectionRule"],
        "prefix_records_per_mate": int(rows[0]["PrefixRecords"]),
        "input_pairs": int(rows[0]["PrefixRecords"]),
        "mates_synchronized": mates_synchronized,
        "normalized_pair_id_sha256": id_digest.hexdigest(),
        "raw_fastq_committed": False,
        "complete_archive_md5_verified": False,
        "archive_identity_boundary": (
            "ENA MD5 and byte counts identify the complete remote files; "
            "this streaming prefix run does not recompute full-file MD5."
        ),
        "mates": summaries,
    }
    summary_path = output_dir / "subset-summary.json"
    summary_path.write_text(
        json.dumps(payload, indent=2, sort_keys=False) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(payload, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
