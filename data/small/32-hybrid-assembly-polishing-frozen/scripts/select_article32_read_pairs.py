#!/usr/bin/env python3
"""Select an exact synchronized FASTQ-pair subset for Article 32.

The one-pass sequential sampler keeps no read pool, preserves source order, and
uses a namespaced MT19937 seed. Gzip mtime is zeroed for byte-stable reruns on
the locked Python/zlib stack.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import io
import json
import platform
import random
import zlib
from pathlib import Path
from typing import TextIO


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--r1", type=Path, required=True)
    parser.add_argument("--r2", type=Path, required=True)
    parser.add_argument("--output-r1", type=Path, required=True)
    parser.add_argument("--output-r2", type=Path, required=True)
    parser.add_argument("--summary", type=Path, required=True)
    parser.add_argument("--total-pairs", type=int, required=True)
    parser.add_argument("--target-pairs", type=int, default=10_000_000)
    parser.add_argument("--seed", type=int, default=20260732)
    return parser.parse_args()


def sha256(path: Path, *, decompress: bool = False) -> str:
    digest = hashlib.sha256()
    opener = gzip.open if decompress else Path.open
    with opener(path, "rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_record(handle: TextIO, path: Path, index: int) -> tuple[str, str, str, str]:
    header = handle.readline()
    if not header:
        raise EOFError
    sequence = handle.readline()
    plus = handle.readline()
    quality = handle.readline()
    if not (sequence and plus and quality):
        raise ValueError(f"Truncated FASTQ record {index} in {path}")
    if not header.startswith("@") or not plus.startswith("+"):
        raise ValueError(f"Malformed FASTQ record {index} in {path}")
    if len(sequence.rstrip("\r\n")) != len(quality.rstrip("\r\n")):
        raise ValueError(f"Sequence/quality mismatch at record {index} in {path}")
    return header, sequence, plus, quality


def normalized_id(header: str) -> str:
    token = header[1:].split(None, 1)[0]
    if token.endswith(("/1", "/2")):
        token = token[:-2]
    return token


def gzip_text_writer(path: Path) -> tuple[io.TextIOWrapper, io.BufferedWriter]:
    raw = path.open("wb")
    zipped = gzip.GzipFile(filename="", mode="wb", fileobj=raw, compresslevel=6, mtime=0)
    return io.TextIOWrapper(zipped, encoding="utf-8", newline=""), raw


def main() -> int:
    args = parse_args()
    for source in (args.r1, args.r2):
        if not source.is_file() or source.stat().st_size == 0:
            raise SystemExit(f"Missing FASTQ archive: {source}")
    if not 0 < args.target_pairs <= args.total_pairs:
        raise SystemExit("target-pairs must be in 1..total-pairs")
    for output in (args.output_r1, args.output_r2, args.summary):
        if output.exists():
            raise SystemExit(f"Refusing to overwrite existing output: {output}")
        output.parent.mkdir(parents=True, exist_ok=True)

    rng = random.Random()
    namespace = f"{args.seed}:MOCK1:ERR9765746"
    rng.seed(namespace, version=2)
    remaining_total = args.total_pairs
    remaining_target = args.target_pairs
    selected = 0
    bases = [0, 0]
    minimum = [None, None]
    maximum = [0, 0]
    pair_ids = hashlib.sha256()

    out1, raw1 = gzip_text_writer(args.output_r1)
    out2, raw2 = gzip_text_writer(args.output_r2)
    try:
        with gzip.open(args.r1, "rt", encoding="utf-8", newline="") as in1, gzip.open(
            args.r2, "rt", encoding="utf-8", newline=""
        ) as in2:
            for index in range(1, args.total_pairs + 1):
                try:
                    rec1 = read_record(in1, args.r1, index)
                    rec2 = read_record(in2, args.r2, index)
                except EOFError as exc:
                    raise ValueError(
                        f"FASTQ ended before declared pair {args.total_pairs}: {index}"
                    ) from exc
                id1, id2 = normalized_id(rec1[0]), normalized_id(rec2[0])
                if id1 != id2:
                    raise ValueError(f"Mate ID mismatch at pair {index}: {id1} != {id2}")

                if rng.random() < (remaining_target / remaining_total):
                    out1.writelines(rec1)
                    out2.writelines(rec2)
                    lengths = (
                        len(rec1[1].rstrip("\r\n")),
                        len(rec2[1].rstrip("\r\n")),
                    )
                    for mate, length in enumerate(lengths):
                        bases[mate] += length
                        minimum[mate] = length if minimum[mate] is None else min(minimum[mate], length)
                        maximum[mate] = max(maximum[mate], length)
                    pair_ids.update(id1.encode())
                    pair_ids.update(b"\n")
                    selected += 1
                    remaining_target -= 1
                remaining_total -= 1

            if in1.readline() or in2.readline():
                raise ValueError("FASTQ contains more pairs than the declared total")
    finally:
        out1.close()
        out2.close()
        raw1.close()
        raw2.close()

    if selected != args.target_pairs or remaining_target or remaining_total:
        raise RuntimeError("Exact sequential sample invariant failed")

    payload = {
        "sample": "MOCK1",
        "run_accession": "ERR9765746",
        "seed": args.seed,
        "rng_namespace": namespace,
        "random_generator": "Python random.Random MT19937; string seed version=2",
        "selection_algorithm": "one-pass exact sequential sampling without replacement",
        "python_version": platform.python_version(),
        "python_implementation": platform.python_implementation(),
        "zlib_build_version": zlib.ZLIB_VERSION,
        "zlib_runtime_version": zlib.ZLIB_RUNTIME_VERSION,
        "source_pairs": args.total_pairs,
        "selected_pairs": selected,
        "selected_fraction": selected / args.total_pairs,
        "selected_pair_id_sha256": pair_ids.hexdigest(),
        "mates": {},
    }
    for mate, path in enumerate((args.output_r1, args.output_r2)):
        payload["mates"][f"R{mate + 1}"] = {
            "bases": bases[mate],
            "minimum_length": minimum[mate],
            "maximum_length": maximum[mate],
            "mean_length": bases[mate] / selected,
            "compressed_sha256": sha256(path),
            "uncompressed_sha256": sha256(path, decompress=True),
            "bytes": path.stat().st_size,
        }
    args.summary.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

