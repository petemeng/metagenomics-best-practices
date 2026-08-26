#!/usr/bin/env python3
"""Select an exact, synchronized FASTQ-pair subset with a fixed seed.

The algorithm performs one-pass sequential sampling without replacement. It
stores no read pool, selects exactly ``target_pairs`` from ``total_pairs``, and
keeps selected records in their original run order. Gzip output uses mtime=0
so reruns with the same Python/zlib stack produce stable bytes.
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
    parser.add_argument("--sample", required=True)
    parser.add_argument("--run-accession", required=True)
    parser.add_argument("--total-pairs", type=int, required=True)
    parser.add_argument("--target-pairs", type=int, required=True)
    parser.add_argument("--seed", type=int, default=20260730)
    return parser.parse_args()


def sha256(path: Path, decompress: bool = False) -> str:
    digest = hashlib.sha256()
    opener = gzip.open if decompress else Path.open
    if decompress:
        handle = opener(path, "rb")
    else:
        handle = opener(path, "rb")
    with handle:
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
    if token.endswith("/1") or token.endswith("/2"):
        token = token[:-2]
    return token


def gzip_text_writer(path: Path) -> tuple[io.TextIOWrapper, io.BufferedWriter]:
    raw = path.open("wb")
    zipped = gzip.GzipFile(filename="", mode="wb", fileobj=raw, compresslevel=6, mtime=0)
    text = io.TextIOWrapper(zipped, encoding="utf-8", newline="")
    return text, raw


def main() -> int:
    args = parse_args()
    for source in (args.r1, args.r2):
        if not source.is_file() or source.stat().st_size == 0:
            raise SystemExit(f"Missing FASTQ archive: {source}")
    if args.target_pairs <= 0 or args.target_pairs > args.total_pairs:
        raise SystemExit("target_pairs must be in 1..total_pairs")
    for output in (args.output_r1, args.output_r2, args.summary):
        if output.exists():
            raise SystemExit(f"Refusing to overwrite existing output: {output}")
        output.parent.mkdir(parents=True, exist_ok=True)

    rng = random.Random()
    rng.seed(f"{args.seed}:{args.sample}", version=2)
    remaining_total = args.total_pairs
    remaining_target = args.target_pairs
    selected = 0
    selected_bases = [0, 0]
    selected_min = [None, None]
    selected_max = [0, 0]
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
                id1 = normalized_id(rec1[0])
                id2 = normalized_id(rec2[0])
                if id1 != id2:
                    raise ValueError(f"Mate ID mismatch at pair {index}: {id1} != {id2}")

                choose = rng.random() < (remaining_target / remaining_total)
                if choose:
                    out1.writelines(rec1)
                    out2.writelines(rec2)
                    lengths = (
                        len(rec1[1].rstrip("\r\n")),
                        len(rec2[1].rstrip("\r\n")),
                    )
                    for mate in (0, 1):
                        selected_bases[mate] += lengths[mate]
                        selected_min[mate] = (
                            lengths[mate]
                            if selected_min[mate] is None
                            else min(selected_min[mate], lengths[mate])
                        )
                        selected_max[mate] = max(selected_max[mate], lengths[mate])
                    pair_ids.update(id1.encode("utf-8"))
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

    if selected != args.target_pairs or remaining_target != 0 or remaining_total != 0:
        raise RuntimeError("Exact sequential sample invariant failed")

    payload = {
        "sample": args.sample,
        "run_accession": args.run_accession,
        "seed": args.seed,
        "rng_namespace": f"{args.seed}:{args.sample}",
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
        "mates": {
            "R1": {
                "bases": selected_bases[0],
                "minimum_length": selected_min[0],
                "maximum_length": selected_max[0],
                "mean_length": selected_bases[0] / selected,
                "compressed_sha256": sha256(args.output_r1),
                "uncompressed_sha256": sha256(args.output_r1, decompress=True),
                "bytes": args.output_r1.stat().st_size,
            },
            "R2": {
                "bases": selected_bases[1],
                "minimum_length": selected_min[1],
                "maximum_length": selected_max[1],
                "mean_length": selected_bases[1] / selected,
                "compressed_sha256": sha256(args.output_r2),
                "uncompressed_sha256": sha256(args.output_r2, decompress=True),
                "bytes": args.output_r2.stat().st_size,
            },
        },
    }
    args.summary.write_text(
        json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    print(json.dumps(payload, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
