#!/usr/bin/env python3
"""Convert an immutable official VFDB FASTA snapshot into an ABRicate database."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import re
import subprocess
from pathlib import Path


HEADER = re.compile(
    r"^(?P<vfg>\w+)(?:\(\w+\|(?P<accession>\w+)(?:\.\d+)?\))?\s+(?P<description>.*)$"
)
GENE = re.compile(r"^\((?P<gene>[^)]+)\)\s*(?P<rest>.*)$")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--database-name", required=True)
    parser.add_argument("--makeblastdb", type=Path, required=True)
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_fasta(path: Path) -> list[tuple[str, str]]:
    opener = gzip.open if path.suffix == ".gz" else open
    rows: list[tuple[str, str]] = []
    header: str | None = None
    seq: list[str] = []
    # VFDB descriptions still contain a small number of legacy Latin-1 bytes.
    with opener(path, "rt", encoding="latin-1") as handle:
        for raw in handle:
            line = raw.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    rows.append((header, "".join(seq)))
                header, seq = line[1:], []
            else:
                seq.append(line.strip().upper())
    if header is not None:
        rows.append((header, "".join(seq)))
    return rows


def main() -> int:
    args = parse_args()
    source = args.input.resolve()
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    parsed = []
    seen_vfg: set[str] = set()
    for header, sequence in read_fasta(source):
        match = HEADER.fullmatch(header)
        if match is None:
            raise ValueError(f"Unrecognized VFDB header: {header}")
        if match["vfg"] in seen_vfg:
            raise ValueError(f"Duplicate VFG identifier: {match['vfg']}")
        seen_vfg.add(match["vfg"])
        description = match["description"]
        gene_match = GENE.match(description)
        gene = gene_match["gene"] if gene_match else match["vfg"]
        clean = re.sub(r"[^AGTC]", "N", sequence)
        parsed.append((gene, match["accession"] or "", match["vfg"], description, clean))
    if not parsed:
        raise ValueError("VFDB input contains no sequences")
    parsed.sort(key=lambda row: (row[0].casefold(), row[1], row[2]))
    fasta = output / "sequences"
    with fasta.open("w", encoding="utf-8") as handle:
        for gene, accession, vfg, description, sequence in parsed:
            handle.write(f">{args.database_name}~~~{gene}~~~{accession}~~~ {vfg} {description}\n")
            for start in range(0, len(sequence), 60):
                handle.write(sequence[start:start + 60] + "\n")
    log = output / "makeblastdb.log"
    subprocess.run([
        str(args.makeblastdb), "-in", str(fasta), "-title", args.database_name,
        "-dbtype", "nucl", "-hash_index", "-logfile", str(log),
    ], check=True)
    manifest = {
        "database": args.database_name,
        "source": str(source),
        "source_bytes": source.stat().st_size,
        "source_sha256": sha256(source),
        "sequences": len(parsed),
        "formatted_bytes": fasta.stat().st_size,
        "formatted_sha256": sha256(fasta),
        "header_contract": "database~~~gene~~~accession~~~VFG description",
    }
    (output / "build-manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(manifest, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
