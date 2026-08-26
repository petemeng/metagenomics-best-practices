#!/usr/bin/env python3
"""Run and freeze representative HUMAnN 3.9 renormalization branches."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import shlex
import subprocess
from pathlib import Path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--environment-prefix", required=True, type=Path)
    parser.add_argument("--article19-dir", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()

    prefix = args.environment_prefix.resolve()
    source_dir = args.article19_dir.resolve()
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    logs_dir = output_dir / "logs"
    logs_dir.mkdir(parents=True, exist_ok=True)

    renorm = prefix / "bin" / "humann_renorm_table"
    humann = prefix / "bin" / "humann"
    input_table = source_dir / "pathabundance-rpk.tsv"
    if not renorm.is_file() or not humann.is_file() or not input_table.is_file():
        raise SystemExit("HUMAnN executables or Article 19 pathway table are missing.")

    env = os.environ.copy()
    env["PYTHONHASHSEED"] = "0"
    version = subprocess.run(
        [str(humann), "--version"],
        check=True,
        capture_output=True,
        text=True,
        env=env,
    ).stdout.strip()
    if "3.9" not in version:
        raise SystemExit(f"Expected HUMAnN 3.9, observed: {version}")

    commands: list[list[str]] = []
    outputs: list[Path] = []
    for mode in ("community", "levelwise"):
        for special in ("y", "n"):
            output = output_dir / f"pathabundance-relab-{mode}-special-{special}.tsv"
            log = logs_dir / f"pathabundance-relab-{mode}-special-{special}.log"
            command = [
                str(renorm),
                "--input", str(input_table),
                "--output", str(output),
                "--units", "relab",
                "--mode", mode,
                "--special", special,
            ]
            result = subprocess.run(
                command,
                check=False,
                capture_output=True,
                text=True,
                env=env,
            )
            log.write_text(
                f"command: {shlex.join(command)}\n"
                f"exit_code: {result.returncode}\n"
                f"stdout:\n{result.stdout}\n"
                f"stderr:\n{result.stderr}\n",
                encoding="utf-8",
            )
            if result.returncode != 0:
                raise SystemExit(f"HUMAnN renormalization failed; see {log}")
            commands.append(command)
            outputs.append(output)

    (output_dir / "tool-versions.tsv").write_text(
        "tool\tversion\nHUMAnN\t3.9\nhumann_renorm_table\t3.9\n",
        encoding="utf-8",
    )
    (output_dir / "commands.sh").write_text(
        "#!/usr/bin/env bash\nset -euo pipefail\nexport PYTHONHASHSEED=0\n"
        + "\n".join(shlex.join(command) for command in commands)
        + "\n",
        encoding="utf-8",
    )

    summary = {
        "status": "completed",
        "humann_version": "3.9",
        "input_file": "data/small/19-humann3-frozen/pathabundance-rpk.tsv",
        "input_sha256": sha256(input_table),
        "units": "relab",
        "modes": ["community", "levelwise"],
        "special_choices": ["y", "n"],
        "actual_branches": 4,
        "pythonhashseed": 0,
        "network_access": False,
    }
    (output_dir / "run-summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    checksum_file = output_dir / "file-checksums.sha256"
    payloads = sorted(
        path for path in output_dir.rglob("*")
        if path.is_file() and path != checksum_file
    )
    checksum_file.write_text(
        "".join(f"{sha256(path)}  {path.relative_to(output_dir)}\n" for path in payloads),
        encoding="utf-8",
    )
    print(json.dumps(summary, ensure_ascii=False, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
