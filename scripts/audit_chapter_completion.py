#!/usr/bin/env python3
"""Audit all 77 QMDs against the publication contract and current QA gate."""

from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from pathlib import Path
from typing import Any

import yaml


SHARED_ROOT = Path(__file__).resolve().parents[2]
if str(SHARED_ROOT) not in sys.path:
    sys.path.insert(0, str(SHARED_ROOT))

from tutorial_automation.manifest import load_manifest


REQUIRED_SECTIONS = (
    "这一步对应论文里的哪张图",
    "理论",
    "准备工作",
    "可复制代码",
    "审计与升级",
    "出版级美化",
    "常见坑",
    "这段 Methods 怎么写",
    "换成你自己的数据怎么做",
    "参考",
)
PROHIBITED_PATTERNS = (
    "作者代码通常长这样",
    "本篇可独立跑通",
    "这体现全系列",
    "即本文",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, default=Path("."))
    parser.add_argument("--manifest", type=Path, default=Path("tutorial.yaml"))
    parser.add_argument("--qa-report", type=Path, default=Path("qa_report.json"))
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--tsv", type=Path)
    parser.add_argument("--require-complete", action="store_true")
    return parser.parse_args()


def frontmatter(text: str) -> dict[str, Any]:
    if not text.startswith("---\n"):
        return {}
    end = text.find("\n---\n", 4)
    if end < 0:
        return {}
    payload = yaml.safe_load(text[4:end])
    return payload if isinstance(payload, dict) else {}


def prose_size(text: str) -> int:
    without_frontmatter = re.sub(r"\A---\n.*?\n---\n", "", text, count=1, flags=re.S)
    without_comments = re.sub(r"<!--.*?-->", "", without_frontmatter, flags=re.S)
    without_code = re.sub(r"```.*?```", "", without_comments, flags=re.S)
    return len(re.sub(r"\s+", "", without_code))


def main() -> int:
    args = parse_args()
    root = args.project_root.resolve()
    manifest_path = (root / args.manifest).resolve()
    qa_path = (root / args.qa_report).resolve()
    normalized = load_manifest(manifest_path)
    raw_manifest = yaml.safe_load(manifest_path.read_text(encoding="utf-8"))
    chapters = raw_manifest["series"]["chapters"]
    validated = {
        int(number)
        for number in raw_manifest["series"].get("validated_articles", [])
    }

    qa: dict[str, Any] = {}
    if qa_path.is_file():
        qa = json.loads(qa_path.read_text(encoding="utf-8"))
    qa_current = (
        qa.get("status") == "passed"
        and qa.get("manifest_hash") == normalized["manifest_hash"]
    )

    rows: list[dict[str, Any]] = []
    for contract in chapters:
        number = int(contract["number"])
        path = root / contract["file"]
        text = path.read_text(encoding="utf-8") if path.is_file() else ""
        metadata = frontmatter(text)
        execute = metadata.get("execute", {})
        missing_sections = [
            section for section in REQUIRED_SECTIONS if section not in text
        ]
        prohibited = [
            pattern for pattern in PROHIBITED_PATTERNS if pattern in text
        ]
        checks = {
            "file_exists": path.is_file(),
            "title_matches": contract["title"] in str(metadata.get("title", "")),
            "not_draft": metadata.get("draft") is not True,
            "execution_declared": (
                isinstance(execute, dict)
                and isinstance(execute.get("eval"), bool)
            ),
            "sections_complete": not missing_sections,
            "substantive_length": prose_size(text) >= 3500,
            "real_input_declared": (
                "data/" in text
                or "curatedMetagenomicData" in text
                or "ExperimentHub" in text
            ),
            "figure_path_declared": "figures/" in text,
            "determinism_declared": (
                "set.seed(" in text
                or "random.seed(" in text
                or "np.random.seed(" in text
                or "--seed" in text
                or "rng_seed" in text
                or "不含随机过程" in text
            ),
            "methods_present": "Methods" in text,
            "references_present": (
                "[@" in text
                or "doi.org/" in text
                or "{#refs}" in text
            ),
            "wechat_contract": isinstance(metadata.get("wechat"), dict),
            "no_prohibited_text": not prohibited,
            "listed_as_validated": number in validated,
            "qa_manifest_current": qa_current,
        }
        ready = all(checks.values())
        rows.append(
            {
                "number": number,
                "file": contract["file"],
                "title": contract["title"],
                "ready": ready,
                "draft": metadata.get("draft") is True,
                "eval": execute.get("eval") if isinstance(execute, dict) else None,
                "prose_characters": prose_size(text),
                "missing_sections": missing_sections,
                "prohibited_patterns": prohibited,
                "failed_checks": [
                    key for key, passed in checks.items() if not passed
                ],
                "checks": checks,
            }
        )

    ready_count = sum(bool(row["ready"]) for row in rows)
    formal_count = sum(not bool(row["draft"]) for row in rows)
    status = "complete" if ready_count == len(rows) == 77 else "in_progress"
    report = {
        "status": status,
        "manifest_hash": normalized["manifest_hash"],
        "qa_status": qa.get("status"),
        "qa_manifest_hash": qa.get("manifest_hash"),
        "qa_manifest_current": qa_current,
        "summary": {
            "total_articles": len(rows),
            "formal_articles": formal_count,
            "validated_articles": len(validated),
            "ready_articles": ready_count,
            "remaining_articles": len(rows) - ready_count,
        },
        "articles": rows,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(report, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )

    if args.tsv:
        args.tsv.parent.mkdir(parents=True, exist_ok=True)
        with args.tsv.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=(
                    "number",
                    "file",
                    "ready",
                    "draft",
                    "eval",
                    "prose_characters",
                    "failed_checks",
                ),
                delimiter="\t",
            )
            writer.writeheader()
            for row in rows:
                writer.writerow(
                    {
                        "number": row["number"],
                        "file": row["file"],
                        "ready": row["ready"],
                        "draft": row["draft"],
                        "eval": row["eval"],
                        "prose_characters": row["prose_characters"],
                        "failed_checks": ",".join(row["failed_checks"]),
                    }
                )

    print(json.dumps(report["summary"], ensure_ascii=False, indent=2))
    if args.require_complete and status != "complete":
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
