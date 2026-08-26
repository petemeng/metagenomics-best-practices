#!/usr/bin/env python3
"""Create missing draft chapters and regenerate the ordered Quarto book contract."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any

import yaml


STUB = """---
title: "{display_title}"
draft: true
execute:
  eval: false
---

<!--
Planned chapter. Do not publish until its public data, lineage card, literature,
commands, outputs, figures, and assertions are recorded in tutorial.yaml and
qa_report.json.status is passed.
-->
"""


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--project-root",
        type=Path,
        default=Path(__file__).resolve().parents[1],
    )
    parser.add_argument(
        "--write-quarto",
        action="store_true",
        help="Regenerate _quarto.yml from tutorial.yaml.",
    )
    return parser.parse_args()


def chapter_map(chapters: list[dict[str, Any]]) -> dict[int, dict[str, Any]]:
    return {int(chapter["number"]): chapter for chapter in chapters}


def quarto_payload(manifest: dict[str, Any]) -> dict[str, Any]:
    series = manifest["series"]
    chapters = chapter_map(series["chapters"])
    book_chapters: list[Any] = [chapters[1]["file"]]

    for part in series["parts"]:
        numbers = [int(number) for number in part["chapters"] if int(number) != 1]
        if not numbers:
            continue
        book_chapters.append(
            {
                "part": part["title"],
                "chapters": [chapters[number]["file"] for number in numbers],
            }
        )

    return {
        "project": {
            "type": "book",
            "output-dir": "_site",
            "execute-dir": "project",
        },
        "book": {
            "title": "宏基因组分析最佳实践 · Shotgun",
            "subtitle": "从真实公开数据到可投稿证据",
            "author": "Songlab",
            "site-url": "https://petemeng.github.io/metagenomics-best-practices/",
            "repo-url": "https://github.com/petemeng/metagenomics-best-practices",
            "repo-actions": ["edit", "issue"],
            "page-navigation": True,
            "search": True,
            "chapters": book_chapters,
        },
        "lang": "zh",
        "execute": {
            "freeze": "auto",
            "warning": False,
            "message": False,
        },
        "format": {
            "html": {
                "theme": {
                    "light": ["cosmo", "styles.scss"],
                    "dark": ["darkly", "styles.scss"],
                },
                "toc": True,
                "toc-depth": 3,
                "code-copy": True,
                "code-tools": True,
                "code-overflow": "wrap",
                "anchor-sections": True,
                "smooth-scroll": True,
            }
        },
        "bibliography": "references.bib",
    }


def main() -> int:
    args = parse_args()
    root = args.project_root.resolve()
    manifest_path = root / "tutorial.yaml"
    manifest = yaml.safe_load(manifest_path.read_text(encoding="utf-8"))
    chapters = manifest["series"]["chapters"]
    created: list[Path] = []

    for chapter in chapters:
        path = root / chapter["file"]
        if int(chapter["number"]) == 1 or path.exists():
            continue
        path.parent.mkdir(parents=True, exist_ok=True)
        display_title = f"第 {int(chapter['number']):02d} 篇 · {chapter['title']}"
        path.write_text(STUB.format(display_title=display_title), encoding="utf-8")
        created.append(path)

    if args.write_quarto:
        quarto_path = root / "_quarto.yml"
        quarto_path.write_text(
            yaml.safe_dump(
                quarto_payload(manifest),
                allow_unicode=True,
                sort_keys=False,
                width=1000,
            ),
            encoding="utf-8",
        )

    print(f"created={len(created)}")
    for path in created:
        print(path.relative_to(root))
    if args.write_quarto:
        print("wrote=_quarto.yml")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
