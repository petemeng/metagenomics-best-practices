#!/usr/bin/env python3
"""Reusable fail-closed checks for the Article 42--44 frozen tutorials."""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import re
from dataclasses import dataclass, field
from pathlib import Path

from PIL import Image


REQUIRED_SECTIONS = (
    "这一步对应论文里的哪张图", "理论", "准备工作", "可复制代码", "审计与升级",
    "出版级美化", "常见坑", "这段 Methods 怎么写", "换成你自己的数据怎么做", "参考",
)


def read_tsv(path: Path) -> list[dict[str, str]]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def as_bool(value: object) -> bool:
    return str(value).strip().lower() in {"true", "t", "1", "yes", "pass"}


@dataclass
class Audit:
    rows: list[dict[str, object]] = field(default_factory=list)

    def add(self, category: str, check: str, status: bool, detail: object) -> None:
        self.rows.append({
            "Category": category, "CheckID": check, "Status": "PASS" if status else "FAIL",
            "Detail": detail if isinstance(detail, str) else json.dumps(detail, ensure_ascii=False, sort_keys=True),
        })

    @property
    def passed(self) -> int:
        return sum(row["Status"] == "PASS" for row in self.rows)

    @property
    def failed(self) -> int:
        return sum(row["Status"] == "FAIL" for row in self.rows)


def audit_checksums(frozen: Path, audit: Audit) -> None:
    manifest = frozen / "file-checksums.sha256"
    audit.add("Checksum", "manifest-exists", manifest.is_file(), str(manifest))
    expected = {}
    for line in manifest.read_text(encoding="utf-8").splitlines():
        digest, relative = line.split("  ", 1)
        expected[relative] = digest
        path = frozen / relative
        observed = sha256(path) if path.is_file() else "MISSING"
        audit.add("Checksum", relative, observed == digest, observed)
    actual = {
        path.relative_to(frozen).as_posix() for path in frozen.rglob("*")
        if path.is_file() and path.name != manifest.name
    }
    audit.add("Checksum", "manifest-complete", set(expected) == actual, {"listed": len(expected), "actual": len(actual)})


def audit_chapter(
    chapter: Path,
    audit: Audit,
    *,
    article: int,
    figure_stems: tuple[str, ...],
    tokens: tuple[str, ...],
) -> None:
    text = chapter.read_text(encoding="utf-8")
    expectations = {
        "draft-false": "draft: false" in text,
        "eval-false": re.search(r"execute:\s*\n\s+eval:\s+false", text) is not None,
        "freeze-auto": "freeze: auto" in text,
        "bibliography": "bibliography: ../references.bib" in text,
        "seed": str(20260700 + article) in text,
        "inline-theme": all(token in text for token in ("pal_pub <-", "theme_pub <-", "save_pub <-")),
        "resource-contract": all(token in text for token in ("RAM", "CPU", "磁盘", "耗时")),
        "methods-template": "Methods template" in text,
        "results-template": "Results template" in text,
        "real-study": all(token in text for token in ("PRJEB52977", "ERR9765746", "ERR9765747")),
        "frozen-input": f"data/small/{article}-" in text and "-frozen" in text,
        "figures": all(f"../figures/{stem}.png" in text for stem in figure_stems),
        "vector-raster-export": all(token in text for token in (".pdf", ".png", ".tiff", "compression = \"lzw\"")),
        "required-tokens": all(token in text for token in tokens),
        "no-source-theme": 'source("R/theme_pub.R")' not in text and "source('R/theme_pub.R')" not in text,
        "no-placeholders": re.search(r"__[A-Z0-9_]+__|TODO|TBD|NNN|vX", text) is None,
        "no-meta-prose": not any(token in text for token in ("本篇可独立", "本文可独立", "全系列约定", "接口只学一次", "作者代码通常长这样", "（即本文）", "无头服务器")),
    }
    for section in REQUIRED_SECTIONS:
        expectations[f"section-{section}"] = re.search(rf"(?m)^##\s+{re.escape(section)}", text) is not None
    for check, status in expectations.items():
        audit.add("Chapter", check, status, check)


def audit_figures(figure_dir: Path, audit: Audit, stems: tuple[str, ...]) -> None:
    for stem in stems:
        for extension in ("pdf", "png", "tiff"):
            path = figure_dir / f"{stem}.{extension}"
            exists = path.is_file() and path.stat().st_size > 0
            audit.add("Figure", f"{stem}-{extension}-exists", exists, str(path))
            if not exists or extension == "pdf":
                continue
            with Image.open(path) as image:
                dpi = tuple(float(value) for value in image.info.get("dpi", (0, 0)))
                width, height = image.size
                audit.add("Figure", f"{stem}-{extension}-dimensions", width >= 1400 and height >= 900, f"{width}x{height}")
                audit.add("Figure", f"{stem}-{extension}-dpi", min(dpi) >= 340, dpi)
                audit.add("Figure", f"{stem}-{extension}-mode", image.mode in {"RGB", "RGBA"}, image.mode)
                if extension == "tiff":
                    compression = image.info.get("compression", "")
                    audit.add("Figure", f"{stem}-tiff-lzw", str(compression).lower() in {"tiff_lzw", "lzw"}, compression)


def finish(
    *, article: int, audit: Audit, output: Path, payload: dict[str, object]
) -> int:
    output.mkdir(parents=True, exist_ok=True)
    write_tsv(output / "validation-checks.tsv", audit.rows)
    result = {
        "article": article, "status": "passed" if audit.failed == 0 else "failed",
        "checks_passed": audit.passed, "checks_failed": audit.failed, **payload,
    }
    (output / "validation-summary.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    (output / "validation.log").write_text(
        f"Article {article} validation: {result['status'].upper()}\nPASS={audit.passed}\nFAIL={audit.failed}\n"
        + "\n".join(f"{row['Status']}\t{row['Category']}\t{row['CheckID']}\t{row['Detail']}" for row in audit.rows) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0 if audit.failed == 0 else 1
