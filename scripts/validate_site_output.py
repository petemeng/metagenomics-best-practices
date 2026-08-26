#!/usr/bin/env python3
"""Validate the rendered Quarto site without network access.

The checker treats ``tutorial.yaml`` as the page contract, confirms that every
declared QMD has a rendered HTML counterpart, and resolves local page/assets
exactly as a browser would.  It intentionally does not probe external URLs:
those are versioned and audited by the per-article source manifests.
"""

from __future__ import annotations

import argparse
import json
import re
from collections import Counter
from html.parser import HTMLParser
from pathlib import Path
from typing import Any
from urllib.parse import unquote, urlsplit

import yaml


MIN_HTML_BYTES = 10_000
EXPECTED_COUNT = 77
IGNORED_SCHEMES = {"data", "http", "https", "javascript", "mailto", "tel"}
FORBIDDEN_RENDERED_PATTERNS = {
    "personal_workspace_path": re.compile(r"/workspace"),
    "redaction_placeholder": re.compile(r"\[REDACTED_LOCAL_PATH\]"),
    "draft_frontmatter": re.compile(r"\bdraft\s*:\s*true\b", re.IGNORECASE),
    "planned_chapter_placeholder": re.compile(r"\bPlanned chapter\b", re.IGNORECASE),
    "missing_cjk_glyphs": re.compile(r"□□"),
}


class ReferenceParser(HTMLParser):
    """Collect local-resource candidates from one rendered page."""

    def __init__(self) -> None:
        super().__init__(convert_charrefs=True)
        self.references: list[tuple[str, str]] = []

    def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        values = dict(attrs)
        if tag in {"img", "script", "source", "video", "audio"}:
            for attr in ("src", "poster"):
                if values.get(attr):
                    self.references.append((f"{tag}:{attr}", values[attr] or ""))
        if tag == "link" and values.get("href"):
            self.references.append(("link:href", values["href"] or ""))
        if tag == "a" and values.get("href"):
            self.references.append(("a:href", values["href"] or ""))
        for attr in ("srcset",):
            if values.get(attr):
                for candidate in (values[attr] or "").split(","):
                    url = candidate.strip().split(" ", 1)[0]
                    if url:
                        self.references.append((f"{tag}:{attr}", url))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--contract", type=Path, default=Path("tutorial.yaml"))
    parser.add_argument("--site-dir", type=Path, default=Path("_site"))
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("qa/site-output-validation.json"),
    )
    return parser.parse_args()


def qmd_to_html(qmd: str) -> Path:
    source = Path(qmd)
    return Path("index.html") if source == Path("index.qmd") else source.with_suffix(".html")


def article_sources(contract: dict[str, Any]) -> list[str]:
    series = contract.get("series")
    if not isinstance(series, dict):
        raise ValueError("tutorial.yaml must contain a series mapping")
    articles = series.get("chapters")
    if not isinstance(articles, list):
        raise ValueError("tutorial.yaml must contain a series.chapters list")
    files: list[str] = []
    for article in articles:
        if not isinstance(article, dict) or not isinstance(article.get("file"), str):
            raise ValueError("every article entry must define a string file path")
        files.append(article["file"])
    return files


def local_target(site_dir: Path, page: Path, raw_url: str) -> Path | None:
    split = urlsplit(raw_url.strip())
    if not raw_url.strip() or split.scheme.lower() in IGNORED_SCHEMES or split.netloc:
        return None
    path_text = unquote(split.path)
    if not path_text:
        return None
    if path_text.startswith("/"):
        candidate = site_dir / path_text.lstrip("/")
    else:
        candidate = page.parent / path_text
    return candidate.resolve(strict=False)


def resolve_existing_target(candidate: Path) -> Path | None:
    if candidate.is_file():
        return candidate
    if candidate.is_dir() and (candidate / "index.html").is_file():
        return candidate / "index.html"
    return None


def main() -> int:
    args = parse_args()
    root = Path.cwd().resolve()
    contract_path = (root / args.contract).resolve()
    site_dir = (root / args.site_dir).resolve()
    output_path = (root / args.output).resolve()

    with contract_path.open("r", encoding="utf-8") as handle:
        contract = yaml.safe_load(handle)
    if not isinstance(contract, dict):
        raise ValueError("tutorial.yaml root must be a mapping")

    qmd_files = article_sources(contract)
    expected_pages = [(site_dir / qmd_to_html(qmd)).resolve() for qmd in qmd_files]
    errors: list[dict[str, Any]] = []
    warnings: list[dict[str, Any]] = []
    checks = Counter()
    reference_counts = Counter()
    local_targets: set[Path] = set()

    checks["article_count"] += 1
    if len(qmd_files) != EXPECTED_COUNT:
        errors.append(
            {"check": "article_count", "expected": EXPECTED_COUNT, "observed": len(qmd_files)}
        )

    duplicate_sources = sorted(path for path, count in Counter(qmd_files).items() if count > 1)
    checks["unique_article_sources"] += 1
    if duplicate_sources:
        errors.append({"check": "unique_article_sources", "paths": duplicate_sources})

    existing_pages: list[Path] = []
    for qmd, page in zip(qmd_files, expected_pages, strict=True):
        checks["expected_html_exists"] += 1
        if not page.is_file():
            errors.append({"check": "expected_html_exists", "source": qmd, "page": str(page)})
            continue
        existing_pages.append(page)
        checks["minimum_html_size"] += 1
        if page.stat().st_size < MIN_HTML_BYTES:
            errors.append(
                {
                    "check": "minimum_html_size",
                    "page": str(page),
                    "bytes": page.stat().st_size,
                    "minimum": MIN_HTML_BYTES,
                }
            )

        rendered = page.read_text(encoding="utf-8", errors="replace")
        for label, pattern in FORBIDDEN_RENDERED_PATTERNS.items():
            checks[f"forbidden_rendered_text:{label}"] += 1
            match = pattern.search(rendered)
            if match:
                errors.append(
                    {
                        "check": f"forbidden_rendered_text:{label}",
                        "page": str(page),
                        "match": match.group(0),
                    }
                )

        parser = ReferenceParser()
        parser.feed(rendered)
        for kind, raw_url in parser.references:
            reference_counts[kind] += 1
            candidate = local_target(site_dir, page, raw_url)
            if candidate is None:
                continue
            # A bare fragment has already been excluded by local_target.  For
            # all remaining local references, both navigation and assets must
            # resolve inside the rendered site.
            checks["local_reference_resolves"] += 1
            resolved = resolve_existing_target(candidate)
            if resolved is None:
                errors.append(
                    {
                        "check": "local_reference_resolves",
                        "page": str(page),
                        "kind": kind,
                        "url": raw_url,
                        "candidate": str(candidate),
                    }
                )
            else:
                local_targets.add(resolved)

    actual_pages = sorted(site_dir.rglob("*.html")) if site_dir.is_dir() else []
    expected_set = set(expected_pages)
    unexpected_pages = sorted(str(path) for path in set(actual_pages) - expected_set)
    missing_from_inventory = sorted(str(path) for path in expected_set - set(actual_pages))
    checks["rendered_html_inventory"] += 1
    if unexpected_pages:
        warnings.append({"check": "rendered_html_inventory", "unexpected": unexpected_pages})
    if missing_from_inventory:
        errors.append({"check": "rendered_html_inventory", "missing": missing_from_inventory})

    report = {
        "status": "passed" if not errors else "failed",
        "contract": str(contract_path),
        "site_dir": str(site_dir),
        "expected_article_count": EXPECTED_COUNT,
        "declared_article_count": len(qmd_files),
        "rendered_page_count": len(actual_pages),
        "validated_page_count": len(existing_pages),
        "check_count": sum(checks.values()),
        "checks_by_type": dict(sorted(checks.items())),
        "references_by_type": dict(sorted(reference_counts.items())),
        "unique_resolved_local_targets": len(local_targets),
        "errors": errors,
        "warnings": warnings,
    }
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, ensure_ascii=False, indent=2))
    return 0 if not errors else 1


if __name__ == "__main__":
    raise SystemExit(main())
