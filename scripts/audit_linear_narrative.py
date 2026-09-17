#!/usr/bin/env python3
"""Audit the 77-chapter linear-narrative refactor.

The baseline records content fingerprints before structural movement.  The
post-refactor audit then proves that prose/code lines, code fences, images,
citations, callouts, and YAML metadata were not lost while also enforcing the
new reader-facing structure.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
from collections import Counter
from pathlib import Path

import yaml


GENERIC_H2 = re.compile(
    r"^##\s+(?:这一步对应论文里的哪张图|理论[：:]|准备工作(?:\s|\{|$)|"
    r"可复制代码(?:\s|[：:]|\{|$)|.*审计与升级|出版级美化|常见坑(?:方框)?(?:\s|[：:]|\{|$)|"
    r"这段 Methods 怎么写|换成你自己的数据怎么做)",
    re.M,
)
NUMBERED_HEADING = re.compile(
    r"^#{2,3}\s+(?:第\s*\d+\s*步(?:\s|[：:｜|])|"
    r"步骤\s*\d+(?:\s|[：:｜|])|"
    r"\d+(?:\.\d+)+(?:[.、：:]\s*|\s+)|"
    r"\d+(?:[、：:]|\.(?!\d))\s*)",
    re.M,
)
META_PHRASES = (
    "本篇可独立跑通",
    "这体现全系列",
    "接口只学一次",
    "作者代码通常长这样",
    "（即本文）",
)
PRIVATE_TOKENS = ("/media/desk16/", "/home/tly9658/", "tly9658")
STRUCTURAL_PREFIX = re.compile(
    r"^(?:第\s*\d+\s*步\s*[：:｜|]?\s*|"
    r"步骤\s*\d+\s*[：:｜|]?\s*|"
    r"\d+(?:\.\d+)+(?:[.、：:]\s*|\s+)|"
    r"\d+(?:[、：:]|\.(?!\d))\s*|"
    r"图\s*\d+\s*[：:]\s*)"
)
LEGACY_SLOT_TITLE = re.compile(
    r"^(?:这一步对应论文里的哪张图|理论[：:]|准备工作|可复制代码|"
    r".*审计与升级|出版级美化|常见坑(?:方框)?|这段 Methods 怎么写|"
    r"换成你自己的数据怎么做|小结|总结|参考(?:文献)?)"
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, default=Path("."))
    parser.add_argument("--manifest", default="tutorial.yaml")
    parser.add_argument("--baseline", type=Path, required=True)
    parser.add_argument("--write-baseline", action="store_true")
    parser.add_argument("--original-dir", type=Path)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--strict", action="store_true")
    return parser.parse_args()


def split_front(text: str) -> tuple[str, str]:
    if not text.startswith("---\n"):
        return "", text
    end = text.find("\n---\n", 4)
    if end < 0:
        return "", text
    return text[: end + 5], text[end + 5 :]


def sha(value: str) -> str:
    return hashlib.sha256(value.encode("utf-8")).hexdigest()


def content_lines(body: str) -> Counter[str]:
    """Fingerprint every non-heading content line, preserving multiplicity."""
    result: Counter[str] = Counter()
    in_fence = False
    fence = ""
    for raw in body.splitlines():
        stripped = raw.strip()
        marker = stripped[:3]
        if marker in {"```", "~~~"}:
            if not in_fence:
                in_fence = True
                fence = marker
            elif marker == fence:
                in_fence = False
                fence = ""
            result[sha(stripped)] += 1
            continue
        if not stripped:
            continue
        if not in_fence and re.match(r"^#{1,6}\s+", raw):
            continue
        result[sha(stripped)] += 1
    return result


def fence_languages(body: str) -> Counter[str]:
    langs: Counter[str] = Counter()
    for line in body.splitlines():
        match = re.match(r"^\s*(?:```|~~~)\{?([A-Za-z0-9_+-]*)", line)
        if match and match.group(1):
            langs[match.group(1).lower()] += 1
    return langs


def semantic_headings(body: str, *, baseline: bool) -> Counter[str]:
    result: Counter[str] = Counter()
    in_fence = False
    fence = ""
    for line in body.splitlines():
        stripped = line.lstrip()
        if stripped.startswith("```") or stripped.startswith("~~~"):
            marker = stripped[:3]
            if not in_fence:
                in_fence = True
                fence = marker
            elif marker == fence:
                in_fence = False
                fence = ""
            continue
        if in_fence:
            continue
        match = re.match(r"^#{2,3}\s+(.+?)\s*$", line)
        if not match:
            continue
        title = re.sub(r"\s*\{#[^}]+\}\s*$", "", match.group(1)).strip()
        if baseline and line.startswith("## ") and LEGACY_SLOT_TITLE.search(title):
            continue
        title = STRUCTURAL_PREFIX.sub("", title, count=1).strip()
        title = re.sub(r"^审计与升级\s*[：:]\s*", "", title).strip()
        if title:
            result[title] += 1
    return result


def metrics(text: str) -> dict:
    front, body = split_front(text)
    lines = content_lines(body)
    images = sorted(
        set(
            re.findall(
                r"(?:!\[[^\]]*\]\(|\bsrc=[\"'])([^)\"']+\.(?:png|jpe?g|svg|pdf|webp))",
                body,
                flags=re.I,
            )
        )
    )
    citations = sorted(set(re.findall(r"(?<![\w.-])@([A-Za-z0-9_:.+-]+)", body)))
    anchors = sorted(set(re.findall(r"\{#([A-Za-z0-9_.:-]+)", body)) - {"sec-pitfalls"})
    return {
        "front_sha256": sha(front.rstrip()),
        "content_line_hashes": sorted(lines.items()),
        "fences": dict(sorted(fence_languages(body).items())),
        "images": images,
        "citations": citations,
        "callouts": len(re.findall(r"^:::\s*\{\.callout-", body, flags=re.M)),
        "anchors": anchors,
        "headings": dict(sorted(semantic_headings(body, baseline=False).items())),
    }


def load_text(root: Path, relative: str, original_dir: Path | None) -> str:
    if original_dir:
        candidate = original_dir / relative
        if candidate.is_file():
            return candidate.read_text(encoding="utf-8")
    return (root / relative).read_text(encoding="utf-8")


def baseline_payload(root: Path, manifest: dict, original_dir: Path | None) -> dict:
    articles: dict[str, dict] = {}
    for chapter in manifest["series"]["chapters"]:
        number = str(int(chapter["number"]))
        relative = chapter["file"]
        text = load_text(root, relative, original_dir)
        item = metrics(text)
        _, body = split_front(text)
        item["headings"] = dict(sorted(semantic_headings(body, baseline=True).items()))
        articles[number] = {"file": relative, **item}
    return {"schema": 1, "articles": articles}


def subset_counter(expected_items: list[list], observed_items: list[list]) -> tuple[bool, int]:
    expected = Counter(dict(expected_items))
    observed = Counter(dict(observed_items))
    missing = sum(max(count - observed[key], 0) for key, count in expected.items())
    return missing == 0, missing


def style_checks(number: int, text: str) -> list[dict]:
    _, body = split_front(text)
    checks: list[tuple[str, bool, str]] = [
        ("no-generic-slot-h2", GENERIC_H2.search(body) is None, "legacy nine-slot H2 absent"),
        ("no-manual-heading-number", NUMBERED_HEADING.search(body) is None, "H2/H3 numbering absent"),
        ("key-takeaways", len(re.findall(r"^##\s+Key Takeaways(?:\s|\{|$)", body, re.M)) == 1, "one Key Takeaways"),
        ("methods-appendix", "## 附录 A：Methods / Results 模板" in body, "Appendix A"),
        ("own-data-appendix", "## 附录 B：换成你自己的数据" in body, "Appendix B"),
        ("references", re.search(r"^##\s+参考文献(?:\s|\{|$)", body, re.M) is not None, "references H2"),
        ("no-author-meta-copy", not any(token in body for token in META_PHRASES), "reader-facing prose only"),
        ("no-private-path", not any(token in body for token in PRIVATE_TOKENS), "no workstation identity"),
    ]
    if number != 71:
        checks.extend(
            [
                ("setup-anchor", "#sec-setup" in body or "#sec-preparation" in body, "website setup anchor"),
                ("code-anchor", "#sec-code" in body, "first executable step anchor"),
                ("theory-localized", re.search(r"^###\s+.+\{#sec-theory\}", body, re.M) is not None, "theory is inside a step"),
                ("audit-localized", re.search(r"^###\s+.+\{#sec-audit\}", body, re.M) is not None, "audit is inside a step"),
            ]
        )
    if number != 75:
        checks.append(("figure-appendix", "## 附录 C：" in body, "Appendix C"))
    return [
        {"name": name, "pass": passed, "detail": detail}
        for name, passed, detail in checks
    ]


def main() -> None:
    args = parse_args()
    root = args.project_root.resolve()
    manifest = yaml.safe_load((root / args.manifest).read_text(encoding="utf-8"))

    if manifest.get('reader_delivery', {}).get('version') == 2:
        from audit_reader_delivery import audit
        payload = audit(root)
        args.report.parent.mkdir(parents=True, exist_ok=True)
        args.report.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + '\n')
        print(json.dumps({k: v for k, v in payload.items() if k != 'reports'}, ensure_ascii=False))
        if args.strict and payload['errors']:
            raise SystemExit(1)
        return

    if args.write_baseline:
        payload = baseline_payload(root, manifest, args.original_dir)
        args.baseline.parent.mkdir(parents=True, exist_ok=True)
        args.baseline.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

    baseline = json.loads(args.baseline.read_text(encoding="utf-8"))
    reports: list[dict] = []
    for chapter in manifest["series"]["chapters"]:
        number = int(chapter["number"])
        relative = chapter["file"]
        text = (root / relative).read_text(encoding="utf-8")
        observed = metrics(text)
        expected = baseline["articles"][str(number)]
        content_ok, missing_lines = subset_counter(
            expected["content_line_hashes"], observed["content_line_hashes"]
        )
        preservation = [
            {"name": "front-matter", "pass": observed["front_sha256"] == expected["front_sha256"], "detail": "unchanged YAML"},
            {"name": "content-lines", "pass": content_ok, "detail": f"missing={missing_lines}"},
            {"name": "code-fences", "pass": observed["fences"] == expected["fences"], "detail": observed["fences"]},
            {"name": "images", "pass": observed["images"] == expected["images"], "detail": f"{len(observed['images'])} paths"},
            {"name": "citations", "pass": observed["citations"] == expected["citations"], "detail": f"{len(observed['citations'])} keys"},
            {"name": "callouts", "pass": observed["callouts"] == expected["callouts"], "detail": observed["callouts"]},
            {"name": "anchors", "pass": set(expected["anchors"]) <= set(observed["anchors"]), "detail": f"{len(observed['anchors'])} retained"},
            {
                "name": "heading-semantics",
                "pass": all(observed["headings"].get(title, 0) >= count for title, count in expected["headings"].items()),
                "detail": f"{len(expected['headings'])} original titles retained",
            },
        ]
        checks = preservation + style_checks(number, text)
        reports.append(
            {
                "article": number,
                "file": relative,
                "pass": all(item["pass"] for item in checks),
                "checks": checks,
            }
        )

    failures = [
        {"article": item["article"], "name": check["name"], "detail": check["detail"]}
        for item in reports
        for check in item["checks"]
        if not check["pass"]
    ]
    payload = {
        "schema": 1,
        "articles": len(reports),
        "passed_articles": sum(item["pass"] for item in reports),
        "failed_checks": len(failures),
        "failures": failures,
        "reports": reports,
    }
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({key: payload[key] for key in ("articles", "passed_articles", "failed_checks")}, ensure_ascii=False))
    if args.strict and failures:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
