#!/usr/bin/env python3
"""Audit rendered website pages and the offline WeChat review bundle."""

from __future__ import annotations

import argparse
import json
import re
from html.parser import HTMLParser
from pathlib import Path
from urllib.parse import unquote, urlsplit

import yaml
from lxml import html as lxml_html


PRIVATE = ("/media/desk16/", "/home/tly9658/", "tly9658", "file://")
LEGACY_HEADINGS = (
    "这一步对应论文里的哪张图",
    "理论：",
    "准备工作",
    "可复制代码",
    "审计与升级",
    "出版级美化",
    "常见坑",
    "这段 Methods 怎么写",
    "换成你自己的数据怎么做",
    "真实论文为",
    "附录 A：Methods / Results 模板",
)


class PublicHTMLParser(HTMLParser):
    def __init__(self) -> None:
        super().__init__(convert_charrefs=True)
        self.headings: list[str] = []
        self.images: list[str] = []
        self.ids: set[str] = set()
        self._heading_tag: str | None = None
        self._heading_parts: list[str] = []

    def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        values = dict(attrs)
        if values.get("id"):
            self.ids.add(str(values["id"]))
        if tag in {"h2", "h3"}:
            self._heading_tag = tag
            self._heading_parts = []
        if tag == "img":
            self.images.append(str(values.get("src") or ""))

    def handle_startendtag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        self.handle_starttag(tag, attrs)

    def handle_data(self, data: str) -> None:
        if self._heading_tag:
            self._heading_parts.append(data)

    def handle_endtag(self, tag: str) -> None:
        if tag == self._heading_tag:
            self.headings.append(" ".join("".join(self._heading_parts).split()))
            self._heading_tag = None
            self._heading_parts = []


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, default=Path("."))
    parser.add_argument("--manifest", default="tutorial.yaml")
    parser.add_argument("--site-dir", default="_site")
    parser.add_argument("--bundle-dir", default="rendered/wechat-linear-review")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--strict", action="store_true")
    return parser.parse_args()


def local_asset(page: Path, root: Path, src: str) -> Path | None:
    parsed = urlsplit(src)
    if parsed.scheme or src.startswith("//") or src.startswith("data:"):
        return None
    path = unquote(parsed.path)
    if not path:
        return None
    if path.startswith("/"):
        return root / path.lstrip("/")
    return (page.parent / path).resolve()


def has_rendered_recommendation(text: str) -> bool:
    """Accept substantive topic-specific callouts on either public surface."""
    document = lxml_html.fromstring(text)
    candidates = document.xpath(
        '//*[contains(@class,"callout-tip") or contains(@class,"callout-important") '
        'or contains(@style,"border-left:4px solid #5f9a7d") '
        'or contains(@style,"border-left:4px solid #8c6ca8")]')
    return any(len(re.sub(r'\s+', '', ' '.join(
        p.text_content() for p in node.xpath('.//p|.//li')))) >= 20 for node in candidates)


def audit_html(page: Path, asset_root: Path, *, wechat: bool) -> list[str]:
    issues: list[str] = []
    if not page.is_file():
        return [f"missing HTML: {page}"]
    text = page.read_text(encoding="utf-8", errors="replace")
    for token in PRIVATE:
        if token in text:
            issues.append(f"private token {token!r} in {page}")
    parser = PublicHTMLParser()
    parser.feed(text)
    headings = [
        re.sub(r"^\d+(?:\.\d+)*\s+", "", heading).strip()
        for heading in parser.headings
    ]
    for heading in headings:
        if any(heading.startswith(token) for token in LEGACY_HEADINGS):
            issues.append(f"legacy heading {heading!r} in {page}")
    # Recommendations may be an early callout; do not require an English
    # slot heading or forbid a preparation ID that contains real input data.
    if not has_rendered_recommendation(text):
        issues.append(f"early decision recommendation missing in {page}")
    for source in parser.images:
        src = source.strip()
        if not src:
            issues.append(f"image without src in {page}")
            continue
        asset = local_asset(page, asset_root, src)
        if asset is not None and not asset.is_file():
            issues.append(f"missing image {src!r} referenced by {page}")
    return issues


def main() -> None:
    args = parse_args()
    root = args.project_root.resolve()
    site = (root / args.site_dir).resolve()
    bundle = (root / args.bundle_dir).resolve()
    manifest = yaml.safe_load((root / args.manifest).read_text(encoding="utf-8"))
    chapters = manifest["series"]["chapters"]
    public_prefix = manifest["publication"]["wechat"]["title_prefix"]

    issues: list[str] = []
    report_path = bundle / "report.json"
    if not report_path.is_file():
        issues.append(f"missing bundle report: {report_path}")
        report = {}
    else:
        report = json.loads(report_path.read_text(encoding="utf-8"))
        expected = {
            "status": "passed",
            "build_mode": "draft_ready",
            "upload_authorized": False,
            "publish_called": False,
            "mass_send_called": False,
            "item_count": 77,
        }
        for key, value in expected.items():
            if report.get(key) != value:
                issues.append(f"bundle report {key}={report.get(key)!r}; expected {value!r}")
        if report.get("errors"):
            issues.append(f"bundle report errors: {report['errors']}")

    expected_titles: set[str] = set()
    embedded_images = 0
    for chapter in chapters:
        number = int(chapter["number"])
        relative = Path(chapter["file"])
        site_page = site / ("index.html" if relative.name == "index.qmd" else relative.with_suffix(".html"))
        issues.extend(audit_html(site_page, site, wechat=False))

        article_dir = bundle / f"{number:02d}"
        article_html = article_dir / "article.html"
        draft_json = article_dir / "draft.json"
        cover = article_dir / "cover.jpg"
        for required in (article_html, draft_json, cover, article_dir / "report.json"):
            if not required.is_file() or required.stat().st_size == 0:
                issues.append(f"missing or empty bundle artifact: {required}")
        issues.extend(audit_html(article_html, article_dir, wechat=True))

        title_topic = str(chapter.get("wechat_title", chapter["title"])).strip()
        expected_title = f"{public_prefix}｜{number}. {title_topic}"
        expected_titles.add(expected_title)
        if draft_json.is_file():
            draft_text = draft_json.read_text(encoding="utf-8", errors="replace")
            for token in PRIVATE:
                if token in draft_text:
                    issues.append(f"private token {token!r} in {draft_json}")
            draft = json.loads(draft_text)
            if draft.get("title") != expected_title:
                issues.append(f"unexpected title in {draft_json}: {draft.get('title')!r}")

        item_report = article_dir / "report.json"
        if item_report.is_file():
            item = json.loads(item_report.read_text(encoding="utf-8"))
            if item.get("errors"):
                issues.append(f"item errors for {number:02d}: {item['errors']}")
            if item.get("build_mode") != "draft_ready":
                issues.append(f"item {number:02d} is not draft_ready")
            # Evidence chapters have no bootstrap to remove; technical
            # chapters may retain real data under the same setup section ID.
            # Validate the delivered boundary rather than a deletion count.
            if draft_json.is_file():
                public = json.loads(draft_json.read_text())["content"]
                if 'class="wechat-omit"' in public:
                    issues.append(f"item {number:02d} retained website-only material")
                if re.search(r'(theme_pub|save_pub)\s*(?:&lt;|<)-\s*function', public):
                    issues.append(f"item {number:02d} retained generic plot helpers")
            embedded_images += int(item.get("embedded_image_count", 0))

    report_titles = {str(item.get("title")) for item in report.get("items", [])}
    if report_titles != expected_titles:
        issues.append("bundle report title set does not match the 77 manifest titles")
    if embedded_images != int(report.get("embedded_image_count", -1)):
        issues.append(
            f"embedded image total mismatch: items={embedded_images}, report={report.get('embedded_image_count')}"
        )

    payload = {
        "status": "passed" if not issues else "failed",
        "articles": len(chapters),
        "site_pages": len(chapters),
        "wechat_items": len(chapters),
        "embedded_images": embedded_images,
        "issues": issues,
        "upload_authorized": False,
        "publish_called": False,
        "mass_send_called": False,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload, ensure_ascii=False, indent=2))
    if args.strict and issues:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
