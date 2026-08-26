#!/usr/bin/env python3
"""Derive a local, mobile-first WeChat preview from the validated Quarto page."""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import subprocess
from copy import deepcopy
from datetime import datetime, timezone
from html import escape
from pathlib import Path
from urllib.parse import unquote, urlsplit

import yaml
from lxml import etree
from lxml import html as lxml_html


ARTICLE_CSS = """
:root {
  color-scheme: light;
  --ink: #172033;
  --muted: #667085;
  --line: #e5e7eb;
  --blue: #0072b2;
  --orange: #d55e00;
  --paper: #ffffff;
  --soft: #f5f8fc;
}
* { box-sizing: border-box; }
body {
  margin: 0;
  background: #eef2f6;
  color: var(--ink);
  font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", "PingFang SC",
    "Hiragino Sans GB", "Microsoft YaHei", sans-serif;
  font-size: 17px;
  line-height: 1.82;
}
.wechat-article {
  width: min(100%, 760px);
  margin: 0 auto;
  padding: 34px 24px 64px;
  background: var(--paper);
}
.article-kicker {
  color: var(--blue);
  font-size: 13px;
  font-weight: 700;
  letter-spacing: .08em;
}
h1 {
  margin: 10px 0 14px;
  font-size: clamp(28px, 7vw, 42px);
  line-height: 1.24;
  letter-spacing: -.02em;
}
.article-deck { margin: 0 0 22px; color: var(--muted); }
.fact-grid {
  display: grid;
  grid-template-columns: repeat(2, minmax(0, 1fr));
  gap: 1px;
  margin: 22px 0 30px;
  overflow: hidden;
  border: 1px solid var(--line);
  border-radius: 14px;
  background: var(--line);
}
.fact-grid div { padding: 12px 14px; background: var(--soft); }
.fact-grid strong { display: block; font-size: 13px; color: var(--muted); }
.fact-grid span { font-weight: 700; }
h2 {
  margin: 46px 0 16px;
  padding-left: 12px;
  border-left: 4px solid var(--blue);
  font-size: 25px;
  line-height: 1.35;
}
h3 { margin: 30px 0 12px; font-size: 20px; line-height: 1.45; }
p { margin: 13px 0; }
a { color: var(--blue); text-decoration: none; }
strong { color: #101828; }
blockquote {
  margin: 20px 0;
  padding: 14px 17px;
  border-left: 4px solid var(--orange);
  background: #fff7ed;
}
figure { margin: 25px 0; }
img { display: block; width: 100%; height: auto; margin: 0 auto; }
figcaption { margin-top: 8px; color: var(--muted); font-size: 14px; text-align: center; }
.quarto-layout-row { display: block !important; }
.quarto-layout-cell { width: 100% !important; margin-bottom: 26px; }
.callout {
  margin: 22px 0;
  overflow: hidden;
  border: 1px solid #cfe4f3;
  border-radius: 12px;
  background: #f2f8fc;
}
.callout-header { padding: 11px 15px; font-weight: 700; color: #075985; }
.callout-body-container { padding: 0 15px 13px; }
.callout-caution { border-color: #fed7aa; background: #fff7ed; }
.callout-caution .callout-header { color: #9a3412; }
.screen-reader-only { display: none; }
details {
  margin: 20px 0;
  padding: 12px 15px;
  border: 1px solid var(--line);
  border-radius: 10px;
}
summary { cursor: pointer; font-weight: 700; }
table {
  display: block;
  width: 100%;
  margin: 18px 0;
  overflow-x: auto;
  border-collapse: collapse;
  font-size: 14px;
  line-height: 1.55;
}
th, td { min-width: 130px; padding: 9px 10px; border: 1px solid var(--line); }
th { background: var(--soft); text-align: left; }
pre {
  overflow-x: auto;
  padding: 14px;
  border-radius: 10px;
  background: #111827;
  color: #f9fafb;
  font-size: 13px;
  line-height: 1.55;
}
code {
  padding: .08em .28em;
  border-radius: 4px;
  background: #f2f4f7;
  color: #9a3412;
  font-family: ui-monospace, SFMono-Regular, Menlo, Consolas, monospace;
}
pre code { padding: 0; background: transparent; color: inherit; }
.resource-card, .next-card {
  margin: 22px 0;
  padding: 16px 18px;
  border-radius: 12px;
  background: var(--soft);
}
.resource-card ul { margin-bottom: 8px; padding-left: 1.25em; }
.resource-link {
  display: inline-block;
  margin-top: 8px;
  padding: 8px 14px;
  border-radius: 999px;
  background: var(--blue);
  color: white;
  font-weight: 700;
}
.next-card strong { color: var(--blue); }
.references { font-size: 14px; color: #475467; }
@media (max-width: 520px) {
  body { background: white; font-size: 16px; }
  .wechat-article { padding: 24px 17px 48px; }
  .fact-grid { grid-template-columns: 1fr; }
  h2 { font-size: 22px; }
}
""".strip()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", default=".")
    parser.add_argument("--manifest", default="tutorial.yaml")
    parser.add_argument("--qa-report", default="qa_report.json")
    parser.add_argument("--article-number", type=int, default=1)
    parser.add_argument("--site-html")
    parser.add_argument("--quiet", action="store_true")
    parser.add_argument(
        "--output-dir",
    )
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def chapter(manifest: dict, number: int) -> dict:
    for item in manifest["series"]["chapters"]:
        if int(item["number"]) == number:
            return item
    raise ValueError(f"Chapter {number} is absent from tutorial.yaml")


def qmd_frontmatter(path: Path) -> dict:
    text = path.read_text(encoding="utf-8")
    if not text.startswith("---\n"):
        raise ValueError(f"QMD frontmatter is missing: {path}")
    end = text.find("\n---\n", 4)
    if end < 0:
        raise ValueError(f"QMD frontmatter is not closed: {path}")
    metadata = yaml.safe_load(text[4:end])
    if not isinstance(metadata, dict):
        raise ValueError(f"QMD frontmatter is invalid: {path}")
    return metadata


def drop_matches(root: etree._Element, xpath: str) -> int:
    matches = list(root.xpath(xpath))
    for node in matches:
        node.drop_tree()
    return len(matches)


def replace_setup(
    section: etree._Element,
    repo_url: str,
    resources: list[dict],
) -> etree._Element:
    resource_items = "\n".join(
        "<li><strong>{label}：</strong>{value}</li>".format(
            label=escape(str(item["label"])),
            value=escape(str(item["value"])),
        )
        for item in resources
    )
    replacement = lxml_html.fragment_fromstring(
        f"""
        <section id="wechat-repro" class="level2">
          <h2>数据与复现入口</h2>
          <div class="resource-card">
            <ul>
              {resource_items}
            </ul>
            <a class="resource-link" href="{repo_url}">完整代码与环境文件</a>
          </div>
        </section>
        """
    )
    parent = section.getparent()
    parent.replace(section, replacement)
    return replacement


def serialize_children(root: etree._Element) -> str:
    return "\n".join(
        etree.tostring(child, encoding="unicode", method="html")
        for child in root
    )


def set_text_only(element: etree._Element, value: str) -> None:
    """Replace a heading without retaining Quarto numbering-node tails."""
    for child in list(element):
        element.remove(child)
    element.text = value


def append_style(element: etree._Element, declaration: str) -> None:
    existing = element.get("style", "").strip().rstrip(";")
    element.set(
        "style",
        f"{existing};{declaration}" if existing else declaration,
    )


def inline_wechat_styles(article: etree._Element) -> None:
    """Apply conservative inline CSS supported by the WeChat editor."""
    rules = [
        ('.//header[contains(@class,"article-header")]', "margin:0 0 26px 0;color:#172033;"),
        ('.//*[contains(@class,"article-kicker")]', "color:#0072b2;font-size:13px;font-weight:700;letter-spacing:1px;"),
        ('.//h1', "margin:10px 0 14px;font-size:30px;line-height:1.35;color:#172033;font-weight:700;"),
        ('.//*[contains(@class,"article-deck")]', "margin:0 0 20px;color:#667085;font-size:17px;line-height:1.75;"),
        ('.//*[contains(@class,"fact-grid")]', "margin:20px 0 28px;border:1px solid #e5e7eb;border-radius:12px;background:#f5f8fc;overflow:hidden;"),
        ('.//*[contains(@class,"fact-grid")]/div', "padding:11px 14px;border-bottom:1px solid #e5e7eb;background:#f5f8fc;"),
        ('.//*[contains(@class,"fact-grid")]//strong', "display:block;color:#667085;font-size:13px;font-weight:700;"),
        ('.//*[contains(@class,"fact-grid")]//span', "display:block;color:#172033;font-size:16px;font-weight:700;"),
        ('.//h2', "margin:40px 0 16px;padding-left:12px;border-left:4px solid #0072b2;color:#172033;font-size:23px;line-height:1.4;font-weight:700;"),
        ('.//h3', "margin:28px 0 12px;color:#172033;font-size:19px;line-height:1.5;font-weight:700;"),
        ('.//p[not(contains(@class,"article-deck"))]', "margin:12px 0;color:#172033;font-size:16px;line-height:1.82;"),
        ('.//a', "color:#0072b2;text-decoration:none;"),
        ('.//blockquote', "margin:20px 0;padding:14px 17px;border-left:4px solid #d55e00;background:#fff7ed;color:#172033;"),
        ('.//figure', "margin:24px 0;"),
        ('.//img', "display:block;width:100%;height:auto;margin:0 auto;"),
        ('.//figcaption', "margin-top:8px;color:#667085;font-size:14px;line-height:1.6;text-align:center;"),
        ('.//*[contains(@class,"quarto-layout-row")]', "display:block;width:100%;"),
        ('.//*[contains(@class,"quarto-layout-cell")]', "display:block;width:100%;margin:0 0 24px;"),
        ('.//*[contains(@class,"callout")]', "margin:20px 0;padding:0;border:1px solid #cfe4f3;border-radius:10px;background:#f2f8fc;overflow:hidden;"),
        ('.//*[contains(@class,"callout-header")]', "padding:10px 14px;color:#075985;font-size:16px;font-weight:700;"),
        ('.//*[contains(@class,"callout-body-container")]', "padding:0 14px 12px;"),
        ('.//*[contains(@class,"callout-caution")]', "border-color:#fed7aa;background:#fff7ed;"),
        ('.//*[contains(@class,"wechat-details")]', "display:block;margin:20px 0;padding:12px 15px;border:1px solid #e5e7eb;border-radius:10px;background:#fafafa;"),
        ('.//*[contains(@class,"wechat-details-title")]', "margin:0 0 10px;color:#172033;font-size:16px;font-weight:700;"),
        ('.//table', "width:100%;margin:18px 0;border-collapse:collapse;color:#172033;font-size:13px;line-height:1.55;"),
        ('.//th', "padding:8px;border:1px solid #d0d5dd;background:#f5f8fc;text-align:left;vertical-align:top;"),
        ('.//td', "padding:8px;border:1px solid #d0d5dd;text-align:left;vertical-align:top;"),
        ('.//code', "padding:1px 4px;border-radius:4px;background:#f2f4f7;color:#9a3412;font-family:Menlo,Consolas,monospace;font-size:0.9em;"),
        ('.//ul | .//ol', "padding-left:1.35em;color:#172033;line-height:1.8;"),
        ('.//li', "margin:5px 0;"),
        ('.//*[contains(@class,"resource-card") or contains(@class,"next-card")]', "margin:20px 0;padding:15px 17px;border-radius:10px;background:#f5f8fc;color:#172033;"),
        ('.//*[contains(@class,"resource-link")]', "display:inline-block;margin-top:8px;padding:7px 13px;border-radius:999px;background:#0072b2;color:#ffffff;font-weight:700;"),
        ('.//*[contains(@class,"references")]', "color:#475467;font-size:14px;line-height:1.7;"),
        ('.//*[contains(@class,"screen-reader-only")]', "display:none;"),
    ]
    for xpath, declaration in rules:
        for element in article.xpath(xpath):
            append_style(element, declaration)


def prepare_wechat_dom(article: etree._Element) -> None:
    """Expand unsupported widgets, clean links, and inline presentation."""
    for details in article.xpath('.//details'):
        details.tag = "section"
        details.set(
            "class",
            (details.get("class", "") + " wechat-details").strip(),
        )
        for summary in details.xpath('./summary'):
            summary.tag = "p"
            summary.set(
                "class",
                (summary.get("class", "") + " wechat-details-title").strip(),
            )

    drop_matches(
        article,
        './/*[contains(concat(" ", normalize-space(@class), " "), " callout-icon-container ")]',
    )
    for link in article.xpath('.//a[@href]'):
        href = link.get("href", "").strip()
        if not href.startswith(("#", "https://", "http://", "mailto:")):
            link.attrib.pop("href", None)
            link.tag = "span"
    for element in article.iter():
        for attribute in list(element.attrib):
            if attribute.startswith("data-"):
                element.attrib.pop(attribute, None)
    inline_wechat_styles(article)


def display_path(path: Path, root: Path) -> str:
    """Return a stable report path without assuming the output is under root."""
    try:
        return str(path.relative_to(root))
    except ValueError:
        return str(path)


def figure_source(root: Path, src: str) -> tuple[Path, Path]:
    """Resolve a rendered figure URL while preserving nested figure folders."""
    clean_src = unquote(urlsplit(src).path).replace("\\", "/")
    marker = "figures/"
    if marker not in clean_src:
        raise ValueError(f"Unexpected image remained in WeChat body: {src}")
    relative = Path(clean_src.split(marker, 1)[1])
    if relative.is_absolute() or ".." in relative.parts:
        raise ValueError(f"Unsafe figure path in WeChat body: {src}")
    return root / "figures" / relative, relative


def rewrite_public_headings(article: etree._Element, title: str) -> None:
    """Make standalone WeChat section labels reader-facing and topic-aware."""
    headings = {
        "sec-target": f"{title}：复现目标",
        "sec-audit": "证据边界与稳健性检查",
        "sec-polish": "让结果经得住审稿：作图与导出",
        "sec-publication": "让结果经得住审稿：作图与导出",
        "sec-pub": "让结果经得住审稿：作图与导出",
        "sec-style": "让结果经得住审稿：作图与导出",
        "sec-figure": "让结果经得住审稿：作图与导出",
        "sec-figures": "让结果经得住审稿：作图与导出",
        "sec-pitfalls": "哪些错误会改变结论",
        "sec-methods": "论文 Methods 应报告什么",
        "sec-own-data": "迁移到自己的数据",
        "sec-yourdata": "迁移到自己的数据",
    }
    for section_id, heading in headings.items():
        matches = article.xpath(f'.//section[@id="{section_id}"]/h2')
        if matches:
            set_text_only(matches[0], heading)


def main() -> int:
    args = parse_args()
    root = Path(args.project_root).resolve()
    manifest_path = (root / args.manifest).resolve()
    qa_path = (root / args.qa_report).resolve()

    manifest = yaml.safe_load(manifest_path.read_text(encoding="utf-8"))
    article_number = int(args.article_number)
    article_contract = chapter(manifest, article_number)
    qmd_path = (root / article_contract["file"]).resolve()
    metadata = qmd_frontmatter(qmd_path)
    wechat = metadata.get("wechat")
    if not isinstance(wechat, dict):
        raise SystemExit(f"Article {article_number:02d} has no wechat frontmatter contract.")

    if args.site_html:
        site_path = (root / args.site_html).resolve()
    else:
        site_relative = (
            Path("_site/index.html")
            if article_number == 1
            else Path("_site") / Path(article_contract["file"]).with_suffix(".html")
        )
        site_path = (root / site_relative).resolve()

    if args.output_dir:
        output_dir = (root / args.output_dir).resolve()
    else:
        output_slug = (
            "01-intro"
            if article_number == 1
            else f"{article_number:02d}-{qmd_path.stem.removeprefix(f'{article_number:02d}-')}"
        )
        output_dir = (root / "rendered/wechat_bundle" / output_slug).resolve()

    qa_report = json.loads(qa_path.read_text(encoding="utf-8"))
    if qa_report.get("status") != "passed":
        raise SystemExit("QA gate is not passed; refusing to create a WeChat preview.")
    if not site_path.is_file():
        raise SystemExit(f"Rendered source page does not exist: {site_path}")

    title = article_contract["title"]
    repo_url = manifest["publication"]["repo"]["public_url"]

    document = lxml_html.parse(str(site_path))
    source_main = document.xpath('//main[@id="quarto-document-content"]')
    if len(source_main) != 1:
        raise SystemExit("Expected one Quarto article main element.")
    article = deepcopy(source_main[0])

    drop_matches(article, './/header[@id="title-block-header"]')
    setup = article.xpath(
        './/section[@id="sec-setup" or @id="sec-preparation"]'
    )
    if len(setup) != 1:
        raise SystemExit("The source article is missing a preparation section.")
    resources = wechat.get("resources", [])
    if not isinstance(resources, list) or not resources:
        raise SystemExit("wechat.resources must contain at least one item.")
    intentionally_removed_images = len(setup[0].xpath('.//img[@src]'))
    replace_setup(setup[0], repo_url, resources)
    for section_id in wechat.get("strip_sections", []):
        stripped = article.xpath(f'.//section[@id="{section_id}"]')
        intentionally_removed_images += sum(
            len(section.xpath('.//img[@src]')) for section in stripped
        )
        drop_matches(article, f'.//section[@id="{section_id}"]')

    code_blocks_removed = drop_matches(article, './/div[contains(concat(" ", normalize-space(@class), " "), " cell-code ")]')
    code_blocks_removed += drop_matches(article, './/div[contains(concat(" ", normalize-space(@class), " "), " sourceCode ")]')
    drop_matches(article, './/button')
    drop_matches(article, './/div[contains(concat(" ", normalize-space(@class), " "), " code-copy-outer-scaffold ")][not(*) and not(normalize-space())]')
    drop_matches(article, './/div[contains(concat(" ", normalize-space(@class), " "), " cell ")][not(*) and not(normalize-space())]')
    drop_matches(article, './/span[contains(concat(" ", normalize-space(@class), " "), " header-section-number ")]')

    code_heading = article.xpath('.//section[@id="sec-code"]/h2')
    if code_heading and wechat.get("code_heading"):
        set_text_only(code_heading[0], str(wechat["code_heading"]))

    refs = article.xpath('.//section[@id="sec-refs"]')
    if refs:
        refs[0].set("class", (refs[0].get("class", "") + " references").strip())

    assets_dir = output_dir / "assets"
    assets_dir.mkdir(parents=True, exist_ok=True)
    image_records: list[dict] = []
    for image in article.xpath(".//img[@src]"):
        src = image.get("src", "")
        try:
            source_image, figure_relative = figure_source(root, src)
        except ValueError as exc:
            raise SystemExit(str(exc)) from exc
        if not source_image.is_file():
            raise SystemExit(f"Image does not exist: {source_image}")
        target_image = assets_dir / figure_relative
        target_image.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source_image, target_image)
        image.set("src", f"assets/{figure_relative.as_posix()}")
        image.set(
            "alt",
            str(
                wechat.get("image_alt", {}).get(
                    figure_relative.as_posix(),
                    wechat.get("image_alt", {}).get(
                        target_image.name,
                        target_image.stem,
                    ),
                )
            ),
        )
        image_records.append(
            {
                "path": display_path(target_image, root),
                "size_bytes": target_image.stat().st_size,
                "sha256": sha256(target_image),
            }
        )

    if article_number < int(manifest["series"]["total_articles"]):
        next_article = chapter(manifest, article_number + 1)
        next_card = lxml_html.fragment_fromstring(
            f"""
            <div class="next-card">
              <strong>下一篇预告</strong>
              <p>{escape(str(next_article["title"]))}</p>
            </div>
            """
        )
        if refs:
            refs[0].addprevious(next_card)
        else:
            article.append(next_card)

    for element in article.iter():
        element.attrib.pop("data-number", None)
        element.attrib.pop("data-anchor-id", None)
        element.attrib.pop("style", None)

    rewrite_public_headings(article, str(title))

    facts = wechat.get("facts", [])
    if not isinstance(facts, list) or len(facts) != 4:
        raise SystemExit("wechat.facts must contain exactly four items.")
    fact_html = "\n".join(
        "<div><strong>{label}</strong><span>{value}</span></div>".format(
            label=escape(str(item["label"])),
            value=escape(str(item["value"])),
        )
        for item in facts
    )
    header = f"""
    <header class="article-header">
      <div class="article-kicker">宏基因组分析最佳实践 · Shotgun</div>
      <h1>{escape(str(title))}</h1>
      <p class="article-deck">{escape(str(wechat["deck"]))}</p>
      <div class="fact-grid">
        {fact_html}
      </div>
    </header>
    """.strip()
    article.insert(0, lxml_html.fragment_fromstring(header))
    prepare_wechat_dom(article)
    body = serialize_children(article)
    standalone = f"""<!doctype html>
<html lang="zh-CN">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>{escape(str(title))}</title>
  <style>{ARTICLE_CSS}</style>
</head>
<body>
  <article class="wechat-article">
    {body}
  </article>
</body>
</html>
"""

    article_html = output_dir / "article.html"
    body_html = output_dir / "article-body.html"
    article_md = output_dir / "article.md"
    article_html.write_text(standalone, encoding="utf-8")
    body_html.write_text(body, encoding="utf-8")

    pandoc = shutil.which("pandoc")
    if pandoc:
        completed = subprocess.run(
            [pandoc, "--from=html", "--to=gfm", "--wrap=none", str(body_html), "-o", str(article_md)],
            check=False,
            capture_output=True,
            text=True,
        )
        if completed.returncode != 0:
            raise SystemExit(f"Pandoc failed while creating article.md: {completed.stderr}")

    rendered_text = article_html.read_text(encoding="utf-8")
    expected_images = int(wechat.get("expected_images", 0))
    expected_retained_images = int(
        wechat.get("expected_retained_images", expected_images)
    )
    checks = {
        "qa_gate_passed": qa_report.get("status") == "passed",
        "has_title": title in rendered_text,
        "has_expected_images": (
            rendered_text.count("<img ") == expected_retained_images
        ),
        "has_no_quarto_source_code": "sourceCode" not in rendered_text and "cell-code" not in rendered_text,
        "has_reproduction_link": repo_url in rendered_text,
        "live_publish_disabled": True,
    }
    status = "passed" if all(checks.values()) else "failed"
    report = {
        "status": status,
        "article_number": article_number,
        "generated_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "source_qmd": str(qmd_path.relative_to(root)),
        "source_html": str(site_path.relative_to(root)),
        "qa_report": str(qa_path.relative_to(root)),
        "qa_run_key": qa_report.get("run_key"),
        "article": display_path(article_html, root),
        "markdown": display_path(article_md, root) if article_md.exists() else None,
        "code_blocks_removed": code_blocks_removed,
        "source_images_expected": expected_images,
        "images_intentionally_removed": intentionally_removed_images,
        "images_retained": len(image_records),
        "images": image_records,
        "checks": checks,
    }
    report_path = output_dir / "report.json"
    report_path.write_text(
        json.dumps(report, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    if not args.quiet:
        print(json.dumps(report, ensure_ascii=False, indent=2))
    return 0 if status == "passed" else 1


if __name__ == "__main__":
    raise SystemExit(main())
