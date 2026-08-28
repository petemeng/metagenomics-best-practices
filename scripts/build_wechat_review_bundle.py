#!/usr/bin/env python3
"""Build 16S-style WeChat review articles from the validated Quarto site.

The website remains the complete computational source.  This offline builder
removes generic environment bootstrap material, keeps decision-changing code
and all meaningful figures, applies the same reader-facing styles as the 16S
series, optimizes article images, and creates deterministic figure-only covers.

Uploading images, creating/replacing drafts, publishing, and mass sending are
separate operations and are never performed by this script.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
from copy import deepcopy
from datetime import datetime, timezone
from pathlib import Path
from typing import Any
from urllib.parse import unquote, urlsplit

import yaml
from lxml import etree, html
from PIL import Image, ImageOps


MAX_TITLE_CHARS = 64
MAX_DIGEST_CHARS = 120
MAX_THUMB_BYTES = 64 * 1024
MAX_ARTICLE_IMAGE_BYTES = 950 * 1024
MAX_ARTICLE_IMAGE_WIDTH = 1600

ROOT_STYLE = (
    "font-family:-apple-system,BlinkMacSystemFont,'PingFang SC',"
    "'Hiragino Sans GB','Noto Sans CJK SC',sans-serif;"
    "color:#2d2d2d;font-size:16px;line-height:1.85;"
    "word-break:break-word;"
)
STYLES = {
    "h2": (
        "margin:38px 0 18px;padding:10px 16px;border-left:4px solid #7c9970;"
        "background:#eef3ea;color:#203124;font-size:22px;line-height:1.45;"
    ),
    "h3": "margin:30px 0 14px;color:#314735;font-size:19px;line-height:1.5;",
    "h4": "margin:24px 0 12px;color:#4b5f4e;font-size:17px;line-height:1.5;",
    "p": "margin:0 0 18px;color:#2d2d2d;font-size:16px;line-height:1.85;",
    "ul": "margin:0 0 18px;padding-left:1.35em;color:#2d2d2d;line-height:1.85;",
    "ol": "margin:0 0 18px;padding-left:1.35em;color:#2d2d2d;line-height:1.85;",
    "li": "margin:0 0 10px;",
    "pre": (
        "margin:20px 0;padding:16px;background:#1f2a20;border-radius:12px;"
        "white-space:pre-wrap;word-break:break-all;color:#eef5ea;"
        "font-family:Menlo,Consolas,monospace;font-size:12px;line-height:1.7;"
    ),
    "code": (
        "padding:2px 5px;background:#f1ede4;border-radius:5px;"
        "font-family:Menlo,Consolas,monospace;font-size:0.9em;color:#8a5a1f;"
    ),
    "a": "color:#8a6428;text-decoration:none;border-bottom:1px solid #cfb17b;",
    "blockquote": (
        "margin:22px 0;padding:16px 18px 2px;background:#f6f1e7;"
        "border-left:4px solid #d7a35b;border-radius:10px;"
    ),
    "table": (
        "width:100%;margin:18px 0;border-collapse:collapse;table-layout:auto;"
        "font-size:13px;line-height:1.55;"
    ),
    "th": "padding:8px 6px;border:1px solid #d8d8d8;background:#eef3ea;text-align:left;",
    "td": "padding:8px 6px;border:1px solid #dedede;vertical-align:top;",
    "figure": "margin:22px 0;text-align:center;",
    "figcaption": (
        "margin:8px 8px 20px;color:#666;font-size:13px;"
        "line-height:1.65;text-align:left;"
    ),
    "img": (
        "display:block;width:100%;max-width:100%;height:auto;"
        "margin:18px auto;border-radius:10px;"
    ),
    "hr": "margin:28px auto;width:72px;height:1px;border:0;background:#d8ccb7;",
}
CALLOUT_STYLES = {
    "caution": "background:#fff4ed;border-left:4px solid #d9734e;",
    "warning": "background:#fff8e8;border-left:4px solid #d7a35b;",
    "important": "background:#f4eef8;border-left:4px solid #8c6ca8;",
    "tip": "background:#edf6f1;border-left:4px solid #5f9a7d;",
    "note": "background:#eef4f7;border-left:4px solid #5f879a;",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", default=".")
    parser.add_argument("--manifest", default="tutorial.yaml")
    parser.add_argument("--qa-report", default="qa_report.json")
    parser.add_argument("--site-dir", default="_site")
    parser.add_argument(
        "--output-dir",
        default="rendered/wechat_review_01_77_16s_style",
    )
    parser.add_argument(
        "--article-number",
        type=int,
        action="append",
        help="Build only this chapter; repeat to select multiple chapters.",
    )
    parser.add_argument("--author", default="Peter")
    parser.add_argument(
        "--review-url",
        help="Public source URL; defaults to publication.repo.public_url.",
    )
    parser.add_argument(
        "--review-only",
        action="store_true",
        help=(
            "Allow local style review when the current QA report is not passed. "
            "The output is marked non-draft-ready and must not be uploaded."
        ),
    )
    return parser.parse_args()


def utc_now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def truncate(text: str, limit: int) -> str:
    normalized = re.sub(r"\s+", " ", text).strip()
    if len(normalized) <= limit:
        return normalized
    return normalized[: max(1, limit - 1)].rstrip() + "…"


def class_tokens(element: etree._Element) -> set[str]:
    return set((element.get("class") or "").split())


def normalized_text(element: etree._Element) -> str:
    return re.sub(r"\s+", " ", "".join(element.itertext())).strip()


def site_html_path(site_dir: Path, qmd_path: str) -> Path:
    qmd = Path(qmd_path)
    if qmd.name == "index.qmd":
        return site_dir / "index.html"
    return site_dir / qmd.with_suffix(".html")


def wechat_draft_titles(manifest: dict[str, Any]) -> dict[int, str]:
    publication = manifest.get("publication", {}).get("wechat", {})
    prefix = str(publication.get("title_prefix", "")).strip()
    max_chars = int(publication.get("title_max_chars", 0))
    if prefix != "宏基因组最佳实践":
        raise RuntimeError("publication.wechat.title_prefix must equal 宏基因组最佳实践")
    if max_chars != MAX_TITLE_CHARS:
        raise RuntimeError(f"publication.wechat.title_max_chars must equal {MAX_TITLE_CHARS}")

    series = manifest.get("series", {})
    chapters = series.get("chapters", [])
    total = int(series.get("total_articles", 0))
    numbers = [int(item.get("number", 0)) for item in chapters]
    if numbers != list(range(1, total + 1)):
        raise RuntimeError("series.chapters must contain consecutive public numbers")

    titles: dict[int, str] = {}
    for item in chapters:
        number = int(item["number"])
        topic = str(item.get("wechat_title", item.get("title", ""))).strip()
        if not topic:
            raise RuntimeError(f"Article {number} has an empty WeChat topic title")
        title = f"{prefix}｜{number}. {topic}"
        if len(title) > max_chars:
            raise RuntimeError(
                f"Article {number} title has {len(title)} characters; maximum is {max_chars}"
            )
        titles[number] = title
    if len(set(titles.values())) != total:
        raise RuntimeError("WeChat public titles must be unique")
    return titles


def create_cover(representative_image: Path, output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    canvas = Image.new("RGB", (900, 383), "white")
    with Image.open(representative_image) as opened:
        figure = ImageOps.exif_transpose(opened)
        if figure.mode in {"RGBA", "LA"} or "transparency" in figure.info:
            rgba = figure.convert("RGBA")
            background = Image.new("RGBA", rgba.size, "white")
            background.alpha_composite(rgba)
            figure = background.convert("RGB")
        else:
            figure = figure.convert("RGB")
        figure = ImageOps.contain(figure, canvas.size, Image.Resampling.LANCZOS)
        canvas.paste(
            figure,
            ((canvas.width - figure.width) // 2, (canvas.height - figure.height) // 2),
        )
    for quality in (82, 74, 68, 62, 56, 50, 44, 38, 34):
        canvas.save(
            output,
            "JPEG",
            quality=quality,
            optimize=True,
            progressive=True,
            subsampling=2,
        )
        if output.stat().st_size <= MAX_THUMB_BYTES:
            return
    raise RuntimeError(f"Cover remains above {MAX_THUMB_BYTES} bytes: {output}")


def optimize_article_image(source: Path, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    with Image.open(source) as opened:
        image = ImageOps.exif_transpose(opened)
        if image.mode in {"RGBA", "LA"} or "transparency" in image.info:
            rgba = image.convert("RGBA")
            background = Image.new("RGBA", rgba.size, "white")
            background.alpha_composite(rgba)
            image = background.convert("RGB")
        else:
            image = image.convert("RGB")
        if image.width > MAX_ARTICLE_IMAGE_WIDTH:
            height = round(image.height * MAX_ARTICLE_IMAGE_WIDTH / image.width)
            image = image.resize(
                (MAX_ARTICLE_IMAGE_WIDTH, height),
                Image.Resampling.LANCZOS,
            )
        quality = 90
        while True:
            image.save(
                destination,
                "JPEG",
                quality=quality,
                optimize=True,
                progressive=True,
                subsampling=2,
            )
            if destination.stat().st_size <= MAX_ARTICLE_IMAGE_BYTES:
                return
            if quality > 64:
                quality -= 8
                continue
            new_width = max(900, round(image.width * 0.82))
            if new_width >= image.width:
                raise RuntimeError(f"Could not optimize article image: {source}")
            new_height = round(image.height * new_width / image.width)
            image = image.resize(
                (new_width, new_height),
                Image.Resampling.LANCZOS,
            )
            quality = 82


def remove_unwanted(main: etree._Element) -> None:
    selectors = [
        ".//script",
        ".//style",
        ".//button",
        ".//nav",
        '//*[@id="title-block-header"]',
        './/*[contains(concat(" ", normalize-space(@class), " "), " code-copy-button ")]',
        './/*[contains(concat(" ", normalize-space(@class), " "), " code-annotation-gutter ")]',
        './/*[contains(concat(" ", normalize-space(@class), " "), " anchorjs-link ")]',
        './/*[contains(concat(" ", normalize-space(@class), " "), " header-section-number ")]',
        './/*[contains(concat(" ", normalize-space(@class), " "), " quarto-title-meta ")]',
    ]
    seen: set[etree._Element] = set()
    for selector in selectors:
        for element in main.xpath(selector):
            if element in seen:
                continue
            seen.add(element)
            parent = element.getparent()
            if parent is not None:
                # ``drop_tree`` preserves the node tail.  Macro-metagenomics
                # pages place the reader-facing heading text in the tail of
                # Quarto's section-number span; ``parent.remove`` would erase
                # both the number and the actual heading.
                element.drop_tree()


def remove_wechat_bootstrap(main: etree._Element) -> int:
    """Remove generic website-only setup while retaining scientific decisions."""
    removed = 0
    for section in list(main.xpath(".//section")):
        headings = section.xpath("./h2[1]")
        section_id = section.get("id") or ""
        heading = normalized_text(headings[0]) if headings else ""
        if section_id not in {"sec-setup", "sec-preparation"} and heading != "准备工作":
            continue
        parent = section.getparent()
        if parent is not None:
            parent.remove(section)
            removed += 1

    for details in list(main.xpath(".//details")):
        summaries = details.xpath("./summary[1]")
        if not summaries:
            continue
        summary = normalized_text(summaries[0])
        if not (
            summary.startswith("展开：")
            and (
                "安装依赖" in summary
                or "定义作图函数" in summary
                or "出版级函数" in summary
            )
        ):
            continue
        parent = details.getparent()
        if parent is not None:
            parent.remove(details)
            removed += 1
    return removed


def remove_explicit_wechat_omissions(main: etree._Element) -> int:
    selector = (
        './/*[contains(concat(" ", normalize-space(@class), " "), '
        '" wechat-omit ")]'
    )
    elements = list(main.xpath(selector))
    marked = set(elements)
    removed = 0
    for element in elements:
        if any(ancestor in marked for ancestor in element.iterancestors()):
            continue
        parent = element.getparent()
        if parent is not None:
            parent.remove(element)
            removed += 1
    return removed


def prune_global_bibliography(main: etree._Element) -> tuple[int, int]:
    """Keep only references actually cited by the current chapter.

    The metagenomics book currently renders the complete repository-wide
    bibliography on every page.  The 16S articles carry chapter-scoped
    references, so reproducing that reader-facing behavior requires pruning
    uncited CSL entries before the HTML is packaged.
    """
    cited: set[str] = set()
    for link in main.xpath('.//a[@href]'):
        href = link.get("href") or ""
        fragment = urlsplit(href).fragment
        if fragment.startswith("ref-"):
            cited.add(fragment)
            # Relative document-fragment links are not useful in WeChat; the
            # citation label remains visible and the full DOI stays below.
            link.attrib.pop("href", None)

    entries = list(main.xpath('.//*[@id and starts-with(@id,"ref-")]'))
    removed = 0
    kept = 0
    for entry in entries:
        if entry.get("id") in cited:
            kept += 1
            continue
        entry.drop_tree()
        removed += 1

    for container in list(main.xpath('.//*[@id="refs"]')):
        if container.xpath('.//*[@id and starts-with(@id,"ref-")]'):
            continue
        section = next(
            (
                ancestor
                for ancestor in container.iterancestors()
                if ancestor.tag == "section"
            ),
            None,
        )
        target = section if section is not None else container
        target.drop_tree()
    return removed, kept


INSTALL_CALL = re.compile(
    r"^\s*(?:install\.packages|BiocManager::install|"
    r"remotes::install_github|pak::pkg_install)\s*\("
)


def strip_leading_install_calls(code_text: str) -> tuple[str, int]:
    lines = code_text.splitlines()
    prefix: list[str] = []
    index = 0
    while index < len(lines) and (
        not lines[index].strip() or lines[index].lstrip().startswith("#|")
    ):
        prefix.append(lines[index])
        index += 1

    removed = 0
    while index < len(lines) and INSTALL_CALL.match(lines[index]):
        depth = lines[index].count("(") - lines[index].count(")")
        index += 1
        while index < len(lines) and depth > 0:
            depth += lines[index].count("(") - lines[index].count(")")
            index += 1
        removed += 1
        while index < len(lines) and not lines[index].strip():
            index += 1

    if not removed:
        return code_text, 0
    return "\n".join(prefix + lines[index:]).rstrip(), removed


def flatten_code(main: etree._Element) -> int:
    stripped_install_calls = 0
    for pre in main.xpath(".//pre"):
        code_text = "".join(pre.itertext()).rstrip()
        code_text, removed = strip_leading_install_calls(code_text)
        stripped_install_calls += removed
        for child in list(pre):
            pre.remove(child)
        pre.text = code_text
        pre.set("style", STYLES["pre"])
    return stripped_install_calls


def transform_special_blocks(main: etree._Element) -> None:
    for details in main.xpath(".//details"):
        details.tag = "section"
        details.set(
            "style",
            "margin:22px 0;padding:16px 18px;background:#faf8f2;"
            "border:1px solid #e6dfd1;border-radius:10px;",
        )
    for summary in main.xpath(".//summary"):
        summary.tag = "p"
        summary.set(
            "style",
            "margin:0 0 14px;color:#314735;font-size:17px;"
            "line-height:1.6;font-weight:bold;",
        )
    for element in main.xpath(
        './/*[contains(concat(" ", normalize-space(@class), " "), " callout ")]'
    ):
        tokens = class_tokens(element)
        kind = next(
            (name for name in CALLOUT_STYLES if f"callout-{name}" in tokens),
            "note",
        )
        element.tag = "section"
        element.set(
            "style",
            "margin:22px 0;padding:16px 18px 4px;border-radius:10px;"
            + CALLOUT_STYLES[kind],
        )


def apply_inline_styles(main: etree._Element) -> None:
    for tag, style in STYLES.items():
        for element in main.xpath(f".//{tag}"):
            if (
                tag == "code"
                and element.getparent() is not None
                and element.getparent().tag == "pre"
            ):
                continue
            element.set("style", style)
    for element in main.xpath(".//strong"):
        element.set("style", "color:#203124;font-weight:700;")
    for element in main.xpath(".//em"):
        element.set("style", "color:#555;font-style:italic;")
    for element in main.xpath(".//div"):
        if not element.get("style"):
            element.tag = "section"


def local_image_source(source_html: Path, src: str) -> Path | None:
    if src.startswith(("http://", "https://", "data:")):
        return None
    clean = unquote(urlsplit(src).path)
    return (source_html.parent / clean).resolve()


def resolve_and_optimize_images(
    main: etree._Element,
    source_html: Path,
    article_dir: Path,
) -> list[dict[str, Any]]:
    records: list[dict[str, Any]] = []
    cache: dict[Path, tuple[str, Path]] = {}
    for image_element in main.xpath(".//img"):
        src = image_element.get("src") or ""
        source = local_image_source(source_html, src)
        if source is None:
            continue
        if not source.exists():
            raise FileNotFoundError(
                f"Missing rendered image {src} referenced by {source_html}"
            )
        if source not in cache:
            relative = f"images/image-{len(cache) + 1:02d}.jpg"
            destination = article_dir / relative
            optimize_article_image(source, destination)
            cache[source] = (relative, destination)
        relative, destination = cache[source]
        image_element.set("src", relative)
        image_element.set("style", STYLES["img"])
        with Image.open(destination) as opened:
            width, height = opened.size
        record = {
            "source_path": str(source),
            "local_path": str(destination.resolve()),
            "relative_src": relative,
            "sha256": sha256(destination),
            "size_bytes": destination.stat().st_size,
            "width": width,
            "height": height,
        }
        if not any(item["local_path"] == record["local_path"] for item in records):
            records.append(record)
    return records


def strip_unsupported_attributes(main: etree._Element) -> None:
    allowed = {"style", "href", "src", "alt", "title", "colspan", "rowspan"}
    for element in main.iter():
        for key in list(element.attrib):
            if key not in allowed:
                del element.attrib[key]
        if element.tag == "a":
            href = element.get("href") or ""
            if href.startswith(("#", "javascript:")):
                element.attrib.pop("href", None)


def sanitize_article(
    source_html: Path,
    article_dir: Path,
) -> tuple[str, list[dict[str, Any]], int, int, int, int, int]:
    document = html.parse(str(source_html)).getroot()
    mains = document.xpath(
        '//main[contains(concat(" ",normalize-space(@class)," ")," content ")]'
        ' | //main[@id="quarto-document-content"] | //main'
    )
    if not mains:
        raise RuntimeError(f"No article content found in {source_html}")
    main = deepcopy(mains[0])
    main.tag = "section"
    remove_unwanted(main)
    removed_references, retained_references = prune_global_bibliography(main)
    removed_omissions = remove_explicit_wechat_omissions(main)
    removed_bootstrap = remove_wechat_bootstrap(main)
    stripped_install_calls = flatten_code(main)
    transform_special_blocks(main)
    apply_inline_styles(main)
    images = resolve_and_optimize_images(main, source_html, article_dir)
    strip_unsupported_attributes(main)
    main.set("style", ROOT_STYLE)
    return (
        etree.tostring(main, encoding="unicode", method="html"),
        images,
        removed_bootstrap,
        removed_omissions,
        stripped_install_calls,
        removed_references,
        retained_references,
    )


def local_preview_html(title: str, content: str) -> str:
    return (
        "<!doctype html><html><head><meta charset='utf-8'>"
        "<meta name='viewport' content='width=device-width,initial-scale=1'>"
        f"<title>{title}</title></head>"
        "<body style='margin:0 auto;padding:24px;max-width:760px;'>"
        f"{content}</body></html>"
    )


def digest_for(chapter: dict[str, Any]) -> str:
    number = int(chapter["number"])
    topic = str(chapter.get("wechat_title", chapter["title"])).strip()
    lead = topic if topic.endswith(("。", "！", "？", "!", "?")) else f"{topic}。"
    text = f"{lead}真实数据、关键分析步骤、结果解释与发表级重绘图。"
    if number == 1:
        text = (
            f"{lead}从测量对象、分析层级和证据边界出发，理解 shotgun DNA "
            "能回答什么、不能回答什么。"
        )
    if number == 77:
        text = (
            f"{lead}把数据、代码、环境、图表和 MAG 提交整理为可核查的发表闭环。"
        )
    return truncate(text, MAX_DIGEST_CHARS)


def representative_cover_source(images: list[dict[str, Any]]) -> tuple[Path, str]:
    """Prefer an original analysis figure over a reproduced paper anchor."""
    reproduced_markers = {"original", "paper", "anchor", "source-figure"}
    for image in images:
        name = Path(str(image["source_path"])).name.lower()
        if not any(marker in name for marker in reproduced_markers):
            return Path(image["local_path"]), "first_original_analysis_figure"
    return Path(images[0]["local_path"]), "first_available_figure"


def selected_chapters(
    chapters: list[dict[str, Any]],
    requested: list[int] | None,
) -> list[dict[str, Any]]:
    if not requested:
        return chapters
    wanted = list(dict.fromkeys(requested))
    known = {int(item["number"]): item for item in chapters}
    missing = [number for number in wanted if number not in known]
    if missing:
        raise RuntimeError(f"Unknown article numbers: {missing}")
    return [known[number] for number in wanted]


def validate_item(
    item: dict[str, Any],
    payload: dict[str, Any],
    draft_ready: bool,
) -> list[str]:
    errors: list[str] = []
    chapter_id = item["chapter_id"]
    expected_prefix = f"宏基因组最佳实践｜{int(chapter_id)}. "
    content = str(payload["content"])
    fragment = html.fragment_fromstring(content, create_parent="section")
    if not payload["title"].startswith(expected_prefix):
        errors.append(f"{chapter_id}: title order prefix is missing")
    if len(payload["title"]) > MAX_TITLE_CHARS:
        errors.append(f"{chapter_id}: title is too long")
    if len(payload["digest"]) > MAX_DIGEST_CHARS:
        errors.append(f"{chapter_id}: digest is too long")
    if item["cover_size_bytes"] > MAX_THUMB_BYTES:
        errors.append(f"{chapter_id}: cover exceeds {MAX_THUMB_BYTES} bytes")
    if item["cover_width"] != 900 or item["cover_height"] != 383:
        errors.append(f"{chapter_id}: cover dimensions are not 900 x 383")
    if item["html_chars"] < 3000:
        errors.append(f"{chapter_id}: article content is unexpectedly short")
    if re.search(r"<(script|style|button|nav)\b", content, flags=re.I):
        errors.append(f"{chapter_id}: unsupported HTML remains")
    if re.search(
        r"<sub[^>]*></sub>\s*~\{|~\{(?:r|bash|python|text)\}",
        content,
        flags=re.I,
    ):
        errors.append(f"{chapter_id}: malformed executable-code fence remains")
    empty_headings = [
        heading
        for heading in fragment.xpath(".//h2 | .//h3 | .//h4")
        if not normalized_text(heading)
    ]
    if empty_headings:
        errors.append(f"{chapter_id}: empty reader-facing heading remains")
    retired_classes = {
        "article-header",
        "article-kicker",
        "fact-grid",
        "resource-card",
        "next-card",
    }
    if any(class_tokens(element) & retired_classes for element in fragment.iter()):
        errors.append(f"{chapter_id}: retired macro-only layout component remains")
    for link in fragment.xpath('.//a[@href]'):
        if not (link.get("href") or "").startswith(("https://", "http://", "mailto:")):
            errors.append(f"{chapter_id}: relative or fragment link remains")
            break
    if len(fragment.xpath('.//img[@src]')) != item["embedded_image_count"]:
        errors.append(f"{chapter_id}: embedded image ledger does not match HTML")
    if re.search(
        r"审阅草稿|开放审阅|GitHub Draft PR|草稿箱继续查看|"
        r"header-section-number|data-local-image|"
        r"宏基因组最佳实践[（(]\d{1,2}/\d{1,2}[）)]|"
        r"第\s*\d{1,2}\s*/\s*77\s*篇",
        content,
        flags=re.I,
    ):
        errors.append(f"{chapter_id}: internal review or numbering metadata remains")
    if re.search(
        r"<h2[^>]*>\s*准备工作\s*</h2>|"
        r"展开：[^<]{0,80}(?:安装依赖|定义作图函数|出版级函数)|"
        r"整仓库(?:运行时|使用者|用户|的一次性|的验收器)|"
        r"只复制(?:本页|本文)|单篇复现|独立运行以上",
        content,
        flags=re.I,
    ):
        errors.append(f"{chapter_id}: website-only bootstrap or maintainer prose remains")
    for code_block in re.findall(r"<pre\b[^>]*>.*?</pre>", content, flags=re.I | re.S):
        if re.search(
            r"(?:install\.packages|BiocManager::install|"
            r"remotes::install_github|pak::pkg_install)\s*\(",
            code_block,
        ):
            errors.append(f"{chapter_id}: package-install command remains in code")
            break
    if "/pull/" in str(payload.get("content_source_url", "")):
        errors.append(f"{chapter_id}: content_source_url points to a pull request")
    if payload.get("thumb_media_id") is not None:
        errors.append(f"{chapter_id}: local payload unexpectedly contains a media id")
    for image_record in item["embedded_images"]:
        if image_record["size_bytes"] > MAX_ARTICLE_IMAGE_BYTES:
            errors.append(f"{chapter_id}: article image exceeds upload budget")
    if draft_ready and item["build_mode"] != "draft_ready":
        errors.append(f"{chapter_id}: draft-ready build mode was not recorded")
    return errors


def build(args: argparse.Namespace) -> dict[str, Any]:
    project = Path(args.project_root).resolve()
    manifest_path = (project / args.manifest).resolve()
    qa_path = (project / args.qa_report).resolve()
    site_dir = (project / args.site_dir).resolve()
    output_dir = (project / args.output_dir).resolve()

    manifest = yaml.safe_load(manifest_path.read_text(encoding="utf-8"))
    qa_report = json.loads(qa_path.read_text(encoding="utf-8"))
    qa_passed = qa_report.get("status") == "passed"
    if not qa_passed and not args.review_only:
        raise RuntimeError(
            "qa_report.json.status must be passed before draft-ready generation; "
            "use --review-only solely for local style review"
        )

    chapters = manifest.get("series", {}).get("chapters", [])
    total = int(manifest.get("series", {}).get("total_articles", 0))
    if total != 77 or len(chapters) != total:
        raise RuntimeError(f"Expected 77 manifest chapters, found {len(chapters)}")
    titles = wechat_draft_titles(manifest)
    chosen = selected_chapters(chapters, args.article_number)
    review_url = args.review_url or manifest["publication"]["repo"]["public_url"]
    draft_ready = qa_passed and not args.review_only

    output_dir.mkdir(parents=True, exist_ok=True)
    items: list[dict[str, Any]] = []
    errors: list[str] = []
    for chapter in chosen:
        number = int(chapter["number"])
        qmd_path = str(chapter["file"])
        source_html = site_html_path(site_dir, qmd_path)
        if not source_html.exists():
            raise FileNotFoundError(f"Missing rendered article: {source_html}")
        article_dir = output_dir / f"{number:02d}"
        article_dir.mkdir(parents=True, exist_ok=True)
        (
            content,
            images,
            removed_bootstrap,
            removed_omissions,
            stripped_install_calls,
            removed_references,
            retained_references,
        ) = sanitize_article(source_html=source_html, article_dir=article_dir)
        if not images:
            raise RuntimeError(
                f"Article {number:02d} has no representative figure for its cover"
            )

        cover = article_dir / "cover.jpg"
        cover_source, cover_selection = representative_cover_source(images)
        create_cover(cover_source, cover)
        with Image.open(cover) as opened:
            cover_width, cover_height = opened.size

        payload = {
            "title": titles[number],
            "author": args.author,
            "digest": digest_for(chapter),
            "content": content,
            "content_source_url": review_url,
            "thumb_media_id": None,
            "show_cover_pic": 1,
            "need_open_comment": 0,
            "only_fans_can_comment": 0,
        }
        article_html = article_dir / "article.html"
        draft_json = article_dir / "draft.json"
        article_html.write_text(
            local_preview_html(titles[number], content),
            encoding="utf-8",
        )
        draft_json.write_text(
            json.dumps(payload, ensure_ascii=False, indent=2) + "\n",
            encoding="utf-8",
        )
        item = {
            "chapter_id": f"{number:02d}",
            "title": titles[number],
            "source_qmd": str((project / qmd_path).resolve()),
            "source_html": str(source_html),
            "article_html": str(article_html),
            "draft_json": str(draft_json),
            "cover_image": str(cover),
            "cover_layout": "representative_figure_only",
            "cover_selection": cover_selection,
            "cover_source_image": str(cover_source),
            "cover_size_bytes": cover.stat().st_size,
            "cover_sha256": sha256(cover),
            "cover_width": cover_width,
            "cover_height": cover_height,
            "html_chars": len(content),
            "removed_bootstrap_block_count": removed_bootstrap,
            "removed_wechat_omit_block_count": removed_omissions,
            "stripped_install_call_count": stripped_install_calls,
            "removed_uncited_reference_count": removed_references,
            "retained_cited_reference_count": retained_references,
            "embedded_image_count": len(images),
            "embedded_images": images,
            "build_mode": "draft_ready" if draft_ready else "review_only",
        }
        item_errors = validate_item(item, payload, draft_ready=draft_ready)
        item["errors"] = item_errors
        errors.extend(item_errors)
        (article_dir / "report.json").write_text(
            json.dumps(item, ensure_ascii=False, indent=2) + "\n",
            encoding="utf-8",
        )
        items.append(item)

    if len({item["title"] for item in items}) != len(items):
        errors.append("Selected draft titles are not unique")
    if not args.article_number and len(items) != 77:
        errors.append(f"Expected 77 items, found {len(items)}")

    status = "failed" if errors else ("passed" if draft_ready else "review_only_passed")
    report = {
        "generated_at": utc_now(),
        "status": status,
        "build_mode": "draft_ready" if draft_ready else "review_only",
        "upload_authorized": False,
        "publish_called": False,
        "mass_send_called": False,
        "project_root": str(project),
        "manifest": str(manifest_path),
        "manifest_sha256": sha256(manifest_path),
        "qa_report": str(qa_path),
        "qa_status": qa_report.get("status"),
        "qa_run_key": qa_report.get("run_key"),
        "qa_manifest_hash": qa_report.get("manifest_hash"),
        "qa_generated_at": qa_report.get("generated_at"),
        "author": args.author,
        "title_style": "宏基因组最佳实践｜N. ",
        "body_style_source": "16S reader-facing WeChat series",
        "review_url": review_url,
        "item_count": len(items),
        "selected_chapters": [item["chapter_id"] for item in items],
        "embedded_image_count": sum(item["embedded_image_count"] for item in items),
        "errors": errors,
        "items": items,
    }
    report_path = output_dir / "report.json"
    report_path.write_text(
        json.dumps(report, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    if errors:
        raise RuntimeError("; ".join(errors))
    return report


def main() -> int:
    report = build(parse_args())
    print(
        json.dumps(
            {
                "status": report["status"],
                "build_mode": report["build_mode"],
                "item_count": report["item_count"],
                "embedded_image_count": report["embedded_image_count"],
                "review_url": report["review_url"],
            },
            ensure_ascii=False,
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
