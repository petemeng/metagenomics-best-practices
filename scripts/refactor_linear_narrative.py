#!/usr/bin/env python3
"""Refactor legacy nine-slot tutorials into evidence-led linear narratives.

The transformation is intentionally structural: it preserves prose, code,
figures, citations, chunk labels, and front matter.  Theory, audit notes, and
pitfall callouts are matched to the analysis step they constrain.  Generic
publication slots become appendices, while the website-only preparation
section keeps its historical id so the WeChat builder can omit bootstrap code.
"""

from __future__ import annotations

import argparse
import json
import math
import re
import shutil
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable

import yaml


H2_RE = re.compile(r"^##\s+(.+?)\s*$", re.M)
H3_RE = re.compile(r"^###\s+(.+?)\s*$", re.M)
HEADING_RE = re.compile(r"^(#{3,6})(\s+)", re.M)
ID_RE = re.compile(r"\s*\{#([^}]+)\}\s*$")
NUMBER_PREFIX_RE = re.compile(
    r"^(?:第\s*\d+\s*步\s*[：:｜|]?\s*|"
    r"步骤\s*\d+\s*[：:｜|]?\s*|"
    r"\d+(?:\.\d+)+(?:[.、：:]\s*|\s+)|"
    r"\d+(?:[、：:]|\.(?!\d))\s*|"
    r"图\s*\d+\s*[：:]\s*)"
)

GENERIC_ROLE_PATTERNS = {
    "target": ("sec-target", r"^这一步对应论文里的哪张图"),
    "theory": ("sec-theory", r"^理论[：:]"),
    "setup": ("sec-setup", r"^准备工作"),
    "setup2": ("sec-preparation", r"^准备工作"),
    "code": ("sec-code", r"^可复制代码"),
    "audit": ("sec-audit", r"(?:审计与升级|从.+审计与升级)"),
    "publication": (
        "sec-polish|sec-publication|sec-pub|sec-figure|sec-figures|sec-style",
        r"^出版级美化",
    ),
    "pitfalls": ("sec-pitfalls", r"^常见坑"),
    "methods": ("sec-methods", r"^这段 Methods 怎么写"),
    "own": (
        "sec-own-data|sec-yourdata",
        r"^换成你自己的数据怎么做",
    ),
    "summary": ("", r"^(?:小结|总结)$"),
    "references": ("sec-references|sec-refs", r"^参考(?:文献)?$"),
}

STOP_TOKENS = {
    "什么", "为什么", "怎么", "怎样", "这里", "本章", "本文", "这个",
    "一个", "一种", "进行", "结果", "数据", "分析", "方法", "步骤", "代码",
    "审计", "升级", "理论", "准备", "工作", "可复制", "出版", "常见", "问题",
    "需要", "不能", "不是", "以及", "into", "from", "with", "using", "and",
    "the", "for", "that", "this", "what", "why", "how",
}

GENERIC_LATIN = {
    "and", "the", "for", "from", "with", "into", "using", "this", "that",
    "data", "result", "results", "analysis", "method", "model", "quality",
    "high", "low", "gene", "genome", "sample", "samples", "step", "code",
    "figure", "table", "group", "groups", "human", "primary", "current",
    "lineage", "value", "values", "effect", "effects", "path", "paths",
}

# Small, reviewable exceptions found during pilot reading.  They resolve
# domain synonyms that lexical matching cannot infer (for example Alpha
# diversity is implemented as Hill numbers).  Keeping them article-scoped
# prevents an innocent word such as "quality" from redirecting other chapters.
PLACEMENT_HINTS: dict[int, list[tuple[str, str]]] = {
    13: [
        (r"Phred", r"FASTQ 前缀|raw FastQC"),
        (r"末端下降|trimming|分别截取|分别过滤", r"paired-end.*fastp"),
        (r"duplication|GC.*WARN", r"raw FastQC|clean FASTQ"),
        (r"文件前缀|prefix 当随机", r"FASTQ 前缀"),
        (r"参数敏感性", r"clean FASTQ"),
        (r"clean reads 还不是", r"完整的一次性流程|离线验证"),
    ],
    22: [
        (r"Alpha|有效特征数|丰富度", r"Hill"),
        (r"relative abundance.*rarefy", r"closure"),
    ],
    44: [
        (r"barrnap|tRNAscan|MIMAG high-quality|完整 MIMAG|>90/<5", r"rRNA/tRNA/CDS"),
        (r"contamination.*chimerism", r"CheckM2.*GUNC"),
        (r"Archaea.*mode", r"rRNA/tRNA/CDS"),
    ],
}


@dataclass
class Section:
    title: str
    body: str
    role: str = "other"
    anchor: str | None = None


@dataclass
class Subsection:
    title: str
    body: str
    anchor: str | None = None
    prefix: list[str] = field(default_factory=list)
    suffix: list[str] = field(default_factory=list)
    pitfalls: list[str] = field(default_factory=list)

    @property
    def searchable(self) -> str:
        return f"{self.title}\n{self.body[:5000]}"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, default=Path("."))
    parser.add_argument("--manifest", default="tutorial.yaml")
    parser.add_argument("--article-number", type=int, action="append")
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--write", action="store_true")
    parser.add_argument("--backup-dir", type=Path)
    return parser.parse_args()


def split_front_matter(text: str) -> tuple[str, str]:
    if not text.startswith("---\n"):
        raise ValueError("QMD does not start with YAML front matter")
    end = text.find("\n---\n", 4)
    if end < 0:
        raise ValueError("Unclosed YAML front matter")
    return text[: end + 5], text[end + 5 :]


def heading_parts(title: str) -> tuple[str, str | None]:
    match = ID_RE.search(title)
    anchor = match.group(1) if match else None
    clean = ID_RE.sub("", title).strip()
    return clean, anchor


def clean_title(title: str) -> str:
    clean, _ = heading_parts(title)
    # Remove exactly one structural prefix.  Repeating the substitution can
    # destroy a real leading quantity: "5.3 18 道发布门" must become
    # "18 道发布门", not "道发布门".
    clean = NUMBER_PREFIX_RE.sub("", clean, count=1).strip()
    return re.sub(r"^审计与升级\s*[：:]\s*", "", clean).strip()


def with_anchor(title: str, anchor: str | None) -> str:
    return f"{title} {{#{anchor}}}" if anchor else title


def role_for(title: str, anchor: str | None) -> str:
    for role, (anchor_pattern, title_pattern) in GENERIC_ROLE_PATTERNS.items():
        if anchor and anchor_pattern and re.fullmatch(anchor_pattern, anchor):
            return "setup" if role == "setup2" else role
        if re.search(title_pattern, title, flags=re.I):
            return "setup" if role == "setup2" else role
    return "other"


def parse_h2_sections(body: str) -> tuple[str, list[Section]]:
    matches = list(H2_RE.finditer(body))
    if not matches:
        return body, []
    intro = body[: matches[0].start()].rstrip()
    sections: list[Section] = []
    for index, match in enumerate(matches):
        end = matches[index + 1].start() if index + 1 < len(matches) else len(body)
        title, anchor = heading_parts(match.group(1))
        section_body = body[match.end() : end].strip("\n")
        sections.append(
            Section(
                title=title,
                body=section_body,
                role=role_for(title, anchor),
                anchor=anchor,
            )
        )
    return intro, sections


def split_h3(body: str) -> tuple[str, list[Subsection]]:
    matches = list(H3_RE.finditer(body))
    if not matches:
        return body.strip(), []
    intro = body[: matches[0].start()].strip()
    subsections: list[Subsection] = []
    for index, match in enumerate(matches):
        end = matches[index + 1].start() if index + 1 < len(matches) else len(body)
        title, anchor = heading_parts(match.group(1))
        subsections.append(
            Subsection(
                title=title,
                body=body[match.end() : end].strip(),
                anchor=anchor,
            )
        )
    return intro, subsections


def rewrite_headings(body: str, *, promote: bool = False) -> str:
    """Clean visible heading prefixes without touching fenced code comments."""
    lines: list[str] = []
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
            lines.append(line)
            continue
        if not in_fence:
            match = re.match(r"^(#{2,6})\s+(.+?)\s*$", line)
            if match:
                title, anchor = heading_parts(match.group(2))
                level = len(match.group(1)) - (1 if promote and len(match.group(1)) >= 3 else 0)
                line = qmd_heading(level, clean_title(title), anchor)
        lines.append(line)
    return "\n".join(lines)


def promote_nested_headings(body: str) -> str:
    return rewrite_headings(body, promote=True)


def tokens(text: str) -> Counter[str]:
    text = re.sub(r"```.*?```", " ", text, flags=re.S)
    latin = [token.lower() for token in re.findall(r"[A-Za-z][A-Za-z0-9_.+-]{2,}", text)]
    cjk_chunks = re.findall(r"[\u3400-\u9fff]{2,}", text)
    cjk: list[str] = []
    for chunk in cjk_chunks:
        cjk.extend(chunk[index : index + 2] for index in range(len(chunk) - 1))
        cjk.extend(chunk[index : index + 3] for index in range(len(chunk) - 2))
    result = Counter(token for token in latin + cjk if token not in STOP_TOKENS)
    return result


def similarity(left: str, right: str) -> float:
    a = tokens(left)
    b = tokens(right)
    if not a or not b:
        return 0.0
    shared = set(a) & set(b)
    numerator = sum(min(a[token], 3) * min(b[token], 3) for token in shared)
    denominator = math.sqrt(
        sum(min(value, 3) ** 2 for value in a.values())
        * sum(min(value, 3) ** 2 for value in b.values())
    )
    return numerator / denominator if denominator else 0.0


def latin_terms(text: str) -> list[str]:
    return [
        token.lower()
        for token in re.findall(r"[A-Za-z][A-Za-z0-9_.+-]{2,}", text)
        if token.lower() not in GENERIC_LATIN
    ]


def chinese_phrases(text: str) -> list[str]:
    """Return informative title fragments, not every noisy character n-gram."""
    text = re.sub(r"[`*'\"“”‘’（）()<>]", " ", text)
    chunks = re.split(r"[：:，,；;。！？!?+×/|→←=]", text)
    phrases: list[str] = []
    for chunk in chunks:
        chunk = re.sub(r"[^\u3400-\u9fff]", "", chunk)
        chunk = re.sub(
            r"^(?:先|再|把|用|从|到|与|和|的|一个|一种|当前|本例|本章)+",
            "",
            chunk,
        )
        if len(chunk) >= 4:
            phrases.append(chunk)
    return phrases


def placement_score(query_title: str, query_body: str, step: Subsection) -> float:
    """Score where an explanatory block belongs.

    Titles carry most of the signal.  Exact tool/statistic names are weighted
    much more heavily than generic prose; the beginning of the step body is a
    secondary signal because commands often contain the decisive name.
    """
    step_title = step.title.lower()
    step_body = step.body[:5000].lower()
    score = 4.0 * similarity(query_title, step.title)
    score += 0.8 * similarity(query_title, step.body[:2500])
    score += 0.35 * similarity(query_body[:1200], step.title)

    identifier_terms = {
        token.lower()
        for token in re.findall(r"[A-Za-z][A-Za-z0-9_.+-]{2,}", query_title)
        if bool(re.search(r"[0-9_.+-]", token))
        or (sum(char.isupper() for char in token) >= 2)
    }
    for term in set(latin_terms(query_title)):
        is_identifier = term in identifier_terms
        weight = 7.0 if is_identifier else 3.5
        if re.search(rf"(?<![a-z0-9]){re.escape(term)}(?![a-z0-9])", step_title):
            score += weight
        elif re.search(rf"(?<![a-z0-9]){re.escape(term)}(?![a-z0-9])", step_body):
            score += weight * 0.45

    for phrase in chinese_phrases(query_title):
        if phrase in re.sub(r"[^\u3400-\u9fff]", "", step.title):
            score += 4.0
        elif phrase in re.sub(r"[^\u3400-\u9fff]", "", step.body[:5000]):
            score += 1.2
    return score


def best_step(
    title: str,
    body: str,
    steps: list[Subsection],
    *,
    article_number: int,
    block_index: int,
    block_count: int,
) -> tuple[int, float, str]:
    query = f"{title}\n{body[:1200]}"
    for query_pattern, step_pattern in PLACEMENT_HINTS.get(article_number, []):
        if not re.search(query_pattern, query, flags=re.I | re.S):
            continue
        matches = [
            index
            for index, step in enumerate(steps)
            if re.search(step_pattern, step.title, flags=re.I)
        ]
        if matches:
            return matches[0], 99.0, "manual_hint"

    scores = [placement_score(title, body, step) for step in steps]
    index = max(range(len(scores)), key=scores.__getitem__)
    best = scores[index]

    # When a title has no usable lexical anchor, preserve narrative order.
    # This is safer than attaching a conceptual paragraph to an arbitrary
    # late-stage plotting or export step because of vocabulary in its body.
    if best < 0.42:
        expected = round((block_index + 0.5) * len(steps) / max(block_count, 1) - 0.5)
        index = min(max(expected, 0), len(steps) - 1)
        return index, best, "order_fallback"
    return index, best, "lexical"


def parse_callouts(body: str) -> tuple[list[str], str]:
    pattern = re.compile(
        r"(?ms)^:::\s*\{\.callout-caution[^}]*\}.*?^:::\s*$"
    )
    callouts = [match.group(0).strip() for match in pattern.finditer(body)]
    remainder = pattern.sub("", body).strip()
    return callouts, remainder


def contains_executable_fence(text: str) -> bool:
    """Whether moving the block can change notebook execution order."""
    return re.search(r"(?m)^\s*(?:```|~~~)\{(?:r|python|bash|sh|julia)\b", text) is not None


def callout_title(callout: str) -> str:
    match = re.search(r'title="([^"]+)"', callout)
    title = match.group(1) if match else "结果解释边界"
    return re.sub(r"^(?:坑|误读)\s*\d+\s*[：:]\s*", "", title)


def setup_heading(section: Section, subsections: list[Subsection]) -> str:
    corpus = " ".join([section.body[:1200], *(item.title for item in subsections)]).lower()
    if any(token in corpus for token in ("ram", "算力", "磁盘", "耗时", "资源门槛")):
        return "先锁定输入、版本和算力边界"
    if any(token in corpus for token in ("wsl", "conda", "mamba", "安装", "环境")):
        return "先确认系统环境与版本来源"
    return "先确认数据来源与运行边界"


def target_heading(topic: str, subsections: list[Subsection]) -> str:
    if subsections:
        return subsections[0].title
    short_topic = topic.split("：", 1)[0].strip()
    return f"真实论文为“{short_topic}”提供了什么证据起点"


def qmd_heading(level: int, title: str, anchor: str | None = None) -> str:
    return f"{'#' * level} {with_anchor(title, anchor)}"


def render_subsection(
    subsection: Subsection,
    *,
    level: int,
    anchor: str | None = None,
) -> str:
    parts = [qmd_heading(level, subsection.title, anchor or subsection.anchor)]
    if subsection.prefix:
        parts.extend(subsection.prefix)
    if subsection.body:
        body = promote_nested_headings(subsection.body) if level == 2 else subsection.body
        parts.append(body)
    if subsection.suffix:
        parts.extend(subsection.suffix)
    if subsection.pitfalls:
        parts.extend(subsection.pitfalls)
    return "\n\n".join(part.strip() for part in parts if part.strip())


def yaml_metadata(front_matter: str) -> dict:
    return yaml.safe_load(front_matter.split("---", 2)[1]) or {}


def topic_from_metadata(metadata: dict) -> str:
    title = str(metadata.get("title", "宏基因组分析"))
    title = re.sub(r"^第\s*\d+\s*篇\s*[·・]\s*", "", title)
    return title


def make_takeaways(metadata: dict, existing: Section | None) -> str:
    if existing and existing.body.strip():
        return qmd_heading(2, "Key Takeaways", existing.anchor) + "\n\n" + existing.body.strip()
    wechat = metadata.get("wechat") or {}
    deck = str(wechat.get("deck", "")).strip()
    facts = wechat.get("facts") or []
    lines: list[str] = []
    if deck:
        lines.append(deck)
    for item in facts[:4]:
        if isinstance(item, dict):
            label = str(item.get("label", "结果")).strip()
            value = str(item.get("value", "")).strip()
            if value:
                lines.append(f"- **{label}：**{value}")
    if not lines:
        lines.append("本章的结论应与实际展示的数据、分析单位和不确定性范围一起阅读。")
    return qmd_heading(2, "Key Takeaways") + "\n\n" + "\n\n".join(lines)


def transform(text: str, article_number: int) -> tuple[str, dict]:
    front, body = split_front_matter(text)
    metadata = yaml_metadata(front)
    topic = topic_from_metadata(metadata)
    intro, sections = parse_h2_sections(body)
    by_role: dict[str, list[Section]] = {}
    for section in sections:
        by_role.setdefault(section.role, []).append(section)

    if not by_role.get("code"):
        return text, {"article": article_number, "status": "skipped", "reason": "no code section"}

    report: dict = {
        "article": article_number,
        "topic": topic,
        "status": "planned",
        "theory_matches": [],
        "audit_matches": [],
        "pitfall_matches": [],
    }

    output: list[str] = [intro.strip()]

    # Evidence anchor: unwrap its first descriptive subheading when available.
    for section in by_role.get("target", []):
        section_intro, target_subsections = split_h3(section.body)
        if target_subsections:
            first = target_subsections[0]
            first.body = "\n\n".join(
                part for part in (section_intro, first.body) if part.strip()
            )
            for index, subsection in enumerate(target_subsections):
                output.append(
                    render_subsection(
                        subsection,
                        level=2,
                        anchor=section.anchor if index == 0 else None,
                    )
                )
        else:
            output.append(
                qmd_heading(2, target_heading(topic, target_subsections), section.anchor)
                + "\n\n"
                + section.body.strip()
            )

    # Website preparation remains visible and reproducible; the historical id
    # lets the WeChat builder remove generic environment bootstrap material.
    for section in by_role.get("setup", []):
        section_intro, setup_subsections = split_h3(section.body)
        parts = [qmd_heading(2, setup_heading(section, setup_subsections), section.anchor)]
        if section_intro:
            parts.append(section_intro)
        for subsection in setup_subsections:
            parts.append(render_subsection(subsection, level=3))
        output.append("\n\n".join(parts))

    code_section = by_role["code"][0]
    code_intro, steps = split_h3(code_section.body)
    if not steps:
        steps = [Subsection(title=code_section.title.split("：", 1)[-1].strip(), body=code_section.body)]
        code_intro = ""
    if code_intro:
        steps[0].body = "\n\n".join(part for part in (code_intro, steps[0].body) if part.strip())

    theory_blocks: list[Subsection] = []
    for section in by_role.get("theory", []):
        theory_intro, subsections = split_h3(section.body)
        if theory_intro:
            if subsections:
                subsections[0].body = "\n\n".join((theory_intro, subsections[0].body))
            else:
                subsections = [Subsection(title=section.title.split("：", 1)[-1].strip(), body=theory_intro)]
        theory_blocks.extend(subsections)
    for index, block in enumerate(theory_blocks):
        if contains_executable_fence(block.body):
            step_index, score, method = 0, 100.0, "execution_order"
        else:
            step_index, score, method = best_step(
                block.title,
                block.body,
                steps,
                article_number=article_number,
                block_index=index,
                block_count=len(theory_blocks),
            )
        anchor = by_role["theory"][0].anchor if index == 0 else block.anchor
        steps[step_index].prefix.append(render_subsection(block, level=3, anchor=anchor))
        report["theory_matches"].append(
            {
                "section": block.title,
                "step": steps[step_index].title,
                "score": round(score, 4),
                "method": method,
            }
        )

    audit_blocks: list[Subsection] = []
    for section in by_role.get("audit", []):
        audit_intro, subsections = split_h3(section.body)
        if audit_intro:
            if subsections:
                subsections[0].body = "\n\n".join((audit_intro, subsections[0].body))
            else:
                subsections = [Subsection(title=section.title, body=audit_intro)]
        audit_blocks.extend(subsections)
    for index, block in enumerate(audit_blocks):
        if contains_executable_fence(block.body):
            step_index, score, method = len(steps) - 1, 100.0, "execution_order"
        else:
            step_index, score, method = best_step(
                block.title,
                block.body,
                steps,
                article_number=article_number,
                block_index=index,
                block_count=len(audit_blocks),
            )
        anchor = by_role["audit"][0].anchor if index == 0 else block.anchor
        steps[step_index].suffix.append(render_subsection(block, level=3, anchor=anchor))
        report["audit_matches"].append(
            {
                "section": block.title,
                "step": steps[step_index].title,
                "score": round(score, 4),
                "method": method,
            }
        )

    pitfall_blocks: list[str] = []
    pitfall_remainder: list[str] = []
    for section in by_role.get("pitfalls", []):
        callouts, remainder = parse_callouts(section.body)
        pitfall_blocks.extend(callouts)
        if remainder:
            pitfall_remainder.append(remainder)
    for index, callout in enumerate(pitfall_blocks):
        title = callout_title(callout)
        if contains_executable_fence(callout):
            step_index, score, method = len(steps) - 1, 100.0, "execution_order"
        else:
            step_index, score, method = best_step(
                title,
                callout,
                steps,
                article_number=article_number,
                block_index=index,
                block_count=len(pitfall_blocks),
            )
        steps[step_index].pitfalls.append(callout)
        report["pitfall_matches"].append(
            {
                "section": callout_title(callout),
                "step": steps[step_index].title,
                "score": round(score, 4),
                "method": method,
            }
        )
    if pitfall_remainder:
        steps[-1].pitfalls.extend(pitfall_remainder)

    for index, step in enumerate(steps):
        output.append(
            render_subsection(
                step,
                level=2,
                anchor=code_section.anchor if index == 0 else None,
            )
        )

    # Article 75 teaches figure construction, so its visual-design sections
    # remain in the main narrative.  Elsewhere they become Appendix C.
    publication_appendix: list[str] = []
    for section in by_role.get("publication", []):
        publication_intro, publication_subsections = split_h3(section.body)
        if article_number == 75:
            if publication_intro:
                output.append(qmd_heading(2, "让图形承担明确的证据任务", section.anchor) + "\n\n" + publication_intro)
            for index, subsection in enumerate(publication_subsections):
                output.append(
                    render_subsection(
                        subsection,
                        level=2,
                        anchor=section.anchor if index == 0 and not publication_intro else None,
                    )
                )
        else:
            parts = [qmd_heading(2, "附录 C：图形与导出规范", section.anchor)]
            if publication_intro:
                parts.append(publication_intro)
            for subsection in publication_subsections:
                parts.append(render_subsection(subsection, level=3))
            publication_appendix.append("\n\n".join(parts))

    summary = by_role.get("summary", [None])[0]
    output.append(make_takeaways(metadata, summary))

    for section in by_role.get("methods", []):
        output.append(
            qmd_heading(2, "附录 A：Methods / Results 模板", section.anchor)
            + "\n\n"
            + section.body.strip()
        )
    for section in by_role.get("own", []):
        output.append(
            qmd_heading(2, "附录 B：换成你自己的数据", section.anchor)
            + "\n\n"
            + section.body.strip()
        )
    output.extend(publication_appendix)

    # Preserve any non-template sections before references.
    for section in by_role.get("other", []):
        output.append(qmd_heading(2, section.title, section.anchor) + "\n\n" + section.body.strip())

    for section in by_role.get("references", []):
        output.append(qmd_heading(2, "参考文献", section.anchor) + "\n\n" + section.body.strip())

    transformed = front.rstrip() + "\n\n" + "\n\n".join(
        part.strip() for part in output if part and part.strip()
    ) + "\n"
    transformed_front, transformed_body = split_front_matter(transformed)
    transformed = transformed_front.rstrip() + "\n\n" + rewrite_headings(transformed_body).strip() + "\n"
    report["status"] = "changed" if transformed != text else "unchanged"
    report["old_h2"] = len(H2_RE.findall(body))
    report["new_h2"] = len(H2_RE.findall(split_front_matter(transformed)[1]))
    return transformed, report


def main() -> None:
    args = parse_args()
    root = args.project_root.resolve()
    manifest = yaml.safe_load((root / args.manifest).read_text(encoding="utf-8"))
    requested = set(args.article_number or [])
    reports: list[dict] = []

    if args.write and args.backup_dir:
        args.backup_dir.mkdir(parents=True, exist_ok=True)

    for chapter in manifest["series"]["chapters"]:
        number = int(chapter["number"])
        if requested and number not in requested:
            continue
        path = root / chapter["file"]
        original = path.read_text(encoding="utf-8")
        if number == 71:
            reports.append({"article": 71, "status": "reference", "reason": "linear baseline"})
            continue
        transformed, report = transform(original, number)
        report["file"] = str(path.relative_to(root))
        reports.append(report)
        if args.write and transformed != original:
            if args.backup_dir:
                backup = args.backup_dir / chapter["file"]
                backup.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(path, backup)
            path.write_text(transformed, encoding="utf-8")

    payload = {
        "articles": len(reports),
        "write": args.write,
        "changed": sum(item.get("status") == "changed" for item in reports),
        "reports": reports,
    }
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    print(json.dumps({key: payload[key] for key in ("articles", "write", "changed")}, ensure_ascii=False))


if __name__ == "__main__":
    main()
