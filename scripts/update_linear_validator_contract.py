#!/usr/bin/env python3
"""Replace legacy nine-slot heading assertions in article validators.

Scientific/result checks stay untouched.  Only tuple/list literals that are
clearly the old structural contract (at least five legacy slot names) are
rewritten to assert the new linear narrative and appendices.
"""

from __future__ import annotations

import ast
import json
import re
from pathlib import Path


LEGACY = (
    "对应论文里的哪张图",
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

CORE = (
    "先确定这一点",
    "#sec-theory",
    "#sec-code",
    "#sec-audit",
    ".callout-caution",
    "#sec-methods",
    "#sec-own-data",
    "参考文献",
)


def line_offsets(text: str) -> list[int]:
    offsets = [0]
    for line in text.splitlines(keepends=True):
        offsets.append(offsets[-1] + len(line))
    return offsets


def literal_strings(node: ast.AST) -> list[str]:
    if not isinstance(node, (ast.Tuple, ast.List)):
        return []
    values: list[str] = []
    for item in node.elts:
        if isinstance(item, ast.Constant) and isinstance(item.value, str):
            values.append(item.value)
    return values


def replacement(node: ast.Tuple | ast.List, indent: str, article: int | None) -> str:
    values = list(CORE)
    values.append("#sec-publication" if article == 75 else "图中应该保留哪些信息")
    opening, closing = ("(", ")") if isinstance(node, ast.Tuple) else ("[", "]")
    body = "\n".join(f'{indent}    "{value}",' for value in values)
    return f"{opening}\n{body}\n{indent}{closing}"


def update(path: Path) -> int:
    if path.name == "validate_article71_sem.py":
        return 0
    text = path.read_text(encoding="utf-8")
    tree = ast.parse(text)
    offsets = line_offsets(text)
    article_match = re.search(r"validate_article(\d+)", path.name)
    article = int(article_match.group(1)) if article_match else None
    edits: list[tuple[int, int, str]] = []
    lines = text.splitlines()
    for node in ast.walk(tree):
        values = literal_strings(node)
        if sum(any(term in value for term in LEGACY) for value in values) < 5:
            continue
        assert isinstance(node, (ast.Tuple, ast.List))
        start = offsets[node.lineno - 1] + node.col_offset
        end = offsets[node.end_lineno - 1] + node.end_col_offset
        indent = re.match(r"^\s*", lines[node.lineno - 1]).group(0)
        edits.append((start, end, replacement(node, indent, article)))

    for start, end, value in sorted(edits, reverse=True):
        text = text[:start] + value + text[end:]
    if edits:
        path.write_text(text, encoding="utf-8")
    return len(edits)


def main() -> None:
    root = Path(__file__).resolve().parents[1]
    if any(re.search(r'^reader-mode:', p.read_text(), re.M)
           for p in [root/'index.qmd']+list((root/'chapters').glob('*.qmd'))):
        raise SystemExit('Historical migration already applied; retain current semantic reader gates.')
    reports = []
    for path in sorted((root / "scripts").glob("validate_article*.py")):
        count = update(path)
        if count:
            reports.append({"file": path.relative_to(root).as_posix(), "literals": count})
    print(json.dumps({"files": len(reports), "literals": sum(item["literals"] for item in reports), "reports": reports}, ensure_ascii=False))


if __name__ == "__main__":
    main()
