"""Render Pandoc math as local images: WeChat does not execute MathJax.

Only trusted, reviewed article formulas are accepted. TeX shell escape is
disabled, dangerous file/definition commands are rejected, and failures stop
the bundle instead of silently stripping symbols. Generated files are local
build artifacts, not article source.
"""
from __future__ import annotations

import hashlib
from pathlib import Path
import re
import subprocess
import tempfile

from lxml import etree
from PIL import Image

FORBIDDEN = re.compile(
    r"\\(?:input|include|openin|openout|read|write|immediate|special|usepackage|"
    r"documentclass|catcode|csname|def|edef|gdef|xdef|newcommand|renewcommand|"
    r"let|loop|repeat|href|url|includegraphics)\b", re.I
)


def render_formula(tex: str, display: bool, cache: Path) -> Path:
    if not tex.strip() or len(tex) > 10000 or FORBIDDEN.search(tex):
        raise ValueError("Empty, oversized or unsafe math expression")
    key = hashlib.sha256((str(display) + tex).encode()).hexdigest()
    destination = cache / (key + ".png")
    if destination.exists():
        return destination
    cache.mkdir(parents=True, exist_ok=True)
    document = (
        "\\documentclass[border=2pt]{standalone}\n"
        "\\usepackage{amsmath,amssymb}\n"
        "\\begin{document}\n$"
        + ("\\displaystyle " if display else "") + tex
        + "$\n\\end{document}\n"
    )
    with tempfile.TemporaryDirectory(prefix="wechat-math-") as tmp:
        work = Path(tmp)
        (work / "formula.tex").write_text(document, encoding="utf-8")
        for command in (
            ["pdflatex", "-no-shell-escape", "-halt-on-error", "-interaction=batchmode", "formula.tex"],
            ["pdftoppm", "-singlefile", "-png", "-r", "240", "formula.pdf", "formula"],
        ):
            result = subprocess.run(command, cwd=work, capture_output=True, timeout=30)
            if result.returncode:
                raise RuntimeError(f"Math rendering failed for expression {tex!r}")
        with Image.open(work / "formula.png") as raster:
            raster.convert("RGB").save(destination)
    return destination


def replace_math(main: etree._Element, cache: Path) -> int:
    count = 0
    for node in list(main.xpath('.//*[contains(concat(" ",normalize-space(@class)," ")," math ")]')):
        text = ''.join(node.itertext()).strip()
        display = "display" in (node.get("class") or "").split()
        opening, closing = (r"\[", r"\]") if display else (r"\(", r"\)")
        if not (text.startswith(opening) and text.endswith(closing)):
            # Citeproc can already render a title symbol as native HTML.
            descendants = list(node.iterdescendants())
            if (descendants and text and not re.search(r"[\\$^_{}]", text)
                    and all(child.tag in {"em", "span", "sub", "sup"} for child in descendants)):
                node.attrib.pop("class", None)
                continue
            raise ValueError("Unexpected rendered math delimiters")
        tex = text[len(opening):-len(closing)]
        path = render_formula(tex, display, cache)
        with Image.open(path) as raster:
            width, height = raster.size
        node.tag = "img"
        for child in list(node):
            node.remove(child)
        node.text = None
        node.attrib.clear()
        node.set("src", str(path))
        node.set("alt", "Equation: " + tex)
        node.set("data-math-style", (
            f"display:{'block' if display else 'inline-block'};"
            f"width:{round(width / 2)}px;max-width:100%;height:auto;"
            f"margin:{'16px auto' if display else '0 2px'};"
            "vertical-align:middle;border-radius:0;"
        ))
        count += 1
    return count
