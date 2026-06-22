#!/usr/bin/env -S uv run --quiet --script
# /// script
# requires-python = ">=3.10"
# dependencies = [
#     "reportlab>=4.0",
#     "markdown>=3.5",
#     "pillow>=10.0",
#     "cairosvg>=2.7",
# ]
# ///
"""Render a markdown file to a stakeholder-styled PDF.

Usage:
    md2pdf.py SRC.md DST.pdf [--title "..."] [--date YYYY-MM-DD]

Layout: NOAA GML stakeholder-report look — Times serif body, navy
(#003366) headings with bottom rules, light-blue (#ccd6e0 / #f0f4f8)
table borders + header tint, Menlo-style code spans on a #f4f4f4 wash.
Header banner is the GML Earth-mosaic + NOAA seagull with a light
cyan band between, drawn on every page.

The CSS this is modeled on lives at miller-ff/methodology.html — keep
the two in rough visual sync.

Markdown features supported:
- ATX headings (# H1, ## H2, ### H3)
- Bold (**text**), italic (*text*), strike (~~text~~), code (`text`)
- Unicode subscripts/superscripts (auto-wrapped in <sub>/<super>)
- Pipe tables with header row + ---- separator
- Bulleted lists (- item) and numbered lists (1. item) with hanging
  indented continuation lines
- Horizontal rules (---)
- Images: ![alt](path) on their own line, scaled to the content width
  (path resolved relative to the source .md; optional trailing text -> caption)
- A leading meta-block of "**Label:** value" lines after the H1 gets
  promoted to a compact title-block

The first H1 in the document becomes the PDF title (and footer text).
The first "**Date:**" meta line becomes the footer date — falls back
to today's date if none present.

Author resolution order (used for PDF metadata + the right side of
the footer):
  1. --author CLI flag
  2. `<!-- author: ... -->` HTML comment anywhere in the doc
     (the comment itself is stripped from the rendered output)
  3. `**Author:**`, `**By:**`, or `**From:**` meta line
  4. Falls back to no author (footer right shows just the date)

Dependencies (auto-installed via uv on first run):
  reportlab, markdown, pillow, cairosvg

To rebuild a doc:
    uv run docs/bin/md2pdf.py docs/foo.md docs/foo.pdf

Or, after `chmod +x docs/bin/md2pdf.py`:
    docs/bin/md2pdf.py docs/foo.md docs/foo.pdf
"""
from __future__ import annotations

import argparse
import datetime as dt
import re
from pathlib import Path

from reportlab.lib import colors
from reportlab.lib.enums import TA_LEFT
from reportlab.lib.pagesizes import LETTER
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.units import inch
from reportlab.platypus import (
    HRFlowable,
    Image,
    KeepTogether,
    Paragraph,
    SimpleDocTemplate,
    Spacer,
    Table,
    TableStyle,
)

# ---------------------------------------------------------------------------
# Paths: assets live next to the script, so the renderer works from any cwd.
# ---------------------------------------------------------------------------

SCRIPT_DIR = Path(__file__).resolve().parent
ASSETS = SCRIPT_DIR.parent / "assets"
GML_LOGO = ASSETS / "GML-logo.jpg"
NOAA_SVG = ASSETS / "noaa-logo.svg"
# Cached high-res PNG render of the NOAA SVG. Lives next to the SVG so
# repeat invocations skip the cairosvg pass.
NOAA_PNG_HIRES = ASSETS / ".noaa-logo-hires.png"


def ensure_noaa_hires() -> Path:
    """Render the NOAA SVG to a 600px PNG once and cache it."""
    if NOAA_PNG_HIRES.exists() and NOAA_PNG_HIRES.stat().st_size > 0:
        return NOAA_PNG_HIRES
    import cairosvg  # type: ignore[import-untyped] -- third-party, no type stubs published

    cairosvg.svg2png(
        url=str(NOAA_SVG),
        write_to=str(NOAA_PNG_HIRES),
        output_width=600,
        output_height=600,
    )
    return NOAA_PNG_HIRES


# ---------------------------------------------------------------------------
# Palette — lifted from miller-ff/methodology.html.
# ---------------------------------------------------------------------------

NAVY = colors.HexColor("#003366")
BODY_TEXT = colors.HexColor("#222222")
BORDER_LIGHT = colors.HexColor("#ccd6e0")
HEADER_BG = colors.HexColor("#f0f4f8")
CODE_BG = colors.HexColor("#f4f4f4")
META_GREY = colors.HexColor("#444444")
BAND_COLOR = colors.HexColor("#c8e3ef")


# ---------------------------------------------------------------------------
# Inline markdown -> reportlab Paragraph mini-HTML.
# ---------------------------------------------------------------------------

_RE_CODE = re.compile(r"`([^`]+)`")
_RE_BOLD = re.compile(r"\*\*([^*]+?)\*\*")
_RE_STRIKE = re.compile(r"~~([^~]+?)~~")
# Italic: single asterisk not adjacent to another asterisk. The look-arounds
# prevent it from chewing on bold markers that haven't been replaced yet.
_RE_ITALIC = re.compile(r"(?<!\*)\*([^*\s][^*]*?)\*(?!\*)")

# Unicode subscripts / superscripts -> ASCII char (then wrap in <sub>/<super>).
# Helvetica/Times don't carry the Unicode-block glyphs and would otherwise
# render them as solid black tofu boxes.
_SUB_MAP = dict(zip(
    "₀₁₂₃₄₅₆₇₈₉₊₋₌₍₎ₐₑₕᵢⱼₖₗₘₙₒₚᵣₛₜᵤᵥₓ",
    list("0123456789+-=()") + list("aehijklmnoprstuvx"),
))
_SUP_MAP = dict(zip(
    "⁰¹²³⁴⁵⁶⁷⁸⁹⁺⁻⁼⁽⁾ⁿⁱ",
    list("0123456789+-=()") + list("ni"),
))


def _translate_unicode_sub_sup(text: str) -> str:
    out: list[str] = []
    i = 0
    while i < len(text):
        c = text[i]
        if c in _SUB_MAP:
            j = i
            run = []
            while j < len(text) and text[j] in _SUB_MAP:
                run.append(_SUB_MAP[text[j]])
                j += 1
            out.append("<sub>" + "".join(run) + "</sub>")
            i = j
        elif c in _SUP_MAP:
            j = i
            run = []
            while j < len(text) and text[j] in _SUP_MAP:
                run.append(_SUP_MAP[text[j]])
                j += 1
            out.append("<super>" + "".join(run) + "</super>")
            i = j
        else:
            out.append(c)
            i += 1
    return "".join(out)


def _escape_for_para(text: str) -> str:
    """Escape <, >, & for reportlab's intra-Paragraph parser, then substitute
    our markdown -> mini-HTML translations on the escaped text.
    """
    text = text.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
    text = _translate_unicode_sub_sup(text)
    # Stash code spans so their contents don't get re-parsed for bold/italic.
    code_chunks: list[str] = []

    def _stash_code(m: re.Match[str]) -> str:
        code_chunks.append(m.group(1))
        return f"\x00CODE{len(code_chunks) - 1}\x00"

    text = _RE_CODE.sub(_stash_code, text)
    text = _RE_BOLD.sub(r"<b>\1</b>", text)
    text = _RE_STRIKE.sub(r"<strike>\1</strike>", text)
    text = _RE_ITALIC.sub(r"<i>\1</i>", text)
    for i, chunk in enumerate(code_chunks):
        text = text.replace(
            f"\x00CODE{i}\x00",
            (
                '<font face="Courier" size="9" color="#003366" '
                'backColor="#f4f4f4">'
                + chunk
                + "</font>"
            ),
        )
    return text


# ---------------------------------------------------------------------------
# Block parsing helpers.
# ---------------------------------------------------------------------------


def _is_table_row(line: str) -> bool:
    return line.lstrip().startswith("|") and line.rstrip().endswith("|")


def _split_table_row(line: str) -> list[str]:
    inner = line.strip()
    inner = inner.removeprefix("|")
    inner = inner.removesuffix("|")
    return [cell.strip() for cell in inner.split("|")]


def _is_table_separator(cells: list[str]) -> bool:
    return all(re.fullmatch(r":?-+:?", c.strip()) for c in cells if c.strip())


# ---------------------------------------------------------------------------
# Paragraph / table / heading styles.
# ---------------------------------------------------------------------------


def make_styles() -> dict[str, ParagraphStyle]:
    base = getSampleStyleSheet()["Normal"]
    return {
        # H1 ≈ CSS 1.5em on 10pt body -> ~15pt; navy with thicker rule.
        "title": ParagraphStyle(
            "Title",
            parent=base,
            fontName="Times-Bold",
            fontSize=17,
            leading=21,
            spaceAfter=2,
            textColor=NAVY,
            keepWithNext=1,   # never orphan a heading at a page bottom
        ),
        # H2 ≈ CSS 1.2em -> ~12pt navy, thin rule.
        "h2": ParagraphStyle(
            "H2",
            parent=base,
            fontName="Times-Bold",
            fontSize=13,
            leading=16,
            spaceBefore=16,
            spaceAfter=2,
            textColor=NAVY,
            keepWithNext=1,
        ),
        # H3 ≈ CSS 1.05em -> ~11pt navy, no rule (matches methodology.html).
        "h3": ParagraphStyle(
            "H3",
            parent=base,
            fontName="Times-Bold",
            fontSize=11,
            leading=14,
            spaceBefore=10,
            spaceAfter=2,
            textColor=NAVY,
            keepWithNext=1,
        ),
        # Body — Georgia isn't built-in; Times-Roman is the nearest available
        # transitional serif. CSS line-height 1.6 -> 16pt leading at 10.5pt.
        "body": ParagraphStyle(
            "Body",
            parent=base,
            fontName="Times-Roman",
            fontSize=10.5,
            leading=16,
            spaceAfter=8,
            alignment=TA_LEFT,
            textColor=BODY_TEXT,
        ),
        "bullet": ParagraphStyle(
            "Bullet",
            parent=base,
            fontName="Times-Roman",
            fontSize=10.5,
            leading=15,
            spaceAfter=4,
            leftIndent=20,
            bulletIndent=6,
            alignment=TA_LEFT,
            textColor=BODY_TEXT,
        ),
        "meta": ParagraphStyle(
            "Meta",
            parent=base,
            fontName="Times-Roman",
            fontSize=10,
            leading=13,
            spaceAfter=1,
            textColor=META_GREY,
        ),
        "table_cell": ParagraphStyle(
            "TableCell",
            parent=base,
            fontName="Times-Roman",
            fontSize=9.5,
            leading=12,
            alignment=TA_LEFT,
            textColor=BODY_TEXT,
        ),
        "table_head": ParagraphStyle(
            "TableHead",
            parent=base,
            fontName="Times-Bold",
            fontSize=9.5,
            leading=12,
            alignment=TA_LEFT,
            textColor=NAVY,
        ),
    }


# ---------------------------------------------------------------------------
# Renderer — markdown text in, list of flowables out.
# ---------------------------------------------------------------------------


_RE_IMAGE = re.compile(
    r'^\s*!\[(?P<alt>[^\]]*)\]\((?P<path>[^)\s]+)(?:\s+"[^"]*")?\)(?P<rest>.*)$')


def render_markdown(md_text: str, styles: dict[str, ParagraphStyle],
                    base_dir: Path | None = None):
    flowables: list = []
    heading_ids: set[int] = set()   # ids of heading flowables, to bundle with a following figure
    lines = md_text.splitlines()
    i = 0
    n = len(lines)

    def flush_paragraph(buf: list[str], style_key: str = "body") -> None:
        if not buf:
            return
        text = " ".join(s.strip() for s in buf).strip()
        if not text:
            return
        flowables.append(Paragraph(_escape_for_para(text), styles[style_key]))

    def image_flowable(src_path: str, alt: str, rest: str = "") -> list:
        """![alt](path) -> the flowable parts of a content-width, centered figure
        (aspect preserved, height capped). Returns a flat list (NOT a KeepTogether)
        so the caller can prepend a preceding heading and wrap once — nesting
        KeepTogethers misrenders. Path is relative to the source .md; trailing text
        after the image (e.g. '*(left panel)*') becomes a small caption."""
        usable = 6.5 * inch
        p = Path(src_path)
        if not p.is_absolute():
            p = (base_dir or Path.cwd()) / p
        if not p.exists():
            return [Paragraph(
                f"<i>[missing image: {_escape_for_para(src_path)}]</i>", styles["body"])]
        try:
            from PIL import Image as _PILImage
            with _PILImage.open(p) as _im:
                iw, ih = _im.size
            aspect = (ih / iw) if iw else 0.6
        except Exception:
            aspect = 0.6
        w, h = usable, usable * aspect
        if h > 7.0 * inch:                       # keep tall figures on one page
            h, w = 7.0 * inch, 7.0 * inch / aspect
        img = Image(str(p), width=w, height=h)
        img.hAlign = "CENTER"
        parts = [Spacer(1, 4), img]
        if rest.strip():
            parts.append(Paragraph(_escape_for_para(rest.strip()), styles["body"]))
        parts.append(Spacer(1, 6))
        return parts

    while i < n:
        line = lines[i]
        stripped = line.strip()

        # Blank line
        if not stripped:
            i += 1
            continue

        # Horizontal rule
        if re.fullmatch(r"-{3,}", stripped):
            flowables.append(Spacer(1, 4))
            flowables.append(
                HRFlowable(
                    width="100%",
                    thickness=0.5,
                    color=BORDER_LIGHT,
                    spaceBefore=2,
                    spaceAfter=6,
                )
            )
            i += 1
            continue

        # Image: ![alt](path) on its own line -> scaled, centered figure.
        # If it immediately follows a heading, bundle the two into one
        # KeepTogether so the heading never orphans at a page bottom
        # (keepWithNext alone misses tall KeepTogether figures).
        m_img = _RE_IMAGE.match(line)
        if m_img:
            parts = image_flowable(
                m_img.group("path"), m_img.group("alt"), m_img.group("rest"))
            # Prepend any run of immediately-preceding heading flowables so a
            # section divider (## ) travels onto the same page as its first
            # subsection (### ) and the figure — one FLAT KeepTogether (nesting
            # KeepTogethers misrenders, splitting the heading off the figure).
            while flowables and id(flowables[-1]) in heading_ids:
                parts.insert(0, flowables.pop())
            flowables.append(KeepTogether(parts))
            i += 1
            continue

        # H1 — 2px navy underline (CSS: border-bottom: 2px solid #003366)
        if stripped.startswith("# "):
            hf = KeepTogether(
                [
                    Paragraph(_escape_for_para(stripped[2:]), styles["title"]),
                    HRFlowable(width="100%", thickness=1.6, color=NAVY,
                               spaceBefore=2, spaceAfter=8),
                ]
            )
            flowables.append(hf)
            heading_ids.add(id(hf))
            i += 1
            continue
        # H2 — 1px light-blue underline (CSS: border-bottom: 1px solid #ccd6e0)
        if stripped.startswith("## "):
            hf = KeepTogether(
                [
                    Paragraph(_escape_for_para(stripped[3:]), styles["h2"]),
                    HRFlowable(width="100%", thickness=0.5, color=BORDER_LIGHT,
                               spaceBefore=1, spaceAfter=4),
                ]
            )
            flowables.append(hf)
            heading_ids.add(id(hf))
            i += 1
            continue
        # H3 — plain (no rule).
        if stripped.startswith("### "):
            hf = Paragraph(_escape_for_para(stripped[4:]), styles["h3"])
            flowables.append(hf)
            heading_ids.add(id(hf))
            i += 1
            continue

        # Table — collect contiguous pipe-table lines.
        if _is_table_row(line):
            rows: list[list[str]] = []
            while i < n and _is_table_row(lines[i]):
                rows.append(_split_table_row(lines[i]))
                i += 1
            header: list[str] | None = None
            data_rows: list[list[str]] = []
            if len(rows) >= 2 and _is_table_separator(rows[1]):
                header = rows[0]
                data_rows = rows[2:]
            else:
                data_rows = rows

            table_data: list[list] = []
            if header is not None:
                table_data.append(
                    [Paragraph(_escape_for_para(c), styles["table_head"]) for c in header]
                )
            for r in data_rows:
                table_data.append(
                    [Paragraph(_escape_for_para(c), styles["table_cell"]) for c in r]
                )

            if not table_data:
                continue
            ncols = max(len(r) for r in table_data)
            for r in table_data:
                while len(r) < ncols:
                    r.append(Paragraph("", styles["table_cell"]))

            # Content-aware column widths. Take the max plain-text length per
            # column (stripping markdown markers), allocate proportionally
            # with min/max bounds, renormalize to usable width.
            usable = 6.5 * inch
            raw_rows: list[list[str]] = []
            if header is not None:
                raw_rows.append(header)
            raw_rows.extend(data_rows)
            col_maxes = [0] * ncols
            for r in raw_rows:
                for ci, cell in enumerate(r[:ncols]):
                    plain = re.sub(r"[*`~]", "", cell)
                    col_maxes[ci] = max(col_maxes[ci], len(plain))
            import math
            # sqrt-damping so a single huge cell doesn't starve sibling columns.
            weights = [max(1.0, math.sqrt(m)) for m in col_maxes]
            min_frac = 0.07
            widths = [
                max(min_frac * usable, w / sum(weights) * usable)
                for w in weights
            ]
            widths = [w * (usable / sum(widths)) for w in widths]

            tbl = Table(
                table_data, colWidths=widths, repeatRows=1 if header else 0
            )
            ts = TableStyle(
                [
                    ("VALIGN", (0, 0), (-1, -1), "TOP"),
                    ("LEFTPADDING", (0, 0), (-1, -1), 6),
                    ("RIGHTPADDING", (0, 0), (-1, -1), 6),
                    ("TOPPADDING", (0, 0), (-1, -1), 5),
                    ("BOTTOMPADDING", (0, 0), (-1, -1), 5),
                    ("GRID", (0, 0), (-1, -1), 0.5, BORDER_LIGHT),
                ]
            )
            if header is not None:
                ts.add("BACKGROUND", (0, 0), (-1, 0), HEADER_BG)
            tbl.setStyle(ts)
            flowables.append(Spacer(1, 4))
            flowables.append(tbl)
            flowables.append(Spacer(1, 6))
            continue

        # Bulleted list — contiguous "- item" lines + indented continuations.
        if re.match(r"-\s+", stripped):
            while i < n:
                ln = lines[i]
                lst = ln.strip()
                if re.match(r"-\s+", lst):
                    body = lst[2:] if lst.startswith("- ") else lst[1:].lstrip()
                    j = i + 1
                    while j < n:
                        nxt = lines[j]
                        if not nxt.strip():
                            break
                        if re.match(r"-\s+", nxt.strip()):
                            break
                        if nxt.startswith("  ") or nxt.startswith("\t"):
                            body += " " + nxt.strip()
                            j += 1
                            continue
                        break
                    flowables.append(
                        Paragraph(
                            _escape_for_para(body),
                            styles["bullet"],
                            bulletText="•",
                        )
                    )
                    i = j
                elif not lst:
                    i += 1
                    break
                else:
                    break
            continue

        # Numbered list — same shape, "N. item".
        if re.match(r"\d+\.\s+", stripped):
            while i < n:
                ln = lines[i]
                lst = ln.strip()
                m = re.match(r"(\d+)\.\s+(.*)", lst)
                if m:
                    num, body = m.group(1), m.group(2)
                    j = i + 1
                    while j < n:
                        nxt = lines[j]
                        if not nxt.strip():
                            break
                        if re.match(r"\d+\.\s+", nxt.strip()):
                            break
                        if re.match(r"-\s+", nxt.strip()):
                            break
                        if nxt.startswith("  ") or nxt.startswith("\t"):
                            body += " " + nxt.strip()
                            j += 1
                            continue
                        break
                    flowables.append(
                        Paragraph(
                            _escape_for_para(body),
                            styles["bullet"],
                            bulletText=f"{num}.",
                        )
                    )
                    i = j
                elif not lst:
                    i += 1
                    break
                else:
                    break
            continue

        # Blockquote — contiguous "> ..." lines as a navy-left-barred callout.
        if stripped.startswith(">"):
            qbuf: list[str] = []
            while i < n and lines[i].strip().startswith(">"):
                q = lines[i].strip()[1:]
                qbuf.append(q[1:] if q.startswith(" ") else q)
                i += 1
            qpara = Paragraph(_escape_for_para(" ".join(qbuf).strip()), styles["body"])
            qt = Table([[qpara]], colWidths=[6.5 * inch])
            qt.setStyle(
                TableStyle(
                    [
                        ("LINEBEFORE", (0, 0), (0, 0), 3, NAVY),
                        ("BACKGROUND", (0, 0), (-1, -1), HEADER_BG),
                        ("LEFTPADDING", (0, 0), (-1, -1), 10),
                        ("RIGHTPADDING", (0, 0), (-1, -1), 8),
                        ("TOPPADDING", (0, 0), (-1, -1), 6),
                        ("BOTTOMPADDING", (0, 0), (-1, -1), 6),
                    ]
                )
            )
            flowables.append(Spacer(1, 4))
            flowables.append(qt)
            flowables.append(Spacer(1, 6))
            continue

        # Default: paragraph — greedy collection of non-blank, non-special lines.
        buf: list[str] = []
        while i < n:
            ln = lines[i]
            s = ln.strip()
            if not s:
                break
            if (
                re.match(r"#{1,6}\s", s)
                or s.startswith(">")
                or re.fullmatch(r"-{3,}", s)
                or _is_table_row(ln)
                or re.match(r"-\s+", s)
                or re.match(r"\d+\.\s+", s)
            ):
                break
            buf.append(ln)
            i += 1
        flush_paragraph(buf)

    return flowables


# ---------------------------------------------------------------------------
# Meta-block detection — promote the leading "**Label:** value" paragraph
# into compact one-per-line meta entries.
# ---------------------------------------------------------------------------

# Looks like a single body paragraph that collapsed several "**Label:** ..."
# lines together. We split it back out at each "**...**:" marker.
_META_LABEL_RE = re.compile(r"<b>([A-Z][A-Za-z ]{0,20}):</b>")


def promote_meta_block(flowables: list, styles: dict[str, ParagraphStyle]) -> list:
    """Re-style the leading paragraph if it's a colon-labeled meta block."""
    out: list = []
    consumed = False
    for f in flowables:
        if not consumed and isinstance(f, Paragraph) and f.style.name == "Body":
            xml = f.text  # reportlab-internal post-escape representation
            # Quick heuristic: leading paragraph that starts with a bold-label
            # like "<b>To:</b>" and has at least 2 such labels.
            if xml.lstrip().startswith("<b>") and len(_META_LABEL_RE.findall(xml)) >= 2:
                # Split at each label, keeping the label.
                labels = list(_META_LABEL_RE.finditer(xml))
                pieces: list[tuple[str, str]] = []  # (label, value)
                for idx, m in enumerate(labels):
                    start = m.start()
                    end = labels[idx + 1].start() if idx + 1 < len(labels) else len(xml)
                    seg = xml[start:end].strip()
                    pieces.append(("", seg))

                for _, seg in pieces:
                    # Re-color the bold label navy.
                    seg = re.sub(
                        r"<b>([^<]+)</b>",
                        r'<font color="#003366"><b>\1</b></font>',
                        seg,
                        count=1,
                    )
                    out.append(Paragraph(seg, styles["meta"]))
                out.append(Spacer(1, 6))
                consumed = True
                continue
        out.append(f)
    return out


# ---------------------------------------------------------------------------
# Title / date extraction.
# ---------------------------------------------------------------------------


def extract_title(md_text: str) -> str | None:
    for line in md_text.splitlines():
        if line.startswith("# "):
            return line[2:].strip()
    return None


def extract_date(md_text: str) -> str | None:
    # Capture just the date token — stop at a "·" separator, a "**" marker, or
    # end of line, so a packed "**Status:** … · **Date:** 2026-06-21 · **Purpose:** …"
    # line yields "2026-06-21", not the trailing fields.
    m = re.search(r"\*\*Date:\*\*\s*([^·*\n]+)", md_text)
    if m:
        return m.group(1).strip().rstrip("*").strip()
    return None


_RE_AUTHOR_COMMENT = re.compile(
    r"<!--\s*author:\s*(.+?)\s*-->", re.IGNORECASE | re.DOTALL
)
_RE_AUTHOR_META = re.compile(
    r"\*\*(?:Author|By|From):\*\*\s*(.+)", re.IGNORECASE
)


def extract_author(md_text: str) -> str | None:
    m = _RE_AUTHOR_COMMENT.search(md_text)
    if m:
        return m.group(1).strip()
    m = _RE_AUTHOR_META.search(md_text)
    if m:
        return m.group(1).strip().rstrip("*").strip()
    return None


_RE_FOOTER_TITLE = re.compile(
    r"<!--\s*footer-title:\s*(.+?)\s*-->", re.IGNORECASE | re.DOTALL
)


def extract_footer_title(md_text: str) -> str | None:
    """A short running-title for the footer (the H1 may be too long to fit beside
    the centred page number). Set via `<!-- footer-title: ... -->`."""
    m = _RE_FOOTER_TITLE.search(md_text)
    return m.group(1).strip() if m else None


def strip_html_comments(md_text: str) -> str:
    """Remove `<!-- ... -->` blocks so they don't render as plain text."""
    return re.sub(r"<!--.*?-->", "", md_text, flags=re.DOTALL)


# ---------------------------------------------------------------------------
# Main build.
# ---------------------------------------------------------------------------

# Unicode sub/superscript digits -> ASCII, for canvas-drawn text (header/footer).
# The Times-Roman base-14 font lacks those glyphs, so e.g. "CH₄" would render as
# "CH■" (a solid black box) in the footer; Paragraph flowables are unaffected.
_CANVAS_SUBSUP = {
    **{0x2080 + d: str(d) for d in range(10)},          # ₀–₉
    0x2070: "0", 0x00B9: "1", 0x00B2: "2", 0x00B3: "3",  # ⁰¹²³
    **{0x2074 + (d - 4): str(d) for d in range(4, 10)},  # ⁴–⁹
}


def _canvas_safe(s: str) -> str:
    """Strip Unicode sub/superscripts so canvas.drawString shows digits, not boxes."""
    return s.translate(_CANVAS_SUBSUP)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("src", type=Path, help="Source markdown file")
    ap.add_argument("dst", type=Path, help="Destination PDF path")
    ap.add_argument(
        "--title",
        help="Document title (defaults to first H1 in SRC)",
    )
    ap.add_argument(
        "--date",
        help="Date shown in the footer (defaults to first **Date:** meta line, "
        "then today)",
    )
    ap.add_argument(
        "--author",
        help="Author shown in the footer + PDF metadata (defaults to "
        "<!-- author: ... --> comment, then **Author:** / **By:** / **From:** "
        "meta line, then unset)",
    )
    ap.add_argument(
        "--footer-title",
        help="Short running-title for the footer (defaults to "
        "<!-- footer-title: ... --> comment, then the document title)",
    )
    args = ap.parse_args()

    md_text = args.src.read_text()
    # Extract metadata from the raw text BEFORE stripping comments — the
    # author/footer-title markers live in comments.
    title = args.title or extract_title(md_text) or args.src.stem
    date = args.date or extract_date(md_text) or dt.date.today().isoformat()
    author = args.author or extract_author(md_text)
    footer_title = args.footer_title or extract_footer_title(md_text) or title
    md_text = strip_html_comments(md_text)

    styles = make_styles()
    flowables = render_markdown(md_text, styles, base_dir=args.src.resolve().parent)
    flowables = promote_meta_block(flowables, styles)

    doc = SimpleDocTemplate(
        str(args.dst),
        pagesize=LETTER,
        leftMargin=1.0 * inch,
        rightMargin=1.0 * inch,
        topMargin=1.25 * inch,  # room for the GML / NOAA header band
        bottomMargin=0.85 * inch,
        title=title,
        author=author or "NOAA CarbonTracker (CT-unified)",
    )

    noaa_path = ensure_noaa_hires()

    def _decorate(canvas, doc):
        canvas.saveState()
        page_w, page_h = LETTER

        # --- Header band + logos --------------------------------------------
        # The band is a thin horizontal strip; logos sit on it, extending
        # above and below. Geometry tuned to the NOAA-GML slide-deck header.
        margin = 0.5 * inch
        band_left = margin + 1.25 * inch  # leaves room for the GML logo
        band_right = page_w - margin - 1.10 * inch  # leaves room for NOAA
        band_h = 0.32 * inch
        band_cy = page_h - 0.65 * inch  # vertical center of the band
        canvas.setFillColor(BAND_COLOR)
        canvas.setStrokeColor(BAND_COLOR)
        canvas.rect(
            band_left,
            band_cy - band_h / 2,
            band_right - band_left,
            band_h,
            fill=1,
            stroke=0,
        )

        logo_h = 0.85 * inch
        # GML logo: 750x556 aspect ratio.
        gml_w = logo_h * (750.0 / 556.0)
        canvas.drawImage(
            str(GML_LOGO),
            margin,
            band_cy - logo_h / 2,
            width=gml_w,
            height=logo_h,
            mask="auto",
            preserveAspectRatio=True,
        )
        # NOAA logo: square.
        noaa_h = logo_h
        canvas.drawImage(
            str(noaa_path),
            page_w - margin - noaa_h,
            band_cy - noaa_h / 2,
            width=noaa_h,
            height=noaa_h,
            mask="auto",
            preserveAspectRatio=True,
        )

        # --- Footer ---------------------------------------------------------
        canvas.setFont("Times-Roman", 8.5)
        canvas.setFillColor(NAVY)
        page_str = f"— {doc.page} —"
        canvas.drawCentredString(page_w / 2.0, 0.5 * inch, page_str)
        canvas.setFillColor(META_GREY)
        # Truncate the running title so a long one can't collide with the
        # centred page number (the overlap reads as a strikethrough).
        foot_title = _canvas_safe(footer_title)
        max_w = (page_w / 2.0
                 - canvas.stringWidth(page_str, "Times-Roman", 8.5) / 2.0
                 - 8) - 1.0 * inch
        if canvas.stringWidth(foot_title, "Times-Roman", 8.5) > max_w:
            while foot_title and canvas.stringWidth(
                    foot_title + "…", "Times-Roman", 8.5) > max_w:
                foot_title = foot_title[:-1]
            foot_title = foot_title.rstrip(" -—") + "…"
        canvas.drawString(1.0 * inch, 0.5 * inch, foot_title)
        # Author · date on the right, or just date if no author was found.
        right_text = f"{author} · {date}" if author else date
        canvas.drawRightString(page_w - 1.0 * inch, 0.5 * inch, _canvas_safe(right_text))
        canvas.restoreState()

    doc.build(flowables, onFirstPage=_decorate, onLaterPages=_decorate)
    size_kb = args.dst.stat().st_size / 1024
    print(
        f"wrote {args.dst} ({size_kb:.1f} KB, title={title!r}, "
        f"date={date!r}, author={author!r})"
    )


if __name__ == "__main__":
    main()
