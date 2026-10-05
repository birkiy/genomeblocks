"""Tiny SVG toolkit for the docs' design diagrams.

Every diagram is laid out by hand on explicit coordinates; this module only
keeps the markup consistent. Colours never appear here: shapes carry classes
(``bx green``, ``ln navy``, ...) that ``_sass/custom/custom.scss`` maps onto
the site's design tokens, so the diagrams follow the theme.
"""
from __future__ import annotations

import math
from html import escape


class Diagram:
    def __init__(self, name: str, w: int, h: int, label: str):
        self.name, self.w, self.h, self.label = name, w, h, label
        self.parts: list[str] = []

    # ── primitives ───────────────────────────────────────────────────────
    def raw(self, s: str):
        self.parts.append(s)

    def rect(self, x, y, w, h, cls="bx", r=8, title=None):
        t = f"<title>{escape(title)}</title>" if title else ""
        self.parts.append(f'<rect class="{cls}" x="{x:g}" y="{y:g}" width="{w:g}" height="{h:g}" rx="{r:g}">{t}</rect>')

    def text(self, x, y, s, cls="t", anchor="middle", weight=None):
        a = f' text-anchor="{anchor}"' if anchor != "start" else ""
        wt = f' font-weight="{weight}"' if weight else ""
        self.parts.append(f'<text class="{cls}" x="{x:g}" y="{y:g}"{a}{wt}>{escape(s)}</text>')

    def rich(self, x, y, spans, cls="s", anchor="start"):
        """One line of text with mixed classes: spans = [(text, cls|None), ...]."""
        a = f' text-anchor="{anchor}"' if anchor != "start" else ""
        inner = "".join(f'<tspan class="{c}">{escape(t)}</tspan>' if c else escape(t) for t, c in spans)
        self.parts.append(f'<text class="{cls}" x="{x:g}" y="{y:g}"{a}>{inner}</text>')

    def line(self, pts, cls="ln"):
        d = " ".join(f"{x:g},{y:g}" for x, y in pts)
        self.parts.append(f'<polyline class="{cls}" points="{d}"/>')

    def path(self, d, cls="ln"):
        self.parts.append(f'<path class="{cls}" d="{d}"/>')

    def circle(self, x, y, r, cls="dot"):
        self.parts.append(f'<circle class="{cls}" cx="{x:g}" cy="{y:g}" r="{r:g}"/>')

    def head(self, x, y, angle, cls="hd", size=7):
        """Arrowhead with its tip at (x, y), pointing along ``angle`` (radians)."""
        a1, a2 = angle + math.radians(152), angle - math.radians(152)
        p = [(x, y), (x + size * math.cos(a1), y + size * math.sin(a1)),
             (x + size * math.cos(a2), y + size * math.sin(a2))]
        d = " ".join(f"{px:.1f},{py:.1f}" for px, py in p)
        self.parts.append(f'<polygon class="{cls}" points="{d}"/>')

    # ── composites ───────────────────────────────────────────────────────
    def arrow(self, pts, kind="", label=None, at=None, anchor="middle", both=False, dash=False):
        """Polyline with an arrowhead at the end (and the start if ``both``).

        ``label`` is drawn at ``at`` (x, y); default is the midpoint of the
        longest segment, nudged above it.
        """
        lcls = "ln" + (f" {kind}" if kind else "") + (" dash" if dash else "")
        hcls = "hd" + (f" {kind}" if kind else "")
        (x1, y1), (x2, y2) = pts[-2], pts[-1]
        ang = math.atan2(y2 - y1, x2 - x1)
        # stop the line a little short so it doesn't poke through the head
        short = list(pts[:-1]) + [(x2 - 4 * math.cos(ang), y2 - 4 * math.sin(ang))]
        if both:
            (a1, b1), (a2, b2) = pts[0], pts[1]
            ang0 = math.atan2(b1 - b2, a1 - a2)
            short[0] = (a1 - 4 * math.cos(ang0), b1 - 4 * math.sin(ang0))
        self.line(short, lcls)
        self.head(x2, y2, ang, hcls)
        if both:
            self.head(pts[0][0], pts[0][1], ang0, hcls)
        if label:
            if at is None:
                segs = list(zip(pts[:-1], pts[1:]))
                (sx, sy), (ex, ey) = max(segs, key=lambda s: math.dist(*s))
                at = ((sx + ex) / 2, (sy + ey) / 2 - 6)
            self.text(at[0], at[1], label, "lbl", anchor)

    def node(self, x, y, w, h, title, sub=None, kind="", mono=False, sub2=None, title_cls=None):
        """A box with a centred title and up to two sub-lines."""
        self.rect(x, y, w, h, "bx" + (f" {kind}" if kind else ""))
        lines = [l for l in (sub, sub2) if l]
        cy = y + h / 2 + 4.5 - 7.5 * len(lines)
        self.text(x + w / 2, cy, title, title_cls or ("tm" if mono else "t"))
        for i, l in enumerate(lines):
            self.text(x + w / 2, cy + 16 + 14 * i, l, "s")

    def chip(self, x, y, w, h, s, kind="", cls="m"):
        self.rect(x, y, w, h, "chip" + (f" {kind}" if kind else ""), r=h / 2)
        self.text(x + w / 2, y + h / 2 + 4, s, cls)

    def cap(self, x, y, s, anchor="start"):
        self.text(x, y, s.upper(), "cap", anchor)

    def interval(self, x1, x2, y, h=10, kind="fm", title=None):
        self.rect(x1, y - h / 2, max(x2 - x1, 2), h, f"iv {kind}", r=2, title=title)

    # ── output ───────────────────────────────────────────────────────────
    def svg(self) -> str:
        body = "\n".join(self.parts)
        return (f'<svg class="gbd" viewBox="0 0 {self.w} {self.h}" role="img" '
                f'aria-label="{escape(self.label)}" xmlns="http://www.w3.org/2000/svg">\n'
                f"{body}\n</svg>\n")
