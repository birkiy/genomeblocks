"""Tiny HTML helpers for notebook reprs (no dependencies, theme-neutral)."""
from __future__ import annotations

from html import escape

_CSS = ("font-family:ui-monospace,Menlo,monospace;font-size:12px;border-collapse:collapse;"
        "margin:2px 0 6px")
_TH = "text-align:right;padding:2px 10px;border-bottom:1px solid #8884;font-weight:600"
_TD = "text-align:right;padding:2px 10px"


def table_html(title, head, rows, gap_after=None, note=None):
    h = "".join(f'<th style="{_TH}">{escape(str(c))}</th>' for c in head)
    out = [f'<div><b style="font-family:system-ui,sans-serif;font-size:13px">{escape(title)}</b>',
           f'<table style="{_CSS}"><thead><tr>{h}</tr></thead><tbody>']
    for k, r in enumerate(rows):
        out.append("<tr>" + "".join(f'<td style="{_TD}">{escape(str(c))}</td>' for c in r) + "</tr>")
        if gap_after is not None and k == gap_after:
            out.append("<tr>" + "".join(f'<td style="{_TD};opacity:.5">⋮</td>' for _ in head) + "</tr>")
    out.append("</tbody></table>")
    if note:
        out.append(f'<div style="font-family:system-ui,sans-serif;font-size:12px;opacity:.75">'
                   f'{escape(note)}</div>')
    out.append("</div>")
    return "".join(out)


def kv_html(title, items, note=None):
    rows = "".join(f'<tr><td style="{_TD};text-align:left;opacity:.7">{escape(str(k))}</td>'
                   f'<td style="{_TD};text-align:left">{escape(str(v))}</td></tr>' for k, v in items)
    n = (f'<div style="font-family:system-ui,sans-serif;font-size:12px;opacity:.75">{escape(note)}</div>'
         if note else "")
    return (f'<div><b style="font-family:system-ui,sans-serif;font-size:13px">{escape(title)}</b>'
            f'<table style="{_CSS}">{rows}</table>{n}</div>')
