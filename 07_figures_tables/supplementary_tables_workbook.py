#!/usr/bin/env python
"""supplementary_tables_workbook.py (2026-10-14): Supplementary Tables 1-20 as ONE Excel workbook, one tab per table part, built from the same
CSV sources and column specifications as the LaTeX tables (SPECS in supplementary_tables_specifications.py) under the Nature Genetics numbering of
manuscript/supplementary_tables_numbering.py. Each tab carries the renumbered title and legend (LaTeX converted to plain text) above the header
row; values are written as numbers at full precision with display formats taken from the LaTeX spec (thousands separators, fixed
decimals, scientific notation for p values). A Contents tab lists every tab with a link. Writes manuscript/Supplementary_Tables.xlsx."""
import importlib.util, pathlib, re, sys
from openpyxl import Workbook
from openpyxl.styles import Alignment, Border, Font, PatternFill, Side
from openpyxl.utils import get_column_letter

REV = pathlib.Path("${ADSP_ROOT}/revision")
CSVDIR = REV / "tables/supp_v3"; OUT = REV / "manuscript/Supplementary_Tables.xlsx"

def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path); m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m); return m
T = load("t2t", REV / "scripts/supplementary_tables_specifications.py"); NG = load("ng", REV / "manuscript/supplementary_tables_numbering.py")

REP = [(r"\times", "×"), (r"\geq", "≥"), (r"\leq", "≤"), (r"\hat\beta", "β̂"), (r"\hat{\beta}", "β̂"), (r"\beta", "β"), (r"\lambda", "λ"), (r"\alpha", "α"),
       (r"\chi", "χ"), (r"\pm", "±"), (r"\%", "%"), (r"\&", "&"), (r"\_", "_"), (r"\#", "#"), (r"\,", ""), (r"\allowbreak{}", ""), (r"\newline", " "),
       (r"\\", " "), ("{,}", ","), ("~", " "), ("--", "–"), (r"\ ", " ")]
def delatex(s):
    for _ in range(4): s = re.sub(r"\\(textbf|textit|texttt|emph|mathrm|textrm|text|mathit|shortstack)\{([^{}]*)\}", r"\2", s)
    s = re.sub(r"\\ref\{[^}]*\}", "", s)
    for a, b in REP: s = s.replace(a, b)
    s = re.sub(r"\^\{([^}]*)\}", r"^\1", s); s = re.sub(r"_\{([^}]*)\}", r"_\1", s)
    s = s.replace("$", "").replace("{", "").replace("}", "")
    return re.sub(r"\s+", " ", s).strip()

def split_bold(cap):
    """(title, legend) from a caption that starts with \\textbf{...}; brace-matched so nested \\textit{} inside the title survives"""
    i = cap.find(r"\textbf{")
    if i < 0: return delatex(cap), ""
    j = i + len(r"\textbf{"); depth = 1; k = j
    while k < len(cap) and depth: depth += {"{": 1, "}": -1}.get(cap[k], 0); k += 1
    return delatex(cap[j:k - 1]), delatex(cap[k:])

def fmt_of(f):
    if f is T.f_int: return "num", "#,##0"
    if f is T.f_sci: return "num", "0.00E+00"
    if getattr(f, "__qualname__", "").startswith("f_num"): return "num", "0." + "0" * f.__closure__[0].cell_contents if f.__closure__[0].cell_contents else "0"
    if f in (T.f_family, T.f_source, T.f_tier, T.f_verdict, T.f_novelty, T.f_bool, T.f_test): return "map", None
    return "text", None

def colinfo(c):
    """(header, csv column or None, formatter or row callable, kind, number format); derived columns are callables on the row"""
    h, spec = c
    if isinstance(spec, tuple): name, f = spec; kind, nf = fmt_of(f); return h, name, f, kind, nf
    return h, None, spec, "row", None

def value(v, kind, f):
    if T.is_na(v): return None
    if kind == "num":
        try: x = float(v); return int(x) if x == int(x) and abs(x) < 1e15 and "e" not in str(v).lower() else x
        except ValueError: return str(v)
    if kind == "map": return delatex(f(v))
    return str(v)

HEAD = Font(bold=True); FILL = PatternFill("solid", fgColor="E6E6E6"); LINE = Border(bottom=Side(style="thin"))
wb = Workbook(); toc = wb.active; toc.title = "Contents"
toc["A1"] = "Supplementary Tables 1–20"; toc["A1"].font = Font(bold=True, size=13)
toc["A2"] = "One table per tab, named ST<number><part>; the title and legend sit above each table. The Supplementary Information PDF holds the Supplementary Figures and Notes."
toc.append([]); toc.append(["Table", "Tab", "Title", "Legend"]); [setattr(c, "font", HEAD) for c in toc[4]]
entries = []
for key, sp in T.SPECS.items():
    m = re.fullmatch(r"S(\d+)([a-z]?)", key)
    if m is None or int(m[1]) not in NG.NUMBER: print("skip", key); continue
    path = CSVDIR / sp["csv"]
    if not path.exists(): print(f"{key}: MISSING {path.name}"); continue
    entries.append((NG.NUMBER[int(m[1])], m[2], key, sp, path))
entries.sort(key=lambda e: (e[0], e[1]))
for new, part, key, sp, path in entries:
    tab = f"ST{new}{part}"; ws = wb.create_sheet(tab)
    title, legend = split_bold(NG.renumber(sp["caption"]))
    rows = T.read_csv(path); cols = sp["cols"]; ncol = len(cols)
    ws["A1"] = title; ws["A1"].font = Font(bold=True, size=12)
    ws["A2"] = legend; ws["A2"].alignment = Alignment(wrap_text=True, vertical="top")
    ws.merge_cells(start_row=2, start_column=1, end_row=2, end_column=max(ncol, 6))
    ws.row_dimensions[2].height = 15 * max(1, min(12, len(legend) // 110 + 1))
    hdr = [delatex(h) for h, _ in cols]; ws.append([]); ws.append(hdr)
    for c in ws[4]: c.font = HEAD; c.fill = FILL; c.border = LINE; c.alignment = Alignment(wrap_text=True, vertical="bottom")
    if any("SKAT–burden omnibus test" in h for h in hdr): ws.row_dimensions[4].height = 60
    info = [colinfo(c) for c in cols]
    for r in rows: ws.append([delatex(f(r)) if kind == "row" else value(r.get(name), kind, f) for (_, name, f, kind, _) in info])
    for j, (_, _, _, kind, nf) in enumerate(info, 1):
        if nf:
            for c in ws.iter_rows(min_row=5, max_row=4 + len(rows), min_col=j, max_col=j):
                if isinstance(c[0].value, (int, float)): c[0].number_format = nf
        width = max([len(hdr[j - 1])] + [len(str(v)) for v in (ws.cell(row=i, column=j).value for i in range(5, min(4 + len(rows), 300) + 1)) if v is not None])
        ws.column_dimensions[get_column_letter(j)].width = min(max(8, width + 2), 60)
        if "SKAT–burden omnibus test" in hdr[j - 1]:
            ws.column_dimensions[get_column_letter(j)].width = min(ws.column_dimensions[get_column_letter(j)].width, 24)
    if sp.get("note"): ws.cell(row=6 + len(rows), column=1, value="Note: " + delatex(sp["note"]))
    ws.freeze_panes = "A5"
    toc.append([f"Supplementary Table {new}{part}", tab, title.split(". ", 1)[1] if ". " in title else title, legend])
    toc.cell(row=toc.max_row, column=2).hyperlink = f"#'{tab}'!A1"; toc.cell(row=toc.max_row, column=2).font = Font(color="0563C1", underline="single")
    print(f"{key} -> {tab}: {len(rows)} rows x {ncol} cols")
for col, w in zip("ABCD", (26, 10, 90, 120)): toc.column_dimensions[col].width = w
for r in toc.iter_rows(min_row=5): r[3].alignment = Alignment(wrap_text=True, vertical="top"); r[2].alignment = Alignment(wrap_text=True, vertical="top")
toc.freeze_panes = "A5"
wb.save(OUT); print("wrote", OUT, "tabs:", len(wb.sheetnames))
