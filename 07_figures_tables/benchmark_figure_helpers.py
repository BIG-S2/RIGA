#!/usr/bin/env python
"""
benchmark_figure_helpers.py (2026-10-09) -- print-size re-render of Figure 3 (fig_plasmode_benchmark_5e8_both).
Same tables and data path as supplementary_figure4_benchmark_ascertainment.py (imported), but the figure is drawn at its final width (7.0 in)
with 7-8 pt type, so nothing is scaled down at \\linewidth.  --layout stacked (default): panels a and b one above the
other; --layout wide: side by side with short two-line titles.  Outputs revision/figures/<out>.png/.pdf.
"""
from __future__ import annotations
import argparse, sys
from pathlib import Path
sys.path.insert(0, "${ADSP_ROOT}/revision/scripts")
import supplementary_figure4_benchmark_ascertainment as F
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
from matplotlib.lines import Line2D

BETTER_SPAN = 0.17
BETTER_Y = 0.05        # vertical position of the "better" arrow, as a fraction of the log y-range
BAND_X = 0.012         # x position of the band label, as a fraction of the x-range
BAND_TXT = "within 2$\\times$ of nominal"
FS = dict(tick=7.0, label=7.5, ann=7.0, band=6.5, better=6.5, title=8.0, sub=6.5, legend=6.5)
TABP = F.TABP


def style(ax, ylim):
    ax.set_facecolor(F.SURF)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(F.INK3); ax.spines[s].set_linewidth(0.6)
    ax.tick_params(colors=F.INK2, length=2.5, width=0.6, labelsize=FS["tick"], pad=2)
    ax.grid(True, color=F.GRID, lw=0.5, zorder=0); ax.set_axisbelow(True)
    ax.set_yscale("log"); ax.set_ylim(*ylim)
    ax.axhspan(0.5, 2.0, color="#eaf1fa", zorder=1, lw=0)
    ax.axhline(1.0, color=F.INK2, lw=0.9, ls=(0, (3, 2)), zorder=2)
    if ylim[1] / ylim[0] < 30:
        ticks = [t for t in (0.3, 0.5, 1, 2, 3, 5) if ylim[0] <= t <= ylim[1]]
        ax.set_yticks(ticks); ax.set_yticklabels([f"{t:g}" for t in ticks]); ax.yaxis.set_minor_locator(mticker.NullLocator())


def draw(ax, d, xcol, secol, ycol, locol, hicol, xlim, ylim, xlabel, title, sub, offsets, title_pt):
    style(ax, ylim)
    for _, r in d.iterrows():
        if not (np.isfinite(r[xcol]) and np.isfinite(r[ycol])):
            continue
        ax.plot([r[xcol], r[xcol]], [max(r[locol], ylim[0] * 1.02), min(r[hicol], ylim[1] * 0.98)], "-", color=r.colour, lw=1.1, alpha=0.5, zorder=3)
        ax.plot([r[xcol] - 1.96 * r[secol], r[xcol] + 1.96 * r[secol]], [r[ycol], r[ycol]], "-", color=r.colour, lw=1.1, alpha=0.5, zorder=3)
        ax.plot([r[xcol]], [r[ycol]], r.marker, ms=(12.0 if r.marker == "*" else 7.0), color=r.colour, mec="white", mew=0.9,
                zorder=(6 if r.key == "oracle" else 5), alpha=(0.9 if r.key == "oracle" else 1.0))
        dx, dy, ha = offsets.get(r.key, (11, 8, "left"))
        arrow = dict(arrowstyle="-", color=F.INK3, lw=0.6, shrinkA=0, shrinkB=3) if abs(dy) >= 14 or abs(dx) >= 24 else None
        ax.annotate(r.short, xy=(r[xcol], r[ycol]), xytext=(dx, dy), textcoords="offset points", fontsize=FS["ann"], color=F.INK, ha=ha, va="center", zorder=6, arrowprops=arrow)
    ax.set_xlim(*xlim)
    ax.set_xlabel(xlabel, fontsize=FS["label"], color=F.INK2, labelpad=3)
    ax.set_ylabel("type-I error, relative to nominal", fontsize=FS["label"], color=F.INK2, labelpad=3)
    ax.text(xlim[0] + BAND_X * (xlim[1] - xlim[0]), 0.58, BAND_TXT, fontsize=FS["band"], color="#5b7fae", ha="left", va="center", zorder=4)
    yb = ylim[0] * (ylim[1] / ylim[0]) ** BETTER_Y
    ax.annotate("better", xy=(xlim[1] - 0.02 * (xlim[1] - xlim[0]), yb), xytext=(xlim[1] - BETTER_SPAN * (xlim[1] - xlim[0]), yb), fontsize=FS["better"], color=F.INK3, va="center", arrowprops=dict(arrowstyle="->", color=F.INK3, lw=0.8))
    # title above the subtitle block, both anchored in points above the axes so the layout does not depend on panel height
    ax.annotate(sub, xy=(0, 1), xycoords="axes fraction", xytext=(0, title_pt[1]), textcoords="offset points", fontsize=FS["sub"], color=F.INK2, va="top", ha="left", linespacing=1.4, annotation_clip=False)
    ax.annotate(title, xy=(0, 1), xycoords="axes fraction", xytext=(0, title_pt[0]), textcoords="offset points", fontsize=FS["title"], fontweight="bold", color=F.INK, va="bottom", ha="left", annotation_clip=False)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--layout", default="stacked", choices=["stacked", "wide"])
    ap.add_argument("--out", default="fig_plasmode_benchmark_5e8_both_v3")
    ap.add_argument("--tag", default="big"); ap.add_argument("--burden-tag", default="main_maf001_genome")
    ap.add_argument("--fpr-table", default=str(TABP / "fpr_bignull_figure_table.csv"))
    ap.add_argument("--burden-fpr-table", default=str(TABP / "fpr_permuted_null_figure_table_maf001_bignull_genome10k.csv"))
    args = ap.parse_args()
    F.BURDEN_MASK = "coding_maf001"; F.BURDEN_FPR_TABLE = Path(args.burden_fpr_table); F.BURDEN_ALPHA = 5e-8; F.LEVEL = 1; F.ALPHA = 5e-8
    F.MASK_DESC[F.BURDEN_MASK] = "genome-wide coding variants, MAF < 0.1%"
    d = F.load(args.tag, Path(args.fpr_table), args.burden_tag)
    gkey = "asc" if args.tag.startswith("asc") else ("main" if args.tag.startswith("main") else "clean")
    gen = F.GEN_NOTE[gkey]
    off_a = dict(F.OFF_A, **F.OFF_TAG.get(gkey, {}).get("a", {})); off_b = dict(F.OFF_B, **F.OFF_TAG.get(gkey, {}).get("b", {}))
    ds = d[np.isfinite(d.power) & np.isfinite(d.ratio)]
    db = d[np.isfinite(d.bpower) & np.isfinite(d.bratio) & ~d.key.str.startswith("regenie")]
    n_null = int(ds.n_null.max()); bn_null = int(np.nan_to_num(db.bn_null.max(), nan=0))
    reg = d[d.key.str.startswith("regenie") & np.isfinite(d.power)]
    reg_note = (f"external tools on {int(reg.n_rep.max())} of them" if len(reg) else "no external tool on these replicates")
    c_lvl = float(ds.c.dropna().iloc[0]) if ds.c.notna().any() else np.nan
    or_txt = f"per-allele OR {np.exp(c_lvl / np.sqrt(2 * 0.3 * 0.7)):.2f} at MAF 0.30" if np.isfinite(c_lvl) else ""
    a_txt = F.ATXT[F.ALPHA]; b_txt = F.ATXT[F.BURDEN_ALPHA]
    ytop = float(max(500.0, 4.0 * np.nanmax(ds.hi)))
    ytop_b = float(max(6.0, 3.0 * np.nanmax(db.bhi))) if np.isfinite(np.nanmax(db.bhi)) else 6.0
    wide = args.layout == "wide"
    if wide:
        # the Springer template prints figures at 5.17 in; design at that width (W = 5.2) so the type is 6-7 pt on the page
        global BETTER_SPAN; BETTER_SPAN = 0.30
        FS.update(tick=6.0, label=6.5, ann=6.2, band=5.6, better=5.8, title=7.0, sub=5.6, legend=5.8)
        W, H = 5.2, 3.5; fig = plt.figure(figsize=(W, H), dpi=400)
        axa = fig.add_axes([0.5 / W, 1.25 / H, 1.95 / W, 1.72 / H]); axb = fig.add_axes([3.0 / W, 1.25 / H, 1.95 / W, 1.72 / H])
        short = {"oracle": "Oracle", "pooled_project": "Pooled", "regenie_pooled_project": "REGENIE", "regenie_gmmat_project": "GMMAT",
                 "regenie_genesis_project": "GENESIS", "riga_projectkey": "RIGA 1:1", "strat_project": "Stratified"}
        d["short"] = d.key.map(short); ds = d[np.isfinite(d.power) & np.isfinite(d.ratio)]; db = d[np.isfinite(d.bpower) & np.isfinite(d.bratio) & ~d.key.str.startswith("regenie")]
        off_a = {"pooled_project": (7, 9, "left"), "regenie_pooled_project": (7, -3, "left"), "regenie_gmmat_project": (-7, 15, "right"),
                 "regenie_genesis_project": (-7, -15, "right"), "riga_projectkey": (7, 8, "left"), "oracle": (9, 9, "left"), "strat_project": (-6, 13, "right")}
        off_b = {"riga_projectkey": (-6, 12, "right"), "pooled_project": (7, 8, "left"), "oracle": (9, 10, "left"), "strat_project": (9, -12, "left")}
        ta = f"a  Single variants on real genotypes\n    type-I error at ${a_txt}$ against power"
        sa = (f"{n_null:,} null tests per design, {int(ds.n_rep_null.max()):,} null replicates\npower on {int(ds.n_rep.max())} replicates (external tools on {int(reg.n_rep.max()) if len(reg) else 0})\n"
              "bars: exact 95% CI (error), 1.96 x replicate SE (power)")
        tb = f"b  Coding rare-variant burden\n    type-I error at ${b_txt}$ against power"
        sb = (f"coding variants with MAF < 0.1%; {bn_null:,} null gene tests\nper design, {int(np.nan_to_num(db.bn_rep_null.max(), nan=0)):,} null replicates; power on {int(np.nan_to_num(db.bn_rep.max(), nan=0))}\n"
              "replicates; REGENIE has no burden test")
        xla = f"power at $\\alpha = {a_txt}$ (causal SNPs, OR 1.20)"; xlb = f"power at $\\alpha = {b_txt}$ (causal genes)"
        tpt = (33, 30); ncol = 1; lb = (0.5, 0.0)
    else:
        W, H = 7.0, 6.9; fig = plt.figure(figsize=(W, H), dpi=400)   # panels 1.9 in tall so the float (figure + caption) fits one page
        axa = fig.add_axes([0.62 / W, 4.45 / H, 6.2 / W, 1.9 / H]); axb = fig.add_axes([0.62 / W, 1.5 / H, 6.2 / W, 1.9 / H])
        ta = f"a  Single variants on real genotypes: error control at ${a_txt}$ against power"
        sa = (f"{n_null:,} null tests per design on {int(ds.n_rep_null.max())} permuted-null replicates; power on {int(ds.n_rep.max())} replicates, {reg_note}.\n"
              f"Bars: exact 95% CI for the type-I error and 1.96 x replicate SE for power.\nGenerator: {gen}.")
        off_a = dict(off_a, oracle=(18, -9, "left"))
        tb = f"b  Rare-variant burden on real genotypes (coding, MAF < 0.1%): error control at ${b_txt}$ against power"
        sb = (f"{bn_null:,} null gene tests per design on {int(np.nan_to_num(db.bn_rep_null.max(), nan=0))} permuted-null replicates; power on {int(np.nan_to_num(db.bn_rep.max(), nan=0))} replicates\n"
              "with per-gene log-odds proportional to |log10 mean MAF|.  REGENIE has no burden test.  Same conventions as panel a.")
        xla = f"power at $\\alpha = {a_txt}$ (clean causal SNPs, {or_txt})"; xlb = f"power at $\\alpha = {b_txt}$ (causal genes, all effect levels)"
        tpt = (36, 32); ncol = 1; lb = (0.5, 0.0)
    fig.patch.set_facecolor("white")
    draw(axa, ds, "power", "power_se", "ratio", "lo", "hi", (0.05, 0.80), (0.3, ytop), xla, ta, sa, off_a, tpt)
    draw(axb, db, "bpower", "bpower_se", "bratio", "blo", "bhi", (0.25, 0.95), (0.3, ytop_b), xlb, tb, sb, off_b, tpt)
    handles = []
    for _, r in d.iterrows():
        if not (np.isfinite(r.power) or np.isfinite(r.bpower)):
            continue
        n = f", {r.n_subjects:,.0f} subjects" if np.isfinite(r.n_subjects) else ""
        handles.append(Line2D([], [], marker=r.marker, ls="", ms=(10.0 if r.marker == "*" else 6.0), color=r.colour, mec="white", mew=0.8, label=f"{r.short}: {r.full}{n}"))
    fig.legend(handles=handles, loc="lower center", ncol=ncol, frameon=False, fontsize=FS["legend"], handletextpad=0.5, columnspacing=2.0, labelspacing=0.4, bbox_to_anchor=lb)
    out = F.REV / "figures" / f"{args.out}.png"
    fig.savefig(out, bbox_inches="tight", facecolor="white", pad_inches=0.04)
    fig.savefig(out.with_suffix(".pdf"), bbox_inches="tight", facecolor="white", pad_inches=0.04)
    plt.close(fig); print("wrote", out)


if __name__ == "__main__":
    main()
