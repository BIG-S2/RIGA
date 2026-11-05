#!/usr/bin/env python
"""
supplementary_figure4_benchmark_ascertainment.py -- the real-genotype counterpart of fig_simulation_benchmark.png: power against type-I error,
both at the same nominal alpha, for the project-keyed designs on the unrelated ADSP pool (27,306 subjects, 281k SNPs).

Panel a  single variants: power at nominal 1e-5 on the clean causal SNPs (per-SD log-odds 0.09 = per-allele OR 1.15 at
         MAF 0.30) against the false-positive rate at 1e-5 on the 20 permuted-null replicates (Y* = the observed labels
         permuted inside cohort x project x centre cells; 5.5M null tests per design).
Panel b  rare-variant burden (chr22 QC7 gene mask, MAC 51-500): power at nominal 1e-3 on the causal genes (all three
         effect levels) against the false-positive rate at 1e-3 on the null genes of the same 20 replicates.
Designs  the generator's own cohort x project x centre cells (oracle, N0); pooled + centre + project (collapsed) + PCs
         (N5); REGENIE Firth + LOCO + centre + project (R5; single variants only, 3 replicates); RIGA 1:1 on the project
         key (M3); project x centre cells as a factor + PCs on every mixed-cell subject (S1); the conditional score
         test on the project x centre strata (C5).  Centre-only / cohort-only / 1:k variants: S79-S84 + tables.
Tags     main (within-cell generator, clean + t-gated causal set, REGENIE on 3 replicates; the main figure),
         asc_clean (ascertainment generator, paired seeds, no REGENIE; S86), clean (within-cell, clean causal set only),
         asc_clean_rg (asc_clean + REGENIE on its own replicates, if plasmode_power_regenie_asc.sbatch is ever run).
Style    fig_simulation_main.py (surface, ink, band within 2x of nominal, exact Clopper-Pearson bars, labels with the
         observed false-positive counts).
Usage    python supplementary_figure4_benchmark_ascertainment.py [--tag main] [--out fig_plasmode_benchmark]
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

REV = Path("${ADSP_ROOT}/revision")
TABP = REV / "tables" / "plasmode_power"
INK, INK2, INK3 = "#0b0b0b", "#52514e", "#8a8985"
GRID, SURF = "#e6e6e3", "#fcfcfb"
BLUE, GREEN, PURPLE, ORANGE = "#2a78d6", "#1baf7a", "#4a3aa7", "#eb6834"
BURDEN_MASK = "gene_mac51_500"      # overridden by --burden-mask
BURDEN_FPR_TABLE = TABP / "fpr_permuted_null_figure_table_maf001.csv"   # burden null FPR rows used with --burden-tag (override: --burden-fpr-table)
MASK_DESC = {"gene_mac51_500": "chr22 genes, MAC 51-500", "gene_mac6_50": "chr22 genes, MAC 6-50", "coding_mac6_50": "chr22 coding variants, MAC 6-50",
             "coding_maf001": "chr22 coding variants, MAF < 0.1%", "gene_maf001": "chr22 genes, all variants MAF < 0.1%",
             "coding_mac6_54": "chr22 coding variants, MAF < 0.1% without MAC 1-5", "gene_mac6_54": "chr22 genes, MAF < 0.1% without MAC 1-5"}
LEVEL = 1                       # effect level index for the single-variant panel (overridden by --level)
ALPHA = 1e-5                    # single-variant alpha (overridden by --alpha)
ATAG = {1e-3: "a1e3", 1e-4: "a1e4", 1e-5: "a1e5", 1e-6: "a1e6", 5e-8: "a5e8"}
BURDEN_ALPHA = 1e-3                                   # panel b alpha (--burden-alpha); the FPR tables label burden panels "burden <mask>, <BLAB>"
BLAB = {1e-3: "1e-3", 1e-4: "1e-4", 1e-5: "1e-5", 5e-8: "5e-8"}
ATXT = {1e-3: "10^{-3}", 1e-4: "10^{-4}", 1e-5: "10^{-5}", 1e-6: "10^{-6}", 5e-8: "5\\times10^{-8}"}
ORACLE = "saturated_cell"       # the generator's own cells (N0); present in every tag
# design -> (short label, full legend name, colour, marker)
DESIGNS = {
    "oracle":                 ("Oracle", "the generator's own cohort x project x center cells", INK3, "*"),
    "pooled_project":         ("Pooled (covariates)", "pooled logistic model, center + project (collapsed) + PC1-20", PURPLE, "^"),
    "regenie_pooled_project": ("REGENIE + center + project", "REGENIE Firth + LOCO, center + project (collapsed) + PC1-20", GREEN, "D"),
    "regenie_gmmat_project":  ("GMMAT + center + project", "GMMAT glmmkin + score test, sparse GRM, center + project (collapsed) + PC1-20", ORANGE, "v"),
    "regenie_genesis_project": ("GENESIS + center + project", "GENESIS fitNullModel + Score.SPA, sparse GRM, center + project (collapsed) + PC1-20", ORANGE, "h"),
    "riga_projectkey":        ("Matched 1:1", "RIGA one-to-one nearest-neighbor pairs on the project key", BLUE, "o"),
    "strat_project":          ("RIGA", "RIGA stratified design: project x center cells as a factor + PC1-20, every mixed-cell subject", BLUE, "s"),
}
GEN_NOTE = {
    "asc": "ascertainment generator (single-label projects keep the subjects whose latent status matches the label)",
    "clean": "within-cell generator (single-label projects keep their labels, so they carry no genetic signal)",
    "main": "within-cell generator (single-label projects keep their labels, so they carry no genetic signal)",
}
# label offsets in points: key -> (dx, dy, ha); a leader line is drawn when the label is displaced far from the marker
OFF_A = {"pooled_project": (11, 10, "left"), "regenie_pooled_project": (11, 8, "left"), "regenie_gmmat_project": (-11, 16, "right"), "regenie_genesis_project": (-11, -14, "right"),
         "riga_projectkey": (11, 8, "left"), "oracle": (16, -20, "left"), "strat_project": (16, 20, "left")}
OFF_B = {"riga_projectkey": (-11, 16, "right"), "pooled_project": (11, 8, "left"), "oracle": (16, 16, "left"),
         "strat_project": (16, -18, "left")}
OFF_TAG = {  # per-tag overrides
    "asc": {"a": {"pooled_project": (-11, 8, "right"), "oracle": (16, -20, "left"), "strat_project": (-16, 22, "right")},
            "b": {}},
}


def test_of(arm: str) -> str:
    return "firth" if arm.startswith("regenie") else ("chi2" if arm.endswith("_clr") else "spa")


def load(tag: str, fpr_table: Path, burden_tag: str = "") -> pd.DataFrame:
    pw = pd.read_csv(TABP / tag / "power_snp.csv")
    bt = burden_tag or tag
    bp = pd.read_csv(TABP / bt / "power_burden.csv") if (TABP / bt / "power_burden.csv").exists() else pd.DataFrame(columns=["arm", "test", "mask", "threshold", "alpha_tag", "level", "power", "se_rep", "n_rep"])
    arms = pd.read_csv(TABP / tag / "arms.csv").set_index("arm")
    fpr = pd.read_csv(fpr_table)
    fpr_b = pd.read_csv(fpr_table) if not burden_tag else pd.concat([fpr, pd.read_csv(BURDEN_FPR_TABLE)], ignore_index=True)   # --burden-fpr-table

    def g(df, c, cast=float, default=np.nan):
        return cast(df[c].iloc[0]) if len(df) else default

    rows = []
    for key, (short, full, col, mk) in DESIGNS.items():
        arm = ORACLE if key == "oracle" else key
        test = test_of(arm)
        d = pw[(pw.arm == arm) & (pw.test == test) & (pw.causal_class == "clean") & (pw.maf_bin == "all") &
               (pw.threshold == "nominal") & (pw.alpha_tag == ATAG[ALPHA]) & (pw.level == LEVEL)]
        f = fpr[(fpr.arm == arm) & (fpr.panel == f"single variants, {ALPHA:g}")]
        b = bp[(bp.arm == arm) & (bp.test == test) & (bp["mask"] == BURDEN_MASK) & (bp.threshold == "nominal") &
               (bp.alpha_tag == ATAG[BURDEN_ALPHA]) & (bp.level == -1)]
        fb = fpr_b[(fpr_b.arm == arm) & (fpr_b.panel == f"burden {BURDEN_MASK}, {BLAB[BURDEN_ALPHA]}")]
        n_subj = float(arms.n_subjects.get(arm, np.nan)) if arm in arms.index else np.nan
        if not np.isfinite(n_subj) and arm.startswith("regenie"):
            n_subj = 27306.0                                    # every external tool analyses the whole pool
            ru = REV.parent / "outputs" / "plasmode" / "power" / ("regenie_asc_clean" if tag.startswith("asc") else "regenie") / "units_plasmode_regenie.csv"  # noqa
            if ru.exists() and len(d):
                u = pd.read_csv(ru)
                u = u[u.arm == arm.replace("regenie_", "")]
                n_subj = float(u.n_subjects.mean()) if len(u) else np.nan
        rows.append(dict(key=key, short=short, full=full, colour=col, marker=mk, arm=arm, n_subjects=n_subj,
                         power=g(d, "power"), power_se=g(d, "se_rep"), n_rep=g(d, "n_rep", int, 0),
                         ratio=g(f, "fpr_1e3_eval_x"), lo=g(f, "fpr_1e3_eval_x_lo"), hi=g(f, "fpr_1e3_eval_x_hi"),
                         c=g(d, "c"), n_fp=g(f, "k" + ATAG[ALPHA][1:], lambda v: int(float(v)) if pd.notna(v) else 0, 0), n_null=g(f, "n", int, 0), n_rep_null=g(f, "n_units", int, 0),
                         bpower=g(b, "power"), bpower_se=g(b, "se_rep"), bn_rep=g(b, "n_rep", int, 0),
                         bratio=g(fb, "fpr_1e3_eval_x"), blo=g(fb, "fpr_1e3_eval_x_lo"), bhi=g(fb, "fpr_1e3_eval_x_hi"),
                         bn_fp=g(fb, "k" + ATAG[BURDEN_ALPHA][1:], lambda v: int(float(v)) if pd.notna(v) else 0, 0) if ("k" + ATAG[BURDEN_ALPHA][1:]) in fb.columns else 0, bn_null=g(fb, "n", int, 0), bn_rep_null=g(fb, "n_units", int, 0)))   # n_units = replicates in the FPR tables
    return pd.DataFrame(rows)


def style(ax, ylim):
    ax.set_facecolor(SURF)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(INK3)
        ax.spines[s].set_linewidth(0.6)
    ax.tick_params(colors=INK2, length=2.5, width=0.6, labelsize=7, pad=2)
    ax.grid(True, color=GRID, lw=0.5, zorder=0)
    ax.set_axisbelow(True)
    ax.set_yscale("log")
    ax.set_ylim(*ylim)
    ax.axhspan(0.5, 2.0, color="#eaf1fa", zorder=1, lw=0)
    ax.axhline(1.0, color=INK2, lw=0.9, ls=(0, (3, 2)), zorder=2)
    if ylim[1] / ylim[0] < 30:
        ticks = [t for t in (0.3, 0.5, 1, 2, 3, 5) if ylim[0] <= t <= ylim[1]]
        ax.set_yticks(ticks)
        ax.set_yticklabels([f"{t:g}" for t in ticks])
        ax.yaxis.set_minor_locator(mticker.NullLocator())


def draw(ax, d, xcol, secol, ycol, locol, hicol, kcol, xlim, ylim, xlabel, title, sub, offsets):
    style(ax, ylim)
    for _, r in d.iterrows():
        if not (np.isfinite(r[xcol]) and np.isfinite(r[ycol])):
            continue
        ax.plot([r[xcol], r[xcol]], [max(r[locol], ylim[0] * 1.02), min(r[hicol], ylim[1] * 0.98)], "-", color=r.colour, lw=1.1, alpha=0.5, zorder=3)
        ax.plot([r[xcol] - 1.96 * r[secol], r[xcol] + 1.96 * r[secol]], [r[ycol], r[ycol]], "-", color=r.colour, lw=1.1, alpha=0.5, zorder=3)
        ax.plot([r[xcol]], [r[ycol]], r.marker, ms=(12.0 if r.marker == "*" else 7.0), color=r.colour, mec="white", mew=0.9,
                zorder=(6 if r.key == "oracle" else 5), alpha=(0.9 if r.key == "oracle" else 1.0))
        dx, dy, ha = offsets.get(r.key, (11, 8, "left"))
        arrow = dict(arrowstyle="-", color=INK3, lw=0.6, shrinkA=0, shrinkB=3) if abs(dy) >= 14 or abs(dx) >= 24 else None
        ax.annotate(r.short, xy=(r[xcol], r[ycol]), xytext=(dx, dy), textcoords="offset points",
                    fontsize=6.6, color=INK, ha=ha, va="center", zorder=6, arrowprops=arrow)
    ax.set_xlim(*xlim)
    ax.set_xlabel(xlabel, fontsize=8, color=INK2, labelpad=4)
    ax.set_ylabel("type-I error, relative to nominal", fontsize=8, color=INK2, labelpad=4)
    ax.text(xlim[0] + 0.012 * (xlim[1] - xlim[0]), 0.58, "within 2$\\times$ of nominal", fontsize=6.2, color="#5b7fae", ha="left", va="center", zorder=4)
    yb = ylim[0] * (ylim[1] / ylim[0]) ** 0.05
    ax.annotate("better", xy=(xlim[1] - 0.02 * (xlim[1] - xlim[0]), yb), xytext=(xlim[1] - 0.17 * (xlim[1] - xlim[0]), yb),
                fontsize=6.4, color=INK3, va="center", arrowprops=dict(arrowstyle="->", color=INK3, lw=0.8))
    t1 = ax.text(0.0, 1.135, title, transform=ax.transAxes, fontsize=8.4, fontweight="bold", color=INK, va="bottom")
    t2 = ax.text(0.0, 1.12, sub, transform=ax.transAxes, fontsize=6.2, color=INK2, va="top", linespacing=1.5)


def main() -> None:
    global BURDEN_MASK, LEVEL, ALPHA, BURDEN_FPR_TABLE, BURDEN_ALPHA
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", default="main", help="main | main_rg | asc_clean | clean | asc_clean_rg (see the docstring)")
    ap.add_argument("--out", default="fig_plasmode_benchmark")
    ap.add_argument("--burden-mask", default=BURDEN_MASK, help="burden mask for panel b (power_burden.csv / FPR table name)")
    ap.add_argument("--fpr-table", default=str(TABP / "fpr_permuted_null_figure_table.csv"), help="permuted-null FPR table (plot_plasmode_null_fpr.py)")
    ap.add_argument("--alpha", type=float, default=1e-5, help="single-variant alpha: 1e-5 (default) or 5e-8 / 1e-6 / 1e-4 / 1e-3")
    ap.add_argument("--level", type=int, default=1, help="effect level index for the single-variant panel")
    ap.add_argument("--burden-tag", default="", help="take panel b from another tag's tables (e.g. main_maf001 for tag big)")
    ap.add_argument("--burden-fpr-table", default=str(BURDEN_FPR_TABLE), help="burden FPR table used with --burden-tag (default: the 20-replicate maf001 table)")
    ap.add_argument("--burden-alpha", type=float, default=1e-3, help="panel b alpha: 1e-3 (default), 1e-4, 1e-5 or 5e-8 (needs matching FPR-table rows)")
    ap.add_argument("--burden-desc", default="", help="panel b description of the burden units (default: MASK_DESC of the mask, chr22 wording)")
    args = ap.parse_args()
    BURDEN_MASK = args.burden_mask
    BURDEN_FPR_TABLE = Path(args.burden_fpr_table)
    BURDEN_ALPHA = args.burden_alpha
    if args.burden_desc:
        MASK_DESC[BURDEN_MASK] = args.burden_desc
    LEVEL = args.level
    ALPHA = args.alpha
    d = load(args.tag, Path(args.fpr_table), args.burden_tag)
    gkey = "asc" if args.tag.startswith("asc") else ("main" if args.tag.startswith("main") else "clean")
    gen = GEN_NOTE[gkey]
    off_a = dict(OFF_A, **OFF_TAG.get(gkey, {}).get("a", {}))
    off_b = dict(OFF_B, **OFF_TAG.get(gkey, {}).get("b", {}))
    ds = d[np.isfinite(d.power) & np.isfinite(d.ratio)]
    db = d[np.isfinite(d.bpower) & np.isfinite(d.bratio) & ~d.key.str.startswith("regenie")]
    n_null = int(ds.n_null.max())
    bn_null = int(np.nan_to_num(db.bn_null.max(), nan=0))
    reg = d[d.key.str.startswith("regenie") & np.isfinite(d.power)]
    reg_note = (f"external tools on {int(reg.n_rep.max())} of them" if len(reg) else "no external tool on these replicates")
    c_lvl = float(ds.c.dropna().iloc[0]) if ds.c.notna().any() else np.nan
    or_txt = f"per-allele OR {np.exp(c_lvl / np.sqrt(2 * 0.3 * 0.7)):.2f} at MAF 0.30" if np.isfinite(c_lvl) else ""
    a_txt = {1e-3: "10^{-3}", 1e-4: "10^{-4}", 1e-5: "10^{-5}", 1e-6: "10^{-6}", 5e-8: "5\\times10^{-8}"}[ALPHA]

    fig, axes = plt.subplots(1, 2, figsize=(13.4, 5.4), dpi=400)
    fig.patch.set_facecolor("white")
    ytop = float(max(500.0, 4.0 * np.nanmax(ds.hi)))          # the 5e-8 run reaches 8e4x for pooled + PCs
    draw(axes[0], ds, "power", "power_se", "ratio", "lo", "hi", "n_fp", (0.05, 0.80), (0.3, ytop),
         f"power at $\\alpha = {a_txt}$ (clean causal SNPs, {or_txt})",
         f"a  Single variants on real genotypes: error control at ${a_txt}$ against power",
         f"{n_null:,} null tests per design on {int(ds.n_rep_null.max())} permuted-null replicates; power on {int(ds.n_rep.max())} replicates, {reg_note}.\n"
         f"Bars: exact 95% CI for the type-I error and 1.96 x replicate SE for power.\nGenerator: {gen}.", off_a)
    ytop_b = float(max(6.0, 3.0 * np.nanmax(db.bhi))) if np.isfinite(np.nanmax(db.bhi)) else 6.0   # widen when a burden arm exceeds the narrow range (5e-8 panels)
    draw(axes[1], db, "bpower", "bpower_se", "bratio", "blo", "bhi", "bn_fp", (0.25, 0.95), (0.3, ytop_b),
         f"power at $\\alpha = {ATXT[BURDEN_ALPHA]}$ (causal genes, all effect levels)",
         f"b  Rare-variant burden on real genotypes ({MASK_DESC.get(BURDEN_MASK, BURDEN_MASK)}): error control at ${ATXT[BURDEN_ALPHA]}$ against power",
         f"{bn_null:,} null gene tests per design on {int(np.nan_to_num(db.bn_rep_null.max(), nan=0))} permuted-null replicates; power on {int(np.nan_to_num(db.bn_rep.max(), nan=0))} replicates\n"
         "with per-gene log-odds proportional to |log10 mean MAF|.  REGENIE has no burden test.\n"
         "Same conventions as panel a" + ("; note the narrower y-axis." if ytop_b <= 6 else "."), off_b)
    handles = []
    for _, r in d.iterrows():
        if not (np.isfinite(r.power) or np.isfinite(r.bpower)):
            continue
        n = f", {r.n_subjects:,.0f} subjects" if np.isfinite(r.n_subjects) else ""
        handles.append(Line2D([], [], marker=r.marker, ls="", ms=(10.0 if r.marker == "*" else 6.0), color=r.colour, mec="white", mew=0.8,
                              label=f"{r.short}: {r.full}{n}"))
    fig.legend(handles=handles, loc="lower center", ncol=2, frameon=False, fontsize=6.6, handletextpad=0.5, columnspacing=2.0,
               labelspacing=0.45, bbox_to_anchor=(0.5, -0.005))
    fig.subplots_adjust(left=0.055, right=0.985, top=0.855, bottom=0.215, wspace=0.24)
    out = REV / "figures" / f"{args.out}.png"
    fig.savefig(out, bbox_inches="tight", facecolor="white")
    fig.savefig(out.with_suffix(".pdf"), bbox_inches="tight", facecolor="white")
    plt.close(fig)
    d.to_csv(TABP / f"{args.out}_table.csv", index=False)
    print("wrote", out)
    print(d[["key", "arm", "n_subjects", "n_rep", "power", "power_se", "ratio", "n_fp", "bpower", "bpower_se", "bratio", "bn_fp"]].round(3).to_string(index=False))


if __name__ == "__main__":
    main()
