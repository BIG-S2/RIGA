#!/usr/bin/env python
"""supplementary_figure7_noncoding_manhattan_qq.py (2026-10-13): as v1 but with the primary statistic SB (STAAR-O without ACAT-V) on sets with
cMAC >= 10, from tables/review/genesets_primary_{noncoding,ncRNA}.csv.gz (gene positions merged from the v1 export).
v1 docstring: Manhattan and Q-Q summary figures (STAARpipelineSummary layout) for the D_S gene-centric results, from the CSVs written
by export_gene_centric_results.R.  Families: noncoding (7 categories, one Manhattan panel each, one Q-Q with all
categories), ncRNA (one category), SKAT-MC (one p per gene).  Two threshold lines per panel: Bonferroni over the tests
of the family (solid) and the STAARpipelineSummary default (dashed).  Outputs S87/S88/S89 + two tables."""
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import beta, chi2

REV = Path("${ADSP_ROOT}/revision")
TAB, FIG = REV / "tables", REV / "figures"
CATS = ["UTR", "upstream", "downstream", "promoter_CAGE", "promoter_DHS", "enhancer_CAGE", "enhancer_DHS"]
PALETTE = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#17becf"]
INK, MUTED = "#1b1b1b", "#6b6b6b"
CHR2 = ["#25567a", "#9dbbd0"]
DEFAULT_ALPHA = {"noncoding": 3.57e-7, "ncRNA": 2.5e-6, "SKAT-MC": 2.5e-6}
SKATMC_NOTE = ""


def gc_lambda(p):
    return float(np.median(chi2.isf(p, 1)) / chi2.isf(0.5, 1))


def bh_q(p):
    n = p.size
    o = np.argsort(p)
    r = np.empty(n)
    r[o] = np.minimum.accumulate((p[o] * n / np.arange(1, n + 1))[::-1])[::-1]
    return np.minimum(r, 1.0)


def genome_axis(d):
    mx = d.groupby("chr").gene_start.max().sort_index()
    off = (mx.shift(1).fillna(0) * 1.02).cumsum()
    d = d.copy()
    d["gx"] = d.gene_start + d.chr.map(off)
    return d, off + mx / 2, float((off + mx).max())


def manhattan_panel(ax, d, gmax, mid, ymax, alpha_bonf, alpha_default, title):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.scatter(d.gx, -np.log10(d.p), s=5, c=[CHR2[c % 2] for c in d.chr], linewidths=0, rasterized=True)
    ax.axhline(-np.log10(alpha_bonf), color="#c1272d", lw=1.0)
    ax.axhline(-np.log10(alpha_default), color=MUTED, lw=0.9, ls=(0, (4, 3)))
    ax.set_xlim(0, gmax)
    ax.set_ylim(0, ymax)
    ax.set_xticks(mid.values)
    ax.set_xticklabels([str(c) for c in mid.index], fontsize=7, color=MUTED)
    ax.tick_params(axis="y", labelsize=7.5, colors=MUTED)
    ax.set_ylabel("-log10 p", fontsize=8, color=INK)
    ax.text(0.003, 0.96, title, transform=ax.transAxes, fontsize=8.5, fontweight="bold", color=INK, va="top")
    hits = d[d.p < alpha_bonf].sort_values("p")
    top_only = hits.empty
    if top_only:
        hits = d.nsmallest(1, "p")
    hits = hits.head(40).sort_values("gx")
    last_gx, stack = -1e12, 0
    for _, r in hits.iterrows():
        stack = stack + 1 if (r.gx - last_gx) < 0.012 * gmax else 0      # stagger labels of neighbouring genes
        last_gx = r.gx
        ax.annotate(r.gene, (r.gx, -np.log10(r.p)), xytext=(4, 3 + 9 * stack), textcoords="offset points", fontsize=7,
                    color=(MUTED if top_only else INK), style=("italic" if top_only else "normal"))


def qq_panel(ax, groups, title):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    nmax = max(len(p) for _, p in groups)
    k = np.arange(1, nmax + 1)
    ex = -np.log10(k / (nmax + 1))
    ax.fill_between(ex, -np.log10(beta.ppf(0.975, k, nmax - k + 1)), -np.log10(beta.ppf(0.025, k, nmax - k + 1)), color="#e8eef3", lw=0)
    xm = 0.0
    for (name, p), col in zip(groups, PALETTE):
        p = np.sort(p)
        e = -np.log10(np.arange(1, p.size + 1) / (p.size + 1))
        ax.scatter(e, -np.log10(p), s=6, c=col, linewidths=0, rasterized=True, label=f"{name}  (lambda {gc_lambda(p):.2f}, n {p.size:,})")
        xm = max(xm, float(e.max()))
    ax.plot([0, xm], [0, xm], color=MUTED, lw=0.9, ls=(0, (4, 3)))
    ax.set_xlabel("expected -log10 p", fontsize=8, color=INK)
    ax.set_ylabel("observed -log10 p", fontsize=8, color=INK)
    ax.tick_params(labelsize=7.5, colors=MUTED)
    ax.legend(frameon=False, fontsize=7, loc="upper left")
    ax.set_title(title, fontsize=8.5, fontweight="bold", color=INK, loc="left")


def family_tables(d, fam):
    n = len(d)
    alpha_bonf = 0.05 / n
    d = d.assign(q=bh_q(d.p.to_numpy()), family=fam)
    d["pass_bonferroni"] = d.p < alpha_bonf
    d["pass_default"] = d.p < DEFAULT_ALPHA[fam]
    d["pass_fdr05"] = d.q < 0.05
    summary = dict(family=fam, n_tests=n, n_categories=d.category.nunique(), lambda_all=round(gc_lambda(d.p.to_numpy()), 3),
                   bonferroni=alpha_bonf, default_alpha=DEFAULT_ALPHA[fam], n_bonferroni=int(d.pass_bonferroni.sum()),
                   n_default=int(d.pass_default.sum()), n_fdr05=int(d.pass_fdr05.sum()), n_below_1e5=int((d.p < 1e-5).sum()))
    return d, alpha_bonf, summary


def figure(fam, d, out, title):
    cats = [c for c in (CATS if fam == "noncoding" else [d.category.iloc[0]]) if c in set(d.category)]
    d, alpha_bonf, summary = family_tables(d, fam)
    dx, mid, gmax = genome_axis(d)
    ymax = float(np.ceil(-np.log10(d.p.min()) + 0.8))
    n_panels = len(cats)
    fig = plt.figure(figsize=(15, 1.55 * n_panels + 2.0), dpi=180)
    gs = fig.add_gridspec(n_panels, 2, width_ratios=[3.2, 1.0], hspace=0.28, wspace=0.16, top=0.94, bottom=0.08)
    axes = [fig.add_subplot(gs[i, 0]) for i in range(n_panels)]
    for ax, c in zip(axes, cats):
        manhattan_panel(ax, dx[dx.category == c], gmax, mid, ymax, alpha_bonf, DEFAULT_ALPHA[fam], c)
    for ax in axes[:-1]:
        ax.set_xticklabels([])
    axes[-1].set_xlabel("chromosome", fontsize=8, color=INK)
    axq = fig.add_subplot(gs[:, 1])
    qq_panel(axq, [(c, d[d.category == c].p.to_numpy()) for c in cats], "Q-Q, all categories")
    fig.suptitle(title, x=0.01, ha="left", fontsize=10, fontweight="bold", color=INK, y=1.01, va="bottom")   # above the panels, clear of the Q-Q heading
    foot = (f"Solid line: Bonferroni over the {summary['n_tests']:,} tests of this family (p < {alpha_bonf:.2e}); dashed line: "
            f"STAARpipelineSummary default (p < {DEFAULT_ALPHA[fam]:.2e}). Labels: genes below the Bonferroni line (grey italic: the panel's "
            f"top gene when none passes). {summary['n_fdr05']} tests at Benjamini-Hochberg q < 0.05. Q-Q band: 95% pointwise. "
            "")
    if fam == "SKAT-MC" and SKATMC_NOTE:
        foot += SKATMC_NOTE                                   # a supplied note describes its own design (replaces the D_S sentence)
    else:
        foot += "D_S: project x center cells as strata, all subjects, QC_label7; SKAT–burden omnibus test = Cauchy combination of all SKAT and burden p values (STAAR-O without ACAT-V), sets with cMAC >= 10."
    if fam == "SKAT-MC" and not SKATMC_NOTE:
        foot += (" CALIBRATION CHECK (2026-09-26): with the outcome permuted inside every project x center cell the same windows give "
                 "lambda 2.02, RP9 3.5e-33 and DTX1 4.7e-29, so these p-values do not reflect the phenotype; the SKAT-MC test is not "
                 "calibrated on this design and its results are not interpretable as they stand.")
    fig.text(0.01, -0.012, foot, fontsize=7, color=MUTED, wrap=True, va="top")
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    print("wrote", out)
    return d, summary


def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--skatmc-csv", default=str(TAB / "skatmc_projcell_all.csv"))
    ap.add_argument("--skatmc-out", default="supp_fig_S89_skatmc_manhattan_qq_projcell.png")
    ap.add_argument("--skatmc-title", default="Supplementary Figure S89. SKAT-MC multicategorical test (Control < MCI < AD), D_S, one test per gene")
    ap.add_argument("--skatmc-note", default="", help="extra footnote text for the SKAT-MC figure (replaces the D_S calibration note)")
    ap.add_argument("--skatmc-only", action="store_true")
    ap.add_argument("--tables-suffix", default="")
    args = ap.parse_args()
    global SKATMC_NOTE
    SKATMC_NOTE = args.skatmc_note
    v1 = pd.read_csv(TAB / "gene_centric_noncoding_projcell_all.csv")
    pos = v1[["gene", "category", "chr", "gene_start"]].drop_duplicates(["gene", "category"])
    parts = []
    for fam in ("noncoding", "ncRNA"):
        t = pd.read_csv(TAB / "review" / f"genesets_primary_{fam}.csv.gz"); t = t[t.tested == True]
        t = t.rename(columns={"Gene name": "gene", "Category": "category", "#SNV": "n_variants"})[["gene", "category", "n_variants", "cMAC", "SB"]].rename(columns={"SB": "p"})
        parts.append(t.merge(pos, on=["gene", "category"], how="left"))
    g = pd.concat(parts, ignore_index=True); g = g[g.gene_start.notna() & g.chr.between(1, 22)]
    runs = [("noncoding", g[g.category.isin(CATS)], "supp_fig_S87_noncoding_manhattan_qq_projcell_SB.png",
             "Supplementary Figure 7. Gene-centric noncoding rare-variant sets with cMAC >= 10, SKAT–burden omnibus test, D_S (stratified on project x center), 7 categories"),
            ("ncRNA", g[g.category == "ncRNA"], "supp_fig_S88_ncRNA_manhattan_qq_projcell_SB.png",
             "Supplementary Figure 8. Gene-centric ncRNA rare-variant sets with cMAC >= 10, SKAT–burden omnibus test, D_S")]
    tables, summaries = [], []
    for fam, d, fname, title in runs:
        if d.empty:
            print("no rows for", fam)
            continue
        dd, s = figure(fam, d, FIG / fname, title)
        tables.append(dd)
        summaries.append(s)
    summ = pd.DataFrame(summaries)
    summ.to_csv(TAB / f"gene_centric_summary_projcell_SB{args.tables_suffix}.csv", index=False)
    hits = pd.concat(tables)
    hits = hits[hits.p < 1e-5].sort_values(["family", "p"])
    hits.to_csv(TAB / f"gene_centric_hits_projcell_SB{args.tables_suffix}.csv", index=False)
    print(summ.to_string(index=False))
    print(hits[hits.pass_bonferroni][["family", "gene", "category", "chr", "n_variants", "p", "q"]].to_string(index=False))


if __name__ == "__main__":
    main()
