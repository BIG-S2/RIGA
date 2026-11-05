#!/usr/bin/env python
"""supplementary_figure9_burden_power.py (2026-10-13) -- supplementary figure: coding burden power against type-I error at 5e-8 under the
designs, from the benchmark's internal burden test (10,000 permuted-null replicates; power on the sparse-spike replicates).
This panel was Figure 3b before the review. The internal burden test uses the saddlepoint approximation and is used only to
compare designs; the discovery statistic itself is evaluated with STAAR on the same null (Figure 3b)."""
import sys
sys.path.insert(0, "${ADSP_ROOT}/revision/scripts")
import supplementary_figure4_benchmark_ascertainment as F
import benchmark_figure_helpers as V
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
TABP = F.TABP
def main():
    F.BURDEN_MASK = "coding_maf001"; F.BURDEN_FPR_TABLE = TABP / "fpr_permuted_null_figure_table_maf001_bignull_genome10k.csv"
    F.BURDEN_ALPHA = 5e-8; F.LEVEL = 1; F.ALPHA = 5e-8; F.MASK_DESC[F.BURDEN_MASK] = "genome-wide coding variants, MAF < 0.1%"
    d = F.load("big", TABP / "fpr_bignull_figure_table.csv", "main_maf001_genome")
    V.FS.update(tick=6.5, label=7.0, ann=6.5, band=6.0, better=6.2, title=7.5, sub=6.0, legend=6.2); V.BETTER_SPAN = 0.30
    short = {"oracle": "Oracle", "pooled_project": "Pooled (covariates)", "riga_projectkey": "Matched 1:1", "strat_project": "RIGA"}
    db = d[np.isfinite(d.bpower) & np.isfinite(d.bratio) & ~d.key.str.startswith("regenie")].copy(); db["short"] = db.key.map(short)
    off_b = {"riga_projectkey": (-6, 12, "right"), "pooled_project": (7, 8, "left"), "oracle": (9, 10, "left"), "strat_project": (9, -12, "left")}
    a_txt = F.ATXT[F.ALPHA]; ytop = float(max(6.0, 3.0 * np.nanmax(db.bhi)))
    W, H = 3.4, 2.9
    fig = plt.figure(figsize=(W, H), dpi=400); fig.patch.set_facecolor("white")
    ax = fig.add_axes([0.55 / W, 0.75 / H, 2.6 / W, 1.55 / H])
    sb = f"{int(np.nan_to_num(db.bn_null.max(), nan=0)):,} null gene tests, {int(np.nan_to_num(db.bn_rep_null.max(), nan=0)):,} null replicates;\npower on {int(np.nan_to_num(db.bn_rep.max(), nan=0))} replicates"
    V.BETTER_Y = 0.93; V.BAND_X = 0.62; V.BAND_TXT = "within 2$\\times$"
    V.draw(ax, db, "bpower", "bpower_se", "bratio", "blo", "bhi", (0.25, 0.95), (0.3, ytop), f"power at $\\alpha = {a_txt}$ (causal genes)",
           f"Coding rare-variant burden (benchmark's internal test):\ntype-I error at ${a_txt}$ against power", sb, off_b, (25, 22))
    handles = [Line2D([], [], marker=r.marker, ls="", ms=(8.0 if r.marker == "*" else 4.8), color=r.colour, mec="white", mew=0.7, label=short[r.key]) for _, r in db.iterrows()]
    fig.legend(handles=handles, loc="lower center", ncol=4, frameon=False, fontsize=V.FS["legend"], handletextpad=0.3, columnspacing=1.0, bbox_to_anchor=(0.5, 0.0))
    out = F.REV / "figures" / "fig_burden_power_supp.png"
    fig.savefig(out, bbox_inches="tight", facecolor="white", pad_inches=0.04); fig.savefig(out.with_suffix(".pdf"), bbox_inches="tight", facecolor="white", pad_inches=0.04)
    plt.close(fig); print("wrote", out); print(db[["key", "bpower", "bpower_se", "bratio", "blo", "bhi", "bn_null", "bn_rep_null", "bn_rep"]].round(3).to_string(index=False))
if __name__ == "__main__":
    main()
