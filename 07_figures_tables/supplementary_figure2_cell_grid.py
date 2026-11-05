#!/usr/bin/env python3
"""supplementary_figure2_cell_grid.py (2026-10-05) -- example grid of the sequencing project x center
cells of the phenotyped ADSP pool (Supplementary Table 1f): cells with both cases and controls
colored by case fraction, single-label cells greyed. Full grid + a compact version for the
Figure 2 schematic. Reads revision/tables/supp_v3/S1f_project_center_cells.csv."""
import pandas as pd, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch, Rectangle
from matplotlib.colors import Normalize
from matplotlib import cm

ROOT = "${ADSP_ROOT}"
d = pd.read_csv(f"{ROOT}/revision/tables/supp_v3/S1f_project_center_cells.csv")
d = d.rename(columns={"Sequencing project": "project", "Sequencing center": "center", "Retained (both labels present)": "retained"})
d["frac"] = d["Cases"] / d["N"]
short = {"Kinghorn Centre for Clinical Genomics (KCCG), Sydney": "KCCG"}
d["center"] = d["center"].map(lambda s: short.get(s, s))
proj_order = d.groupby("project")["N"].sum().sort_values(ascending=False).index.tolist()
cent_order = d.groupby("center")["N"].sum().sort_values(ascending=False).index.tolist()
cmap = matplotlib.colormaps["Blues"]; norm = Normalize(0, 1)
GREY = "#c8c8c8"

def draw(projects, fname, title, cell_h=0.42, cell_w=0.95, fontsize=6.5):
    sub = d[d.project.isin(projects)]
    nrow, ncol = len(projects), len(cent_order)
    fig, ax = plt.subplots(figsize=(1.9 + cell_w * ncol, 1.6 + cell_h * nrow))
    for i, p in enumerate(projects):
        for j, c in enumerate(cent_order):
            r = sub[(sub.project == p) & (sub.center == c)]
            if r.empty: continue
            r = r.iloc[0]
            y = nrow - 1 - i
            if r.retained:
                col = cmap(0.15 + 0.75 * norm(r.frac)); txt = f"{int(r.N):,}\n{100*r.frac:.0f}% cases"
                ax.add_patch(Rectangle((j, y), 1, 1, facecolor=col, edgecolor="white", lw=1.2))
                ax.text(j + 0.5, y + 0.5, txt, ha="center", va="center", fontsize=fontsize,
                        color="white" if r.frac > 0.55 else "#102a43")
            else:
                lab = "all controls" if r.Cases == 0 else "all cases"
                ax.add_patch(Rectangle((j, y), 1, 1, facecolor=GREY, edgecolor="white", lw=1.2, hatch="///"))
                ax.text(j + 0.5, y + 0.5, f"{int(r.N):,}\n{lab}", ha="center", va="center", fontsize=fontsize, color="#444444")
    ax.set_xlim(0, ncol); ax.set_ylim(0, nrow)
    ax.set_xticks(np.arange(ncol) + 0.5); ax.set_xticklabels(cent_order, rotation=45, ha="left", fontsize=8)
    ax.xaxis.tick_top()
    ax.set_yticks(np.arange(nrow) + 0.5); ax.set_yticklabels(projects[::-1], fontsize=7.5)
    ax.tick_params(length=0)
    for s in ax.spines.values(): s.set_visible(False)
    ax.set_title(title, fontsize=10, loc="left", pad=36)
    sm = cm.ScalarMappable(norm=norm, cmap=cmap); sm.set_array([])
    cb = fig.colorbar(sm, ax=ax, fraction=0.03, pad=0.02, shrink=0.5 if nrow > 20 else 0.8)
    cb.set_label("case fraction (cells with both labels)", fontsize=8); cb.ax.tick_params(labelsize=7)
    ax.legend(handles=[Patch(facecolor=GREY, hatch="///", edgecolor="white", label="single-label cell: no within-cell contrast, set aside (15 cells, 7,235 subjects)"),
                       Patch(facecolor=cmap(0.55), edgecolor="white", label="informative cell: both cases and controls, stratum of the design (64 cells, 30,894 subjects)")],
              loc="upper left", bbox_to_anchor=(0, -0.01), fontsize=7.5, frameon=False)
    fig.savefig(fname, dpi=300, bbox_inches="tight", facecolor="white"); plt.close(fig)
    print("wrote", fname)

n_cells = len(d); n_ret = int(d.retained.sum())
print(f"{n_cells} cells, {n_ret} informative, {n_cells - n_ret} single-label; subjects retained {int(d.loc[d.retained,'N'].sum()):,} of {int(d.N.sum()):,}")
draw(proj_order, f"{ROOT}/revision/figures/fig2_cell_grid_full.png",
     f"Sequencing project x center cells of the phenotyped pool ({n_cells} cells; rows = projects, columns = centers)")
# compact version: the 16 largest projects plus the single-label projects, for the schematic
single = d.loc[~d.retained, "project"].unique().tolist()
compact = [p for p in proj_order[:16]] + [p for p in single if p not in proj_order[:16]][:6]
draw(compact, f"{ROOT}/revision/figures/fig2_cell_grid_compact.png",
     "Example for the Figure 2 schematic: 16 largest projects and six single-label projects", cell_h=0.5, fontsize=7)
