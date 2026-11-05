"""Publication-ready figures (Nature Genetics styling) for the 6-scenario
simulation benchmark.

Outputs (to --out-dir, default outputs/power_pure_only):
  main_figure.{png,pdf,svg}      4-panel composite figure (a-d)
  fig_calibration.{png,pdf}      lambda_GC + tail FPR side-by-side (panel a+b)
  fig_fpr_levels.{png,pdf}       FPR/alpha at four alpha levels
  fig_power.{png,pdf}            Power at GW significance
  cache_metrics.parquet          Aggregated long-form metrics (for fast re-runs)

Styling choices:
  - Wong colorblind-safe palette
  - Sans-serif (Arial / DejaVu Sans fallback)
  - Larger publication fonts, minimal grid, no top/right spines
  - Panel labels (a, b, c, d) in bold uppercase at top-left
  - Vector PDF as primary deliverable
"""
from __future__ import annotations
import argparse
import glob
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
import numpy as np
import pandas as pd

# -- Nature Genetics typography ---------------------------------------------
# Prefer Arial, fall back to Helvetica, then DejaVu Sans (matplotlib default).
PREFERRED_FONTS = ["Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"]
_available = {f.name for f in font_manager.fontManager.ttflist}
FONT_FAMILY = next((f for f in PREFERRED_FONTS if f in _available), "DejaVu Sans")

plt.rcParams.update({
    "font.family":           FONT_FAMILY,
    "font.size":             11,
    "axes.titlesize":        13,
    "axes.titleweight":      "bold",
    "axes.labelsize":        12,
    "axes.labelweight":      "regular",
    "axes.linewidth":        1.0,
    "axes.spines.top":       False,
    "axes.spines.right":     False,
    "xtick.labelsize":       10,
    "ytick.labelsize":       10,
    "xtick.major.width":     1.0,
    "ytick.major.width":     1.0,
    "xtick.major.size":      4,
    "ytick.major.size":      4,
    "legend.fontsize":       10,
    "legend.title_fontsize": 11,
    "legend.frameon":        False,
    "figure.titlesize":      14,
    "figure.dpi":            100,
    "savefig.dpi":           400,
    "savefig.bbox":          "tight",
    "savefig.transparent":   False,
    "pdf.fonttype":          42,   # editable text in vector PDFs
    "ps.fonttype":           42,
})

# -- Scenario / method ordering ---------------------------------------------
LEGACY_KEY_REMAP = {
    "2pure_same_cov":   "2pure_homo",
    "2pure_hetero_cov": "2pure_hetero",
}
DROP_KEYS = {"2SamePop_hetero_cov"}

SCENARIOS = [
    "2pure_homo",
    "2pure_hetero",
    "1pure_1admixed_homo",
    "1pure_1admixed_hetero",
    "2admixed_homo",
    "2admixed_hetero",
]
SCEN_LABEL = {
    "2pure_homo":            "2 pure\nhomog.",
    "2pure_hetero":          "2 pure\nheterog.",
    "1pure_1admixed_homo":   "1 pure +\n1 admixed\nhomog.",
    "1pure_1admixed_hetero": "1 pure +\n1 admixed\nheterog.",
    "2admixed_homo":         "2 admixed\nhomog.",
    "2admixed_hetero":       "2 admixed\nheterog.",
}

METHODS = [
    "combined_pcs",
    "combined_population",
    "combined_study_pcs",
    "combined_population_pcs",
    "matched_population",
]
METHOD_LABEL = {
    "combined_pcs":            "Pooled + PCs",
    "combined_population":     "Pooled + Pop",
    "combined_study_pcs":      "Pooled + Study + PCs",
    "combined_population_pcs": "Pooled + Pop + PCs",
    "matched_population":      "Matched 1:1 + Pop",
}
# Wong colorblind-safe palette; RIGA highlighted in vermillion for emphasis.
METHOD_COLOR = {
    "combined_pcs":            "#0072B2",   # blue
    "combined_population":     "#009E73",   # bluish green
    "combined_study_pcs":      "#56B4E9",   # sky blue
    "combined_population_pcs": "#E69F00",   # orange
    "matched_population":      "#D55E00",   # vermillion (RIGA)
}
METHOD_MARKER = {
    "combined_pcs":            "o",
    "combined_population":     "s",
    "combined_study_pcs":      "^",
    "combined_population_pcs": "D",
    "matched_population":      "*",
}
ALPHA_INFO = [
    (5e-8, "5e8", r"$\alpha = 5\!\times\!10^{-8}$"),
    (5e-7, "5e7", r"$\alpha = 5\!\times\!10^{-7}$"),
    (1e-5, "1e5", r"$\alpha = 1\!\times\!10^{-5}$"),
    (1e-3, "1e3", r"$\alpha = 1\!\times\!10^{-3}$"),
]


# -- Data loading -----------------------------------------------------------
def _load(run_dir: str, cache: Path) -> pd.DataFrame:
    """Load (or cache) the aggregated long-form metrics dataframe."""
    if cache.exists():
        print(f"Loading cached metrics: {cache}", flush=True)
        df = pd.read_parquet(cache)
    else:
        files = sorted(glob.glob(f"{run_dir}/power_run_*_metrics.csv"))
        if not files:
            raise FileNotFoundError(run_dir)
        print(f"Loading {len(files)} CSV files from {run_dir} ...", flush=True)
        df = pd.concat([pd.read_csv(f) for f in files], ignore_index=True)
        df = df[~df["scenario_key"].isin(DROP_KEYS)].copy()
        df["scenario_key"] = df["scenario_key"].replace(LEGACY_KEY_REMAP)
        cache.parent.mkdir(parents=True, exist_ok=True)
        df.to_parquet(cache, index=False)
        print(f"Cached metrics -> {cache}", flush=True)
    df = df[df["scenario_key"].isin(SCENARIOS)]
    df = df[df["method"].isin(METHODS)]
    return df


# -- Helpers ----------------------------------------------------------------
def _x_positions(jitter: float = 0.16):
    scen_x = np.arange(len(SCENARIOS))
    width = jitter * 2
    method_off = np.linspace(-width, width, len(METHODS))
    return scen_x, method_off


def _fmt(ratio: float) -> str:
    if ratio < 10:
        return f"{ratio:.2f}"
    if ratio < 100:
        return f"{ratio:.1f}"
    return f"{ratio:.0f}"


def _set_scenario_axis(ax, rotate: int = 0):
    scen_x = np.arange(len(SCENARIOS))
    ax.set_xticks(scen_x)
    ax.set_xticklabels([SCEN_LABEL[s] for s in SCENARIOS],
                       rotation=rotate, ha="center", linespacing=1.1)
    ax.set_xlim(scen_x[0] - 0.5, scen_x[-1] + 0.5)


def _annotate(ax, x, y, txt, dy=8, fontsize=8):
    ax.annotate(txt, (x, y), textcoords="offset points",
                xytext=(0, dy), ha="center",
                fontsize=fontsize, color="#222222")


def _panel_label(ax, label, x=-0.13, y=1.04, fontsize=15):
    ax.text(x, y, label, transform=ax.transAxes,
            fontsize=fontsize, fontweight="bold", va="bottom", ha="left")


def _legend_handles():
    """Build legend handles (one per method + nominal line)."""
    handles = []
    for m in METHODS:
        h = plt.Line2D([0], [0], marker=METHOD_MARKER[m], linestyle="",
                       markersize=10, color=METHOD_COLOR[m],
                       markeredgecolor="black", markeredgewidth=0.6,
                       label=METHOD_LABEL[m])
        handles.append(h)
    nominal = plt.Line2D([0], [0], color="#666666", linestyle="--",
                         linewidth=1.2, label="Nominal (calibrated)")
    return handles + [nominal]


# -- Individual panels ------------------------------------------------------
def panel_lambda(ax, df):
    scen_x, method_off = _x_positions()
    for j, m in enumerate(METHODS):
        for i, sc in enumerate(SCENARIOS):
            sub = df[(df["scenario_key"] == sc) & (df["method"] == m)]
            mean = sub["lambda_gc_null"].mean()
            sd   = sub["lambda_gc_null"].std()
            x = scen_x[i] + method_off[j]
            ax.errorbar(x, mean, yerr=sd, fmt=METHOD_MARKER[m],
                        color=METHOD_COLOR[m], markersize=9,
                        capsize=2.5, capthick=0.9, elinewidth=0.9,
                        markeredgecolor="black", markeredgewidth=0.5,
                        zorder=3)
    ax.axhline(1.0, color="#666666", linestyle="--", linewidth=1.2, zorder=1)
    _set_scenario_axis(ax)
    ax.set_ylabel(r"$\lambda_{GC}^{\,\mathrm{null}}$")
    ax.set_title("Genomic-control inflation under the null", pad=8)
    ax.grid(True, axis="y", alpha=0.18, linewidth=0.6)
    ax.set_axisbelow(True)


def panel_tail_fpr(ax, df, alpha=5e-8, alpha_tag="5e8",
                   alpha_label=r"$\alpha\!=\!5\!\times\!10^{-8}$"):
    scen_x, method_off = _x_positions()
    for j, m in enumerate(METHODS):
        for i, sc in enumerate(SCENARIOS):
            sub = df[(df["scenario_key"] == sc) & (df["method"] == m)]
            ratio = sub[f"fpr_{alpha_tag}"].mean() / alpha
            x = scen_x[i] + method_off[j]
            ax.scatter(x, ratio, marker=METHOD_MARKER[m], s=110,
                       color=METHOD_COLOR[m], edgecolor="black",
                       linewidth=0.5, zorder=3)
    ax.axhline(1.0, color="#666666", linestyle="--", linewidth=1.2, zorder=1)
    _set_scenario_axis(ax)
    ax.set_ylabel(r"FPR / $\alpha$  (log scale)")
    ax.set_title(f"Tail false-positive rate at {alpha_label}", pad=8)
    ax.set_yscale("log")
    ax.set_ylim(0.3, 2500)
    ax.grid(True, axis="y", which="both", alpha=0.18, linewidth=0.6)
    ax.set_axisbelow(True)


def panel_fpr_alpha_grid(ax_list, df):
    """Four panels of FPR/alpha at increasingly stringent alpha levels."""
    scen_x, method_off = _x_positions()
    for ax_idx, (ax, (alpha, atag, label)) in enumerate(zip(ax_list, ALPHA_INFO)):
        col = f"fpr_{atag}"
        for j, m in enumerate(METHODS):
            for i, sc in enumerate(SCENARIOS):
                sub = df[(df["scenario_key"] == sc) & (df["method"] == m)]
                ratio = sub[col].mean() / alpha
                x = scen_x[i] + method_off[j]
                ax.scatter(x, ratio, marker=METHOD_MARKER[m], s=90,
                           color=METHOD_COLOR[m], edgecolor="black",
                           linewidth=0.5, zorder=3)
        ax.axhline(1.0, color="#666666", linestyle="--", linewidth=1.2, zorder=1)
        _set_scenario_axis(ax)
        ax.set_title(label, pad=6, fontsize=12)
        ax.set_yscale("log")
        ax.set_ylim(0.3, 2500)
        ax.grid(True, axis="y", which="both", alpha=0.18, linewidth=0.6)
        ax.set_axisbelow(True)
        if ax_idx == 0:
            ax.set_ylabel(r"FPR / $\alpha$  (log scale)")


def panel_power(ax, df):
    scen_x, method_off = _x_positions()
    for j, m in enumerate(METHODS):
        for i, sc in enumerate(SCENARIOS):
            sub = df[(df["scenario_key"] == sc) & (df["method"] == m)]
            mean = sub["power_5e8"].mean()
            sd = sub["power_5e8"].std()
            x = scen_x[i] + method_off[j]
            ax.errorbar(x, mean, yerr=sd, fmt=METHOD_MARKER[m],
                        color=METHOD_COLOR[m], markersize=9,
                        capsize=2.5, capthick=0.9, elinewidth=0.9,
                        markeredgecolor="black", markeredgewidth=0.5,
                        zorder=3)
    _set_scenario_axis(ax)
    ax.set_ylabel(r"Power at $\alpha\!=\!5\!\times\!10^{-8}$")
    ax.set_title("Statistical power", pad=8)
    ax.set_ylim(0, 1.0)
    ax.grid(True, axis="y", alpha=0.18, linewidth=0.6)
    ax.set_axisbelow(True)


# -- Composite figures ------------------------------------------------------
def main_figure(df, out_path: Path):
    """Nature-Genetics-style 4-panel composite (a-d)."""
    fig = plt.figure(figsize=(13, 9.5), constrained_layout=False)
    gs = fig.add_gridspec(
        nrows=2, ncols=2,
        left=0.07, right=0.78, top=0.94, bottom=0.10,
        hspace=0.55, wspace=0.28,
    )
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[1, 0])
    ax_d = fig.add_subplot(gs[1, 1])

    panel_lambda(ax_a, df)
    panel_tail_fpr(ax_b, df)
    panel_tail_fpr(ax_c, df, alpha=5e-7, alpha_tag="5e7",
                   alpha_label=r"$\alpha\!=\!5\!\times\!10^{-7}$")
    panel_power(ax_d, df)

    _panel_label(ax_a, "a")
    _panel_label(ax_b, "b")
    _panel_label(ax_c, "c")
    _panel_label(ax_d, "d")

    handles = _legend_handles()
    fig.legend(
        handles=handles,
        labels=[h.get_label() for h in handles],
        loc="center left",
        bbox_to_anchor=(0.79, 0.5),
        frameon=False, title="Method",
        handletextpad=0.6, borderaxespad=0.0,
    )

    for ext in ("png", "pdf", "svg"):
        fig.savefig(out_path.with_suffix(f".{ext}"))
    plt.close(fig)


def fig_calibration(df, out_path: Path):
    # Bump every font size by +3 for this figure only (scoped via rc_context).
    _bump = {
        "font.size":             14,   # 11 + 3
        "axes.titlesize":        16,   # 13 + 3
        "axes.labelsize":        15,   # 12 + 3
        "xtick.labelsize":       13,   # 10 + 3
        "ytick.labelsize":       13,   # 10 + 3
        "legend.fontsize":       13,   # 10 + 3
        "legend.title_fontsize": 14,   # 11 + 3
    }
    with matplotlib.rc_context(_bump):
        # Wider canvas so the larger (+3 pt) multi-line x-axis labels don't collide.
        fig, axes = plt.subplots(1, 2, figsize=(16, 5.8))
        plt.subplots_adjust(left=0.06, right=0.80, top=0.90, bottom=0.22, wspace=0.24)
        panel_lambda(axes[0], df)
        panel_tail_fpr(axes[1], df)
        _panel_label(axes[0], "a", fontsize=18)   # 15 + 3
        _panel_label(axes[1], "b", fontsize=18)   # 15 + 3
        handles = _legend_handles()
        fig.legend(handles=handles, labels=[h.get_label() for h in handles],
                   loc="center left", bbox_to_anchor=(0.81, 0.5),
                   frameon=False, title="Method")
        for ext in ("png", "pdf"):
            fig.savefig(out_path.with_suffix(f".{ext}"))
        plt.close(fig)


def fig_fpr_levels(df, out_path: Path):
    fig, axes = plt.subplots(1, 4, figsize=(20, 5.5), sharey=True)
    plt.subplots_adjust(left=0.05, right=0.85, top=0.90, bottom=0.18, wspace=0.10)
    panel_fpr_alpha_grid(axes, df)
    for ax, lab in zip(axes, "abcd"):
        _panel_label(ax, lab, x=-0.05)
    handles = _legend_handles()
    fig.legend(handles=handles, labels=[h.get_label() for h in handles],
               loc="center left", bbox_to_anchor=(0.86, 0.5),
               frameon=False, title="Method")
    fig.suptitle(r"False-positive rate inflation across significance thresholds",
                 y=0.99, fontsize=14)
    for ext in ("png", "pdf"):
        fig.savefig(out_path.with_suffix(f".{ext}"))
    plt.close(fig)


def fig_power(df, out_path: Path):
    fig, ax = plt.subplots(figsize=(9, 5.5))
    plt.subplots_adjust(left=0.10, right=0.70, top=0.90, bottom=0.20)
    panel_power(ax, df)
    handles = _legend_handles()
    fig.legend(handles=handles, labels=[h.get_label() for h in handles],
               loc="center left", bbox_to_anchor=(0.71, 0.5),
               frameon=False, title="Method")
    for ext in ("png", "pdf"):
        fig.savefig(out_path.with_suffix(f".{ext}"))
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run-dir", default="${ADSP_ROOT}/outputs/power_runs_v5")
    ap.add_argument("--out-dir", default="${ADSP_ROOT}/outputs/power_pure_only")
    args = ap.parse_args()
    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    cache = out / "cache_metrics.parquet"

    df = _load(args.run_dir, cache)
    n_reps = df["replicate"].nunique()
    n_scen = df["scenario_key"].nunique()
    n_meth = df["method"].nunique()
    n_cfg  = df["config_key"].nunique()
    print(f"Loaded {n_reps} reps x {n_scen} scenarios x {n_meth} methods x {n_cfg} configs "
          f"= {len(df):,} rows", flush=True)
    print(f"  font.family = {FONT_FAMILY}", flush=True)

    main_figure(df,     out / "main_figure")
    print(f"  -> {out / 'main_figure.pdf'}", flush=True)
    fig_calibration(df, out / "fig_calibration")
    print(f"  -> {out / 'fig_calibration.pdf'}", flush=True)
    fig_fpr_levels(df,  out / "fig_fpr_levels")
    print(f"  -> {out / 'fig_fpr_levels.pdf'}", flush=True)
    fig_power(df,       out / "fig_power")
    print(f"  -> {out / 'fig_power.pdf'}", flush=True)


if __name__ == "__main__":
    main()
