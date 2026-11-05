#!/usr/bin/env python
"""supplementary_figure5_stratified_designs.py -- the two new stratified designs against the reference designs on the genome-wide coding masks
(2026-10-02): (a) coding-mask false-positive ratio vs threshold (permuted null, 1,000 replicates), (b) synonymous units,
(c) coding power at nominal thresholds (60 within-cell sparse-spike replicates).  Input: strata_designs_vs_reference.csv."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

T = "${ADSP_ROOT}/revision/tables/plasmode_power/strata_designs_vs_reference.csv"
OUT = "${ADSP_ROOT}/revision/figures/fig_strata_designs_vs_reference.png"
ARMS = [("saturated_cell", "Oracle: generating cells", "#4D4D4D", "o"),
        ("strat_project", "RIGA: project x center strata (+ PC1-20)", "#0072B2", "s"),
        ("strat_pcl", "S3 strata: project x center x LPA (+ PC1-20)", "#009E73", "D"),
        ("riga_projectkey", "Matched 1:1, project key", "#56B4E9", "^"),
        ("pooled_project", "Pooled: center + project + PC1-20", "#D55E00", "v"),
        ("strat_lpa", "S2 strata: LPA x sex x age band (+ PC1-20)", "#CC79A7", "P")]
d = pd.read_csv(T)
nul = d[d.kind == "null_fpr_ratio"]
pw = d[(d.kind == "power") & (d.level == -1)]
fig, axes = plt.subplots(1, 3, figsize=(15, 4.6))
for ax, mask, title in [(axes[0], "coding_maf001", "a  Coding units (26,768 genes): false-positive ratio"), (axes[1], "syn_maf001", "b  Synonymous units (18,544 genes): false-positive ratio")]:
    for arm, lab, col, mk in ARMS:
        t = nul[(nul.arm == arm) & (nul["mask"] == mask)].sort_values("alpha", ascending=False)
        t = t[t.alpha >= 1e-6]
        y = t.ratio.to_numpy(); lo = t.lo.to_numpy(); hi = t.hi.to_numpy()
        ax.errorbar(t.alpha, np.maximum(y, 0.05), yerr=[np.maximum(y - lo, 0), np.maximum(hi - y, 0)], color=col, marker=mk, ms=6, lw=1.6, capsize=2, label=lab)
    ax.axhline(1, color="0.5", lw=0.8, ls="--")
    ax.set_xscale("log"); ax.set_yscale("log"); ax.invert_xaxis()
    ax.set_xlabel("nominal threshold"); ax.set_ylabel("false positives / expected"); ax.set_title(title, loc="left", fontsize=10)
    ax.grid(True, which="major", color="0.9", lw=0.6)
ax = axes[2]
for arm, lab, col, mk in ARMS:
    t = pw[(pw.arm == arm) & (pw["mask"] == "coding_maf001")].sort_values("alpha", ascending=False)
    ax.errorbar(t.alpha, t.power, yerr=t.se, color=col, marker=mk, ms=6, lw=1.6, capsize=2, label=lab)
ax.set_xscale("log"); ax.invert_xaxis(); ax.set_ylim(0.3, 0.9)
ax.set_xlabel("nominal threshold"); ax.set_ylabel("power (spiked coding units, all effect levels)")
ax.set_title("c  Coding units: power (60 replicates)", loc="left", fontsize=10); ax.grid(True, color="0.9", lw=0.6)
axes[0].legend(fontsize=7.5, loc="upper left", frameon=False)
fig.suptitle("Stratified designs with simulated outcomes on observed ADSP genotypes\n(genome-wide MAF < 0.1% coding SNVs, 1,000 permutations)", fontsize=11)
fig.tight_layout()
fig.savefig(OUT, dpi=200)
print("wrote", OUT)
