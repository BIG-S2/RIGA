#!/usr/bin/env python
"""supplementary_figure3_rare_variant_tools.py (2026-10-13; the SKAT + burden omnibus of this study, from fpr_rv_native_tests.csv, replaces the
benchmark's saddlepoint burden row) -- rare-variant TOOL comparison on the plasmode permuted null (2026-10-03): false-positive ratio vs threshold for
eight tests under three designs (oracle cells, project x center strata, pooled covariates), on all coding genes (top row) and on
genes with cMAC >= 10 (bottom row).  Inputs: the long tables of rare_variant_tools_false_positives.py (all genes, cMAC >= 10)."""
import argparse
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

TESTS = [("SB", "SKAT + burden omnibus (this study)", "#000000", "o"),
         ("saige_burden", "SAIGE-GENE+ burden", "#56B4E9", "s"),
         ("saige_skato", "SAIGE-GENE+ SKAT-O", "#0072B2", "D"),
         ("regenie_burden", "REGENIE Firth burden", "#009E73", "^"),
         ("regenie_skato", "REGENIE SKAT-O", "#E69F00", "v"),
         ("staar_B_1_25", "STAAR burden (1,25)", "#CC79A7", "P"),
         ("staar_O", "STAAR-O (no annotation)", "#D55E00", "X"),
         ("smmat_O", "SMMAT SKAT-O", "#999999", "h")]
DESIGNS = [("saturated_cell", "Oracle: generating cells"), ("strat_project", "RIGA: project x center strata + PC1-20"),
           ("pooled_project", "Pooled: center + project + PC1-20")]
ap = argparse.ArgumentParser()
ap.add_argument("--all-table", required=True)
ap.add_argument("--cmac-table", required=True)
ap.add_argument("--out", required=True)
a = ap.parse_args()
nat = pd.read_csv("${ADSP_ROOT}/revision/tables/plasmode_power/fpr_rv_native_tests.csv")
nat = nat[nat.test == "SB"].rename(columns={"design": "covset", "test": "tool_test", "boot_lo": "lo", "boot_hi": "hi"}); nat["mask"] = "coding"
def with_sb(t, floor):
    return pd.concat([t[t.tool_test != "harness_burden_spa"], nat[nat.cmac_floor == floor][["covset", "tool_test", "mask", "alpha", "n_rep", "ratio", "lo", "hi"]]], ignore_index=True)
tabs = [("all coding genes (26,768)", with_sb(pd.read_csv(a.all_table), 0)), ("genes with cMAC >= 10 (20,361)", with_sb(pd.read_csv(a.cmac_table), 10))]
fig, axes = plt.subplots(2, 3, figsize=(15.5, 8.6), sharey="row")
for r, (rowlab, t) in enumerate(tabs):
    t = t[(t["mask"] == "coding") & (t.alpha >= 1e-6)]
    for c, (arm, dlab) in enumerate(DESIGNS):
        ax = axes[r, c]
        reps = []
        for key, lab, col, mk in TESTS:
            u = t[(t.covset == arm) & (t.tool_test == key)].sort_values("alpha", ascending=False)
            if u.empty:
                continue
            reps.append(f"{lab.split(' ')[0]} {int(u.n_rep.max())}")
            y = u.ratio.to_numpy(); lo = u.lo.to_numpy(); hi = u.hi.to_numpy()
            ax.errorbar(u.alpha, np.maximum(y, 0.05), yerr=[np.maximum(y - lo, 0), np.maximum(hi - y, 0)], color=col, marker=mk, ms=5.5, lw=1.4, capsize=2, label=lab)
        ax.axhline(1, color="0.5", lw=0.8, ls="--")
        ax.set_xscale("log"); ax.set_yscale("log"); ax.invert_xaxis()
        ax.grid(True, which="major", color="0.9", lw=0.6)
        ax.set_title(f"{'abc'[c] if r == 0 else 'def'[c]}  {dlab}" if r == 0 else f"{'def'[c]}  {dlab}", loc="left", fontsize=9.5)
        if c == 0:
            ax.set_ylabel(f"false positives / expected\n{rowlab}")
        if r == 1:
            ax.set_xlabel("nominal threshold")
        ax.text(0.02, 0.02, "replicates: " + ", ".join(sorted(set(reps))), transform=ax.transAxes, fontsize=6.5, color="0.35", va="bottom")
# legend from the panel that carries every test (the oracle panel has no SAIGE run)
hl = {}
for ax in axes.ravel():
    for h, l in zip(*ax.get_legend_handles_labels()):
        hl.setdefault(l, h)
order = [lab for _, lab, _, _ in TESTS if lab in hl]
axes[0, 0].legend([hl[l] for l in order], order, fontsize=7.5, loc="upper left", frameon=False)
fig.suptitle("Rare-variant tests on the plasmode permuted null: the design decides, the test sets the floor (genome-wide MAF < 0.1% coding SNVs; normal approximation throughout)", fontsize=11)
fig.tight_layout()
fig.savefig(a.out, dpi=200)
print("wrote", a.out)
