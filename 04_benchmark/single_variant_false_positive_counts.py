#!/usr/bin/env python
"""single_variant_false_positive_counts.py (2026-10-13) -- the Figure 3a table (same layout as plasmode_bignull_fpr_table.py)
for the normal-approximation score test (p_chi2, the test the ADSP single-variant scan uses) after removing the benchmark
panel's heterozygote-excess sites (pooled heterozygote excess at P < 1e-4, the direction and threshold of QC_label7; the
pruned panel never had a Hardy-Weinberg filter). External tools (REGENIE/GMMAT/GENESIS npz) get the same site exclusion.
Writes revision/tables/plasmode_power/fpr_bignull_figure_table_hetfilt_chi2.csv (+ _by_replicate.csv)."""
import glob, re
from pathlib import Path
import numpy as np, pandas as pd
from scipy.stats import beta
ROOT = Path("${ADSP_ROOT}")
HARDY = ROOT / "revision/tables/plasmode_power/hardy/all.hardy"
OUT = ROOT / "revision/tables/plasmode_power/fpr_bignull_figure_table_hetfilt_chi2.csv"
ALPHAS = [1e-3, 1e-4, 1e-5, 5e-8]; TAG = {1e-3: "1e3", 1e-4: "1e4", 1e-5: "1e5", 5e-8: "5e8"}
h = pd.read_csv(HARDY, sep="\t"); bad = h[(h["O(HET_A1)"] > h["E(HET_A1)"]) & (h["P"] < 1e-4)]
badset = set(zip(bad["#CHROM"].astype(int), bad["ID"].str.split("_").str[1].astype(int)))
print("heterozygote-excess sites removed:", len(badset), "of", len(h))
def mask_of(z):
    chrpos = list(zip(z["chr"].astype(int), z["pos"].astype(int)))
    return np.array([cp in badset for cp in chrpos])
def cp(k, n, a):
    lo = beta.ppf(0.025, k, n - k + 1) if k > 0 else 0.0; hi = beta.ppf(0.975, k + 1, n - k)
    return k / (n * a), lo / a, hi / a
rows = []; hetmask = None
for f in sorted(glob.glob(str(ROOT / "outputs/plasmode/power/bignull/npz/rep*_seed*.npz"))):
    rep = int(re.search(r"rep(\d+)_", Path(f).name).group(1)); z = np.load(f)
    if hetmask is None: hetmask = mask_of(z)
    keep = ~z["excluded"].astype(bool) & ~hetmask
    for arm in [str(a) for a in z["arms"]]:
        k = f"p_chi2__{arm}"
        if k not in z.files: continue
        p = z[k].astype(np.float64)[keep]; p = p[np.isfinite(p) & (p > 0) & (p <= 1)]
        rows.append(dict(arm=arm, rep=rep, n=int(p.size), **{f"k{TAG[a]}": int((p < a).sum()) for a in ALPHAS}))
for f in sorted(glob.glob(str(ROOT / "outputs/plasmode/power/tools_bignull/npz/regenie_*_rep*.npz"))):
    rep = int(re.search(r"rep(\d+)", Path(f).name).group(1)); z = np.load(f, allow_pickle=True)
    covset = str(z["covset"][0]) if "covset" in z.files else re.search(r"regenie_(.+)_rep", Path(f).name).group(1)
    keep = ~z["excluded"].astype(bool) & ~mask_of(z)
    p = z["p"].astype(np.float64)[keep]; p = p[np.isfinite(p) & (p > 0) & (p <= 1)]
    rows.append(dict(arm=f"regenie_{covset}", rep=rep, n=int(p.size), **{f"k{TAG[a]}": int((p < a).sum()) for a in ALPHAS}))
d = pd.DataFrame(rows); d.to_csv(OUT.with_name(OUT.stem + "_by_replicate.csv"), index=False)
agg = d.groupby("arm").agg(n_units=("rep", "nunique"), n=("n", "sum"), **{f"k{TAG[a]}": (f"k{TAG[a]}", "sum") for a in ALPHAS}).reset_index()
tabs = []
for a in ALPHAS:
    t = agg.copy(); xs = [cp(int(k), int(n), a) for k, n in zip(t[f"k{TAG[a]}"], t.n)]
    t["fpr_1e3_eval_x"] = [x[0] for x in xs]; t["fpr_1e3_eval_x_lo"] = [x[1] for x in xs]; t["fpr_1e3_eval_x_hi"] = [x[2] for x in xs]
    t["n_ok"] = t.n_units; t["n_failed"] = 0; t["panel"] = f"single variants, {a:g}"; tabs.append(t)
res = pd.concat(tabs, ignore_index=True); res.to_csv(OUT, index=False)
print(res[res.panel.str.endswith(("1e-05", "5e-08"))][["panel", "arm", "n_units", "n", "fpr_1e3_eval_x", "fpr_1e3_eval_x_lo", "fpr_1e3_eval_x_hi"]].round(2).to_string(index=False))
