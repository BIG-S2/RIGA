#!/usr/bin/env python
"""single_variant_false_positive_intervals.py (2026-10-13, review point 8) -- replicate-based uncertainty for the single-variant
false-positive ratios of the 1,000-permutation big null. Within a replicate the 281,156 tests are dependent (LD), so the
Clopper-Pearson interval over pooled tests is too narrow; the replicate is the independent unit. For every design and both
p-value types (saddlepoint and the normal-approximation score test that the deployed single-variant scan uses) this writes
the observed/expected ratio with (i) the pooled Clopper-Pearson interval and (ii) a percentile bootstrap interval over
replicates (2,000 resamples). SNPs flagged `excluded` (APOE region and windows around known AD loci) are not counted, as in
the published tables. Writes revision/tables/plasmode_power/bignull_replicate_ci.csv."""
import glob
import numpy as np
import pandas as pd
from scipy.stats import beta

ROOT = "${ADSP_ROOT}/outputs/plasmode/power/bignull/npz"
OUT = "${ADSP_ROOT}/revision/tables/plasmode_power/bignull_replicate_ci.csv"
ARMS = ["saturated_cell", "strat_project", "riga_projectkey", "pooled_project", "pooled_pcs", "riga_pair33"]
ALPHAS = [1e-3, 1e-5, 5e-8]
files = sorted(glob.glob(f"{ROOT}/rep*.npz"))
cnt = {(a, t, al): [] for a in ARMS for t in ("spa", "chi2") for al in ALPHAS}
ntest = {(a, t): [] for a in ARMS for t in ("spa", "chi2")}
for f in files:
    z = np.load(f)
    keep = ~z["excluded"].astype(bool)
    for a in ARMS:
        for t in ("spa", "chi2"):
            k = f"p_{t}__{a}"
            if k not in z.files:
                continue
            p = z[k][keep]; p = p[np.isfinite(p)]
            ntest[(a, t)].append(p.size)
            for al in ALPHAS:
                cnt[(a, t, al)].append(int(np.sum(p < al)))
rng = np.random.default_rng(20261013)
rows = []
for a in ARMS:
    for t in ("spa", "chi2"):
        n = np.asarray(ntest[(a, t)], float)
        if n.size == 0:
            continue
        for al in ALPHAS:
            c = np.asarray(cnt[(a, t, al)], float)
            K, E = c.sum(), (n * al).sum()
            lo = beta.ppf(0.025, K, n.sum() - K + 1) / al if K > 0 else 0.0
            hi = beta.ppf(0.975, K + 1, n.sum() - K) / al
            idx = rng.integers(0, c.size, size=(2000, c.size))
            boot = c[idx].sum(axis=1) / (n[idx].sum(axis=1) * al)
            rows.append(dict(design=a, test=t, alpha=al, n_rep=c.size, observed=int(K), expected=E, ratio=K / E,
                             cp_lo=lo, cp_hi=hi, boot_lo=np.percentile(boot, 2.5), boot_hi=np.percentile(boot, 97.5)))
res = pd.DataFrame(rows)
res.to_csv(OUT, index=False)
print(res.round(3).to_string(index=False))
