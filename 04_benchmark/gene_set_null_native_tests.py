#!/usr/bin/env python
"""gene_set_null_native_tests.py (2026-10-13) -- the STAAR-family tests on the permuted gene-set null (rv_tools genome runs;
100 replicates for the stratified and pooled designs, 50 for the oracle): the primary SB (Cauchy combination of the SKAT and
burden p values = STAAR-O without ACAT-V), STAAR-O, burden and SKAT, with and without the cMAC >= 10 floor. For each,
the observed/expected ratio at 1e-3, 1e-4 and the gene-scan threshold 2.5e-6 with a percentile bootstrap over replicates,
and the share of null scans with any BH FDR-5% discovery (Clopper-Pearson). Read-only on inputs.
Writes revision/tables/plasmode_power/fpr_rv_native_tests.csv."""
import glob
import numpy as np, pandas as pd
from scipy.stats import beta
ROOT = "${ADSP_ROOT}/outputs/plasmode/rv_tools/genome/tools"
OUT = "${ADSP_ROOT}/revision/tables/plasmode_power/fpr_rv_native_tests.csv"
ALPHAS = [1e-3, 1e-4, 1e-5, 2.5e-6, 1e-6]
def cauchy(P):
    P = np.clip(np.asarray(P, float), 1e-300, 1 - 1e-16)
    t = np.where(P < 1e-15, 1.0 / (P * np.pi), np.tan((0.5 - P) * np.pi)).mean(axis=1)
    return np.where(t > 1e15, 1.0 / (t * np.pi), 0.5 - np.arctan(t) / np.pi)
def bh_n(p, q=0.05):
    p = np.sort(p[np.isfinite(p)]); m = len(p); ok = np.nonzero(p <= q * np.arange(1, m + 1) / m)[0]
    return 0 if len(ok) == 0 else ok[-1] + 1
rng = np.random.default_rng(20261013); rows = []
for design in ("saturated_cell", "strat_project", "pooled_project"):
    files = sorted(f for f in glob.glob(f"{ROOT}/{design}_rep*_genes.tsv.gz") if "_syn_" not in f)
    acc = {}
    for f in files:
        d = pd.read_csv(f, sep="\t"); d = d[np.isfinite(d.p_staar_O)]
        tests = {"SB": cauchy(d[["p_staar_S_1_25", "p_staar_B_1_25", "p_staar_S_1_1", "p_staar_B_1_1"]].to_numpy()), "STAAR-O": d.p_staar_O.to_numpy(),
                 "burden": cauchy(d[["p_staar_B_1_25", "p_staar_B_1_1"]].to_numpy()), "SKAT": cauchy(d[["p_staar_S_1_25", "p_staar_S_1_1"]].to_numpy())}
        for floor in (0, 10):
            keep = (d.cMAC >= floor).to_numpy()
            for t, p in tests.items():
                pk = p[keep]; a = acc.setdefault((t, floor), dict(n=[], any=[], k={al: [] for al in ALPHAS}))
                a["n"].append(pk.size); a["any"].append(int(bh_n(pk) > 0))
                for al in ALPHAS: a["k"][al].append(int(np.sum(pk < al)))
    for (t, floor), a in acc.items():
        n = np.asarray(a["n"], float); nrep = n.size; nany = int(np.sum(a["any"]))
        lo_s = beta.ppf(0.025, nany, nrep - nany + 1) if nany > 0 else 0.0; hi_s = beta.ppf(0.975, nany + 1, nrep - nany) if nany < nrep else 1.0
        for al in ALPHAS:
            k = np.asarray(a["k"][al], float); idx = rng.integers(0, nrep, size=(2000, nrep))
            boot = k[idx].sum(axis=1) / (n[idx].sum(axis=1) * al)
            rows.append(dict(design=design, test=t, cmac_floor=floor, n_rep=nrep, n_tests_mean=n.mean(), alpha=al, observed=int(k.sum()),
                             expected=float((n * al).sum()), ratio=float(k.sum() / (n * al).sum()), boot_lo=np.percentile(boot, 2.5), boot_hi=np.percentile(boot, 97.5),
                             share_any_fdr=nany / nrep, share_lo=lo_s, share_hi=hi_s))
res = pd.DataFrame(rows); res.to_csv(OUT, index=False)
print(res[(res.alpha == 2.5e-6)][["design", "test", "cmac_floor", "n_rep", "observed", "ratio", "boot_lo", "boot_hi", "share_any_fdr"]].round(3).to_string(index=False))
