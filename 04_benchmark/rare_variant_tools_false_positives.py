#!/usr/bin/env python
"""
rare_variant_tools_false_positives.py -- pooled false-positive counts of the rare-variant TOOLS on the plasmode permuted null (2026-10-02).
Reads every COMPLETED per-replicate gene table under <root>/{tools,regenie,saige}/<covset>_rep<k>_genes.tsv.gz (a table counts only
when its sibling <stem>_meta.json exists: the R driver checkpoints partial tables), one row per gene (x mask for REGENIE / SAIGE),
p-value columns p_<tool>_<test>.  --harness adds the harness's own burden test (saddlepoint score on the allele-count collapse, key
burden_spa__coding_maf001__<arm> in the genome null npz files) for the same arms and replicates on (a) all tool genes and (b) the
REGENIE-tested genes (REGENIE ignores masks with MAC < 5), given as `npzdir=gene_table.csv` pairs, one per genome part.
Output: one row per covset x tool_test x mask x alpha: n_rep, n_tests, k, ratio (= k / n_tests / alpha), exact Clopper-Pearson 95% CI.
Usage: python rare_variant_tools_false_positives.py --root outputs/plasmode/rv_tools/genome --out <csv> [--harness d1/npz=gt1.csv,d2/npz=gt2.csv] [--max-rep N]
"""
from __future__ import annotations

import argparse
import glob
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import beta

ALPHAS = [1e-3, 1e-4, 1e-5, 2.5e-6, 1e-6]
HARNESS_ARMS = ["pooled_project", "strat_project", "saturated_cell"]


def cp(k: int, n: int, a: float):
    lo = beta.ppf(0.025, k, n - k + 1) if k > 0 else 0.0
    hi = beta.ppf(0.975, k + 1, n - k) if k < n else 1.0
    return k / n / a, lo / a, hi / a


class Tally:
    def __init__(self):
        self.k, self.n, self.reps = {}, {}, {}

    def add(self, key, p: np.ndarray, rep: int):
        ok = np.isfinite(p) & (p >= 0) & (p <= 1)
        self.n[key] = self.n.get(key, 0) + int(ok.sum())
        self.reps.setdefault(key, set()).add(rep)
        kk = self.k.setdefault(key, np.zeros(len(ALPHAS), dtype=np.int64))
        for j, al in enumerate(ALPHAS):
            kk[j] += int((p[ok] < al).sum())


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--harness", default="", help="comma list of npzdir=gene_table.csv (harness null runs, one per genome part)")
    ap.add_argument("--max-rep", type=int, default=0)
    ap.add_argument("--min-cmac", type=float, default=0, help="keep only genes with harness cMAC >= this (tool tables by gene name, harness by unit)")
    ap.add_argument("--gene-table", default="${ADSP_ROOT}/outputs/plasmode/burden_genome_maf001/gene_table_coding_maf001.csv")
    a = ap.parse_args()
    cmac = pd.read_csv(a.gene_table).drop_duplicates("unit").set_index("unit")["cMAC"] if a.min_cmac > 0 else None
    files = [f for f in sorted(glob.glob(str(Path(a.root) / "*" / "*_rep*_genes.tsv.gz"))) if Path(f.replace("_genes.tsv.gz", "_meta.json")).exists()]
    T = Tally()
    genes_all, genes_regenie = set(), set()
    for f in files:
        m = re.search(r"/([^/]+)/([^/]+?)_rep(\d+)_genes\.tsv\.gz$", f)   # family = the directory right above the file
        family, covset, rep = m.group(1), m.group(2), int(m.group(3))
        if a.max_rep and rep > a.max_rep:
            continue
        d = pd.read_csv(f, sep="\t")
        if "mask" not in d.columns:
            d["mask"] = "coding"
        if cmac is not None:
            d = d[d["gene"].astype(str).map(cmac).fillna(0).to_numpy() >= a.min_cmac]
        genes_all.update(d["gene"].astype(str))
        if family == "regenie":
            genes_regenie.update(d.loc[d["mask"] == "coding", "gene"].astype(str))
        for col in [c for c in d.columns if c.startswith("p_")]:
            for mask, g in d.groupby("mask"):
                T.add((covset, col[2:], str(mask)), pd.to_numeric(g[col], errors="coerce").to_numpy(dtype=float), rep)
    print(f"{len(files)} completed gene tables, {len(genes_all):,} genes seen, {len(genes_regenie):,} REGENIE-tested coding genes", flush=True)
    if a.harness:
        allreps = set().union(*T.reps.values()) if T.reps else set()
        for pair in a.harness.split(","):
            npzdir, gtf = pair.split("=")
            gt = pd.read_csv(gtf)
            names = gt["unit"].astype(str)
            big = (gt["cMAC"].to_numpy() >= a.min_cmac) if cmac is not None else np.ones(len(gt), dtype=bool)
            subsets = {"harness_burden_spa": names.isin(genes_all).to_numpy() & big}
            if genes_regenie:
                subsets["harness_burden_spa_regeniegenes"] = names.isin(genes_regenie).to_numpy() & big
            for f in sorted(glob.glob(str(Path(npzdir) / "*.npz"))):
                rep = int(re.search(r"rep(\d+)", Path(f).name).group(1))
                if rep not in allreps:
                    continue
                z = np.load(f, allow_pickle=False)
                for arm in HARNESS_ARMS:
                    key = f"burden_spa__coding_maf001__{arm}"
                    if key not in z.files:
                        continue
                    p = z[key].astype(float)
                    if len(p) != len(gt):
                        raise SystemExit(f"{f}: {len(p)} units != gene table {len(gt)} ({gtf})")
                    for lab, keep in subsets.items():
                        T.add((arm, lab, "coding"), p[keep], rep)
    rows = []
    for key, kk in sorted(T.k.items()):
        for j, al in enumerate(ALPHAS):
            r, lo, hi = cp(int(kk[j]), T.n[key], al)
            rows.append(dict(covset=key[0], tool_test=key[1], mask=key[2], alpha=al, n_rep=len(T.reps[key]), n_tests=T.n[key], k=int(kk[j]), ratio=r, lo=lo, hi=hi))
    out = pd.DataFrame(rows)
    out.to_csv(a.out, index=False)
    print("wrote", a.out)
    for (covset, mask), t in out.groupby(["covset", "mask"]):
        print(f"\n== {covset} / {mask}: ratio to nominal (k)")
        piv = t.pivot(index="tool_test", columns="alpha", values="ratio").round(2)
        kp = t.pivot(index="tool_test", columns="alpha", values="k")
        nn = t.groupby("tool_test").agg(n_rep=("n_rep", "max"), n_tests=("n_tests", "max"))
        nn.columns = pd.MultiIndex.from_tuples([("n", "rep"), ("n", "tests")])
        print(pd.concat({"ratio": piv, "k": kp}, axis=1).join(nn).to_string())


if __name__ == "__main__":
    main()
