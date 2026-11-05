#!/usr/bin/env python
"""
stratified_designs_compare.py -- side-by-side calibration and power of the two NEW stratified designs (strat_lpa = LPA cluster x Sex x
age band strata; strat_pcl = project x center x LPA cluster strata) against the reference designs of the genome-wide burden runs
(2026-10-01).  Null error rates: the strata null run (1,000 permutations, two arms) and the main genome-wide null (same Y*, reference
arms) in the long-table format of plasmode_burden_bignull_fpr_table.py / plasmode_burden_bignull_tail.py; power: the strata power
run (plasmode_burden_power_simple.py), which scored the reference arms on the same replicates.
Usage: python stratified_designs_compare.py --null-strata <long> --null-main <long> --tail-strata <tail> --tail-main <tail> --power <csv> --out <csv>
"""
import argparse

import pandas as pd

ARMS = ["saturated_cell", "strat_project", "riga_projectkey", "pooled_project", "strat_lpa", "strat_pcl"]
LABEL = {"saturated_cell": "N0 oracle (cohort x project x center cells)", "strat_project": "S1 project x center strata + PCs",
         "riga_projectkey": "M3 RIGA project key", "pooled_project": "N5 pooled + center + project + PCs",
         "strat_lpa": "S2 LPA x sex x age-band strata + PCs", "strat_pcl": "S3 project x center x LPA strata + PCs"}
MASKS = ["coding_maf001", "syn_maf001"]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--null-strata", required=True)
    ap.add_argument("--null-main", required=True)
    ap.add_argument("--tail-strata", default="")
    ap.add_argument("--tail-main", default="")
    ap.add_argument("--power", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    nul = pd.concat([pd.read_csv(a.null_strata), pd.read_csv(a.null_main)], ignore_index=True)
    nul = nul[nul.arm.isin(ARMS) & nul["mask"].isin(MASKS)]
    nul = nul.drop_duplicates(subset=["arm", "mask", "alpha"], keep="first")
    if a.tail_strata and a.tail_main:
        tail = pd.concat([pd.read_csv(a.tail_strata), pd.read_csv(a.tail_main)], ignore_index=True)
        tail = tail[tail.arm.isin(ARMS) & tail["mask"].isin(MASKS) & tail.alpha.isin([2.5e-6, 1e-6])]
        tail = tail.drop_duplicates(subset=["arm", "mask", "alpha"], keep="first")
        nul = pd.concat([nul, tail[["arm", "mask", "alpha", "n_rep", "n_tests", "k", "ratio", "lo", "hi"]]], ignore_index=True)
    nul["kind"] = "null_fpr_ratio"
    pw = pd.read_csv(a.power)
    pw = pw[pw.arm.isin(ARMS) & pw["mask"].isin(MASKS)]
    pw["kind"] = "power"
    out = pd.concat([nul, pw], ignore_index=True, sort=False)
    out["design"] = out.arm.map(LABEL)
    out.to_csv(a.out, index=False)
    print("wrote", a.out)
    for mask in MASKS:
        t = nul[nul["mask"] == mask].pivot(index="arm", columns="alpha", values="ratio").reindex(ARMS).round(2)
        k = nul[nul["mask"] == mask].pivot(index="arm", columns="alpha", values="k").reindex(ARMS)
        print(f"\n== NULL {mask}: FPR ratio to nominal (k events)\n{pd.concat({'ratio': t, 'k': k}, axis=1).to_string()}")
        p = pw[(pw["mask"] == mask) & (pw.level == -1)].pivot(index="arm", columns="alpha", values="power").reindex(ARMS).round(3)
        print(f"\n== POWER {mask} (all levels, nominal thresholds)\n{p.to_string()}")


if __name__ == "__main__":
    main()
