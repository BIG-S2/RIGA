#!/usr/bin/env python
"""supplementary_table19_two_tests.py (2026-10-14): Supplementary Table 19 from single_variant_saddlepoint_retest.R (the FDR-5% single variants of
Supplementary Table 2 re-tested with the normal-approximation and the saddlepoint score tests on the deployed null). Prints the
numbers the caption needs (count below the full-scan boundary under each test, and the non-APOE variants' status)."""
import csv, math
from pathlib import Path
REV = Path("${ADSP_ROOT}/revision"); BND = 3.905e-7
r = list(csv.DictReader(open(REV / "tables/review/sv_fdr_spa_retest.csv")))
out = []
for x in r:
    pn = float(x["p_normal"]); ps = float(x["p_spa"]) if x["p_spa"] not in ("", "NA") else float("nan")
    out.append({"CHR": x["CHR"], "POS": x["POS"], "REF": x["REF"], "ALT": x["ALT"], "MAC": int(float(x["MAC"])), "MAF": x["MAF"], "P normal approximation": pn, "P saddlepoint": ps,
                "Ratio": ps / pn if pn > 0 else float("nan"), "Below full-scan boundary (saddlepoint)": "yes" if ps < BND else "no"})
out.sort(key=lambda z: (int(z["CHR"]), int(z["POS"])))
with open(REV / "tables/supp_v3/S19_sv_two_tests.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(out[0].keys())); w.writeheader(); w.writerows(out)
n = len(out); below = sum(z["Below full-scan boundary (saddlepoint)"] == "yes" for z in out)
non = [z for z in out if not (z["CHR"] == "19" and 44.4e6 <= float(z["POS"]) <= 46.4e6)]
print(f"S19: {n} variants; below boundary under the saddlepoint test: {below}; non-APOE: {len(non)}, of which below: {sum(z['Below full-scan boundary (saddlepoint)']=='yes' for z in non)}")
rat = sorted(z["Ratio"] for z in non if not math.isnan(z["Ratio"])); print("non-APOE ratio range: %.2f-%.2f" % (rat[0], rat[-1]))
for z in non: print(z["CHR"], z["POS"], z["MAC"], "%.3g -> %.3g" % (z["P normal approximation"], z["P saddlepoint"]), z["Below full-scan boundary (saddlepoint)"])
