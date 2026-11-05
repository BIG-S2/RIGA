#!/usr/bin/env python
"""null_model_sensitivity_summary.py (2026-10-13, review points 4 and 7): robustness of every reported signal to the null-model
specification, from the second-pass re-tests (targeted_retest_sb_pass.R → tables/review/targeted_sb_<tag>.csv): the full
model, four named ablations (nocls: no latent-class terms; pc20: PC1-20; noAPOEsel: APOE-extreme projects excluded;
noSelected: APOE-extreme and early-onset projects excluded) and leave-one-project-out over every project. Writes
revision/tables/supp_v3/S18_robustness_specifications.csv (one row per signal x named specification, plus the leave-one-out
range) and prints the summary."""
import glob, re
from pathlib import Path
import numpy as np, pandas as pd
REV = Path("${ADSP_ROOT}/revision"); R = REV / "tables/review"
BOUND = {"coding": 2.16e-6, "noncoding": 4.92e-6, "ncRNA": 6.49e-6, "sv": 3.88e-7}   # FDR critical values of the full scans
rows = []
for f in sorted(glob.glob(str(R / "targeted_sb_*.csv"))):
    d = pd.read_csv(f); d["tag"] = re.sub(r"^targeted_sb_(.*)\.csv$", r"\1", Path(f).name); rows.append(d[d.kind == "gene set"])
svrows = []   # 2026-10-14: single variants of the QC_label7 rescan from targeted_sv_pass.R (replaces the six-variant rows of the first pass)
for f in sorted(glob.glob(str(R / "targeted_sv_*.csv"))):
    d = pd.read_csv(f); d["tag"] = re.sub(r"^targeted_sv_(.*)\.csv$", r"\1", Path(f).name); d["kind"] = "single variant"; svrows.append(d)
D = pd.concat(rows + svrows, ignore_index=True)
D["p"] = np.where(D.kind == "gene set", D.SB, D.pvalue); D["fam"] = np.where(D.kind == "gene set", D.family, "sv")
D["signal"] = np.where(D.kind == "gene set", D.gene + " " + D.cat.fillna(""), D.gene + " chr" + D.chr.astype("Int64").astype(str) + ":" + D.pos.astype("Int64").astype(str))
D["below_boundary"] = D.p < D.fam.map(BOUND)
named = ["full", "nocls", "pc20", "noAPOEsel", "noSelected"]
out = []
for sig, g in D.groupby("signal", sort=False):
    fam = g.fam.iloc[0]; rec = {"Signal": sig, "Family": fam, "FDR boundary": BOUND[fam]}
    for t in named:
        v = g[g.tag == t].p; rec[t] = float(v.iloc[0]) if len(v) else np.nan
    loo = g[g.tag.str.startswith("loo_")]
    rec["LOO n"] = int(loo.p.notna().sum()); rec["LOO min p"] = float(loo.p.min()); rec["LOO max p"] = float(loo.p.max()); rec["LOO n below boundary"] = int(loo.below_boundary.sum())
    worst = loo.loc[loo.p.idxmax()] if loo.p.notna().any() else None; rec["LOO worst tag"] = worst.tag if worst is not None else ""
    out.append(rec)
S = pd.DataFrame(out); S.to_csv(REV / "tables/supp_v3/S18_robustness_specifications.csv", index=False)
pd.set_option("display.width", 250); print(S.round(10).to_string(index=False))
