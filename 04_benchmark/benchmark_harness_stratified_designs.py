#!/usr/bin/env python
"""
benchmark_harness_stratified_designs.py -- benchmark_harness_rare_burden.py plus two stratified designs (2026-10-01, user request):
  strat_lpa   S2  strata = LPA ancestry cluster x Sex x age band (284 mixed strata, 27,117 subjects); Age + PC1-20 + strata
                  (Sex is absorbed by the strata)
  strat_pcl   S3  strata = sequencing project x center x LPA cluster (208 mixed strata, 20,361 subjects); Age + Sex + PC1-20 + strata
Both keep every subject of a stratum that contains cases and controls, exactly as strat_project (S1) does for project x center.
The harness files are untouched: the arm registry and build_frames are extended in memory, everything else passes through.
"""
from __future__ import annotations

import sys
from pathlib import Path

SCRIPTS = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPTS))
import benchmark_harness_rare_burden as W  # noqa: E402
import benchmark_harness as P  # noqa: E402

NEW = {"strat_lpa": ("S2", "strata = LPA cluster x Sex x age band + Age + PC1-20, mixed strata"),
       "strat_pcl": ("S3", "strata = project x center x LPA cluster + Age + Sex + PC1-20, mixed strata")}
P.ARM_ORDER = list(P.ARM_ORDER) + [a for a in NEW if a not in P.ARM_ORDER]
P.ARM_LABEL = dict(P.ARM_LABEL, **{a: v[0] for a, v in NEW.items()})
if hasattr(P, "ARM_DESC"):
    P.ARM_DESC = dict(P.ARM_DESC, **{a: v[1] for a, v in NEW.items()})
_orig_build_frames = P.build_frames


def build_frames(pool, y, matched, py_sets, arms):
    frames, clr = _orig_build_frames(pool, y, matched, py_sets, [a for a in arms if a not in NEW])
    if any(a in arms for a in NEW):
        base = P.base_frame(pool, y)
        pcs20 = [c for c in P.H.PCS20 if c in base.columns]
        if "strat_lpa" in arms:
            lab = (pool["lpa_cluster"].astype(str) + "|" + pool["Sex"].astype(str) + "|" + pool["age_band"].astype(str)).to_numpy()
            rows = P.mixed_rows(lab, y)
            d = base.iloc[rows].reset_index(drop=True)
            d["strat_lpa"] = lab[rows]
            frames["strat_lpa"] = (d, rows, ["Age"] + pcs20 + ["strat_lpa"], None)
        if "strat_pcl" in arms:
            lab = (pool["Study"].astype(str) + "|" + pool["center"].astype(str) + "|" + pool["lpa_cluster"].astype(str)).to_numpy()
            rows = P.mixed_rows(lab, y)
            d = base.iloc[rows].reset_index(drop=True)
            d["strat_pcl"] = lab[rows]
            frames["strat_pcl"] = (d, rows, ["Age", "Sex"] + pcs20 + ["strat_pcl"], None)
    return frames, clr


P.build_frames = build_frames

if __name__ == "__main__":
    print(f"[strata wrapper] arms added: {list(NEW)}", flush=True)
    sys.exit(W.main())
