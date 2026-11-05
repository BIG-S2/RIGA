#!/usr/bin/env python
"""
gmmat_genesis_single_variant_driver.py -- run GMMAT and GENESIS (gmmat_genesis_single_variant.R) on one replicate of the large
common-variant simulation and store the p-values in the npz layout the aggregator already reads for REGENIE
(plasmode_power_aggregate.load_null_pools / load_power: <out-root>/npz/regenie_<covset>_rep<k>.npz with p, excluded,
covset), under the covset names gmmat_project and genesis_project (arm names regenie_gmmat_project /
regenie_genesis_project once REG_SETS is patched by power_aggregate.py; the test label "firth" that the
loader attaches is a bookkeeping string: GMMAT is a plain score test, GENESIS a saddlepoint-corrected score test).

Steps: 1. find <native-root>/npz/rep<k>_seed*.npz (a bignull Y* file or a power-harness replicate), read ystar (+ keep
          for ascertainment replicates: dropped subjects get Y = NA);
       2. write the IID,Y csv, run the R driver (--threads);
       3. convert <prefix>_pvals.tsv.gz to the two npz files (p float32 in .bim order, tested, excluded from
          plasmode_null_harness.snp_annotation, seed, replicate, ystar, keep) + one row per tool in
          <out-root>/units_plasmode_tools.csv (n_tested, lambda, FPR x nominal at 1e-3 / 1e-5 on the non-excluded SNPs).
Idempotent per (replicate, tool).  Never edits an existing file.
"""
from __future__ import annotations

import argparse
import glob
import json
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import chi2

ROOT = Path("${ADSP_ROOT}")
SCRIPTS = ROOT / "revision" / "scripts"
sys.path.insert(0, str(SCRIPTS))
import benchmark_null_harness as H  # noqa: E402

TOOLS = {"gmmat": ("p_gmmat", "gmmat_project"), "genesis": ("p_genesis_spa", "genesis_project")}
T0 = time.time()


def log(m: str) -> None:
    print(f"[tools {time.time() - T0:7.1f}s] {m}", flush=True)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--replicate", type=int, required=True)
    ap.add_argument("--native-root", required=True, help="dir with npz/rep<k>_seed*.npz carrying ystar")
    ap.add_argument("--out-root", required=True)
    ap.add_argument("--which", default="gmmat,genesis")
    ap.add_argument("--threads", type=int, default=int(os.environ.get("SLURM_CPUS_PER_TASK", "8")))
    ap.add_argument("--scratch", default=os.environ.get("SIM_SCRATCH", "/dev/shm/plasmode_tools"))
    ap.add_argument("--pool-csv", default=str(H.OUT_ROOT / "plasmode_pool.csv"))
    ap.add_argument("--bfile", default=str(H.OUT_ROOT / "pool_281k"))
    ap.add_argument("--persnp", default=str(H.PERSNP))
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()
    rep = int(args.replicate)
    which = [w for w in args.which.split(",") if w]
    out_root = Path(args.out_root)
    npz_dir = out_root / "npz"
    npz_dir.mkdir(parents=True, exist_ok=True)
    todo = [w for w in which if args.force or not (npz_dir / f"regenie_{TOOLS[w][1]}_rep{rep:03d}.npz").exists()]
    if not todo:
        log(f"replicate {rep}: {which} already done -> skip")
        return 0
    hits = sorted(glob.glob(str(Path(args.native_root) / "npz" / f"rep{rep:03d}_seed*.npz")))
    if not hits:
        raise SystemExit(f"no npz for replicate {rep} under {args.native_root}/npz")
    z = np.load(hits[-1], allow_pickle=False)
    pool = pd.read_csv(args.pool_csv, low_memory=False, usecols=["IID"])
    N = len(pool)
    ystar = z["ystar"].astype(np.int8)
    keep = z["keep"].astype(bool) if "keep" in z.files else np.ones(N, dtype=bool)
    if ystar.size == int(keep.sum()) and ystar.size != N:      # ascertainment layout: ystar stored for the kept subjects only
        full = np.zeros(N, dtype=np.int8)
        full[keep] = ystar
        ystar = full
    if ystar.size != N:
        raise SystemExit(f"ystar {ystar.size} != pool {N}")
    seed = int(z["seed"]) if "seed" in z.files else 0
    scratch = Path(args.scratch)
    scratch.mkdir(parents=True, exist_ok=True)
    ycsv = scratch / f"ystar_rep{rep:03d}.csv"
    pd.DataFrame({"IID": pool.IID.astype(str), "Y": np.where(keep, ystar.astype(float), np.nan)}).to_csv(ycsv, index=False, na_rep="NA")
    log(f"replicate {rep} seed {seed}: {int(ystar[keep].sum())} cases / {int(keep.sum())} analyzed ({int((~keep).sum())} dropped); tools {todo}")
    prefix = scratch / f"tools_rep{rep:03d}"
    cmd = ["Rscript", str(SCRIPTS / "gmmat_genesis_single_variant.R"), "--ystar", str(ycsv), "--out", str(prefix),
           "--which", ",".join(todo), "--threads", str(args.threads)]
    log("Rscript " + " ".join(cmd[2:]))
    rc = subprocess.run(cmd).returncode
    if rc != 0:
        raise SystemExit(f"R driver exit {rc}")
    pv = pd.read_csv(f"{prefix}_pvals.tsv.gz", sep="\t", dtype={"snp": str})
    meta = json.loads(Path(f"{prefix}_meta.json").read_text())
    bim = pd.read_csv(str(args.bfile) + ".bim", sep=r"\s+", header=None, names=["chr", "id", "cm", "pos", "a1", "a2"], dtype={1: str})
    if len(pv) != len(bim) or not np.array_equal(pv.snp.astype(str).values, bim.id.astype(str).values):
        raise SystemExit("p-value table is not in .bim order")
    ann = H.snp_annotation(bim, Path(args.persnp) if args.persnp else None)
    excluded = ann.excluded.to_numpy()
    rows = []
    failed = []
    for w in todo:
        col, covset = TOOLS[w]
        p = pv[col].to_numpy(dtype=np.float64)
        tested = np.isfinite(p)
        if not tested.any():
            log(f"{w}: no p-values (the R fit failed: {meta.get(w, {}).get('error', 'see the R log')}); nothing stored")
            failed.append(w)
            continue
        pp = p[tested & ~excluded]
        n = int(pp.size)
        lam = float(np.median(chi2.isf(pp, 1)) / chi2.isf(0.5, 1)) if n else np.nan
        row = dict(tool=w, covset=covset, replicate=rep, seed=seed, n_tested=int(tested.sum()), n=n, lambda_gc=lam,
                   **{f"n_lt_{t}": int((pp < a).sum()) for a, t in ((1e-3, "1e3"), (1e-4, "1e4"), (1e-5, "1e5"), (1e-6, "1e6"), (5e-8, "5e8"))},
                   converged=meta.get(w, {}).get("converged"), seconds=meta.get(w, {}).get("null_fit_s", np.nan))
        for a, t in ((1e-3, "1e3"), (1e-5, "1e5")):
            row[f"fpr_{t}_x"] = row[f"n_lt_{t}"] / n / a if n else np.nan
        rows.append(row)
        np.savez_compressed(npz_dir / f"regenie_{covset}_rep{rep:03d}.npz", p=p.astype(np.float32), tested=tested, excluded=excluded,
                            covset=np.array([covset]), tool=np.array([w]), seed=np.uint64(seed), replicate=np.int64(rep),
                            ystar=ystar.astype(np.int8), keep=keep, chr=ann.chr.to_numpy().astype(np.int16), pos=ann.pos.to_numpy().astype(np.int64))
        log(f"{w}: tested {int(tested.sum())}, lambda {lam:.3f}, FPR x nominal 1e-3 {row['fpr_1e3_x']:.2f}, 1e-5 {row['fpr_1e5_x']:.1f}")
    if rows:
        H._append_csv(out_root / "units_plasmode_tools.csv", rows)
    (out_root / "meta").mkdir(exist_ok=True)
    shutil.move(f"{prefix}_meta.json", out_root / "meta" / f"tools_rep{rep:03d}_meta.json")   # scratch is /dev/shm: no cross-device rename
    log(f"done in {time.time() - T0:.0f} s" + (f"; FAILED tools: {failed}" if failed else ""))
    return 3 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
