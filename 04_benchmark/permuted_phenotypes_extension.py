#!/usr/bin/env python
"""
permuted_phenotypes_extension.py -- EXTENSION of the permuted-null phenotypes (2026-09-28): 9,000 more Y* for the genome-wide
burden null so that burden error rates can be read at 5e-8 (about 13 expected events per design on the coding mask with
10,000 permutations).  Same scheme as permuted_phenotypes.py (labels permuted inside every cohort x project x centre cell,
plasmode_null_harness.draw_ystar scheme A) with a NEW seed root, SeedSequence(20260928).spawn(N), written to
<out-root>/seed_manifest_bignull.csv in the same format; files <out-root>/ystar/npz/rep<k>_seed<seed>.npz, k = 1..N, with the
keys the harness reads.  The harness run numbers these replicates 1000 + k so they never collide with the first 1,000.
All replicates are written by ONE process (the pool is read once); existing files are skipped.
Usage: python permuted_phenotypes_extension.py [--out-root outputs/plasmode/bignull_ext] [--n-tasks 9000] [--start 1] [--end N]
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

SCRIPTS = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPTS))
import permuted_phenotypes as Y  # noqa: E402

Y.ROOT_SEED = 20260928   # new root: the first 1,000 permutations used 20260923


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-root", default=str(Y.ROOT / "outputs" / "plasmode" / "bignull_ext"))
    ap.add_argument("--pool-csv", default=str(Y.H.OUT_ROOT / "plasmode_pool.csv"))
    ap.add_argument("--n-tasks", type=int, default=9000)
    ap.add_argument("--start", type=int, default=1)
    ap.add_argument("--end", type=int, default=0)
    args = ap.parse_args()
    out_root = Path(args.out_root)
    m = Y.manifest(out_root / "seed_manifest_bignull.csv", args.n_tasks)
    npz_dir = out_root / "ystar" / "npz"
    npz_dir.mkdir(parents=True, exist_ok=True)
    pool = pd.read_csv(args.pool_csv, low_memory=False)
    ad = pool["AD"].to_numpy().astype(np.int8)
    cells = pool["cell_A"].astype(str)
    ad_by_cell = pd.Series(ad).groupby(cells.values).sum().sort_index()
    end = args.end or args.n_tasks
    n_new = 0
    for k in range(args.start, end + 1):
        seed = int(m.loc[m.task == k, "seed"].iloc[0])
        out = npz_dir / f"rep{k:03d}_seed{seed}.npz"
        if out.exists():
            continue
        ystar, _ = Y.H.draw_ystar(pool, "A", seed)
        same = (pd.Series(ystar).groupby(cells.values).sum().sort_index() == ad_by_cell).all()
        if not same:
            raise SystemExit(f"replicate {k}: per-cell case counts not preserved (bug)")
        np.savez_compressed(out, ystar=ystar.astype(np.int8), seed=np.uint64(seed), replicate=np.int64(k),
                            scheme=np.array(["A"]), n_cells=np.int64(cells.nunique()), n_cases=np.int64(int(ystar.sum())))
        n_new += 1
        if n_new % 500 == 0:
            print(f"  {k}: {n_new} written", flush=True)
    print(f"done: {n_new} new files, {len(list(npz_dir.glob('rep*.npz')))} present in {npz_dir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
