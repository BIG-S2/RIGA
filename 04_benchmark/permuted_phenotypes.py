#!/usr/bin/env python
"""
permuted_phenotypes.py -- permuted-null phenotypes for the LARGE common-variant simulation (2026-09-23).

One replicate = the observed AD labels of the 27,306-subject plasmode pool permuted inside each cohort x project x
centre cell (cell_A), exactly as deliverable B's scheme A (plasmode_null_harness.draw_ystar: rng_y =
default_rng(SeedSequence(seed).spawn(2)[0]), one permutation per cell), with NEW seeds
(SeedSequence(20260923).spawn(N_TASKS), written once to <out-root>/seed_manifest_bignull.csv).

Writes <out-root>/ystar/npz/rep<k>_seed<seed>.npz with the keys benchmark_harness.py --ystar-from and
regenie_single_variant.py --native-root read: ystar (int8, pool row order), seed (uint64), replicate, scheme, cells.
Idempotent; never edits an existing file.

Usage: python permuted_phenotypes.py --replicate k [--out-root outputs/plasmode/bignull] [--n-tasks 1000]
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path("${ADSP_ROOT}")
sys.path.insert(0, str(ROOT / "revision" / "scripts"))
import benchmark_null_harness as H  # noqa: E402  (draw_ystar; read-only import)

ROOT_SEED = 20260923


def manifest(path: Path, n: int) -> pd.DataFrame:
    if path.exists():
        m = pd.read_csv(path)
        if len(m) < n:
            raise SystemExit(f"{path} has {len(m)} tasks < {n}")
        return m
    children = np.random.SeedSequence(ROOT_SEED).spawn(n)
    m = pd.DataFrame([dict(task=i + 1, task_type="bignull_A", spawn_key=i,
                           seed=int(c.generate_state(1, dtype=np.uint64)[0])) for i, c in enumerate(children)])
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(".tmp.csv")
    m.to_csv(tmp, index=False)
    tmp.rename(path)
    return m


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--replicate", type=int, required=True)
    ap.add_argument("--out-root", default=str(ROOT / "outputs" / "plasmode" / "bignull"))
    ap.add_argument("--pool-csv", default=str(H.OUT_ROOT / "plasmode_pool.csv"))
    ap.add_argument("--n-tasks", type=int, default=1000)
    args = ap.parse_args()
    out_root = Path(args.out_root)
    m = manifest(out_root / "seed_manifest_bignull.csv", args.n_tasks)
    row = m[m.task == args.replicate]
    if not len(row):
        raise SystemExit(f"replicate {args.replicate} not in the manifest")
    seed = int(row.seed.iloc[0])
    npz_dir = out_root / "ystar" / "npz"
    npz_dir.mkdir(parents=True, exist_ok=True)
    out = npz_dir / f"rep{args.replicate:03d}_seed{seed}.npz"
    if out.exists():
        print(f"{out.name} exists -> skip")
        return 0
    pool = pd.read_csv(args.pool_csv, low_memory=False)
    ystar, _ = H.draw_ystar(pool, "A", seed)
    ad = pool["AD"].to_numpy().astype(np.int8)
    cells = pool["cell_A"].astype(str)
    same = (pd.Series(ystar).groupby(cells.values).sum().sort_index() == pd.Series(ad).groupby(cells.values).sum().sort_index()).all()
    if not same:
        raise SystemExit("per-cell case counts not preserved (bug)")
    np.savez_compressed(out, ystar=ystar.astype(np.int8), seed=np.uint64(seed), replicate=np.int64(args.replicate),
                        scheme=np.array(["A"]), n_cells=np.int64(cells.nunique()), n_cases=np.int64(int(ystar.sum())))
    print(f"wrote {out.name}: {int(ystar.sum())} cases (real {int(ad.sum())}), {cells.nunique()} cells, "
          f"identity with real labels {float((ystar == ad).mean()):.4f}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
