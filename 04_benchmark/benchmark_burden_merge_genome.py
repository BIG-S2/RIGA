#!/usr/bin/env python
"""
benchmark_burden_merge_genome.py -- concatenate the per-chromosome MAF < 0.1% burden channels
(outputs/plasmode/burden_chr<k>_maf001/, k = 1..22; chr22 from plasmode_burden_chr22_maf001.R, the others from
benchmark_burden_test_by_chromosome.R) into ONE directory that the harness reads through --burden-dir:
  burden_<mask>.npy        units x 27,306 float32, pool row order, chromosomes 1..22 in order (written with open_memmap,
                           so peak memory is one chromosome's matrix at a time)
  gene_table_<mask>.csv    the per-chromosome tables with a `chr` column and a genome-wide unit_index; positions unchanged
  masks_manifest.csv, annotation_categories.csv (copied from chr22), burden_build.json (sources and per-chromosome counts)
Positions stay within-chromosome on purpose: the harness compares unit spans by position only, so a causal unit's 1-Mb
linked window also excludes same-coordinate units on other chromosomes (a few conservative exclusions in power runs,
none in null runs, which spike nothing).  Window units are already named chr<k>_win<idx>, so unit names stay unique.
Usage: python benchmark_burden_merge_genome.py [--out outputs/plasmode/burden_genome_maf001] [--chrs 1,2,...,22]
"""
from __future__ import annotations

import argparse
import json
import shutil
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path("${ADSP_ROOT}")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(ROOT / "outputs/plasmode/burden_genome_maf001"))
    ap.add_argument("--chrs", default=",".join(str(c) for c in range(1, 23)))
    args = ap.parse_args()
    chrs = [int(c) for c in args.chrs.split(",")]
    src = {c: ROOT / "outputs/plasmode" / f"burden_chr{c}_maf001" for c in chrs}
    missing = [c for c in chrs if not (src[c] / "masks_manifest.csv").exists()]
    if missing:
        raise SystemExit(f"chromosomes without a finished build: {missing}")
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    ref = pd.read_csv(src[22] / "masks_manifest.csv") if 22 in chrs else pd.read_csv(src[chrs[0]] / "masks_manifest.csv")
    masks = ref["mask"].tolist()
    rows, build = [], {"sources": {str(c): str(src[c]) for c in chrs}, "per_chromosome": {}}
    for m in masks:
        mats = {c: np.load(src[c] / f"burden_{m}.npy", mmap_mode="r") for c in chrs}
        tabs = {c: pd.read_csv(src[c] / f"gene_table_{m}.csv") for c in chrs}
        n_cols = {int(a.shape[1]) for a in mats.values()}
        if len(n_cols) != 1:
            raise SystemExit(f"{m}: sample counts differ across chromosomes: {n_cols}")
        for c in chrs:
            if mats[c].shape[0] != len(tabs[c]):
                raise SystemExit(f"{m} chr{c}: {mats[c].shape[0]} matrix rows vs {len(tabs[c])} table rows")
        n_units = sum(int(a.shape[0]) for a in mats.values())
        N = n_cols.pop()
        dst = out / f"burden_{m}.npy"
        M = np.lib.format.open_memmap(dst, mode="w+", dtype=np.float32, shape=(n_units, N))
        r0 = 0
        parts = []
        for c in chrs:
            a = mats[c]
            M[r0:r0 + a.shape[0]] = a
            t = tabs[c].copy()
            t.insert(1, "chr", c)
            t["unit_index_chr"] = t["unit_index"]
            t["unit_index"] = np.arange(r0, r0 + len(t))
            parts.append(t)
            build["per_chromosome"].setdefault(str(c), {})[m] = dict(n_units=int(len(t)), n_variants=int(t.n_variants.sum()))
            r0 += a.shape[0]
        M.flush()
        # spot check: first and last row of every chromosome block equal the source
        for c in chrs:
            t = parts[chrs.index(c)]
            i0, i1 = int(t.unit_index.iloc[0]), int(t.unit_index.iloc[-1])
            if not (np.array_equal(M[i0], mats[c][0]) and np.array_equal(M[i1], mats[c][-1])):
                raise SystemExit(f"{m} chr{c}: copied rows differ from the source")
        del M
        T = pd.concat(parts, ignore_index=True)
        T.to_csv(out / f"gene_table_{m}.csv", index=False)
        r = ref.set_index("mask").loc[m]
        rows.append(dict(mask=m, unit_type=r.unit_type, mac_lo=int(r.mac_lo), mac_hi=int(r.mac_hi), n_units=int(len(T)),
                         n_variants=int(T.n_variants.sum()), total_cMAC=float(T.cMAC.sum()),
                         npy=str(dst), gene_table=str(out / f"gene_table_{m}.csv")))
        print(f"{m:16s} {len(T):7d} units  {int(T.n_variants.sum()):10d} variants  {n_units * N * 4 / 1e9:5.2f} GB", flush=True)
    pd.DataFrame(rows).to_csv(out / "masks_manifest.csv", index=False)
    ann = src[22] / "annotation_categories.csv" if 22 in chrs else src[chrs[0]] / "annotation_categories.csv"
    if ann.exists():
        shutil.copy(ann, out / "annotation_categories.csv")
    (out / "burden_build.json").write_text(json.dumps(build, indent=1))
    print("wrote", out)


if __name__ == "__main__":
    main()
