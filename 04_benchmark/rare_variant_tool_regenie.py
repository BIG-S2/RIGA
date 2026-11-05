#!/usr/bin/env python
"""
rare_variant_tool_regenie.py -- REGENIE v3.2.9 GENE-BASED tests (burden Firth + SKAT / SKAT-O / ACAT-V / ACAT-O) on the plasmode
permuted null (2026-10-01), on the SAME gene x variant sets as the harness's coding_maf001 mask (plasmode_rv_groups.py:
regenie.setlist / regenie.anno / regenie.maskdef; masks coding = syn+nonsyn, nonsyn, syn) and the SAME Y* of one replicate
(outputs/plasmode/power/bignull/npz/rep<k>_seed*.npz), under one covariate set:
  pooled_project  N5  Age + Sex + PC1-20 + centre + Study_cov (all subjects); step 1 REUSED from the common-variant tool run
                      (outputs/plasmode/power/tools_bignull/logs/pooled_project_rep<k>/step1_pred.list: same Y*, same covariates)
  strat_project   S1  Age + Sex + PC1-20 + project x centre cell (categorical), subjects of mixed cells only; step 1 fitted here
  saturated_cell  N0  Age + Sex + cohort x project x centre cell (categorical), mixed cells only; step 1 fitted here
Step 2: --bt --firth --approx --pThresh 0.05, --build-mask sum (the harness collapses by allele COUNT), --aaf-bins 0.01 and
--vc-maxAAF 0.01 (every exported variant has pool MAF < 0.1%; the bin only has to contain them all), --vc-tests skato,acato.
Outputs: <out-root>/regenie/<covset>_rep<k>_genes.tsv.gz (gene, mask, p_regenie_burden / _skat / _skato / _acatv / _acato),
         <out-root>/regenie/<covset>_rep<k>_meta.json, logs under <out-root>/regenie/logs/<covset>_rep<k>/
"""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path("${ADSP_ROOT}")
SCRIPTS = ROOT / "revision" / "scripts"
sys.path.insert(0, str(SCRIPTS))
import regenie_single_variant as R  # noqa: E402  (build_inputs, find_npz, REGENIE_BIN)
import benchmark_null_harness as H  # noqa: E402

T0 = time.time()


def log(m: str) -> None:
    print(f"[rv-regenie {time.time() - T0:7.1f}s] {m}", flush=True)


def lam(p: np.ndarray) -> float:
    from scipy.stats import chi2
    p = p[np.isfinite(p) & (p > 0) & (p <= 1)]
    return float(np.median(chi2.isf(p, 1)) / chi2.isf(0.5, 1)) if len(p) else float("nan")


def mixed_rows(lab: np.ndarray, y: np.ndarray) -> np.ndarray:
    d = pd.DataFrame({"l": lab, "y": y})
    ct = d.groupby("l").y.agg(["sum", "size"])
    ok = ct.index[(ct["sum"] > 0) & (ct["sum"] < ct["size"])]
    return np.flatnonzero(d.l.isin(ok).to_numpy())


def write_inputs(pool: pd.DataFrame, ystar: np.ndarray, covset: str, workdir: Path) -> dict:
    """pheno / covar files (FID 0 = the exported .fam convention; REGENIE matches samples by FID+IID, order is free)."""
    fam = pd.DataFrame({"FID": "0", "IID": pool.IID.astype(str).values})
    if covset == "pooled_project":
        inp = R.build_inputs(pool, ystar, fam, "pooled_project", workdir)
        inp["rows"] = np.arange(len(pool))
        return inp
    lab = (pool["Study"].astype(str) + "|" + pool["center"].astype(str)).to_numpy() if covset == "strat_project" else pool["cell_A"].astype(str).to_numpy()
    rows = mixed_rows(lab, ystar)
    sub = pool.iloc[rows].reset_index(drop=True)
    inp = R.build_inputs(sub, ystar[rows], fam.iloc[rows].reset_index(drop=True), "pooled_pcs", workdir)   # Age + Sex + PC1-20, no categoricals
    cov = pd.read_csv(inp["covar"], sep="\t", dtype={"FID": str, "IID": str})
    cov["cell"] = pd.Series(lab[rows]).str.replace(r"[^A-Za-z0-9_.-]", "_", regex=True).to_numpy()
    if covset == "saturated_cell":
        cov = cov[["FID", "IID", "Age", "Sex", "cell"]]
        inp["quant"] = ["Age", "Sex"]
    cov.to_csv(inp["covar"], sep="\t", index=False)
    inp["cats"] = ["cell"]
    inp["levels"] = {"cell": int(cov["cell"].nunique())}
    inp["rows"] = rows
    log(f"covset {covset}: {len(rows)} subjects in {inp['levels']['cell']} mixed cells, {int(ystar[rows].sum())} cases")
    return inp


def regenie_common(args, inp: dict, bed: str) -> list:
    keep = inp["quant"] + inp["cats"]
    cmd = [args.regenie, "--bed", bed, "--phenoFile", str(inp["pheno"]), "--covarFile", str(inp["covar"]),
           "--covarColList", ",".join(keep), "--bt", "--threads", str(args.threads)]
    if inp["cats"]:
        cmd += ["--catCovarList", ",".join(inp["cats"]), "--maxCatLevels", "999"]
    return cmd


def run(cmd: list, stdout: Path, what: str) -> None:
    log(f"{what}: " + " ".join(cmd))
    t = time.time()
    rc = subprocess.run(cmd, stdout=open(stdout, "w"), stderr=subprocess.STDOUT).returncode
    if rc != 0:
        raise SystemExit(f"REGENIE {what} exit {rc}; see {stdout}")
    log(f"{what} done in {time.time() - t:.0f}s")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--replicate", type=int, required=True)
    ap.add_argument("--covset", default="pooled_project", choices=["pooled_project", "strat_project", "saturated_cell"])
    ap.add_argument("--bfile", required=True, help="merged coding-variant PLINK prefix (plasmode_rv_merge.sh)")
    ap.add_argument("--groups-dir", required=True, help="directory with regenie.setlist / regenie.anno / regenie.maskdef")
    ap.add_argument("--native-root", default=str(ROOT / "outputs/plasmode/power/bignull"))
    ap.add_argument("--out-root", default=str(ROOT / "outputs/plasmode/rv_tools"))
    ap.add_argument("--step1-reuse", default=str(ROOT / "outputs/plasmode/power/tools_bignull/logs"),
                    help="pooled_project only: <dir>/pooled_project_rep<k>/step1_pred.list is reused when present")
    ap.add_argument("--step1-bfile", default=str(H.OUT_ROOT / "pool_281k"))
    ap.add_argument("--step1-extract", default=str(H.OUT_ROOT / "pool_281k_step1_extract_mac5.txt"))
    ap.add_argument("--pool-csv", default=str(H.OUT_ROOT / "plasmode_pool.csv"))
    ap.add_argument("--manifest", default=str(ROOT / "outputs/plasmode/bignull/seed_manifest_bignull.csv"))
    ap.add_argument("--regenie", default=R.REGENIE_BIN)
    ap.add_argument("--threads", type=int, default=int(os.environ.get("SLURM_CPUS_PER_TASK", "8")))
    ap.add_argument("--scratch", default=os.environ.get("SIM_SCRATCH", "/dev/shm/plasmode_rv_regenie"))
    ap.add_argument("--tag", default="")
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()
    rep = int(args.replicate)
    out_dir = Path(args.out_root) / f"regenie{args.tag}"
    out_dir.mkdir(parents=True, exist_ok=True)
    stem = f"{args.covset}_rep{rep:03d}"
    out_tsv = out_dir / f"{stem}_genes.tsv.gz"
    if out_tsv.exists() and not args.force:
        log(f"{out_tsv.name} exists -> skip")
        return 0
    pool = pd.read_csv(args.pool_csv, low_memory=False)
    src = R.find_npz(Path(args.native_root) / "npz", rep, Path(args.manifest))
    z = np.load(src, allow_pickle=True)
    ystar = z["ystar"].astype(np.int64)
    seed = int(z["seed"]) if "seed" in z.files else -1
    if len(ystar) != len(pool):
        raise SystemExit(f"{src.name}: ystar {len(ystar)} != pool {len(pool)}")
    log(f"replicate {rep}: Y* from {src.name} ({int(ystar.sum())} cases / {len(ystar)}), covset {args.covset}")
    workdir = Path(args.scratch) / stem
    logdir = out_dir / "logs" / stem
    workdir.mkdir(parents=True, exist_ok=True)
    logdir.mkdir(parents=True, exist_ok=True)
    inp = write_inputs(pool, ystar, args.covset, workdir)
    # ---- step 1 (LOCO ridge predictions on the 281k LD-pruned common panel)
    pred = Path(args.step1_reuse) / f"pooled_project_rep{rep:03d}" / "step1_pred.list" if args.covset == "pooled_project" else None
    step1_reused = bool(pred and pred.exists())
    if not step1_reused:
        s1 = logdir / "step1"
        # subset covariate sets: SNPs with MAC >= 5 INSIDE the analyzed subset (plasmode_rv_step1_extract.sbatch); the pool-level list
        # holds SNPs that are monomorphic in the subset and step 1 aborts on them (smoke job 3268958)
        sub_list = Path(args.out_root).parent / f"step1_extract_{args.covset}.txt"
        ext = sub_list if sub_list.exists() else Path(args.step1_extract)
        extract = ["--extract", str(ext)] if args.step1_extract and ext.exists() else []
        log(f"step 1 extract list: {ext}")
        cmd1 = regenie_common(args, inp, args.step1_bfile) + extract + ["--step", "1", "--bsize", "1000", "--loocv", "--lowmem",
                                                                         "--lowmem-prefix", str(workdir / "l0"), "--out", str(s1)]
        run(cmd1, logdir / "step1.stdout", "step 1")
        pred = Path(str(s1) + "_pred.list")
    else:
        log(f"step 1 reused: {pred}")
    # ---- step 2 gene-based
    gd = Path(args.groups_dir)
    s2 = logdir / "step2"
    cmd2 = regenie_common(args, inp, args.bfile) + ["--step", "2", "--pred", str(pred), "--firth", "--approx", "--pThresh", "0.05",
                                                      "--set-list", str(gd / "regenie.setlist"), "--anno-file", str(gd / "regenie.anno"),
                                                      "--mask-def", str(gd / "regenie.maskdef"), "--aaf-bins", "0.01", "--build-mask", "sum",
                                                      "--vc-tests", "skato,acato", "--vc-maxAAF", "0.01", "--bsize", "200", "--out", str(s2)]
    run(cmd2, logdir / "step2.stdout", "step 2")
    res = Path(str(s2) + "_Y.regenie")
    r = pd.read_csv(res, sep=r"\s+", comment="#")
    # ID = <gene>.<mask>.<aaf bin>; gene names may contain dots (FP565260.3) and the bin is "0.01" or "singleton" (REGENIE adds a
    # singleton bin by default); keep the 0.01 bin only, which holds every exported variant (pool MAF < 0.1%)
    parts = r["ID"].astype(str).str.extract(r"^(?P<gene>.+)\.(?P<mask>coding|nonsyn|syn)\.(?P<bin>singleton|all|[0-9.]+)$")
    r = pd.concat([r, parts], axis=1)
    r = r[r["mask"].notna() & (r["bin"] != "singleton")].copy()
    r["test"] = r["TEST"].astype(str).str.replace("ADD-", "", regex=False).str.replace("ADD", "burden", regex=False).str.lower()
    r["p"] = np.power(10.0, -pd.to_numeric(r["LOG10P"], errors="coerce"))
    wide = r.pivot_table(index=["gene", "mask"], columns="test", values="p", aggfunc="first")
    wide.columns = [f"p_regenie_{c}" for c in wide.columns]
    n = r[r["test"] == "burden"].set_index(["gene", "mask"])[["A1FREQ", "N"]].rename(columns={"A1FREQ": "mask_aaf", "N": "n_regenie"})
    wide = wide.join(n, how="left").reset_index()
    wide.to_csv(out_tsv, sep="\t", index=False)
    pc = [c for c in wide.columns if c.startswith("p_regenie_")]
    meta = dict(covset=args.covset, replicate=rep, seed=seed, n_subjects=int(len(inp["rows"])), n_cases=int(ystar[inp["rows"]].sum()),
                levels=inp["levels"], step1_reused=step1_reused, n_rows=int(len(wide)),
                masks={m: int((wide["mask"] == m).sum()) for m in sorted(wide["mask"].unique())},
                lambda_coding={c: lam(wide.loc[wide["mask"] == "coding", c].to_numpy(dtype=float)) for c in pc},
                n_tested_coding={c: int(np.isfinite(wide.loc[wide["mask"] == "coding", c].to_numpy(dtype=float)).sum()) for c in pc},
                build_mask="sum", aaf_bins="0.01", vc_tests="skato,acato", seconds=time.time() - T0)
    (out_dir / f"{stem}_meta.json").write_text(json.dumps(meta, indent=1, default=float))
    log(f"wrote {out_tsv.name}: {len(wide)} gene x mask rows; coding lambda " +
        ", ".join(f"{k.replace('p_regenie_', '')} {v:.3f}" for k, v in meta["lambda_coding"].items()))
    return 0


if __name__ == "__main__":
    sys.exit(main())
