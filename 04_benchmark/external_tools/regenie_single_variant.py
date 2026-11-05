#!/usr/bin/env python
"""
regenie_single_variant.py -- deliverable B, EXTERNAL-TOOL arm of the plasmode harness.

One SLURM task = one (covariate set, replicate).  It re-runs a replicate of
benchmark_null_harness.py through REGENIE v3.2.9 instead of the harness's own
score test, so that prediction B10 (REGENIE Firth with PCs only is calibrated,
1.00-1.10, because the LOCO ridge absorbs the centre / project cell structure
as "polygenic") can be read against the native arms on the SAME Y*, the SAME
27,306 subjects and the SAME 281,156 LD-pruned SNPs.

 covariate sets (4)      the native arm they mirror
   pooled_pcs      N1    Age + Sex + PC1-20
   pooled_center   N2    N1 + centre (11 merged levels, --catCovarList)
   pooled_cohort   N4    N2 + Cohort_cov (levels with n < 50 or monomorphic Y*
                         collapsed to Other, plasmode_null_harness.collapse_levels)
   pooled_project  N5    N2 + Study_cov  (same rule on the GCAD project label)
 replicates (3)          1, 2, 3 of scheme A

Y* is NOT redrawn: it is read from the replicate's npz sidecar written by
benchmark_null_harness.py (outputs/plasmode/native/npz/rep<k>_seed<seed>.npz,
key 'ystar'), so a REGENIE row and the native rows of the same replicate are
the same permutation.  The task therefore REQUIRES that replicate to have been
run natively first.

 step 1  --bt --loocv --bsize 1000 --lowmem (--lowmem-prefix in the scratch)
 step 2  --bt --firth --approx --pThresh 0.05
 LOG10P of the ADD test is matched to the .bim by variant id; p = 10^-LOG10P,
 untested variants (REGENIE's --minMAC) stay NaN.

Outputs (outputs/plasmode/regenie/)
   npz/regenie_<set>_rep<k>.npz       p, log10p, beta, se, chisq, a1freq,
                                      tested, ystar, seed, replicate
   units_plasmode_regenie.csv         one row per (set, replicate) x test, the
                                      metrics() block of the harness (copied
                                      verbatim below) on the same subsets
   bybin_plasmode_regenie.csv         MAF / region / missingness / t-gated bins
   logs/<set>_rep<k>/                 REGENIE step 1 + step 2 logs and .regenie

Never edits an existing file; heavy compute only via plasmode_regenie.sbatch.
"""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import time
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np
import pandas as pd
from scipy.stats import chi2

ROOT = Path("${ADSP_ROOT}")
SCRIPTS = ROOT / "revision" / "scripts"
OUT_ROOT = ROOT / "outputs" / "plasmode"
REGENIE_BIN = os.environ.get("REGENIE_BIN", "${HOME_DIR}/bin/regenie")
sys.path.insert(0, str(SCRIPTS))
import benchmark_null_harness as H   # noqa: E402  (pandas/numpy/statsmodels only; no simulator)

SETS = {
    "pooled_pcs": dict(label="N1", cat=[]),
    "pooled_center": dict(label="N2", cat=["center"]),
    "pooled_cohort": dict(label="N4", cat=["center", "Cohort_cov"]),
    "pooled_project": dict(label="N5", cat=["center", "Study_cov"]),
}
LABEL2SET = {v["label"]: k for k, v in SETS.items()}
T0 = time.time()


def log(msg: str) -> None:
    print(f"[{time.time() - T0:8.1f}s] {msg}", flush=True)


# ----------------------------------------------------------------- metrics --
CHI2_MED = chi2.ppf(0.5, 1)
ALPHAS = ((1e-3, "1e3"), (1e-4, "1e4"), (1e-5, "1e5"), (5e-8, "5e8"))


def metrics(p: np.ndarray) -> Dict:
    """COPY of plasmode_null_harness.metrics (itself a copy of
    arenaA_pathC_repair.metrics L330-350).  Verified against H.metrics at start-up."""
    ok = np.isfinite(p) & (p > 0) & (p <= 1)
    pp = p[ok]
    out: Dict = dict(n=int(pp.size))
    if pp.size == 0:
        out["lambda_gc"] = np.nan
        for _, tag in ALPHAS:
            out[f"n_lt_{tag}"] = 0
            out[f"fpr_{tag}"] = np.nan
            out[f"fpr_{tag}_x"] = np.nan
        out["min_p"] = np.nan
        return out
    out["lambda_gc"] = float(np.median(chi2.isf(pp, 1)) / CHI2_MED)
    for a, tag in ALPHAS:
        k = int((pp < a).sum())
        out[f"n_lt_{tag}"] = k
        out[f"fpr_{tag}"] = k / pp.size
        out[f"fpr_{tag}_x"] = k / pp.size / a
    out["min_p"] = float(pp.min())
    return out


def _assert_metrics_copy() -> None:
    rng = np.random.default_rng(11)
    p = rng.random(50_000) ** 3
    a, b = metrics(p), H.metrics(p)
    bad = [k for k in a if not (a[k] == b[k] or (isinstance(a[k], float) and np.isclose(a[k], b[k], rtol=0, atol=0)))]
    if bad:
        raise SystemExit(f"metrics() copy has drifted from the harness: {bad}")


# ------------------------------------------------------------------- inputs --
def find_npz(native_npz: Path, replicate: int, manifest: Path) -> Path:
    hits = sorted(native_npz.glob(f"rep{replicate:03d}_seed*.npz"))
    if not hits:
        seed = H.read_manifest(manifest, replicate)
        cand = native_npz / f"rep{replicate:03d}_seed{seed}.npz"
        raise SystemExit(f"no native npz for replicate {replicate} (looked for {cand}); "
                         "run plasmode_null.sbatch for that replicate first")
    if len(hits) > 1:
        log(f"WARNING: {len(hits)} npz for replicate {replicate}; using {hits[-1].name}")
    return hits[-1]


def build_inputs(pool: pd.DataFrame, ystar: np.ndarray, fam: pd.DataFrame, covset: str,
                 workdir: Path) -> Dict:
    """pheno / covar files in fam order; returns the file paths and column lists."""
    n = len(pool)
    if len(fam) != n or not np.array_equal(fam.IID.astype(str).values, pool.IID.astype(str).values):
        raise SystemExit("fam order != plasmode_pool.csv order; REGENIE inputs would be misaligned")
    d = pd.DataFrame({"FID": fam.FID.astype(str).values, "IID": fam.IID.astype(str).values})
    pheno = d.copy()
    pheno["Y"] = ystar.astype(np.int64)          # REGENIE --bt: 0 = control, 1 = case
    ph = workdir / "pheno.txt"
    pheno.to_csv(ph, sep="\t", index=False, na_rep="NA")

    cov = d.copy()
    age = pd.to_numeric(pool["Age"], errors="coerce")
    cov["Age"] = age.fillna(age.mean()).to_numpy(dtype=np.float64)
    cov["Sex"] = pd.to_numeric(pool["Sex"], errors="coerce").fillna(0).to_numpy(dtype=np.int64)
    quant = ["Age", "Sex"]
    for c in H.PCS20:
        if c in pool.columns:
            v = pd.to_numeric(pool[c], errors="coerce")
            cov[c] = v.fillna(v.mean()).to_numpy(dtype=np.float64)
            quant.append(c)
    cats: List[str] = []
    for c in SETS[covset]["cat"]:
        if c == "center":
            cov["center"] = pool["center"].astype(str).to_numpy()
        elif c == "Cohort_cov":
            cov["Cohort_cov"] = H.collapse_levels(pool["Cohort"], ystar).astype(str).to_numpy()
        elif c == "Study_cov":
            cov["Study_cov"] = H.collapse_levels(pool["Study"], ystar).astype(str).to_numpy()
        cov[c] = pd.Series(cov[c]).str.replace(r"[^A-Za-z0-9_.-]", "_", regex=True).to_numpy()
        cats.append(c)
    cv = workdir / "covar.txt"
    cov.to_csv(cv, sep="\t", index=False, na_rep="NA")
    lev = {c: int(pd.Series(cov[c]).nunique()) for c in cats}
    log(f"covar set {covset}: {len(quant)} quantitative ({quant[0]},{quant[1]},PC1..PC{len(quant) - 2}), "
        f"categorical {lev if lev else '{}'}; cases {int(ystar.sum())} / {n}")
    return dict(pheno=ph, covar=cv, quant=quant, cats=cats, n=n, levels=lev)


def run_regenie(args, inp: Dict, workdir: Path, logdir: Path) -> Path:
    bfile = str(args.bfile)
    thr = str(args.threads)
    # --covarColList is a KEEP FILTER over the covariate file, so the categorical columns must
    # be listed there as well as in --catCovarList (as in REGENIE's own documented example);
    # --maxCatLevels raises the default cap of 10 (centre has 11 levels, Cohort_cov ~30).
    keep = inp["quant"] + inp["cats"]
    common = ["--bed", bfile, "--phenoFile", str(inp["pheno"]), "--covarFile", str(inp["covar"]),
              "--covarColList", ",".join(keep), "--bt", "--threads", thr]
    if inp["cats"]:
        common += ["--catCovarList", ",".join(inp["cats"]), "--maxCatLevels", "99"]
    s1 = logdir / "step1"
    # Step 1 aborts on any zero-variance SNP ("!! Uh-oh, SNP ... has low variance", job 1765486, 2026-09-20): the 281k
    # set was LD-pruned on all 58,507 subjects and a few SNPs are monomorphic in the 27,306-subject unrelated pool.
    # --extract restricts STEP 1 ONLY to SNPs with pool MAC >= 5 (plink2 --freq counts, pool_281k_step1_extract_mac5.txt);
    # step 2 still tests every SNP (its own --minMAC 5 default skips the untestable ones).
    extract = [] if not args.step1_extract or not Path(args.step1_extract).exists() else ["--extract", str(args.step1_extract)]
    if args.step1_extract and not extract:
        log(f"WARNING: step-1 extract list {args.step1_extract} not found; running step 1 on all SNPs")
    cmd1 = [args.regenie, "--step", "1"] + common + extract + [
        "--bsize", str(args.step1_bsize), "--loocv", "--lowmem",
        "--lowmem-prefix", str(workdir / "l0"), "--out", str(s1)]
    log("step 1: " + " ".join(cmd1))
    t = time.time()
    rc = subprocess.run(cmd1, stdout=open(logdir / "step1.stdout", "w"), stderr=subprocess.STDOUT).returncode
    if rc != 0:
        raise SystemExit(f"REGENIE step 1 exit {rc}; see {logdir / 'step1.stdout'} and {s1}.log")
    log(f"step 1 done in {time.time() - t:.0f}s")
    pred = Path(str(s1) + "_pred.list")
    if not pred.exists():
        raise SystemExit(f"REGENIE step 1 wrote no {pred}")
    s2 = logdir / "step2"
    cmd2 = [args.regenie, "--step", "2"] + common + [
        "--bsize", str(args.step2_bsize), "--firth", "--approx", "--pThresh", "0.05",
        "--pred", str(pred), "--out", str(s2)]
    log("step 2: " + " ".join(cmd2))
    t = time.time()
    rc = subprocess.run(cmd2, stdout=open(logdir / "step2.stdout", "w"), stderr=subprocess.STDOUT).returncode
    if rc != 0:
        raise SystemExit(f"REGENIE step 2 exit {rc}; see {logdir / 'step2.stdout'} and {s2}.log")
    log(f"step 2 done in {time.time() - t:.0f}s")
    res = Path(str(s2) + "_Y.regenie")
    if not res.exists():
        hits = sorted(logdir.glob("step2_*.regenie"))
        if not hits:
            raise SystemExit(f"REGENIE step 2 wrote no .regenie next to {s2}")
        res = hits[0]
    return res


def parse_regenie(res: Path, bim: pd.DataFrame) -> Dict[str, np.ndarray]:
    """LOG10P of the ADD test by variant id, reindexed to the .bim order."""
    r = pd.read_csv(res, sep=r"\s+", comment="#")
    if "TEST" in r.columns:
        r = r[r["TEST"].astype(str).str.upper().isin(("ADD", "ADD-FIRTH", "ADD-SPA"))]
    r = r.drop_duplicates(subset="ID", keep="first").set_index("ID")
    r = r.reindex(bim["id"].astype(str).values)
    log10p = pd.to_numeric(r.get("LOG10P"), errors="coerce").to_numpy(dtype=np.float64)
    p = np.power(10.0, -log10p)
    out = dict(log10p=log10p, p=p, tested=np.isfinite(log10p))
    for src, dst in (("BETA", "beta"), ("SE", "se"), ("CHISQ", "chisq"), ("A1FREQ", "a1freq")):
        out[dst] = (pd.to_numeric(r.get(src), errors="coerce").to_numpy(dtype=np.float64)
                    if src in r.columns else np.full(len(bim), np.nan))
    return out


# ------------------------------------------------------------------- driver --
def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("Outputs")[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--covset", required=True, help="pooled_pcs | pooled_center | pooled_cohort | pooled_project (or N1/N2/N4/N5)")
    ap.add_argument("--replicate", type=int, required=True)
    ap.add_argument("--scheme", default=None, help="default: from the replicate (1-20 -> A)")
    ap.add_argument("--pool-csv", default=str(OUT_ROOT / "plasmode_pool.csv"))
    ap.add_argument("--manifest", default=str(OUT_ROOT / "seed_manifest.csv"))
    ap.add_argument("--bfile", default=str(OUT_ROOT / "pool_281k"))
    ap.add_argument("--native-root", default=str(OUT_ROOT / "native"))
    ap.add_argument("--out-root", default=str(OUT_ROOT / "regenie"))
    ap.add_argument("--persnp", default=str(H.PERSNP))
    ap.add_argument("--scratch", default=os.environ.get("SIM_SCRATCH", "/dev/shm/plasmode_regenie"))
    ap.add_argument("--regenie", default=REGENIE_BIN)
    ap.add_argument("--threads", type=int, default=int(os.environ.get("SLURM_CPUS_PER_TASK", "8")))
    ap.add_argument("--step1-bsize", type=int, default=1000)
    ap.add_argument("--step1-extract", default=str(OUT_ROOT / "pool_281k_step1_extract_mac5.txt"),
                    help="SNP list for REGENIE step 1 only (pool MAC >= 5); '' to disable")
    ap.add_argument("--step2-bsize", type=int, default=400)
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()

    covset = LABEL2SET.get(args.covset, args.covset)
    if covset not in SETS:
        raise SystemExit(f"--covset {args.covset}: expected one of {sorted(SETS)} or {sorted(LABEL2SET)}")
    _assert_metrics_copy()
    rep = int(args.replicate)
    scheme = args.scheme or H.task_scheme(rep)
    out_root = Path(args.out_root)
    npz_dir = out_root / "npz"
    npz_dir.mkdir(parents=True, exist_ok=True)
    units_csv = out_root / "units_plasmode_regenie.csv"
    bybin_csv = out_root / "bybin_plasmode_regenie.csv"
    npz_path = npz_dir / f"regenie_{covset}_rep{rep:03d}.npz"
    if npz_path.exists() and not args.force:
        log(f"{npz_path.name} exists -> skip (use --force)")
        return 0
    log(f"covset={covset} ({SETS[covset]['label']}) replicate={rep} scheme={scheme} host={os.uname().nodename}")

    pool = pd.read_csv(args.pool_csv, low_memory=False)
    src = find_npz(Path(args.native_root) / "npz", rep, Path(args.manifest))
    z = np.load(src, allow_pickle=True)
    ystar = z["ystar"].astype(np.int64)
    seed = int(z["seed"]) if "seed" in z.files else H.read_manifest(Path(args.manifest), rep)
    if len(ystar) != len(pool):
        raise SystemExit(f"{src.name}: ystar {len(ystar)} != pool {len(pool)}")
    log(f"Y* from {src.name}: seed {seed}, {int(ystar.sum())} cases / {len(ystar)} (real AD {int(pool.AD.sum())})")

    fam = pd.read_csv(str(args.bfile) + ".fam", sep=r"\s+", header=None,
                      names=["FID", "IID", "PAT", "MAT", "SEX", "PHENO"], dtype={0: str, 1: str})
    bim = pd.read_csv(str(args.bfile) + ".bim", sep=r"\s+", header=None,
                      names=["chr", "id", "cm", "pos", "a1", "a2"], dtype={1: str})
    workdir = Path(args.scratch) / f"{covset}_rep{rep:03d}"
    logdir = out_root / "logs" / f"{covset}_rep{rep:03d}"
    workdir.mkdir(parents=True, exist_ok=True)
    logdir.mkdir(parents=True, exist_ok=True)

    inp = build_inputs(pool, ystar, fam, covset, workdir)
    res = run_regenie(args, inp, workdir, logdir)
    rr = parse_regenie(res, bim)
    log(f"{res.name}: {int(rr['tested'].sum())} of {len(bim)} variants tested "
        f"({100 * rr['tested'].mean():.2f}%), min p {np.nanmin(rr['p']):.3g}")

    # ---- metrics on the harness's subsets (same annotation, same masks)
    a1 = rr["a1freq"]
    maf = np.where(np.isfinite(a1), np.minimum(a1, 1.0 - a1), np.nan)
    ann = H.snp_annotation(bim, Path(args.persnp) if args.persnp else None)
    masks = H.subset_masks(ann, maf)
    rng = np.random.default_rng(np.random.SeedSequence(int(seed)).spawn(2)[1])
    common = dict(design=scheme, scheme=scheme, scenario="adsp_pool", seed=int(seed), replicate=rep,
                  n_pool=int(len(pool)), n_snps_bed=int(len(bim)), n_excluded_windows=int(ann.excluded.sum()),
                  bfile=Path(args.bfile).name, tool="regenie", is_external=True,
                  driver=Path(__file__).name, arm=covset, arm_label=SETS[covset]["label"],
                  is_matched=False, n_pairs=np.nan, n_subjects=int(len(pool)),
                  n_cases_arm=int(ystar.sum()), covset=covset, ok=True, matcher="none",
                  covariates=";".join(inp["cats"]) + ";Age;Sex;PC1..PC20",
                  n_tested=int(rr["tested"].sum()), regenie_version="3.2.9", error="")
    p = rr["p"]
    unit_rows = [dict(common, test="firth", **metrics(p[masks["all"]]),
                      **{f"{k}_incl": v for k, v in metrics(p[masks["all_incl_windows"]]).items()},
                      **{f"{k}_clean": v for k, v in metrics(p[masks["clean"]]).items()},
                      **H.window_stats(p, ann.window.to_numpy(), masks["all"], rng),
                      minp_apoe_window=(float(np.nanmin(p[masks["apoe_window"]])) if masks["apoe_window"].any() else np.nan),
                      minp_known_windows=(float(np.nanmin(p[masks["known_windows"]])) if masks["known_windows"].any() else np.nan))]
    bybin_rows = []
    for lab, mk in masks.items():
        if lab in ("all", "apoe_window", "known_windows"):
            continue
        btype = "maf" if lab.startswith("maf") else ("batch_class" if "tgated" in lab else "subset")
        bybin_rows.append(dict(common, test="firth", bin_type=btype, bin_label=lab,
                               n_bin=int(mk.sum()), **metrics(p[mk])))
    r0 = unit_rows[0]
    log(f"[{SETS[covset]['label']} {covset} firth] lam={r0['lambda_gc']:.3f} 1e-3={r0['fpr_1e3_x']:.2f}x "
        f"1e-5={r0['fpr_1e5_x']:.1f}x clean lam={r0['lambda_gc_clean']:.3f} "
        f"windows<1e-5={r0.get('n_windows_hit_1e5', np.nan)}")

    for path in (units_csv, bybin_csv):
        if path.exists() and path.stat().st_size > 0:
            d = pd.read_csv(path)
            if {"seed", "arm"} <= set(d.columns):
                keep = ~((d.seed == int(seed)) & (d.arm == covset))
                if (~keep).any():
                    d[keep].to_csv(path, index=False)
    H._append_csv(units_csv, unit_rows)
    H._append_csv(bybin_csv, bybin_rows)
    np.savez_compressed(npz_path, p=p.astype(np.float32), log10p=rr["log10p"].astype(np.float32),
                        beta=rr["beta"].astype(np.float32), se=rr["se"].astype(np.float32),
                        chisq=rr["chisq"].astype(np.float32), a1freq=a1.astype(np.float32),
                        tested=rr["tested"], ystar=ystar.astype(np.int8), maf=maf.astype(np.float32),
                        seed=np.uint64(seed), replicate=np.int64(rep),   # uint64: manifest seeds exceed 2**63 (2026-09-20)
                        chr=ann.chr.to_numpy().astype(np.int16), pos=ann.pos.to_numpy().astype(np.int64),
                        excluded=ann.excluded.to_numpy(), covset=np.array([covset]))
    (logdir / "summary.json").write_text(json.dumps(
        dict(covset=covset, replicate=rep, seed=int(seed), scheme=scheme, levels=inp["levels"],
             n_tested=int(rr["tested"].sum()), lambda_gc=r0["lambda_gc"], lambda_gc_clean=r0["lambda_gc_clean"],
             fpr_1e3_x=r0["fpr_1e3_x"], fpr_1e5_x=r0["fpr_1e5_x"], t_total_s=time.time() - T0),
        indent=1, default=float))
    log(f"wrote {npz_path.name}, +1 unit row, +{len(bybin_rows)} bybin rows; total {time.time() - T0:.0f}s")
    return 0


if __name__ == "__main__":
    sys.exit(main())
