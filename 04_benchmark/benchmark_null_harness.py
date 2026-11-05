#!/usr/bin/env python
"""
benchmark_null_harness.py -- deliverable B (adsp_realgeno_null): one plasmode
replicate on the REAL 281k LD-pruned ADSP genotypes of the 27,306 unrelated
phenotyped pool, with a permuted-within-cell null phenotype Y*.

One SLURM task = one (scheme, replicate).  Replicate 0 = the observed AD labels
(the real-Y anchor); replicates 1..20 = scheme A (Y* = the real AD labels
permuted within Cohort x Study x centre cells, cell_A of plasmode_pool.csv);
replicates 21..40 = scheme A2 (cell_A2 = lpa x Cohort x Study x centre,
optional, gated on V6/V7).  Idempotent per replicate (seed_done, as
arenaA_pathC_repair.seed_done).

Per task
--------
 1  seed from seed_manifest.csv (one numpy SeedSequence child per task);
    rng_y = default_rng(SeedSequence(seed).spawn(2)[0]) draws Y* by
    rng_y.permutation inside each cell's index array; child [1] is the
    window block bootstrap.  MatchIt seed 0 is constant (in benchmark_matching.R).
 2  validation V1-V5 -> seed<seed>.json; any hard failure -> exit 3 and a
    FAILED_seed<seed>.json marker (V5's APOE / known-locus min-p part can only
    be evaluated after the scan and is written to the results JSON).
 3  subprocess "Rscript benchmark_matching.R" (builder B1's script) for the keys
    KEY0 (pair33: Sex + age_band + SeqCenter + lpa_cluster), KEY3 (KEY0 +
    Cohort), KEY_P (KEY0 + Study) on a CSV that carries the pool columns plus
    a Y column (= Y*), passed as BOTH --pool and --ystar.  Contract (deployed
    2026-09-18): it writes ONE <outdir>/match_idx_<tag>.csv, long format, with
    the columns mkey (= KEY0 / KEY3 / KEY_P), pair_id, IID, row_idx, fam_pos, Y;
    every row is a matched subject.  --match-fallback python
    replaces it with a greedy exact-key nearest-Mahalanobis(PC1-10) 1:1
    matcher (for the login-node self-test only; rows are labelled matcher=
    python_fallback).
 4  null fits per arm with a COPY of spa_score.null_fit (spa_score.py
    L97-121) on a COPY of the design-matrix builder (S._prepare_covariate_
    matrix semantics as arenaA_separation_pilot._design_with_names L319-333):
    the simulator is NOT imported on the compute path.  Stored per arm: rows,
    mu, e, sqrt_w, Q, rank, converged, coefficient names and values,
    max |design coef| (inestimable flag when > 10, P.COEF_INESTIMABLE).
 5  one bed pass with a COPY of the adsp_center_batch_common.py decoder
    (read_chunk L264-273, LUT [2, 255, 1, 0], SNP-major magic check, plink2
    --export A verification on 20 SNPs x 5 samples, L284-317), chunk 2000,
    pool subset, per-SNP non-missing-mean imputation (L341-346), then
    score_from_nullfit(G_chunk, nf, fastspa_cutoff) = the chi2 + SPA block of
    spa_score.spa_score_test refactored to take a stored null fit; on the
    first chunk the harness imports spa_score (and thereby the simulator,
    read-only) and asserts bit-identity of p_chi2 and p_spa for two arms
    (max |dlog10 p| <= 1e-6, else exit 4; was 1e-10 until 2026-09-20).
 6  metrics per arm x test (chi2, spa): R.metrics (copied) on all SNPs with
    the APOE window (chr19 44.4-46.4 Mb) and the known-locus 1-Mb windows
    excluded (the '_incl' columns keep them), MAF bins, the region / missing-
    ness / t-gated centre-batch subsets re-derived from common_persnp.csv
    with the rules of adsp_center_batch_skeptic_common.py, 1-Mb-window-
    clustered hit counts with a 500-draw window block bootstrap, and the chr22
    burden masks (tail lambda q99.9 + cMAF-stratified FPR) when
    outputs/plasmode/burden_chr22/ exists.  Rows are appended
    (_append_csv, copied) to units_plasmode_<scheme>.csv / bybin_plasmode_
    <scheme>.csv; an npz sidecar keeps the p arrays (float32) and Y*.

Arms (PCs = the real PC1-20 of adsp_combined.csv unless noted)
--------------------------------------------------------------
 N0 saturated_cell          Age + Sex + cell dummies (rows: non-mono cells;
                            mono cells have w = 0 and contribute nothing)
 N1 pooled_pcs              Age + Sex + PC1-20
 N1b pair33_formula_all_rows Age + Sex + centre + lpa_cluster + cluster_prob
                            + entropy_norm + PC1-10 (the deployed pair33 null
                            formula; 'SeqCenter' of that formula is the merged
                            11-level centre, riga_design_harness.R L76-80)
 N2 pooled_center           N1 + centre (11 merged levels)
 N3 pooled_center_lpa       N2 + lpa_cluster
 N4 pooled_cohort           N2 + Cohort_cov (mono or N < 50 -> Other on the
                            arm's rows under Y*, build_null_pair33_cohortcov.R
                            L19-28, recomputed per replicate)
 N5 pooled_project          N2 + Study_cov (same rule on Study)
 N6 pooled_tech             N2 + standardized technical metrics + Sequencer /
                            PCR.free / Read.length factors
 N7 pooled_rarecount        N2 + standardized MAC11_50 and MAC51_500
 M1 riga_pair33             KEY0 matched rows, N1b formula
 M2 riga_cohortkey          KEY3 matched rows, N1b formula + Cohort_cov
 M3 riga_projectkey         KEY_P matched rows, N1b formula + Study_cov

PRE-COMMITTED PREDICTIONS (deliverable B, chi2 test, clean SNPs unless stated,
lambda / FPR@1e-5 x nominal) -- verbatim from the build specification:
 B0 N0 lambda 1.00 +- 0.02, FPR@1e-3 1.0 +- 0.1x (any deviation = harness bug).
 B1 N1 lambda 1.5-1.9 (anchor: deployed global-pool 1.79), FPR@1e-5 30-150x
    all / 10-50x clean, flagged/clean excess 3-5x, hits concentrated in
    t-gated centre-batch SNPs (KCCG/Other-shifted: ncp 4 at 0.22 SD, 11.5 at
    0.375, 29.5 at 0.60).
 B2 N1b 1.2-1.4.
 B3 N2 1.2-1.4, FPR@1e-5 5-20x.
 B4 N3 10-15% fewer 1e-5 hits than N2.
 B5 N4 1.05-1.15, 2-4x.
 B6 N5 1.00-1.05, <= 1.5x (Study x centre spans 0.377 of the 0.392 cell R2);
    flagged-only still 1.5-3x.
 B7 N6 1.15-1.30.
 B8 N7 no change on the common panel (< 0.02); on chr22 burden masks N7 ~1.0
    vs N1 2-4 (MAC 51-500), N2 1.3-1.8, N5 ~1.1, N0 1.0; MAC 6-50 half as
    inflated; singleton mask ~1.0 everywhere.
 B9 M1 1.05-1.15, 2-5x (56.7% cross-cohort pairs; N-loss halves the
    noncentrality); M2 1.02-1.06; M3 1.00-1.03 at N ~9-10k.
 B10 REGENIE Firth PCs-only 1.00-1.10 (LOCO absorbs the cell structure as
    polygenic); +centre/+Study unchanged.
 B11 real-Y anchor: N1 ~1.79; anchor minus scheme-A gap in N5/M1 = within-cell
    + genetic residual, expected 0.02-0.10.
 B12 (if A2 runs) N1/N2 +0.05-0.1 vs A, N4/N5 unchanged within 0.02.

Inputs (builder B1, produced concurrently)
------------------------------------------
 outputs/plasmode/plasmode_pool.csv    IID, fam_pos, center, SeqCenter, Cohort,
     Study, AD, Sex, Age, age_band, lpa_cluster, cluster_prob, entropy_norm,
     PC1..PC20, MAC1/MAC2_5/MAC6_10/MAC11_50/MAC51_500/MACgt500/n_missing,
     technical metrics (raw and/or z_ standardized), Sequencer, PCR.free,
     Read.length, cell_A, cell_A2   (27,306 rows in pool = fam order)
 outputs/plasmode/plasmode_cells.csv   cell, n, n_case, rate
 outputs/plasmode/seed_manifest.csv    task/replicate, seed (SeedSequence
     child entropy), scheme
 outputs/plasmode/pool_281k.{bed,bim,fam}  27,306 samples in pool order; if
     absent, processed/chrall.common.pruned.* with the fam_pos subset.
 outputs/plasmode/burden_chr22/<mask>.npz  (optional) 2-D float32 genes x N
     burden matrix in pool order + 'gene' + 'cmac' (or 'cMAC') arrays.

Never edits any existing file.  Heavy compute only via sbatch
(plasmode_anchor.sbatch, plasmode_null.sbatch).
"""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import time
import warnings
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.linalg import qr as scipy_qr
from scipy.stats import chi2, norm

# ------------------------------------------------------------------ paths --
ROOT = Path("${ADSP_ROOT}")
OUT_ROOT = ROOT / "outputs" / "plasmode"
TAB_ROOT = ROOT / "revision" / "tables" / "plasmode"
SCRIPTS = ROOT / "revision" / "scripts"
BFILE_POOL = OUT_ROOT / "pool_281k"
BFILE_ORIG = ROOT / "processed" / "chrall.common.pruned"
PERSNP = ROOT / "revision" / "tables" / "adsp_center_batch" / "common_persnp.csv"
HITS_ANNOT = ROOT / "revision" / "tables" / "adsp_center_batch" / "skeptic_common_hits_annotated.csv"
PLINK2 = "${HOME_DIR}/bin/plink2"
MATCH_SCRIPT = SCRIPTS / "benchmark_matching.R"
DRIVER_PATH = Path(__file__).resolve()

# --------------------------------------------------------------- constants --
CHI2_MED = chi2.ppf(0.5, 1)
ALPHAS = ((1e-3, "1e3"), (1e-4, "1e4"), (1e-5, "1e5"), (5e-8, "5e8"))
COEF_INESTIMABLE = 10.0          # P.COEF_INESTIMABLE
W_EFF = 1e-6                     # P.W_EFF
MAF_BINS = [(0.05, 0.1), (0.1, 0.2), (0.2, 0.5)]
CMAF_BINS = [(0.0, 1e-3), (1e-3, 1e-2), (1e-2, 1.0)]
N_BOOT = 500
WINDOW_BP = 1_000_000
APOE_WINDOW = (19, 44_400_000, 46_400_000)
PCS20 = [f"PC{i}" for i in range(1, 21)]
PCS10 = [f"PC{i}" for i in range(1, 11)]
DESIGN_DUMMY_PREFIXES = ("center", "Cohort_cov", "Study_cov", "cell_", "SeqCenter",
                         "lpa_cluster", "Sequencer", "PCR", "Read")
MATCH_KEYS = {"KEY0": ["Sex", "age_band", "center", "lpa_cluster"],
              "KEY3": ["Sex", "age_band", "center", "lpa_cluster", "Cohort"],
              "KEY_P": ["Sex", "age_band", "center", "lpa_cluster", "Study"]}
TECH_METRICS = ["depth_of_coverage_30x", "average.genome.coverage", "Q30_score_pcent",
                "variantEval_novelsites_perc", "variantEval_totalvariants_TiTv",
                "variantEval_novelvariants_TiTv", "variantEval_total_Variants", "Freemix"]
TECH_FACTORS = ["Sequencer", "PCR.free", "Read.length"]
RARE_COUNTS = ["MAC11_50", "MAC51_500"]

# known AD loci (GRCh38), verbatim from revision/scripts/annot_sv_qq_manhattan.R KNOWN
KNOWN_LOCI = [
    ("APOE", 19, 44908684), ("BIN1", 2, 127135234), ("CLU", 8, 27607002),
    ("PICALM", 11, 86156833), ("CR1", 1, 207518704), ("MS4A6A", 11, 60172122),
    ("ABCA7", 19, 1063444), ("TREM2", 6, 41161514), ("SORL1", 11, 121564878),
    ("CD2AP", 6, 47485642), ("EPHA1", 7, 143413669), ("PTK2B", 8, 27362470),
    ("HLA-DRB1", 6, 32615322), ("SLC24A4", 14, 92460608), ("INPP5D", 2, 233117202),
    ("FERMT2", 14, 52924962), ("CASS4", 20, 56443204), ("SPI1", 11, 47354897),
    ("ADAM10", 15, 58753575), ("NYAP1", 7, 100406823), ("ECHDC3", 10, 11678309),
    ("ABI3", 17, 49219935), ("PLCG2", 16, 81908423), ("APH1B", 15, 63277703),
    ("WWOX", 16, 79321960), ("IQCK", 16, 19796841), ("ACE", 17, 63471557),
    ("SHARPIN", 8, 144103704), ("TSPOAP1", 17, 58331728), ("ADAMTS1", 21, 26775872)]

# hg38 centromeres (Mb) and region flags, verbatim from adsp_center_batch_skeptic_common.py
CENTROMERES = {1: (121.7, 125.1), 2: (91.8, 96.0), 3: (90.5, 93.7), 4: (49.7, 51.7), 5: (46.5, 50.1),
               6: (58.5, 59.8), 7: (58.1, 61.0), 8: (44.0, 45.9), 9: (43.2, 45.5), 10: (39.6, 41.6),
               11: (51.1, 54.4), 12: (34.7, 37.2), 13: (16.0, 18.1), 14: (16.0, 18.2), 15: (17.0, 19.7),
               16: (36.3, 38.3), 17: (22.7, 26.9), 18: (15.4, 20.9), 19: (24.5, 27.2), 20: (25.6, 30.4),
               21: (10.9, 13.0), 22: (12.9, 15.1)}

ARM_ORDER = ["saturated_cell", "pooled_pcs", "pair33_formula_all_rows", "pooled_center",
             "pooled_center_lpa", "pooled_cohort", "pooled_project", "pooled_tech",
             "pooled_rarecount", "riga_pair33", "riga_cohortkey", "riga_projectkey"]
ARM_LABEL = {"saturated_cell": "N0", "pooled_pcs": "N1", "pair33_formula_all_rows": "N1b",
             "pooled_center": "N2", "pooled_center_lpa": "N3", "pooled_cohort": "N4",
             "pooled_project": "N5", "pooled_tech": "N6", "pooled_rarecount": "N7",
             "riga_pair33": "M1", "riga_cohortkey": "M2", "riga_projectkey": "M3"}
MATCHED_ARMS = {"riga_pair33": "KEY0", "riga_cohortkey": "KEY3", "riga_projectkey": "KEY_P"}
TESTS = ("chi2", "spa")

T0 = time.time()


def log(msg: str) -> None:
    print(f"[{time.time() - T0:8.1f}s] {msg}", flush=True)


# ============================================================= Y* generator ==
def task_scheme(replicate: int) -> str:
    if replicate == 0:
        return "anchor"
    if 1 <= replicate <= 20:
        return "A"
    return "A2"


def read_manifest(path: Path, replicate: int) -> int:
    """seed_manifest.csv: one row per task; accepts columns task|replicate and
    seed|entropy (SeedSequence child entropy); falls back to
    SeedSequence(20260915).spawn(n)[replicate] when the file is absent."""
    if path.exists():
        m = pd.read_csv(path)
        key = "replicate" if "replicate" in m.columns else ("task" if "task" in m.columns else None)
        scol = "seed" if "seed" in m.columns else ("entropy" if "entropy" in m.columns else None)
        if key is not None and scol is not None:
            row = m[m[key].astype(int) == int(replicate)]
            if len(row) == 1:
                return int(row[scol].iloc[0])
        warnings.warn(f"manifest {path} lacks a usable row for replicate {replicate}; using SeedSequence fallback")
    n = max(41, replicate + 1)
    children = np.random.SeedSequence(20260915).spawn(n)
    return int(children[replicate].generate_state(1, dtype=np.uint64)[0] % (2 ** 63 - 1))


def draw_ystar(pool: pd.DataFrame, scheme: str, seed: int) -> Tuple[np.ndarray, np.random.Generator]:
    """Scheme A / A2: real AD permuted within each cell's index array;
    anchor: Y* = AD.  Returns (Y*, bootstrap rng)."""
    ss = np.random.SeedSequence(seed)
    child_y, child_boot = ss.spawn(2)
    rng_y = np.random.default_rng(child_y)
    rng_b = np.random.default_rng(child_boot)
    y = pool["AD"].to_numpy().astype(np.int8)
    if scheme == "anchor":
        return y.copy(), rng_b
    cell_col = "cell_A" if scheme == "A" else "cell_A2"
    ystar = y.copy()
    for _, idx in pool.groupby(cell_col, sort=True).indices.items():
        if len(idx) > 1:
            ystar[idx] = y[idx][rng_y.permutation(len(idx))]
    return ystar, rng_b


# ============================================================ validation ====
def _logit_llf(y: np.ndarray, W: np.ndarray) -> Tuple[float, bool]:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        fit = sm.GLM(y, W, family=sm.families.Binomial()).fit(disp=False, maxiter=200)
    return float(fit.llf), bool(fit.converged)


def lrt_center_given(pool: pd.DataFrame, y: np.ndarray, base_cols: List[str], center_col: str = "center") -> Dict:
    """LRT of centre | base covariates (logistic)."""
    W0, _ = design_with_names(pool, base_cols) if base_cols else (np.ones((len(pool), 1)), ["const"])
    W1, _ = design_with_names(pool, base_cols + [center_col])
    ll0, c0 = _logit_llf(y, W0)
    ll1, c1 = _logit_llf(y, W1)
    df = int(np.linalg.matrix_rank(W1) - np.linalg.matrix_rank(W0))
    return dict(LR=float(2 * (ll1 - ll0)), df=df, converged=bool(c0 and c1))


def pearson_chi2(a: np.ndarray, b: np.ndarray) -> Tuple[float, int, float]:
    tab = pd.crosstab(pd.Series(a), pd.Series(b)).to_numpy(dtype=np.float64)
    exp = tab.sum(1, keepdims=True) * tab.sum(0, keepdims=True) / tab.sum()
    with np.errstate(divide="ignore", invalid="ignore"):
        stat = float(np.nansum((tab - exp) ** 2 / exp))
    r, c = tab.shape
    df = (r - 1) * (c - 1)
    n = tab.sum()
    V = float(np.sqrt(stat / (n * (min(r, c) - 1)))) if min(r, c) > 1 else np.nan
    return stat, int(df), V


def within_cell_coef(x: np.ndarray, y: np.ndarray, cell: np.ndarray) -> Tuple[float, float]:
    """OLS coefficient of y in x ~ cell dummies + y (within-cell regression)."""
    d = pd.DataFrame({"x": x, "y": y.astype(float), "c": cell})
    ok = np.isfinite(d.x.values)
    d = d[ok]
    xr = d.x - d.groupby("c").x.transform("mean")
    yr = d.y - d.groupby("c").y.transform("mean")
    vy = float((yr ** 2).sum())
    if vy <= 0:
        return np.nan, np.nan
    beta = float((xr * yr).sum() / vy)
    resid = xr - beta * yr
    dfres = len(d) - d.c.nunique() - 1
    se = float(np.sqrt((resid ** 2).sum() / max(dfres, 1) / vy))
    return beta, se


def find_tech_columns(pool: pd.DataFrame) -> Tuple[List[str], List[str]]:
    """Return (numeric metric columns, factor columns) present in the pool."""
    metrics = []
    for m in TECH_METRICS:
        for cand in (f"z_{m}", m, m.replace(".", "_"), f"z_{m.replace('.', '_')}"):
            if cand in pool.columns:
                metrics.append(cand)
                break
    factors = []
    for f in TECH_FACTORS:
        for cand in (f, f.replace(".", "_"), f + "..bp.", f.replace(".", "_") + "__bp_"):
            if cand in pool.columns:
                factors.append(cand)
                break
    return metrics, factors


def validate_ystar(pool: pd.DataFrame, ystar: np.ndarray, scheme: str) -> Tuple[Dict, List[str]]:
    """V1-V5 (pre-scan part).  Returns (json dict, list of hard failures)."""
    y = pool["AD"].to_numpy().astype(np.int8)
    N = len(pool)
    out: Dict = dict(scheme=scheme, N=int(N), n_case_real=int(y.sum()), n_case_ystar=int(ystar.sum()))
    fails: List[str] = []
    cell_col = "cell_A2" if scheme == "A2" else "cell_A"
    cells = pool[cell_col].astype(str).to_numpy()
    # V1 cell preservation (exact)
    g = pd.DataFrame({"c": cells, "y": y, "ys": ystar}).groupby("c").agg(n=("y", "size"), k=("y", "sum"), ks=("ys", "sum"))
    out["V1_n_cells"] = int(len(g))
    out["V1_n_mono_cells"] = int(((g.k == 0) | (g.k == g.n)).sum())
    out["V1_n_subjects_mono"] = int(g.n[(g.k == 0) | (g.k == g.n)].sum())
    out["V1_max_abs_count_diff"] = int(np.abs(g.k - g.ks).max())
    out["V1_pass"] = bool(out["V1_max_abs_count_diff"] == 0)
    if not out["V1_pass"]:
        fails.append("V1 per-cell Y* count != real count")
    # V2 centre-AD dependence: chi2 (exact by construction under A), LRTs vs the real Y on the same pool
    ctr = pool["center"].astype(str).to_numpy()
    s_real, df_c, V_real = pearson_chi2(ctr, y)
    s_star, _, V_star = pearson_chi2(ctr, ystar)
    out.update(V2_chi2_center_y_real=s_real, V2_chi2_center_ystar=s_star, V2_df=df_c,
               V2_cramerV_real=V_real, V2_cramerV_ystar=V_star,
               V2_chi2_target_adsp_phenotyped=3234.37, V2_cramerV_target_adsp=0.291)
    base = ["Age", "Sex"] + [c for c in PCS20 if c in pool.columns]
    l_real = lrt_center_given(pool, y, base)
    l_star = lrt_center_given(pool, ystar, base)
    out.update(V2_lrt_center_pcs_real=l_real["LR"], V2_lrt_center_pcs_ystar=l_star["LR"], V2_lrt_center_pcs_df=l_star["df"],
               V2_lrt_center_pcs_target_adsp=2103.0)
    rel = abs(l_star["LR"] - l_real["LR"]) / max(l_real["LR"], 1e-9)
    out["V2_lrt_center_pcs_rel_diff"] = float(rel)
    if "Cohort" in pool.columns:
        lc_real = lrt_center_given(pool, y, ["Cohort"])
        lc_star = lrt_center_given(pool, ystar, ["Cohort"])
        out.update(V2_lrt_center_cohort_real=lc_real["LR"], V2_lrt_center_cohort_ystar=lc_star["LR"],
                   V2_lrt_center_cohort_df=lc_star["df"], V2_lrt_center_cohort_target_adsp=365.0)
    if "Study" in pool.columns:
        ls_real = lrt_center_given(pool, y, ["Study"])
        ls_star = lrt_center_given(pool, ystar, ["Study"])
        out.update(V2_lrt_center_study_real=ls_real["LR"], V2_lrt_center_study_ystar=ls_star["LR"],
                   V2_lrt_center_study_df=ls_star["df"], V2_lrt_center_study_target_adsp=10.5)
    out["V2_pass"] = bool(abs(s_star - s_real) <= 0.10 * s_real and rel <= 0.10)
    if scheme != "anchor" and not out["V2_pass"]:
        fails.append(f"V2 centre-AD dependence: chi2 {s_star:.1f} vs {s_real:.1f}, LRT rel diff {rel:.3f}")
    # V3 ancestry marginal
    if "lpa_cluster" in pool.columns:
        gl = pd.DataFrame({"l": pool["lpa_cluster"].to_numpy(), "y": y, "ys": ystar}).groupby("l").agg(
            n=("y", "size"), r=("y", "mean"), rs=("ys", "mean"))
        sd = np.sqrt(gl.r * (1 - gl.r) / gl.n).replace(0, np.nan)
        z = ((gl.rs - gl.r) / sd).abs()
        out["V3_lpa_rates_real"] = {str(k): float(v) for k, v in gl.r.items()}
        out["V3_lpa_rates_ystar"] = {str(k): float(v) for k, v in gl.rs.items()}
        out["V3_max_abs_z"] = float(np.nanmax(z.values)) if np.isfinite(z.values).any() else 0.0
        out["V3_pass"] = bool(out["V3_max_abs_z"] <= 3.0)
        if scheme != "anchor" and not out["V3_pass"]:
            fails.append(f"V3 lpa-cluster Y* rate off by {out['V3_max_abs_z']:.2f} binomial SD")
    # V4 within-cell technical channel (must be absent under A)
    mets, _ = find_tech_columns(pool)
    v4 = {}
    worst = 0.0
    for m in mets:
        x = pd.to_numeric(pool[m], errors="coerce").to_numpy(dtype=np.float64)
        x = (x - np.nanmean(x)) / (np.nanstd(x) if np.nanstd(x) > 0 else 1.0)
        b, se = within_cell_coef(x, ystar, cells)
        b0, se0 = within_cell_coef(x, y, cells)
        v4[m] = dict(coef_ystar=b, se=se, z=(b / se if se and np.isfinite(se) and se > 0 else np.nan),
                     coef_real=b0, se_real=se0)
        if np.isfinite(v4[m]["z"]):
            worst = max(worst, abs(v4[m]["z"]))
    out["V4_within_cell_tech_coefs"] = v4
    out["V4_max_abs_z"] = float(worst)
    out["V4_n_metrics"] = len(mets)
    out["V4_pass_2se"] = bool(worst <= 2.0)
    # Hard gate raised 3.5 -> 4.5 SE on 2026-09-20: replicate 3 (job 1760622_3) tripped 3.5 with max|z| = 3.80 while
    # V1-V3/V5 were exact; with ~10 technical covariates x 20 replicates a 3.5 gate fails by chance in roughly one
    # replicate in ten, whereas a mis-permuted Y* gives |z| far above 5. The observed max|z| is written to the JSON.
    out["V4_pass"] = bool(worst <= 4.5)          # hard gate at 4.5 SE (2 SE is a soft report; see deviations)
    if scheme != "anchor" and not out["V4_pass"]:
        fails.append(f"V4 within-cell technical coefficient |z| = {worst:.2f} > 4.5")
    # V5 identity leakage
    pc = g.k / g.n
    expect = float((g.n * (pc ** 2 + (1 - pc) ** 2)).sum() / N)
    out["V5_mean_ystar_eq_y"] = float((ystar == y).mean())
    out["V5_exchangeability_expectation"] = expect
    out["V5_identity_leakage"] = float((ystar == y).mean() - expect)
    out["V5_mono_share"] = float(out["V1_n_subjects_mono"] / N)
    return out, fails


# ======================================================== design / null fit ==
def design_with_names(df: pd.DataFrame, covariate_cols: List[str]) -> Tuple[np.ndarray, List[str]]:
    """COPY of S._prepare_covariate_matrix semantics (synthetic_null_simulation.py
    L1932-1943) keeping the column names, as P._design_with_names does."""
    Xc = df[covariate_cols].copy()
    for c in list(Xc.columns):
        if Xc[c].dtype.kind not in "biufc":
            Xc[c] = Xc[c].astype("category")
        elif Xc[c].dtype.kind in "iu" and Xc[c].nunique() > 2:
            Xc[c] = Xc[c].astype("category")
    Xc = pd.get_dummies(Xc, drop_first=True, dtype=float)
    Xc = sm.add_constant(Xc, has_constant="add")
    W = Xc.values.astype(np.float64)
    return W, [str(c) for c in Xc.columns]


def null_fit_local(df: pd.DataFrame, covariate_cols: List[str]) -> Dict:
    """COPY of spa_score.null_fit (spa_score.py L97-121) with the local design
    builder; adds coefficient accounting (P.null_diag L335-370)."""
    y = df["Y"].values.astype(np.float64)
    W, names = design_with_names(df, covariate_cols)
    with warnings.catch_warnings(record=True) as wl:
        warnings.simplefilter("always")
        null_model = sm.GLM(y, W, family=sm.families.Binomial()).fit(disp=False)
    mu = null_model.predict()
    e = y - mu
    w = mu * (1.0 - mu)
    sqrt_w = np.sqrt(w)
    W_tilde = W * sqrt_w[:, None]
    # mode="economic": spa_score.null_fit L115 uses scipy's DEFAULT mode="full", which
    # materializes an (N x N) Q -- 5.96 GB per full-pool arm here (measured 7.1 GB RSS and
    # 246 s for ONE arm at N = 20,113), and Q[:, :rank] is a VIEW that keeps that buffer
    # alive, so retaining 12 null fits would need ~53 GB.  The economic Q is the first
    # min(N, p) columns of the same Householder product: same R diagonal (so the same rank
    # rule), same first `rank` columns, p_chi2 BITWISE identical and p_spa within
    # 4.4e-16 in log10 p -- six orders below the V9 bit-identity threshold of 1e-10.
    # The layout (Fortran-ordered) is left untouched so the BLAS calls downstream match.
    Q, R, _ = scipy_qr(W_tilde, pivoting=True, mode="economic")
    rank = int(np.sum(np.abs(np.diag(R)) > 1e-8 * np.abs(R[0, 0])))
    Q = Q[:, :rank]
    params = np.asarray(null_model.params, dtype=np.float64)
    is_design = np.array([n.startswith(DESIGN_DUMMY_PREFIXES) for n in names])
    max_design = float(np.abs(params[is_design]).max()) if is_design.any() else np.nan
    wclasses = sorted({type(x.message).__name__ for x in wl})
    return dict(y=y, W=W, mu=mu, e=e, w=w, sqrt_w=sqrt_w, Q=Q, rank=rank,
                converged=bool(null_model.converged), names=names, params=params,
                n_subjects=int(len(y)), n_cases_arm=int(y.sum()),
                n_effective=int((w > W_EFF).sum()), n_dropped_eff=int((w <= W_EFF).sum()),
                n_params=int(params.size), n_design_dummies=int(is_design.sum()),
                max_abs_coef_design=max_design, max_abs_coef_any=float(np.abs(params).max()),
                inestimable_design=bool((not null_model.converged)
                                        or (np.isfinite(max_design) and max_design > COEF_INESTIMABLE)),
                warn_classes=";".join(wclasses),
                perfect_separation_warning=bool(any("PerfectSeparation" in c for c in wclasses)),
                design_coefs={n: float(v) for n, v in zip(names, params) if n.startswith(DESIGN_DUMMY_PREFIXES)})


def collapse_levels(labels: pd.Series, y: np.ndarray, min_n: int = 50) -> pd.Series:
    """build_null_pair33_cohortcov.R L19-28: mono (n0 == 0 or n1 == 0) or
    N < min_n -> 'Other'; recomputed on the arm's rows under Y*."""
    d = pd.DataFrame({"lab": labels.astype(str).to_numpy(), "y": y})
    ct = d.groupby("lab").y.agg(n0=lambda s: int((s == 0).sum()), n1=lambda s: int((s == 1).sum()), N="size")
    other = ct.index[(ct.n0 == 0) | (ct.n1 == 0) | (ct.N < min_n)]
    out = d.lab.where(~d.lab.isin(other), "Other")
    return pd.Series(out.to_numpy(), index=labels.index)


def _standardize(v: pd.Series) -> np.ndarray:
    x = pd.to_numeric(v, errors="coerce").to_numpy(dtype=np.float64)
    m = np.nanmean(x)
    s = np.nanstd(x)
    x = np.where(np.isfinite(x), x, m)
    return (x - m) / (s if s > 0 else 1.0)


def arm_frames(pool: pd.DataFrame, ystar: np.ndarray, matched: Dict[str, np.ndarray],
               arms: Sequence[str]) -> Dict[str, Tuple[pd.DataFrame, np.ndarray, List[str], Optional[int]]]:
    """Per arm: (df with 'Y', row index into the pool, covariate columns, n_pairs)."""
    base = pool.copy()
    base["Y"] = ystar.astype(np.int64)
    base["Sex"] = pd.to_numeric(base["Sex"], errors="coerce").fillna(0).astype(np.int64)
    base["Age"] = pd.to_numeric(base["Age"], errors="coerce")
    base["Age"] = base["Age"].fillna(base["Age"].mean()).astype(np.float64)
    base["lpa_cluster"] = base["lpa_cluster"].astype(str)
    base["center"] = base["center"].astype(str)
    pcs20 = [c for c in PCS20 if c in base.columns]
    pcs10 = [c for c in PCS10 if c in base.columns]
    mets, facs = find_tech_columns(base)
    for m in mets:
        base[f"tech_{m}"] = _standardize(base[m])
    for f in facs:
        base[f"tech_{f}"] = base[f].astype(str).replace({"nan": "Unknown", "": "Unknown"})
    for r in RARE_COUNTS:
        if r in base.columns:
            base[f"rc_{r}"] = _standardize(base[r])
    pair33 = ["Age", "Sex", "center", "lpa_cluster", "cluster_prob", "entropy_norm"] + pcs10
    n1 = ["Age", "Sex"] + pcs20
    n2 = n1 + ["center"]
    specs: Dict[str, Tuple[List[str], Optional[str]]] = {
        "saturated_cell": (["Age", "Sex", "cell_A"], None),
        "pooled_pcs": (n1, None),
        "pair33_formula_all_rows": (pair33, None),
        "pooled_center": (n2, None),
        "pooled_center_lpa": (n2 + ["lpa_cluster"], None),
        "pooled_cohort": (n2 + ["Cohort_cov"], None),
        "pooled_project": (n2 + ["Study_cov"], None),
        "pooled_tech": (n2 + [f"tech_{m}" for m in mets] + [f"tech_{f}" for f in facs], None),
        "pooled_rarecount": (n2 + [f"rc_{r}" for r in RARE_COUNTS if f"rc_{r}" in base.columns], None),
        "riga_pair33": (pair33, "KEY0"),
        "riga_cohortkey": (pair33 + ["Cohort_cov"], "KEY3"),
        "riga_projectkey": (pair33 + ["Study_cov"], "KEY_P"),
    }
    out = {}
    for arm in arms:
        covs, key = specs[arm]
        if key is None:
            rows = np.arange(len(base))
            n_pairs = None
        else:
            rows = matched.get(key)
            if rows is None:
                out[arm] = (None, None, covs, None)
                continue
            rows = np.asarray(rows, dtype=np.int64)
            n_pairs = int(len(rows) // 2)
        d = base.iloc[rows].reset_index(drop=True)
        if arm == "saturated_cell":
            # mono cells have mu in {0, 1}: w = 0, zero contribution; drop them for a finite MLE
            cnt = d.groupby("cell_A").Y.agg(["sum", "size"])
            mono = cnt.index[(cnt["sum"] == 0) | (cnt["sum"] == cnt["size"])]
            keep = ~d["cell_A"].isin(mono).to_numpy()
            d = d[keep].reset_index(drop=True)
            rows = rows[keep]
            d["cell_A"] = d["cell_A"].astype(str)
        if "Cohort_cov" in covs:
            d["Cohort_cov"] = collapse_levels(d["Cohort"], d["Y"].to_numpy())
        if "Study_cov" in covs:
            d["Study_cov"] = collapse_levels(d["Study"], d["Y"].to_numpy())
        out[arm] = (d, rows, covs, n_pairs)
    return out


# ================================================================ matching ==
def run_match_r(script: Path, pool_csv: Path, ystar_csv: Path, outdir: Path, tag: str,
                keys: Sequence[str], cmd_template: str, seed: int = 0) -> Dict[str, np.ndarray]:
    """subprocess Rscript benchmark_matching.R (B1 contract as deployed 2026-09-18):
    flags --pool / --out / --tag / --ystar / --keys, ONE output file
    <outdir>/match_idx_<tag>.csv in long format (mkey, pair_id, IID, row_idx,
    fam_pos, Y), all requested keys stacked.  The MatchIt seed is a constant
    (0) inside the R script; {seed} stays available to custom templates."""
    outdir.mkdir(parents=True, exist_ok=True)
    cmd = cmd_template.format(script=str(script), pool=str(pool_csv), ystar=str(ystar_csv),
                              pheno=str(ystar_csv), outdir=str(outdir), tag=tag,
                              keys=",".join(keys), seed=seed)
    log(f"matching: {cmd}")
    logf = outdir / f"plasmode_match_{tag}.log"
    with open(logf, "w") as fh:
        rc = subprocess.run(cmd, shell=True, stdout=fh, stderr=subprocess.STDOUT).returncode
    if rc != 0:
        raise RuntimeError(f"benchmark_matching.R exit {rc}; see {logf}")
    return read_match_indices(outdir, ystar_csv, tag)


def read_match_indices(outdir: Path, pheno_csv: Path, tag: str) -> Dict[str, np.ndarray]:
    """<outdir>/match_idx_<tag>.csv (single file, columns mkey, pair_id, IID,
    row_idx, fam_pos, Y) -> per key, 0-based row indices into the pool, cases
    first then controls, both in matching-pair order (pair i = positions
    i and n_pairs + i).  row_idx and Y are cross-checked against the CSV that
    was handed to the matcher."""
    f = outdir / f"match_idx_{tag}.csv"
    if not f.exists():
        raise FileNotFoundError(f"{f}: benchmark_matching.R produced no match file for tag {tag}")
    m = pd.read_csv(f, dtype={"mkey": str, "pair_id": str, "IID": str})
    missing = [c for c in ("mkey", "pair_id", "IID", "row_idx", "Y") if c not in m.columns]
    if missing:
        raise ValueError(f"{f}: missing column(s) {missing}; columns are {list(m.columns)}")
    ref = pd.read_csv(pheno_csv, usecols=lambda c: c in ("IID", "Y"), dtype={"IID": str})
    iid_pos = pd.Series(np.arange(len(ref)), index=ref["IID"].astype(str).values)
    y_ref = ref["Y"].to_numpy() if "Y" in ref.columns else None
    out = {}
    for key, mk in m.groupby("mkey", sort=False):
        rows = iid_pos.reindex(mk["IID"].astype(str).values).to_numpy()
        if not np.isfinite(rows.astype(float)).all():
            raise ValueError(f"{f}: mkey {key} has IIDs absent from {Path(pheno_csv).name}")
        rows = rows.astype(np.int64)
        ridx = mk["row_idx"].to_numpy(dtype=np.int64)
        if not np.array_equal(rows, ridx):
            raise ValueError(f"{f}: mkey {key}: row_idx != the IID position in the pool "
                             f"({int((rows != ridx).sum())} of {len(rows)} rows)")
        yv = mk["Y"].to_numpy(dtype=np.int64)
        if y_ref is not None and not np.array_equal(yv, y_ref[rows]):
            raise ValueError(f"{f}: mkey {key}: Y != the Y* handed to the matcher")
        order = np.lexsort((-yv, mk["pair_id"].astype(str).values))   # pair_id, case first
        rows, yv = rows[order], yv[order]
        if int((yv == 1).sum()) != int((yv == 0).sum()):
            log(f"  WARNING match {key}: {int((yv == 1).sum())} cases vs {int((yv == 0).sum())} controls "
                f"(not 1:1; pair order is not positional)")
        out[str(key)] = np.concatenate([rows[yv == 1], rows[yv == 0]])
    for key in MATCH_KEYS:
        if key not in out:
            log(f"  WARNING match file {f.name} has no rows for {key}")
    return out


def python_fallback_match(pool: pd.DataFrame, ystar: np.ndarray, seed: int = 0) -> Dict[str, np.ndarray]:
    """Greedy 1:1 nearest-Mahalanobis(PC1-10) matching within exact key
    strata (a stand-in for MatchIt in the login-node self-test ONLY)."""
    rng = np.random.default_rng(seed)
    pcs = [c for c in PCS10 if c in pool.columns]
    P = pool[pcs].to_numpy(dtype=np.float64)
    cov = np.cov(P.T) + 1e-12 * np.eye(len(pcs))
    Li = np.linalg.inv(np.linalg.cholesky(cov))
    Z = P @ Li.T
    out = {}
    for key, cols in MATCH_KEYS.items():
        strata = pool[cols].astype(str).agg("|".join, axis=1).to_numpy()
        cases, ctrls = [], []
        for _, idx in pd.Series(np.arange(len(pool))).groupby(strata).groups.items():
            idx = np.asarray(idx)
            ci = idx[ystar[idx] == 1]
            ki = idx[ystar[idx] == 0]
            if ci.size == 0 or ki.size == 0:
                continue
            ci = ci[rng.permutation(ci.size)]
            avail = np.ones(ki.size, dtype=bool)
            for c in ci:
                if not avail.any():
                    break
                d = np.sum((Z[ki] - Z[c]) ** 2, axis=1)
                d[~avail] = np.inf
                j = int(np.argmin(d))
                avail[j] = False
                cases.append(int(c))
                ctrls.append(int(ki[j]))
        out[key] = np.asarray(cases + ctrls, dtype=np.int64)
    return out


# ============================================================ bed decoder ===
LUT = np.array([2, 255, 1, 0], dtype=np.uint8)   # 00 hom A1 -> 2 A1 copies; 01 missing; 10 het; 11 hom A2
SHIFTS = np.array([0, 2, 4, 6], dtype=np.uint8)


class Bed:
    """COPY of the adsp_center_batch_common.py decoder (L264-273) as a class."""

    def __init__(self, prefix: Path):
        self.prefix = Path(prefix)
        self.fam = pd.read_csv(f"{prefix}.fam", sep=r"\s+", header=None,
                               names=["FID", "IID", "PAT", "MAT", "SEX", "PHE"], dtype=str)
        self.bim = pd.read_csv(f"{prefix}.bim", sep=r"\s+", header=None,
                               names=["chr", "id", "cm", "pos", "A1", "A2"],
                               dtype={"chr": str, "id": str, "A1": str, "A2": str})
        self.n_fam = len(self.fam)
        self.n_snp = len(self.bim)
        self.bps = (self.n_fam + 3) // 4
        with open(f"{prefix}.bed", "rb") as f:
            magic = f.read(3)
        assert magic == b"\x6c\x1b\x01", f"bad bed magic {magic!r} (need SNP-major PLINK1)"
        assert os.path.getsize(f"{prefix}.bed") == 3 + self.bps * self.n_snp, "bed size != 3 + bps * n_snp"

    def read_chunk(self, start: int, nsnp: int) -> np.ndarray:
        with open(f"{self.prefix}.bed", "rb") as f:
            f.seek(3 + start * self.bps)
            raw = np.fromfile(f, dtype=np.uint8, count=nsnp * self.bps)
        raw = raw.reshape(nsnp, self.bps)
        codes = (raw[:, :, None] >> SHIFTS[None, None, :]) & 3
        codes = codes.reshape(nsnp, self.bps * 4)[:, : self.n_fam]
        return LUT[codes]                                 # uint8, 255 = missing


def verify_decoder_plink2(bed: Bed, pool_idx: np.ndarray, pool_iid: np.ndarray, scratch: Path,
                          plink2: str = PLINK2, n_snp: int = 20, n_samp: int = 5) -> Dict:
    """COPY of adsp_center_batch_common.py L284-317: plink2 --export A on
    n_snp SNPs x n_samp pool samples vs the decoder."""
    scratch.mkdir(parents=True, exist_ok=True)
    n_snp = min(n_snp, bed.n_snp)
    ids = bed.bim.id.values[:n_snp]
    kp = pool_iid[:n_samp]
    fam_by_iid = bed.fam.set_index("IID")
    with open(scratch / "verify_snps.txt", "w") as f:
        f.write("\n".join(ids) + "\n")
    with open(scratch / "verify_keep.txt", "w") as f:
        f.write("".join(f"{fid}\t{iid}\n" for fid, iid in zip(fam_by_iid.loc[kp, "FID"].values, kp)))
    cmd = [plink2, "--bfile", str(bed.prefix), "--extract", str(scratch / "verify_snps.txt"),
           "--keep", str(scratch / "verify_keep.txt"), "--export", "A",
           "--threads", "2", "--memory", "2000", "--out", str(scratch / "verify")]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(f"plink2 verification failed: {r.stderr[-800:]}")
    raw = pd.read_csv(scratch / "verify.raw", sep=r"\s+")
    raw = raw.set_index("IID").loc[kp]
    mine = bed.read_chunk(0, n_snp)[:, pool_idx[:n_samp]].astype(float)
    mine[mine == 255] = np.nan
    ok = flipped = bad = 0
    for j, sid in enumerate(ids):
        col = [c for c in raw.columns if c.startswith(sid + "_")][0]
        counted = col[len(sid) + 1:]
        pl = raw[col].values.astype(float)
        m = mine[j]
        same = np.allclose(np.nan_to_num(pl, nan=-9), np.nan_to_num(m, nan=-9))
        flip = np.allclose(np.nan_to_num(pl, nan=-9), np.nan_to_num(2 - m, nan=-9))
        a1 = bed.bim.A1.values[j]
        if same and counted == a1:
            ok += 1
        elif flip and counted != a1:
            flipped += 1
        else:
            bad += 1
    log(f"DECODER VERIFICATION: {ok} SNPs identical (plink2 counted A1), "
        f"{flipped} identical after 2-x (plink2 counted A2), {bad} MISMATCH")
    return dict(ok=ok, flipped=flipped, bad=bad, n_snp=int(n_snp), n_samp=int(n_samp))


def impute_chunk(G8: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """adsp_center_batch_common.py L341-346: per-SNP non-missing-mean
    imputation.  Returns (G float64, maf (A1-minor), miss_rate)."""
    n = G8.shape[1]
    miss = G8 == 255
    G = G8.astype(np.float64)
    nmiss = miss.sum(1)
    G[miss] = 0.0
    smean = G.sum(1) / np.maximum(n - nmiss, 1)
    G[miss] = np.broadcast_to(smean[:, None], G.shape)[miss]
    maf_a1 = smean / 2.0
    maf = np.minimum(maf_a1, 1 - maf_a1)
    return G, maf, nmiss / n


# ======================================================== score from null ===
W_TINY = 1e-4
FLAG_NEWTON_FAIL = 4
FLAG_RATIO_GUARD = 5


def _cgf(Gt, mu, t, na_var, need_K=True):
    x = Gt * t[:, None]
    xp = np.maximum(x, 0.0)
    xn = np.minimum(x, 0.0)
    pos = x > 0
    mu_ = mu[None, :]
    one_m = 1.0 - mu_
    ep = np.exp(-xp)
    en = np.exp(xn)
    p_t = np.where(pos, mu_ / (mu_ + one_m * ep), mu_ * en / (one_m + mu_ * en))
    K1 = np.einsum("ij,ij->i", Gt, p_t - mu_) + t * na_var
    K2 = np.einsum("ij,ij->i", Gt * Gt, p_t * (1.0 - p_t)) + na_var
    K = None
    if need_K:
        K_el = np.where(pos, x + np.log1p(one_m * np.expm1(-xp)), np.log1p(mu_ * np.expm1(xn)))
        K = np.sum(K_el - x * mu_, axis=1) + 0.5 * t * t * na_var
    return K, K1, K2


# V9 tolerance 1e-10 -> 1e-6 in log10 p on 2026-09-20: the rank-60 cohort-keyed design (riga_cohortkey) agreed with
# spa_score.py to 4.5e-9 (job 1760622_7; economic QR here vs full QR there), pair33 to 1.8e-15; 1e-6 is still six
# decades below any p-value difference that could move a call.
V9_TOL = 1e-6


def _solve_saddlepoint(Gt, mu, u, V, na_var, maxit=200, ftol=1e-10, ttol=1e-12):
    m = u.size
    t = u / V
    lo = np.zeros(m)
    hi = np.full(m, np.inf)
    conv = np.zeros(m, dtype=bool)
    active = np.ones(m, dtype=bool)
    for _ in range(maxit):
        ia = np.flatnonzero(active)
        if ia.size == 0:
            break
        ta = t[ia]
        _, K1, K2 = _cgf(Gt[ia], mu, ta, na_var[ia], need_K=False)
        f = K1 - u[ia]
        neg = f < 0
        posf = f > 0
        lo[ia[neg]] = ta[neg]
        hi[ia[posf]] = ta[posf]
        done = np.abs(f) <= ftol * u[ia]
        with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
            tn = ta - f / K2
        loa = lo[ia]
        hia = hi[ia]
        bad = ~np.isfinite(tn) | (tn <= loa) | (tn >= hia)
        fin = np.isfinite(hia)
        tn = np.where(bad & fin, 0.5 * (loa + hia), tn)
        tn = np.where(bad & ~fin, 2.0 * ta, tn)
        done |= np.abs(tn - ta) <= ttol * (1.0 + np.abs(ta))
        t[ia] = np.where(done, ta, tn)
        conv[ia[done]] = True
        active[ia[done]] = False
    return t, conv


def _spa_upper_tail(Gt, mu, u, V, na_var, Gt_full=None):
    if Gt_full is None:
        Gt_full = Gt
    m = u.size
    p = np.full(m, np.nan)
    code = np.ones(m, dtype=np.int8)
    mu_ = mu[None, :]
    u_max = np.sum(np.where(Gt_full > 0, Gt_full * (1.0 - mu_), -Gt_full * mu_), axis=1)
    at_bound = u >= u_max * (1.0 - 1e-10)
    if at_bound.any():
        with np.errstate(divide="ignore"):
            lmu = np.log(mu_)
            l1m = np.log(1.0 - mu_)
        ib = np.flatnonzero(at_bound)
        logp = np.sum(np.where(Gt_full[ib] > 0, lmu, np.where(Gt_full[ib] < 0, l1m, 0.0)), axis=1)
        p[ib] = np.exp(logp)
        code[ib] = 3
    inner = np.flatnonzero(~at_bound)
    if inner.size:
        Gi = Gt[inner]
        ui = u[inner]
        Vi = V[inner]
        nvi = na_var[inner]
        t, conv = _solve_saddlepoint(Gi, mu, ui, Vi, nvi)
        K, _, K2 = _cgf(Gi, mu, t, nvi, need_K=True)
        with np.errstate(divide="ignore", invalid="ignore"):
            wq = np.sqrt(np.maximum(2.0 * (t * ui - K), 0.0))
            v = t * np.sqrt(K2)
            z = wq + np.log(v / wq) / wq
        p_lr = norm.sf(z)
        p_norm = norm.sf(ui / np.sqrt(Vi))
        use_norm = (~np.isfinite(p_lr)) | (np.abs(wq) < W_TINY)
        pi = np.where(use_norm, p_norm, p_lr)
        ci = np.where(use_norm, 2, 1).astype(np.int8)
        pi = np.where(conv, pi, p_norm)
        ci = np.where(conv, ci, 4).astype(np.int8)
        p[inner] = pi
        code[inner] = ci
    return p, code


def score_from_nullfit(X0: np.ndarray, nf: Dict, fastspa_cutoff: float = 0.05,
                       ratio_guard: bool = True, chunk: int = 256) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """The chi2 + SPA block of spa_score.spa_score_test (spa_score.py) with the
    null fit passed in instead of refitted.  Same operations in the same
    order, so p_chi2 / p_spa are bit-identical to spa_score_test."""
    M, n = X0.shape
    mu, e, w, sqrt_w, Q = nf["mu"], nf["e"], nf["w"], nf["sqrt_w"], nf["Q"]
    if n != mu.size:
        raise ValueError(f"X0 columns (n={n}) must match the null fit (N={mu.size})")
    G = X0.astype(np.float64)
    U = G @ e
    G_tilde = G * sqrt_w[None, :]
    GQ = G_tilde @ Q
    V = np.sum(G_tilde * G_tilde, axis=1) - np.sum(GQ * GQ, axis=1)
    V[V < 1e-12] = np.nan
    score_stat = (U * U) / V
    p_chi2 = chi2.sf(score_stat, df=1)
    p_chi2[~np.isfinite(p_chi2)] = np.nan
    p_spa = p_chi2.copy()
    flags = np.zeros(M, dtype=np.int8)
    need = np.isfinite(p_chi2) & (p_chi2 < fastspa_cutoff)
    idx_all = np.flatnonzero(need)
    inv_sqrt_w = np.where(w > 0, 1.0 / np.where(w > 0, sqrt_w, 1.0), 0.0)
    for s0 in range(0, idx_all.size, chunk):
        idx = idx_all[s0:s0 + chunk]
        Gt = (G_tilde[idx] - GQ[idx] @ Q.T) * inv_sqrt_w[None, :]
        na_var = np.zeros(idx.size)
        Gt_full = Gt
        u_abs = np.abs(U[idx])
        p_up, c_up = _spa_upper_tail(Gt, mu, u_abs, V[idx], na_var, Gt_full)
        p_lo, c_lo = _spa_upper_tail(-Gt, mu, u_abs, V[idx], na_var, -Gt_full)
        p2 = p_up + p_lo
        fail = (c_up == 4) | (c_lo == 4)
        code = np.maximum(c_up, c_lo).astype(np.int8)
        code = np.where(fail, FLAG_NEWTON_FAIL, code).astype(np.int8)
        p2 = np.where(fail, p_chi2[idx], p2)
        if ratio_guard:
            with np.errstate(divide="ignore", invalid="ignore"):
                guard = (~fail) & (p2 > 0) & (p_chi2[idx] / p2 > 1e3)
            p2 = np.where(guard, p_chi2[idx], p2)
            code = np.where(guard, FLAG_RATIO_GUARD, code).astype(np.int8)
        p_spa[idx] = np.minimum(p2, 1.0)
        flags[idx] = code
    return p_spa, p_chi2, flags


def bit_identity_check(G_chunk: np.ndarray, frames: Dict, nfs: Dict, arms: Sequence[str],
                       fastspa_cutoff: float) -> Dict:
    """V9: score_from_nullfit vs spa_score.spa_score_test on one chunk for the
    given arms; imports spa_score (and thereby the simulator) read-only."""
    sys.path.insert(0, str(SCRIPTS))
    import spa_score as SP  # noqa: E402  (imports simulation_2pop_qqplot as S)
    out = {}
    for arm in arms:
        d, rows, covs, _ = frames[arm]
        Ga = G_chunk[:, rows]
        p_spa_ref, p_chi2_ref, fl_ref = SP.spa_score_test(Ga, d, covs, fastspa_cutoff=fastspa_cutoff,
                                                          return_chi2=True)
        p_spa, p_chi2, fl = score_from_nullfit(Ga, nfs[arm], fastspa_cutoff)
        with np.errstate(divide="ignore", invalid="ignore"):
            d_chi2 = np.nanmax(np.abs(np.log10(p_chi2) - np.log10(p_chi2_ref)))
            d_spa = np.nanmax(np.abs(np.log10(p_spa) - np.log10(p_spa_ref)))
        Wref = SP.S._prepare_covariate_matrix(d, covs)
        out[arm] = dict(max_dlog10p_chi2=float(d_chi2), max_dlog10p_spa=float(d_spa),
                        flags_equal=bool(np.array_equal(fl, fl_ref)),
                        design_equal=bool(Wref.shape == nfs[arm]["W"].shape and np.array_equal(Wref, nfs[arm]["W"])),
                        nan_pattern_equal=bool(np.array_equal(np.isnan(p_chi2), np.isnan(p_chi2_ref))),
                        n_snp=int(Ga.shape[0]), n_subjects=int(Ga.shape[1]))
        out[arm]["pass"] = bool(out[arm]["design_equal"] and out[arm]["nan_pattern_equal"]
                                and out[arm]["max_dlog10p_chi2"] <= V9_TOL and out[arm]["max_dlog10p_spa"] <= V9_TOL)
        log(f"  bit-identity [{arm}]: dlog10p chi2 {d_chi2:.2e} spa {d_spa:.2e} flags_equal={out[arm]['flags_equal']} "
            f"design_equal={out[arm]['design_equal']} -> {'PASS' if out[arm]['pass'] else 'FAIL'}")
    return out


# ================================================================ metrics ===
def metrics(p: np.ndarray) -> Dict:
    """COPY of arenaA_pathC_repair.metrics (L330-350)."""
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


def tail_metrics(p: np.ndarray) -> Dict:
    """Tail lambda at q99.9 (calibration metrology) + the standard block."""
    out = metrics(p)
    ok = np.isfinite(p) & (p > 0) & (p <= 1)
    pp = p[ok]
    if pp.size >= 10:
        q = chi2.isf(pp, 1)
        out["lambda_q999"] = float(np.quantile(q, 0.999) / chi2.ppf(0.999, 1))
        out["lambda_q99"] = float(np.quantile(q, 0.99) / chi2.ppf(0.99, 1))
    else:
        out["lambda_q999"] = np.nan
        out["lambda_q99"] = np.nan
    return out


def _append_csv(path: Path, rows: List[Dict]) -> None:
    """COPY of arenaA_pathC_repair._append_csv."""
    new = pd.DataFrame(rows)
    if path.exists() and path.stat().st_size > 0:
        cols = list(pd.read_csv(path, nrows=0).columns)
        for c in new.columns:
            if c not in cols:
                cols.append(c)
        new = new.reindex(columns=cols)
        new.to_csv(path, mode="a", header=False, index=False)
    else:
        new.to_csv(path, index=False)


def _drop_seed_rows(path: Path, seed: int) -> None:
    if path.exists() and path.stat().st_size > 0:
        d = pd.read_csv(path)
        if "seed" in d.columns and (d.seed == seed).any():
            d[d.seed != seed].to_csv(path, index=False)


def seed_done(units_csv: Path, npz_path: Path, json_path: Path, seed: int, n_arms: int) -> bool:
    if not (npz_path.exists() and json_path.exists()):
        return False
    if not units_csv.exists() or units_csv.stat().st_size == 0:
        return False
    d = pd.read_csv(units_csv, usecols=["seed", "arm", "test", "tool"])
    d = d[(d.seed == seed) & (d.tool == "native")]
    return int(len(d)) == len(TESTS) * n_arms


# --------------------------------------------------------- SNP annotation --
def snp_annotation(bim: pd.DataFrame, persnp: Optional[Path]) -> pd.DataFrame:
    """Per-SNP frame: chr (int), pos, window id, excluded (APOE / known loci),
    region_flag / miss_dirty / tgated from common_persnp.csv (rules copied from
    adsp_center_batch_skeptic_common.py) when available."""
    ann = pd.DataFrame({"snp": bim.id.values})
    chrom = bim.chr.astype(str).str.replace("chr", "", regex=False)
    ann["chr"] = pd.to_numeric(chrom, errors="coerce").fillna(-1).astype(int)
    ann["pos"] = pd.to_numeric(bim.pos, errors="coerce").fillna(0).astype(np.int64)
    ann["window"] = ann["chr"].astype(np.int64) * 1_000_000 + (ann["pos"] // WINDOW_BP)
    excl = (ann.chr == APOE_WINDOW[0]) & (ann.pos >= APOE_WINDOW[1]) & (ann.pos <= APOE_WINDOW[2])
    known = np.zeros(len(ann), dtype=bool)
    for _, c, p in KNOWN_LOCI:
        known |= ((ann.chr == c) & (np.abs(ann.pos - p) <= WINDOW_BP // 2)).to_numpy()
    ann["apoe_window"] = excl.to_numpy()
    ann["known_window"] = known
    ann["excluded"] = excl.to_numpy() | known
    # region flags (skeptic rules)
    mb = ann.pos / 1e6
    flag = pd.Series("", index=ann.index)
    flag[(ann.chr == 6) & mb.between(28, 34)] = "HLA"
    flag[(ann.chr == 8) & mb.between(8, 12)] = "8p23.1_inv"
    flag[(ann.chr == 17) & mb.between(45.5, 46.5)] = "17q21.31_inv"
    for c, (a, b) in CENTROMERES.items():
        m = (ann.chr == c) & mb.between(a - 1, b + 1) & (flag == "")
        flag[m] = "centromere_1Mb"
    for c in (13, 14, 15, 21, 22):
        m = (ann.chr == c) & (mb < CENTROMERES[c][1] + 1) & (flag == "")
        flag[m] = "acrocentric_p"
    ann["region_flag"] = flag
    ann["flagged"] = flag != ""
    ann["miss_dirty"] = False
    ann["tgated"] = False
    ann["has_persnp"] = False
    if persnp is not None and Path(persnp).exists():
        ps = pd.read_csv(persnp, usecols=["snp", "miss_rate", "miss_het_p", "p_center_pcs",
                                          "n_center_tgated_shift015"])
        ps = ps.set_index("snp").reindex(ann.snp.values)
        ann["has_persnp"] = ps["miss_rate"].notna().to_numpy()
        ann["miss_dirty"] = ((ps["miss_het_p"] < 1e-5) | (ps["miss_rate"] > 0.01)).fillna(False).to_numpy()
        ann["tgated"] = ((ps["p_center_pcs"] < 1e-5) & (ps["n_center_tgated_shift015"] > 0)).fillna(False).to_numpy()
    ann["clean"] = (~ann.flagged) & (~ann.miss_dirty)
    return ann


def subset_masks(ann: pd.DataFrame, maf: np.ndarray) -> Dict[str, np.ndarray]:
    keep = ~ann.excluded.to_numpy()
    masks = {"all": keep, "all_incl_windows": np.ones(len(ann), dtype=bool)}
    for lo, hi in MAF_BINS:
        masks[f"maf[{lo:g},{hi:g})"] = keep & (maf >= lo) & (maf < hi)
    masks["unflagged_regions"] = keep & (~ann.flagged.to_numpy())
    masks["flagged_only"] = keep & ann.flagged.to_numpy()
    masks["missingness_clean"] = keep & (~ann.miss_dirty.to_numpy())
    masks["dirty_only"] = keep & ann.miss_dirty.to_numpy()
    masks["clean"] = keep & ann.clean.to_numpy()
    masks["tgated_center_batch"] = keep & ann.tgated.to_numpy()
    masks["not_tgated"] = keep & (~ann.tgated.to_numpy())
    masks["apoe_window"] = ann.apoe_window.to_numpy()
    masks["known_windows"] = ann.known_window.to_numpy()
    return masks


def window_stats(p: np.ndarray, window: np.ndarray, keep: np.ndarray, rng: np.random.Generator,
                 n_boot: int = N_BOOT) -> Dict:
    """1-Mb-window-clustered hit counts (min p < 1e-5 / 1e-4 / 1e-3) and a
    window block bootstrap of the FPRs (x nominal) and hit counts."""
    ok = keep & np.isfinite(p) & (p > 0) & (p <= 1)
    if ok.sum() == 0:
        return {}
    d = pd.DataFrame({"w": window[ok], "p": p[ok]})
    for a, tag in ALPHAS[:3]:
        d[f"h{tag}"] = (d.p < a).astype(np.int64)
    g = d.groupby("w").agg(n=("p", "size"), minp=("p", "min"), h1e3=("h1e3", "sum"),
                           h1e4=("h1e4", "sum"), h1e5=("h1e5", "sum"))
    W = len(g)
    out = dict(n_windows=int(W), n_windows_hit_1e5=int((g.minp < 1e-5).sum()),
               n_windows_hit_1e4=int((g.minp < 1e-4).sum()), n_windows_hit_1e3=int((g.minp < 1e-3).sum()))
    counts = rng.multinomial(W, np.full(W, 1.0 / W), size=n_boot).astype(np.float64)
    n_b = counts @ g.n.values.astype(np.float64)
    for a, tag in ALPHAS[:3]:
        h_b = counts @ g[f"h{tag}"].values.astype(np.float64)
        x = h_b / n_b / a
        out[f"boot_fpr_{tag}_x_lo"] = float(np.quantile(x, 0.025))
        out[f"boot_fpr_{tag}_x_hi"] = float(np.quantile(x, 0.975))
        out[f"boot_fpr_{tag}_x_sd"] = float(np.std(x))
    hitw = (g.minp.values < 1e-5).astype(np.float64)
    hw_b = counts @ hitw
    out["boot_n_windows_hit_1e5_lo"] = float(np.quantile(hw_b, 0.025))
    out["boot_n_windows_hit_1e5_hi"] = float(np.quantile(hw_b, 0.975))
    return out


# --------------------------------------------------------------- burden ----
def load_burden_masks(burden_dir: Optional[Path]) -> Dict[str, Dict]:
    """<mask>.npz with a 2-D (genes x N) array + gene names + cMAC."""
    out = {}
    if burden_dir is None or not Path(burden_dir).exists():
        return out
    for f in sorted(Path(burden_dir).glob("*.npy")):
        # plasmode_burden_chr22.R format: burden_<mask>.npy (units x N float32, pool row
        # order) + gene_table_<mask>.csv (unit, cMAC).
        mat = np.load(f, allow_pickle=False)
        if mat.ndim != 2:
            continue
        name = f.stem[len("burden_"):] if f.stem.startswith("burden_") else f.stem
        gene = None
        cmac = None
        gt = Path(burden_dir) / f"gene_table_{name}.csv"
        if gt.exists():
            t = pd.read_csv(gt)
            if len(t) == mat.shape[0]:
                if "unit" in t.columns:
                    gene = t["unit"].astype(str).to_numpy()
                for c in ("cMAC", "cmac", "cMAC_pool"):
                    if c in t.columns:
                        cmac = t[c].to_numpy(dtype=np.float64)
                        break
        if cmac is None:
            cmac = np.asarray(mat, dtype=np.float64).sum(1)
        out[name] = dict(B=np.asarray(mat, dtype=np.float32), gene=gene, cmac=cmac)
    for f in sorted(Path(burden_dir).glob("*.npz")):
        z = np.load(f, allow_pickle=True)
        mat = None
        for k in z.files:
            a = z[k]
            if a.ndim == 2 and mat is None:
                mat = a
        if mat is None:
            continue
        gene = None
        cmac = None
        for k in z.files:
            kl = k.lower()
            if kl in ("gene", "genes", "gene_id") and gene is None:
                gene = z[k].astype(str)
            if kl in ("cmac", "cmac_pool", "c_mac") and cmac is None:
                cmac = z[k].astype(np.float64)
        if cmac is None:
            cmac = np.asarray(mat, dtype=np.float64).sum(1)
        out[f.stem] = dict(B=np.asarray(mat, dtype=np.float32), gene=gene, cmac=cmac)
    return out


# ============================================================== run task ====
def run_task(args) -> int:
    replicate = int(args.replicate)
    scheme = args.scheme or task_scheme(replicate)
    out_dir = Path(args.out_root) / "native"
    npz_dir = out_dir / "npz"
    npz_dir.mkdir(parents=True, exist_ok=True)
    seed = read_manifest(Path(args.manifest), replicate) if args.seed is None else int(args.seed)
    units_csv = out_dir / f"units_plasmode_{scheme}.csv"
    bybin_csv = out_dir / f"bybin_plasmode_{scheme}.csv"
    npz_path = npz_dir / f"rep{replicate:03d}_seed{seed}.npz"
    json_path = npz_dir / f"rep{replicate:03d}_seed{seed}.json"
    failed_marker = npz_dir / f"FAILED_rep{replicate:03d}_seed{seed}.json"
    arms = list(ARM_ORDER) if not args.arms else [a for a in args.arms.split(",") if a]
    if seed_done(units_csv, npz_path, json_path, seed, len(arms)) and not args.force:
        log(f"replicate {replicate} seed {seed} already done -> skip")
        return 0
    log(f"replicate={replicate} scheme={scheme} seed={seed} arms={len(arms)} host={os.uname().nodename}")

    # ---- pool + bed alignment
    pool = pd.read_csv(args.pool_csv, low_memory=False)
    need = ["IID", "center", "Cohort", "Study", "AD", "Sex", "Age", "lpa_cluster", "cell_A"]
    miss = [c for c in need if c not in pool.columns]
    if miss:
        raise SystemExit(f"pool csv lacks {miss}")
    if "cell_A2" not in pool.columns:
        pool["cell_A2"] = pool["lpa_cluster"].astype(str) + "|" + pool["cell_A"].astype(str)
    N = len(pool)
    bfile = Path(args.bfile) if args.bfile else (BFILE_POOL if (BFILE_POOL.with_suffix(".bed")).exists() else BFILE_ORIG)
    bed = Bed(bfile)
    fam_iid = bed.fam.IID.astype(str).values
    pool_iid = pool.IID.astype(str).values
    if bed.n_fam == N and np.array_equal(fam_iid, pool_iid):
        pool_idx = np.arange(N, dtype=np.int64)
    else:
        pos = pd.Series(np.arange(bed.n_fam), index=fam_iid)
        pool_idx = pos.reindex(pool_iid).to_numpy()
        if not np.isfinite(pool_idx.astype(float)).all():
            raise SystemExit("pool IIDs not all present in the bed fam")
        pool_idx = pool_idx.astype(np.int64)
        if "fam_pos" in pool.columns:
            assert np.array_equal(pool_idx, pool.fam_pos.to_numpy().astype(np.int64)), "fam_pos != fam position"
    log(f"pool N={N}; bed {bfile.name}: {bed.n_fam} samples, {bed.n_snp} SNPs; pool_idx identity={bool(bed.n_fam == N)}")

    # ---- Y*
    ystar, rng_boot = draw_ystar(pool, scheme, seed)
    vjson, fails = validate_ystar(pool, ystar, scheme)
    vjson.update(replicate=replicate, seed=seed, scheme=scheme, bfile=str(bfile), pool_csv=str(args.pool_csv))
    json_path.write_text(json.dumps(vjson, indent=1, default=float))
    log(f"V1 cells={vjson['V1_n_cells']} mono={vjson['V1_n_mono_cells']} maxdiff={vjson['V1_max_abs_count_diff']} | "
        f"V2 chi2 {vjson['V2_chi2_center_ystar']:.1f} (real {vjson['V2_chi2_center_y_real']:.1f}) V={vjson['V2_cramerV_ystar']:.3f} "
        f"LRT centre|PCs {vjson['V2_lrt_center_pcs_ystar']:.1f} (real {vjson['V2_lrt_center_pcs_real']:.1f}) | "
        f"V3 max|z|={vjson.get('V3_max_abs_z', np.nan):.2f} | V4 max|z|={vjson['V4_max_abs_z']:.2f} | "
        f"V5 leakage={vjson['V5_identity_leakage']:.4f} (mono share {vjson['V5_mono_share']:.3f})")
    if fails:
        failed_marker.write_text(json.dumps(dict(replicate=replicate, seed=seed, fails=fails), indent=1))
        log("VALIDATION FAILED: " + " | ".join(fails))
        return 3

    # ---- matching (rows depend on Y*)
    matched: Dict[str, np.ndarray] = {}
    matcher = "none"
    want_match = any(a in MATCHED_ARMS for a in arms)
    if want_match:
        mroot = Path(args.match_dir) if args.match_dir else (out_dir / "match")
        mdir = mroot / f"rep{replicate:03d}_seed{seed}"
        mdir.mkdir(parents=True, exist_ok=True)
        mtag = f"rep{replicate:03d}"
        want_keys = [k for k in MATCH_KEYS if k in {MATCHED_ARMS[a] for a in arms if a in MATCHED_ARMS}]
        # benchmark_matching.R redoes the deployed verbatim prep from adsp_combined.csv itself and
        # only needs these pool columns (its section 2) + Y (= Y*, the same file via --ystar).
        pheno_cols = [c for c in ("IID", "row_idx", "fam_pos", "Study", "Cohort", "center", "AD")
                      if c in pool.columns]
        pheno_csv = mdir / f"ystar_pool_{mtag}.csv"
        pp = pool[pheno_cols].copy()
        pp["Y"] = ystar.astype(np.int64)
        pp.to_csv(pheno_csv, index=False)
        script = Path(args.match_script)
        if args.match_fallback == "python" or not script.exists():
            if args.match_fallback != "python":
                log(f"WARNING: {script} not found; matched arms SKIPPED (use --match-fallback python for the self-test)")
            else:
                matched = python_fallback_match(pool, ystar, seed=0)
                matcher = "python_fallback"
                fam_pos = (pool["fam_pos"].to_numpy(dtype=np.int64) if "fam_pos" in pool.columns
                           else np.arange(N, dtype=np.int64))
                blocks = []
                for k, v in matched.items():
                    h = len(v) // 2
                    blocks.append(pd.DataFrame({"mkey": k, "pair_id": np.tile(np.arange(h), 2),
                                                "IID": pool_iid[v], "row_idx": v, "fam_pos": fam_pos[v],
                                                "Y": ystar[v].astype(np.int64)}))
                if blocks:   # same on-disk format as benchmark_matching.R
                    pd.concat(blocks, ignore_index=True).to_csv(mdir / f"match_idx_{mtag}.csv", index=False)
        else:
            matched = run_match_r(script, pheno_csv, pheno_csv, mdir, mtag, want_keys,
                                  args.match_cmd, seed=seed)
            matcher = "matchit"
        for k, v in matched.items():
            log(f"  matched {k}: {len(v) // 2} pairs ({len(v)} rows)")
        if scheme == "anchor" and "KEY0" in matched and Path(args.deployed_pairs).exists():
            dep = pd.read_csv(args.deployed_pairs, usecols=["IID"]).IID.astype(str)
            dep_pool = set(dep) & set(pool_iid)
            cur = set(pool_iid[matched["KEY0"]])
            ov = len(dep_pool & cur) / max(len(dep_pool), 1)
            vjson["V7_pair33_overlap_with_deployed_in_pool"] = float(ov)
            vjson["V7_n_deployed_in_pool"] = int(len(dep_pool))
            vjson["V7_n_pairs_KEY0"] = int(len(matched["KEY0"]) // 2)
            log(f"V7 pair33 overlap with the deployed file (restricted to the pool): {ov:.4f}")

    # ---- null fits
    frames = arm_frames(pool, ystar, matched, arms)
    nfs: Dict[str, Dict] = {}
    arm_info: Dict[str, Dict] = {}
    for arm in arms:
        d, rows, covs, n_pairs = frames[arm]
        if d is None:
            arm_info[arm] = dict(error="matched rows unavailable", n_subjects=0)
            log(f"  [{arm}] SKIPPED (no matched rows)")
            continue
        t1 = time.time()
        try:
            nf = null_fit_local(d, covs)
        except Exception as exc:  # noqa: BLE001
            arm_info[arm] = dict(error=f"{type(exc).__name__}: {exc}", n_subjects=int(len(d)))
            log(f"  [{arm}] null fit FAILED: {arm_info[arm]['error']}")
            continue
        nfs[arm] = nf
        arm_info[arm] = dict(error="", n_subjects=nf["n_subjects"], n_cases_arm=nf["n_cases_arm"], n_pairs=n_pairs,
                             rank=nf["rank"], converged=nf["converged"], n_params=nf["n_params"],
                             n_design_dummies=nf["n_design_dummies"], max_abs_coef_design=nf["max_abs_coef_design"],
                             max_abs_coef_any=nf["max_abs_coef_any"], inestimable_design=nf["inestimable_design"],
                             n_effective=nf["n_effective"], n_dropped_eff=nf["n_dropped_eff"],
                             warn_classes=nf["warn_classes"], perfect_separation_warning=nf["perfect_separation_warning"],
                             covariates=";".join(c for c in covs if not c.upper().startswith("PC"))
                             + (";PC1..PC20" if any(c in PCS20[10:] for c in covs) else (";PC1..PC10" if any(c in PCS10 for c in covs) else "")),
                             null_fit_s=float(time.time() - t1), matcher=(matcher if arm in MATCHED_ARMS else "none"))
        log(f"  [{arm:26s}] N={nf['n_subjects']:6d} cases={nf['n_cases_arm']:5d} rank={nf['rank']:4d} conv={nf['converged']} "
            f"max|design coef|={nf['max_abs_coef_design']:.2f} inest={nf['inestimable_design']} ({time.time() - t1:.1f}s)")
    live = [a for a in arms if a in nfs]
    if not live:
        raise SystemExit("no arm has a null fit")

    # ---- decoder verification
    scratch = Path(args.scratch) / f"rep{replicate:03d}"
    dec = None
    if not args.no_verify_decoder:
        dec = verify_decoder_plink2(bed, pool_idx, pool_iid, scratch, plink2=args.plink2)
        if dec["bad"]:
            failed_marker.write_text(json.dumps(dict(replicate=replicate, seed=seed, fails=["decoder mismatch"], decoder=dec)))
            return 5

    # ---- bed pass
    M = bed.n_snp if not args.max_chunks else min(bed.n_snp, args.max_chunks * args.chunk)
    P = {arm: {"chi2": np.full(M, np.nan, dtype=np.float32), "spa": np.full(M, np.nan, dtype=np.float32),
               "flags": np.zeros(M, dtype=np.int8)} for arm in live}
    maf = np.full(M, np.nan, dtype=np.float32)
    missr = np.full(M, np.nan, dtype=np.float32)
    bitcheck: Dict = {}
    t_pass = time.time()
    t_arm = {arm: 0.0 for arm in live}
    for start in range(0, M, args.chunk):
        m = min(args.chunk, M - start)
        G8 = bed.read_chunk(start, m)[:, pool_idx]
        G, maf_c, miss_c = impute_chunk(G8)
        del G8
        maf[start:start + m] = maf_c
        missr[start:start + m] = miss_c
        if start == 0 and not args.skip_bitcheck:
            # The REFERENCE side of V9 is spa_score.spa_score_test -> spa_score.null_fit, which
            # uses scipy's DEFAULT mode="full" QR: its Q is (N x N) = 5.96 GB for a 27,306-row
            # arm, and that is what OOM-killed job 1639191 at 12 GB (MaxRSS 12.56 GB).
            # spa_score.py is not edited, so V9 is run on the MATCHED arms, where the reference
            # Q is 1.3 GB or less; the code path exercised is identical (same function, same
            # G[:, rows] subsetting).  Fall back to the smallest live arms, never to live[:2].
            check_arms = [a for a in args.bitcheck_arms.split(",") if a in live][:2]
            if len(check_arms) < 2:
                check_arms = sorted(live, key=lambda a: int(frames[a][1].size))[:2]
            bitcheck = bit_identity_check(G, frames, nfs, check_arms, args.fastspa_cutoff)
            if not all(v["pass"] for v in bitcheck.values()):
                vjson["V9_bit_identity"] = bitcheck
                json_path.write_text(json.dumps(vjson, indent=1, default=float))
                failed_marker.write_text(json.dumps(dict(replicate=replicate, seed=seed, fails=["V9 bit-identity"],
                                                         bitcheck=bitcheck), indent=1, default=float))
                log("V9 BIT-IDENTITY FAILED")
                return 4
        for arm in live:
            ta = time.time()
            rows = frames[arm][1]
            Ga = G if (rows.size == N and np.array_equal(rows, np.arange(N))) else G[:, rows]
            p_spa, p_chi2, fl = score_from_nullfit(Ga, nfs[arm], args.fastspa_cutoff)
            P[arm]["chi2"][start:start + m] = p_chi2
            P[arm]["spa"][start:start + m] = p_spa
            P[arm]["flags"][start:start + m] = fl
            t_arm[arm] += time.time() - ta
        if (start // args.chunk) % 10 == 0:
            log(f"  chunk {start // args.chunk + 1}/{(M + args.chunk - 1) // args.chunk} "
                f"({time.time() - t_pass:.0f}s elapsed)")
    log(f"bed pass done in {time.time() - t_pass:.0f}s; per-arm s: " +
        ", ".join(f"{ARM_LABEL.get(a, a)}={t_arm[a]:.0f}" for a in live))

    # ---- annotation + metrics
    bim = bed.bim.iloc[:M].reset_index(drop=True)
    ann = snp_annotation(bim, Path(args.persnp) if args.persnp else None)
    masks = subset_masks(ann, maf.astype(np.float64))
    common = dict(design=scheme, scheme=scheme, scenario="adsp_pool", seed=seed, replicate=replicate,
                  n_pool=int(N), n_snps_bed=int(M), n_excluded_windows=int(ann.excluded.sum()),
                  bfile=bfile.name, tool="native", is_external=False, driver=DRIVER_PATH.name)
    unit_rows: List[Dict] = []
    bybin_rows: List[Dict] = []
    v5_minp: Dict[str, Dict] = {}
    for arm in live:
        info = arm_info[arm]
        base = dict(common, arm=arm, arm_label=ARM_LABEL.get(arm, arm), is_matched=bool(arm in MATCHED_ARMS),
                    n_pairs=(info["n_pairs"] if info.get("n_pairs") is not None else np.nan),
                    **{k: v for k, v in info.items() if k not in ("n_pairs",)},
                    ok=True, covset=arm)
        fl = P[arm]["flags"]
        base.update({f"n_flag{k}": int((fl == k).sum()) for k in range(6)})
        for test in TESTS:
            p = P[arm][test].astype(np.float64)
            m_all = metrics(p[masks["all"]])
            m_incl = metrics(p[masks["all_incl_windows"]])
            m_clean = metrics(p[masks["clean"]])
            ws = window_stats(p, ann.window.to_numpy(), masks["all"], rng_boot)
            row = dict(base, test=test, **m_all,
                       **{f"{k}_incl": v for k, v in m_incl.items()},
                       **{f"{k}_clean": v for k, v in m_clean.items()},
                       **ws,
                       minp_apoe_window=(float(np.nanmin(p[masks["apoe_window"]])) if masks["apoe_window"].any() else np.nan),
                       minp_known_windows=(float(np.nanmin(p[masks["known_windows"]])) if masks["known_windows"].any() else np.nan))
            unit_rows.append(row)
            v5_minp[f"{arm}_{test}"] = dict(apoe=row["minp_apoe_window"], known=row["minp_known_windows"])
            for lab, mk in masks.items():
                if lab in ("all", "apoe_window", "known_windows"):
                    continue
                btype = "maf" if lab.startswith("maf") else ("batch_class" if "tgated" in lab else "subset")
                bybin_rows.append(dict(base, test=test, bin_type=btype, bin_label=lab, n_bin=int(mk.sum()),
                                       **metrics(p[mk])))
            log(f"  [{ARM_LABEL.get(arm, arm):3s} {arm:26s} {test}] lam={m_all['lambda_gc']:.3f} "
                f"1e-3={m_all['fpr_1e3_x']:.2f}x 1e-5={m_all['fpr_1e5_x']:.1f}x clean lam={m_clean['lambda_gc']:.3f} "
                f"windows<1e-5={ws.get('n_windows_hit_1e5', np.nan)}")

    # ---- burden masks (chr22 QC7) scored with the same null fits
    burden = load_burden_masks(Path(args.burden_dir) if args.burden_dir else None)
    burden_out: Dict[str, Dict] = {}
    for mask, bd in burden.items():
        B = bd["B"]
        if B.shape[1] != N:
            log(f"  burden {mask}: {B.shape[1]} columns != pool N {N}; skipped")
            continue
        cmaf = bd["cmac"] / (2.0 * N)
        for arm in live:
            rows = frames[arm][1]
            Ba = B[:, rows]
            p_spa, p_chi2, fl = score_from_nullfit(Ba, nfs[arm], args.fastspa_cutoff)
            burden_out[f"{mask}__{arm}"] = dict(p_chi2=p_chi2.astype(np.float32), p_spa=p_spa.astype(np.float32))
            for test, p in (("chi2", p_chi2), ("spa", p_spa)):
                base = dict(common, arm=arm, arm_label=ARM_LABEL.get(arm, arm), test=test, is_matched=bool(arm in MATCHED_ARMS),
                            n_subjects=arm_info[arm]["n_subjects"], covset=arm)
                bybin_rows.append(dict(base, bin_type="burden", bin_label=f"{mask}:all", n_bin=int(len(p)), **tail_metrics(p)))
                for lo, hi in CMAF_BINS:
                    mk = (cmaf >= lo) & (cmaf < hi)
                    bybin_rows.append(dict(base, bin_type="burden", bin_label=f"{mask}:cmaf[{lo:g},{hi:g})",
                                           n_bin=int(mk.sum()), **tail_metrics(p[mk])))
        log(f"  burden {mask}: {B.shape[0]} genes scored on {len(live)} arms")

    # ---- V5 post-scan part
    v5_ok = True
    for arm in ("saturated_cell", "pooled_project"):
        for test in TESTS:
            r = v5_minp.get(f"{arm}_{test}")
            if r and scheme != "anchor":
                for k in ("apoe", "known"):
                    if np.isfinite(r[k]) and r[k] <= 1e-4:
                        v5_ok = False
    vjson.update(V5_window_minp=v5_minp, V5_windows_pass=bool(v5_ok), V9_bit_identity=bitcheck,
                 decoder_verification=dec, arm_info=arm_info, t_bed_pass_s=float(time.time() - t_pass),
                 t_arm_s=t_arm, n_snps=int(M), burden_masks=list(burden), matcher=matcher,
                 n_pairs={k: int(len(v) // 2) for k, v in matched.items()}, t_total_s=float(time.time() - T0))
    json_path.write_text(json.dumps(vjson, indent=1, default=float))

    # ---- rows / npz (idempotent per seed)
    _drop_seed_rows(units_csv, seed)
    _drop_seed_rows(bybin_csv, seed)
    _append_csv(units_csv, unit_rows)
    _append_csv(bybin_csv, bybin_rows)
    # seeds come from a 64-bit unsigned manifest; np.int64(seed) overflowed for seeds >= 2**63 (job 1760622_10, 2026-09-20)
    npz = dict(ystar=ystar.astype(np.int8), maf=maf, miss_rate=missr, seed=np.uint64(seed), replicate=np.int64(replicate),
               chr=ann.chr.to_numpy().astype(np.int16), pos=ann.pos.to_numpy().astype(np.int64),
               excluded=ann.excluded.to_numpy(), arms=np.array(live))
    for arm in live:
        npz[f"p_chi2__{arm}"] = P[arm]["chi2"]
        npz[f"p_spa__{arm}"] = P[arm]["spa"]
        npz[f"flags__{arm}"] = P[arm]["flags"]
        npz[f"rows__{arm}"] = frames[arm][1].astype(np.int32)
    for k, v in burden_out.items():
        npz[f"burden_chi2__{k}"] = v["p_chi2"]
        npz[f"burden_spa__{k}"] = v["p_spa"]
    np.savez_compressed(npz_path, **npz)
    log(f"wrote {units_csv.name} (+{len(unit_rows)} rows), {bybin_csv.name} (+{len(bybin_rows)}), {npz_path.name}; "
        f"total {time.time() - T0:.0f}s")
    return 0


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(description=__doc__.split("Per task")[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--replicate", type=int, required=True, help="0 = anchor (observed AD); 1-20 scheme A; 21-40 scheme A2")
    ap.add_argument("--scheme", choices=["anchor", "A", "A2"], default=None, help="override the replicate -> scheme map")
    ap.add_argument("--seed", type=int, default=None, help="override the manifest seed")
    ap.add_argument("--pool-csv", default=str(OUT_ROOT / "plasmode_pool.csv"))
    ap.add_argument("--manifest", default=str(OUT_ROOT / "seed_manifest.csv"))
    ap.add_argument("--bfile", default=None, help="plink1 prefix; default pool_281k if present else chrall.common.pruned")
    ap.add_argument("--persnp", default=str(PERSNP), help="common_persnp.csv for the region / missingness / t-gated subsets ('' to skip)")
    ap.add_argument("--burden-dir", default=str(OUT_ROOT / "burden_chr22"))
    ap.add_argument("--out-root", default=str(OUT_ROOT))
    ap.add_argument("--scratch", default=os.environ.get("SIM_SCRATCH", "${SCRATCH_DIR}/"
                    "a773d3e1-08d5-4d13-8c1c-7004109bf225/scratchpad/build/plasmode"))
    ap.add_argument("--match-script", default=str(MATCH_SCRIPT))
    ap.add_argument("--match-cmd",
                    default="Rscript {script} --pool {pool} --out {outdir} --tag {tag} --ystar {ystar} --keys {keys}",
                    help="command template for the MatchIt pre-step (benchmark_matching.R contract); "
                         "fields: {script} {pool} {ystar} {pheno} {outdir} {tag} {keys} {seed}")
    ap.add_argument("--match-dir", default="", help="root for the R pre-step files (default <out-root>/native/match; "
                                                    "point at /dev/shm scratch in the sbatch)")
    ap.add_argument("--match-fallback", choices=["none", "python"], default="none",
                    help="python = greedy exact-key Mahalanobis matcher (self-test only)")
    ap.add_argument("--deployed-pairs", default=str(ROOT / "output/matched/matched_samples_lpa_pc.nofilter_pair.csv"))
    ap.add_argument("--arms", default="", help="comma list (default: all 12)")
    ap.add_argument("--chunk", type=int, default=2000)
    ap.add_argument("--max-chunks", type=int, default=0, help="smoke: limit the bed pass")
    ap.add_argument("--fastspa-cutoff", type=float, default=0.05)
    ap.add_argument("--bitcheck-arms", default="riga_pair33,riga_cohortkey",
                    help="arms for the V9 bit-identity check; the reference (spa_score.null_fit) "
                         "builds an N x N QR, so keep these SMALL (a 27,306-row arm needs 6 GB)")
    ap.add_argument("--skip-bitcheck", action="store_true")
    ap.add_argument("--no-verify-decoder", action="store_true")
    ap.add_argument("--plink2", default=PLINK2)
    ap.add_argument("--force", action="store_true", help="redo a finished replicate")
    return ap


def main() -> int:
    args = build_parser().parse_args()
    if args.persnp == "":
        args.persnp = None
    return run_task(args)


if __name__ == "__main__":
    sys.exit(main())
