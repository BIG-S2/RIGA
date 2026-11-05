#!/usr/bin/env python
"""
benchmark_harness.py -- plasmode POWER replicate on the REAL 281k LD-pruned ADSP
genotypes of the 27,306-subject unrelated pool (deliverable B's pool), with a known
causal architecture spiked into a phenotype that keeps the real Cohort x project x
centre structure.

Generator (one replicate = one SLURM task; seed from seed_manifest_power.csv)
-------------------------------------------------------------------------------
 1  causal single variants: drawn from the pruned set, 10 per (class x MAF bin x
    effect level); classes 'clean' (unflagged region, missingness-clean, NOT
    centre-batch t-gated, outside the APOE / known-locus windows, not in the 1,536-SNP
    genotype-class artifact list) and 'tgated' (the centre/project-differentiated
    class of common_persnp.csv, same other filters); MAF bins [0.05,0.1) [0.1,0.2)
    [0.2,0.5) on the pool allele counts; effect levels c = 0.06 / 0.09 / 0.12 log-odds
    PER SD OF GENOTYPE (beta_j = c / sd(G_j), random sign), so every MAF bin carries the
    same information; causal SNPs are >= 1 Mb apart.
 2  causal burden units: 8 per (mask x level) from the chr22 QC7 masks coding_mac1_5,
    coding_mac6_50, gene_mac6_50, gene_mac51_500 with 20 <= cMAC <= 1000, distinct genes;
    per-allele log-odds gamma = a * |log10(mean per-variant MAF in the unit)| with
    a = 0.1 / 0.2 / 0.3 (the SKAT convention, same sign for every allele of the unit).
 3  eta_i = sum_j beta_j G_ij + sum_u gamma_u B_ui.  --spike snp (default; the 281k scan,
    out subdir native) spikes the SNPs only; --spike burden spikes the burden units only
    and skips the scan (out subdir burden, ~5 min per replicate); both = one replicate.
 4  Y_i ~ Bernoulli(expit(alpha_cell + eta_i)) inside every Cohort x project x centre
    cell (cell_A) that has both cases and controls, alpha_cell solved so that the
    EXPECTED case count equals the observed one; cells with only cases or only
    controls keep their observed labels (they carry the confounding and no signal:
    conservative for the pooled designs) -- OR, with --mono-ascertain r, their subjects
    get a latent status from the same model and only those matching the project's
    label are kept (ascertainment emulation: retained subjects carry the case/control
    enrichment of a recruited sample, matched / stratified designs unchanged).  With
    eta = 0 this is the binomial version of deliverable B's within-cell permutation null.
 5  validation VP1-VP5 (per-cell counts, centre-AD dependence, within-cell technical
    channel: marginal VP4 is a soft report because technically driven causal dosages
    make Y track the metrics by construction, the hard gate VP4b conditions on eta;
    identity leakage; VP6 = realized vs expected information on the causal SNPs,
    post-scan) -> the replicate's json.
 6  matching: MatchIt 1:1 (benchmark_matching.R, keys KEY0 / KEY3 / KEY_P) + a greedy
    1:k nearest-Mahalanobis(PC1-10) matcher inside exact KEY_P strata (k = 1, 2, 3).
 7  arms (below), one bed pass, chr22 burden masks, per-SNP / per-unit p arrays
    (npz), per-causal-SNP rows, per-stratum power rows and internal-null rows.

Arms
----
 unconditional (plasmode_null_harness null fit + score / SPA; test chi2 + spa):
  N0 saturated_cell          Age + Sex + cell dummies, mixed cells (oracle for the null;
                             under a polygenic alternative it lacks the within-cell ancestry adjustment)
  N0p saturated_cell_pcs     Age + Sex + PC1-20 + cell dummies, mixed cells (oracle)
  N1 pooled_pcs              Age + Sex + PC1-20                      all 27,306
  N2 pooled_center           N1 + centre
  N4 pooled_cohort           N2 + Cohort_cov (mono / n<50 -> Other)
  N5 pooled_project          N2 + Study_cov  (mono / n<50 -> Other)
  S1 strat_project           N1 + project x centre cell (full factor, nothing collapsed),
                             rows of the mixed project x centre cells
  M1 riga_pair33             MatchIt 1:1 on KEY0 (centre key), pair33 formula
  M2 riga_cohortkey          MatchIt 1:1 on KEY3 (+Cohort), + Cohort_cov
  M3 riga_projectkey         MatchIt 1:1 on KEY_P (+project), + Study_cov
  P1/P2/P3 riga_projectkey_py1/2/3   greedy 1:k on KEY_P, + Study_cov
 conditional (stratified / conditional-logistic score test at beta = 0 with the
 exact within-stratum variance; test chi2 only, family 'clr'):
  C0 strat_cell_clr          strata = cell_A            (exact oracle)
  C5 strat_project_clr       strata = project x centre
  C2 strat_keyP_clr          strata = KEY_P (Sex x age band x centre x lpa x project)
  C3 riga_projectkey_clr     strata = the MatchIt KEY_P pairs
  CP1/CP2/CP3 riga_projectkey_py{1,2,3}_clr   strata = the greedy 1:k sets

--ystar-from <null npz> replaces steps 1-5 by the stored permuted-null Y* of a
deliverable-B replicate (and its MatchIt KEY_P pairs from rows__riga_projectkey), so the
NEW arms get their own null distribution on the same 20 null replicates
(outputs/plasmode/power/nullarms/).

Outputs (<out-root>/<native|nullarms>/)
  npz/rep<k>_seed<seed>.npz   y (= ystar, the key regenie_single_variant.py reads), ad, eta, p_true, maf, chr, pos, excluded, near_causal,
                              causal_* arrays, per arm p_chi2__/p_spa__/flags__/rows__/
                              strata__, burden_chi2__<mask>__<arm>, burden_spa__...,
                              burden_causal_level__<mask>, burden_gamma__<mask>,
                              burden_linked__<mask>
  npz/rep<k>_seed<seed>.json  validation, causal design, arm_info, timing
  units_plasmode_<tag>_rep<k>.csv   (one file per replicate for every table below)
                              per arm x test: internal-null metrics (non-causal SNPs
                              >= 1 Mb from every causal SNP, APOE / known windows out),
                              causal mean chi2 and power at nominal alphas
  power_plasmode_<tag>.csv    per arm x test x (class, maf_bin, level): n, mean chi2,
                              n_lt_<alpha>
  causal_plasmode_<tag>.csv   per causal SNP x arm x test: p, chi2
  burden_plasmode_<tag>.csv   per mask x arm x test x group (null / causal level):
                              n, tail metrics, n_lt_<alpha>, mean chi2
Never edits an existing file.  Heavy compute only via sbatch.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
import warnings
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
import scipy.sparse as sps
from scipy.optimize import brentq
from scipy.special import expit
from scipy.stats import chi2

ROOT = Path("${ADSP_ROOT}")
SCRIPTS = ROOT / "revision" / "scripts"
sys.path.insert(0, str(SCRIPTS))
import benchmark_null_harness as H  # noqa: E402  (read-only import; nothing in it is modified)

OUT_ROOT = H.OUT_ROOT / "power"
ACOUNT = H.OUT_ROOT / "pool_281k_counts.acount"
ARTIFACT_LIST = H.OUT_ROOT / "pool_281k_step1_excluded.csv"
DRIVER_PATH = Path(__file__).resolve()

# ------------------------------------------------------------ causal design --
SNP_LEVELS = (0.06, 0.09, 0.12)            # log-odds per SD of genotype
SNP_CLASSES = ("clean", "tgated")
N_PER_SNP_STRATUM = {"clean": 10, "tgated": 5}   # per (class x MAF bin x level)
MIN_CAUSAL_SPACING = 1_000_000
BURDEN_MASKS = ("coding_mac1_5", "coding_mac6_50", "gene_mac6_50", "gene_mac51_500")
BURDEN_LEVELS = (0.1, 0.2, 0.3)            # gamma = a * |log10(mean per-variant MAF)|
N_PER_BURDEN_STRATUM = 8
BURDEN_CMAC_RANGE = (20, 1000)
POWER_ALPHAS = ((1e-3, "1e3"), (1e-4, "1e4"), (1e-5, "1e5"), (1e-6, "1e6"), (5e-8, "5e8"))
CHI2_CAP = 300.0                            # float32 p underflows near chi2 ~ 170; cap the statistic

UNCOND_ARMS = ["saturated_cell", "saturated_cell_pcs", "pooled_pcs", "pooled_center", "pooled_cohort", "pooled_project",
               "strat_project", "riga_pair33", "riga_cohortkey", "riga_projectkey",
               "riga_projectkey_py1", "riga_projectkey_py2", "riga_projectkey_py3"]
CLR_ARMS = ["strat_cell_clr", "strat_project_clr", "strat_keyP_clr", "riga_projectkey_clr",
            "riga_projectkey_py1_clr", "riga_projectkey_py2_clr", "riga_projectkey_py3_clr"]
ARM_ORDER = UNCOND_ARMS + CLR_ARMS
NEW_ARMS = ["saturated_cell_pcs", "strat_project", "riga_projectkey_py1", "riga_projectkey_py2", "riga_projectkey_py3"] + CLR_ARMS
ARM_LABEL = dict(H.ARM_LABEL, saturated_cell_pcs="N0p", strat_project="S1", riga_projectkey_py1="P1", riga_projectkey_py2="P2",
                 riga_projectkey_py3="P3", strat_cell_clr="C0", strat_project_clr="C5", strat_keyP_clr="C2",
                 riga_projectkey_clr="C3", riga_projectkey_py1_clr="CP1", riga_projectkey_py2_clr="CP2",
                 riga_projectkey_py3_clr="CP3")
PY_MATCH_K = {"riga_projectkey_py1": 1, "riga_projectkey_py2": 2, "riga_projectkey_py3": 3}
H_ARMS = {"saturated_cell", "pooled_pcs", "pooled_center", "pooled_cohort", "pooled_project",
          "riga_pair33", "riga_cohortkey", "riga_projectkey"}

T0 = time.time()


def log(msg: str) -> None:
    print(f"[{time.time() - T0:8.1f}s] {msg}", flush=True)


# ================================================================== seeds ====
def read_power_manifest(path: Path, replicate: int) -> int:
    if path.exists():
        m = pd.read_csv(path)
        row = m[m["task"].astype(int) == int(replicate)]
        if len(row) == 1:
            return int(row["seed"].iloc[0])
    n = max(40, replicate + 1)
    children = np.random.SeedSequence(20260921).spawn(n)
    return int(children[replicate - 1].generate_state(1, dtype=np.uint64)[0])


def write_power_manifest(path: Path, n: int = 40) -> None:
    if path.exists():
        return
    children = np.random.SeedSequence(20260921).spawn(n)
    rows = [dict(task=i + 1, task_type="power_native", spawn_key=i,
                 seed=int(c.generate_state(1, dtype=np.uint64)[0])) for i, c in enumerate(children)]
    path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(path, index=False)


# ========================================================== causal SNP draw ==
def pool_maf_from_acount(acount: Path, bim: pd.DataFrame) -> np.ndarray:
    a = pd.read_csv(acount, sep="\t", usecols=["ID", "ALT_CTS", "OBS_CT"], dtype={"ID": str})
    a = a.set_index("ID").reindex(bim.id.astype(str).values)
    if a["OBS_CT"].isna().any():
        raise SystemExit(f"{acount.name}: {int(a['OBS_CT'].isna().sum())} bim ids missing")
    alt = a["ALT_CTS"].to_numpy(dtype=np.float64)
    obs = a["OBS_CT"].to_numpy(dtype=np.float64)
    f = np.where(obs > 0, alt / np.maximum(obs, 1), np.nan)
    return np.minimum(f, 1.0 - f)


def draw_causal_snps(ann: pd.DataFrame, maf: np.ndarray, artifact: np.ndarray, rng: np.random.Generator,
                     n_per_class: Dict[str, int] = N_PER_SNP_STRATUM, classes: Sequence[str] = SNP_CLASSES) -> pd.DataFrame:
    """Stratified random draw with >= MIN_CAUSAL_SPACING between causal SNPs."""
    base_ok = (~ann.excluded.to_numpy()) & (~artifact) & np.isfinite(maf) & (maf >= 0.05)
    cls_mask = {
        "clean": base_ok & ann.clean.to_numpy() & (~ann.tgated.to_numpy()),
        "tgated": base_ok & ann.tgated.to_numpy() & (~ann.flagged.to_numpy()) & (~ann.miss_dirty.to_numpy()),
    }
    chrom = ann.chr.to_numpy()
    pos = ann.pos.to_numpy()
    chosen: List[int] = []
    rows = []
    for cls in classes:
        n_per = int(n_per_class[cls])
        for lo, hi in H.MAF_BINS:
            elig = np.flatnonzero(cls_mask[cls] & (maf >= lo) & (maf < hi))
            if elig.size < n_per * len(SNP_LEVELS) * 3:
                log(f"  WARNING: only {elig.size} eligible SNPs for {cls} maf[{lo},{hi})")
            order = elig[rng.permutation(elig.size)]
            need = n_per * len(SNP_LEVELS)
            picked: List[int] = []
            for j in order:
                if len(picked) >= need:
                    break
                if chosen:
                    c_arr = np.asarray(chosen)
                    close = (chrom[c_arr] == chrom[j]) & (np.abs(pos[c_arr] - pos[j]) < MIN_CAUSAL_SPACING)
                    if close.any():
                        continue
                picked.append(int(j))
                chosen.append(int(j))
            if len(picked) < need:
                raise SystemExit(f"could not place {need} causal SNPs in {cls} maf[{lo},{hi}) (got {len(picked)})")
            lev = np.repeat(np.arange(len(SNP_LEVELS)), n_per)
            for j, li in zip(picked, lev):
                rows.append(dict(snp_index=int(j), snp=str(ann.snp.iloc[j]), chr=int(chrom[j]), pos=int(pos[j]),
                                 causal_class=cls, maf_bin=f"maf[{lo:g},{hi:g})", maf_pool=float(maf[j]),
                                 level=int(li), c=float(SNP_LEVELS[li]),
                                 sign=int(rng.choice([-1, 1]))))
    return pd.DataFrame(rows)


def near_causal_mask(ann: pd.DataFrame, causal: pd.DataFrame, window: int = MIN_CAUSAL_SPACING) -> np.ndarray:
    chrom = ann.chr.to_numpy()
    pos = ann.pos.to_numpy()
    near = np.zeros(len(ann), dtype=bool)
    for c, p in zip(causal.chr.to_numpy(), causal.pos.to_numpy()):
        near |= (chrom == c) & (np.abs(pos - p) <= window)
    return near


# ======================================================= burden causal draw ==
def load_gene_tables(burden_dir: Path, masks: Dict[str, Dict]) -> Dict[str, pd.DataFrame]:
    out = {}
    for m in masks:
        f = burden_dir / f"gene_table_{m}.csv"
        if f.exists():
            t = pd.read_csv(f)
            if len(t) == masks[m]["B"].shape[0]:
                out[m] = t
    return out


WIN_BP = 100_000


def unit_spans(tables: Dict[str, pd.DataFrame]) -> Dict[str, pd.DataFrame]:
    """gene_table_<mask>.csv carries a NaN pos_min for many units (builder defect; pos_max
    is intact).  Windows: pos_min = window index x 100 kb by construction.  Genes: the
    minimum pos_min / maximum pos_max over every mask's table of that gene; a gene with
    no finite pos_min anywhere gets pos_max - 200 kb (generous -> more windows flagged as
    linked, i.e. conservative for the burden null).  Returns the tables with the columns
    span_lo / span_hi / span_imputed added."""
    gene_lo: Dict[str, float] = {}
    gene_hi: Dict[str, float] = {}
    for t in tables.values():
        if str(t.unit_type.iloc[0]) == "win100k":
            continue
        for u, lo, hi in zip(t.unit.astype(str), t.pos_min, t.pos_max):
            if np.isfinite(lo):
                gene_lo[u] = min(gene_lo.get(u, np.inf), float(lo))
            if np.isfinite(hi):
                gene_hi[u] = max(gene_hi.get(u, -np.inf), float(hi))
    out = {}
    for m, t in tables.items():
        t = t.copy()
        if str(t.unit_type.iloc[0]) == "win100k":
            widx = t.unit.astype(str).str.extract(r"win(\d+)$")[0].astype(float)
            lo = widx * WIN_BP
            hi = t.pos_max.astype(float)
            imp = t.pos_min.isna()
        else:
            lo = t.unit.astype(str).map(gene_lo).astype(float)
            hi = t.unit.astype(str).map(gene_hi).astype(float)
            hi = hi.where(np.isfinite(hi), t.pos_max.astype(float))
            imp = ~np.isfinite(lo)
            lo = lo.where(np.isfinite(lo), hi - 2 * WIN_BP)
        t["span_lo"] = lo.to_numpy()
        t["span_hi"] = hi.to_numpy()
        t["span_imputed"] = imp.to_numpy()
        out[m] = t
    return out


def draw_causal_burden(tables: Dict[str, pd.DataFrame], N: int, rng: np.random.Generator,
                       n_per: int = N_PER_BURDEN_STRATUM) -> pd.DataFrame:
    used: set = set()
    rows = []
    for m in BURDEN_MASKS:
        if m not in tables:
            log(f"  WARNING: burden mask {m} not available; no causal units drawn from it")
            continue
        t = tables[m]
        elig = t[(t.cMAC >= BURDEN_CMAC_RANGE[0]) & (t.cMAC <= BURDEN_CMAC_RANGE[1]) & (t.n_variants >= 2)]
        elig = elig[~elig.unit.astype(str).isin(used)]
        order = elig.iloc[rng.permutation(len(elig))]
        need = n_per * len(BURDEN_LEVELS)
        if len(order) < need:
            raise SystemExit(f"burden mask {m}: only {len(order)} eligible units for {need} causal draws")
        pick = order.iloc[:need]
        lev = np.repeat(np.arange(len(BURDEN_LEVELS)), n_per)
        for (_, r), li in zip(pick.iterrows(), lev):
            mean_maf = float(r.cMAC) / (2.0 * N * float(r.n_variants))
            gamma = BURDEN_LEVELS[li] * abs(np.log10(max(mean_maf, 1e-7)))
            used.add(str(r.unit))
            rows.append(dict(mask=m, unit_index=int(r.unit_index), unit=str(r.unit), unit_type=str(r.unit_type),
                             n_variants=int(r.n_variants), cMAC=float(r.cMAC), cMAF=float(r.cMAF),
                             pos_min=int(r.span_lo), pos_max=int(r.span_hi), span_imputed=bool(r.span_imputed),
                             mean_var_maf=mean_maf,
                             level=int(li), a=float(BURDEN_LEVELS[li]), gamma=float(gamma)))
    return pd.DataFrame(rows)


def burden_linked_flags(tables: Dict[str, pd.DataFrame], causal_b: pd.DataFrame, window: int = 0) -> Dict[str, np.ndarray]:
    """Per mask: True for units that share a gene with a causal unit (any mask), for window
    units that overlap a causal gene's span, and -- when window > 0 -- for ANY unit whose span
    lies within `window` bp of a causal gene's span (rare haplotypes are shared over hundreds
    of kb, so neighbouring genes' burdens track a causal burden: the dense 96-unit run showed
    the cell-stratified oracle's non-causal gene units at 55-200x FPR@1e-3); the causal units
    themselves included."""
    genes = set(causal_b.unit.astype(str)) if len(causal_b) else set()
    spans = causal_b[["pos_min", "pos_max"]].to_numpy() if len(causal_b) else np.zeros((0, 2))
    out = {}
    for m, t in tables.items():
        linked = t.unit.astype(str).isin(genes).to_numpy()
        if len(spans):
            lo = t.span_lo.to_numpy()[:, None]
            hi = t.span_hi.to_numpy()[:, None]
            w = int(window) if str(t.unit_type.iloc[0]) != "win100k" else max(int(window), 0)
            ov = (lo <= spans[None, :, 1] + w) & (hi >= spans[None, :, 0] - w)
            if str(t.unit_type.iloc[0]) == "win100k" or w > 0:
                linked |= ov.any(axis=1)
        out[m] = linked
    return out


# ============================================================== generator ====
def generate_y(pool: pd.DataFrame, eta: np.ndarray, rng: np.random.Generator) -> Tuple[np.ndarray, np.ndarray, pd.DataFrame]:
    """Y ~ Bernoulli(expit(alpha_cell + eta)) in mixed cells, alpha_cell solved for the
    observed expected count; mono cells keep the observed labels."""
    ad = pool["AD"].to_numpy().astype(np.int8)
    cells = pool["cell_A"].astype(str).to_numpy()
    y = ad.copy()
    p_true = ad.astype(np.float64).copy()
    rows = []
    for c, idx in pd.Series(np.arange(len(pool))).groupby(cells).indices.items():
        idx = np.asarray(idx)
        k = int(ad[idx].sum())
        n = int(idx.size)
        if k == 0 or k == n or n < 2:
            rows.append(dict(cell=c, n=n, k=k, mono=True, alpha=np.nan, expected=float(k), realized=k))
            continue
        e = eta[idx]
        f = lambda a: float(expit(a + e).sum() - k)  # noqa: E731
        alpha = brentq(f, -60.0, 60.0, xtol=1e-10, maxiter=500)
        p = expit(alpha + e)
        yi = (rng.random(n) < p).astype(np.int8)
        y[idx] = yi
        p_true[idx] = p
        rows.append(dict(cell=c, n=n, k=k, mono=False, alpha=float(alpha), expected=float(p.sum()),
                         realized=int(yi.sum()), sd_expected=float(np.sqrt((p * (1 - p)).sum()))))
    return y, p_true, pd.DataFrame(rows)


def within_cell_coef_partial(x: np.ndarray, y: np.ndarray, z: np.ndarray, cell: np.ndarray) -> Tuple[float, float]:
    """OLS coefficient of y in x ~ cell dummies + y + z (z = the burden part of eta)."""
    d = pd.DataFrame({"x": x, "y": y.astype(float), "z": z.astype(float), "c": cell})
    d = d[np.isfinite(d.x.values)]
    R = np.column_stack([d.y - d.groupby("c").y.transform("mean"), d.z - d.groupby("c").z.transform("mean")])
    xr = (d.x - d.groupby("c").x.transform("mean")).to_numpy()
    XtX = R.T @ R
    if np.linalg.matrix_rank(XtX) < 2:
        return np.nan, np.nan
    beta = np.linalg.solve(XtX, R.T @ xr)
    resid = xr - R @ beta
    dfres = len(d) - d.c.nunique() - 2
    s2 = float((resid ** 2).sum() / max(dfres, 1))
    se = float(np.sqrt(s2 * np.linalg.inv(XtX)[0, 0]))
    return float(beta[0]), se


def validate_power(pool: pd.DataFrame, y: np.ndarray, p_true: np.ndarray, cells_df: pd.DataFrame,
                   eta_burden: Optional[np.ndarray] = None) -> Tuple[Dict, List[str]]:
    """eta_burden (historical name) = the generator's full linear predictor eta; VP4b conditions on it."""
    ad = pool["AD"].to_numpy().astype(np.int8)
    N = len(pool)
    out: Dict = dict(N=int(N), n_case_real=int(ad.sum()), n_case_y=int(y.sum()))
    fails: List[str] = []
    mixed = cells_df[~cells_df.mono]
    z = (mixed.realized - mixed.expected) / mixed.sd_expected.replace(0, np.nan)
    out["VP1_n_cells"] = int(len(cells_df))
    out["VP1_n_mono_cells"] = int(cells_df.mono.sum())
    out["VP1_n_subjects_mono"] = int(cells_df.n[cells_df.mono].sum())
    out["VP1_max_abs_z_count"] = float(np.nanmax(np.abs(z.values))) if len(z) else 0.0
    out["VP1_mono_exact"] = bool(np.array_equal(y[pool["cell_A"].astype(str).isin(cells_df.cell[cells_df.mono]).to_numpy()],
                                                ad[pool["cell_A"].astype(str).isin(cells_df.cell[cells_df.mono]).to_numpy()]))
    # gate 4.5 -> 6.0 (2026-09-21): task 1956231_45 hit 4.58 by chance (85 mixed cells x 60 replicates; small cells are
    # skewed binomials); a mis-generated Y (labels permuted across cells) gives |z| far above 10.  The value is recorded.
    out["VP1_pass"] = bool(out["VP1_mono_exact"] and out["VP1_max_abs_z_count"] <= 6.0)
    if not out["VP1_pass"]:
        fails.append(f"VP1 per-cell counts: max |z| {out['VP1_max_abs_z_count']:.2f}, mono exact {out['VP1_mono_exact']}")
    ctr = pool["center"].astype(str).to_numpy()
    s_real, df_c, V_real = H.pearson_chi2(ctr, ad)
    s_y, _, V_y = H.pearson_chi2(ctr, y)
    out.update(VP2_chi2_center_y_real=s_real, VP2_chi2_center_y=s_y, VP2_df=df_c, VP2_cramerV_real=V_real, VP2_cramerV_y=V_y)
    base = ["Age", "Sex"] + [c for c in H.PCS20 if c in pool.columns]
    l_real = H.lrt_center_given(pool, ad, base)
    l_y = H.lrt_center_given(pool, y, base)
    rel = abs(l_y["LR"] - l_real["LR"]) / max(l_real["LR"], 1e-9)
    out.update(VP2_lrt_center_pcs_real=l_real["LR"], VP2_lrt_center_pcs_y=l_y["LR"], VP2_lrt_center_pcs_rel_diff=float(rel))
    out["VP2_pass"] = bool(abs(s_y - s_real) <= 0.20 * s_real and rel <= 0.20)   # structure report; Bernoulli draws give 0.00-0.13 within-cell
    if not out["VP2_pass"]:
        fails.append(f"VP2 centre-AD dependence: chi2 {s_y:.1f} vs {s_real:.1f}, LRT rel diff {rel:.3f}")
    cells = pool["cell_A"].astype(str).to_numpy()
    mets, _ = H.find_tech_columns(pool)
    worst = 0.0
    v4 = {}
    for m in mets:
        x = pd.to_numeric(pool[m], errors="coerce").to_numpy(dtype=np.float64)
        x = (x - np.nanmean(x)) / (np.nanstd(x) if np.nanstd(x) > 0 else 1.0)
        b, se = H.within_cell_coef(x, y, cells)
        zz = b / se if se and np.isfinite(se) and se > 0 else np.nan
        v4[m] = dict(coef=b, se=se, z=zz)
        if np.isfinite(zz):
            worst = max(worst, abs(zz))
    out["VP4_within_cell_tech_coefs"] = v4
    out["VP4_max_abs_z"] = float(worst)
    out["VP4_pass_marginal"] = bool(worst <= 4.5)      # SOFT: a causal set whose dosages are technically driven (the
    # centre-batch t-gated SNP class, rare-allele burdens) makes Y track the sequencing metrics within cells BY
    # CONSTRUCTION (job 1939566 tasks 4/6/7 |z| 8.4 with SNP-only spike-ins; burden smoke 9.9).  The hard gate is the
    # coefficient of Y given the generator's linear predictor eta (VP4b): conditional on eta, Y is Bernoulli noise.
    if eta_burden is not None:
        worst_b = 0.0
        v4b = {}
        for m in mets:
            x = pd.to_numeric(pool[m], errors="coerce").to_numpy(dtype=np.float64)
            x = (x - np.nanmean(x)) / (np.nanstd(x) if np.nanstd(x) > 0 else 1.0)
            b, se = within_cell_coef_partial(x, y, eta_burden, cells)
            zz = b / se if se and np.isfinite(se) and se > 0 else np.nan
            v4b[m] = dict(coef=b, se=se, z=zz)
            if np.isfinite(zz):
                worst_b = max(worst_b, abs(zz))
        out["VP4b_within_cell_tech_coefs_given_eta"] = v4b
        out["VP4b_max_abs_z"] = float(worst_b)
        out["VP4b_pass"] = bool(worst_b <= 4.5)
        if not out["VP4b_pass"]:
            fails.append(f"VP4b within-cell technical coefficient given eta |z| = {worst_b:.2f} > 4.5")
    else:
        out["VP4_pass"] = out["VP4_pass_marginal"]
        if not out["VP4_pass"]:
            fails.append(f"VP4 within-cell technical coefficient |z| = {worst:.2f} > 4.5")
    g = pd.DataFrame({"c": cells, "y": ad}).groupby("c").y.agg(["sum", "size"])
    pc = g["sum"] / g["size"]
    expect = float((g["size"] * (pc ** 2 + (1 - pc) ** 2)).sum() / N)
    out["VP5_mean_y_eq_ad"] = float((y == ad).mean())
    out["VP5_exchangeability_expectation"] = expect
    out["VP5_identity_leakage"] = float((y == ad).mean() - expect)
    return out, fails


# ============================================================ CLR (strata) ==
class CLR:
    """Conditional-logistic / stratified score test at beta = 0 with the exact
    within-stratum (permutation) variance:
        U = sum_s sum_i g_i (y_i - ybar_s)
        V = sum_s [n1_s n0_s / (n_s (n_s - 1))] sum_i (g_i - gbar_s)^2
    Rows of strata without both a case and a control are dropped."""

    def __init__(self, rows: np.ndarray, strata: np.ndarray, y_pool: np.ndarray):
        rows = np.asarray(rows, dtype=np.int64)
        codes, _ = pd.factorize(pd.Series(strata).astype(str), sort=True)
        yy = y_pool[rows].astype(np.float64)
        n_s = np.bincount(codes).astype(np.float64)
        n1 = np.bincount(codes, weights=yy)
        n0 = n_s - n1
        ok_s = (n1 > 0) & (n0 > 0) & (n_s >= 2)
        keep = ok_s[codes]
        self.rows = rows[keep]
        codes = codes[keep]
        codes, _ = pd.factorize(codes, sort=True)
        self.codes = codes
        yy = yy[keep]
        self.y = yy
        n_s = np.bincount(codes).astype(np.float64)
        n1 = np.bincount(codes, weights=yy)
        n0 = n_s - n1
        self.n_strata = int(n_s.size)
        self.n_rows = int(self.rows.size)
        self.n_cases = int(n1.sum())
        self.ybar = n1 / n_s
        self.yc = yy - self.ybar[codes]
        self.kappa = n1 * n0 / (n_s * (n_s - 1.0))
        self.n_s = n_s
        self.inv_n = 1.0 / n_s
        self.S = sps.csr_matrix((np.ones(self.n_rows), (np.arange(self.n_rows), codes)),
                                shape=(self.n_rows, self.n_strata))
        self.ST = self.S.T.tocsr()
        # effective sample size of the stratified design: sum 4 n1 n0 / n
        self.n_effective = float((4.0 * n1 * n0 / n_s).sum())
        self.info_weight = float((n1 * n0 / n_s).sum())          # sum_s n_s pi_s (1 - pi_s)

    def score(self, G: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """G: (m x n_rows) float64 -> (p, stat, V)."""
        U = G @ self.yc
        sums = np.asarray(self.ST @ G.T).T                    # m x S
        gbar = sums * self.inv_n[None, :]
        ss = np.asarray(self.ST @ (G * G).T).T - self.n_s[None, :] * gbar * gbar
        V = ss @ self.kappa
        V = np.where(V > 1e-12, V, np.nan)
        stat = U * U / V
        p = chi2.sf(stat, df=1)
        return p, stat, V


# ============================================================ 1:k matcher ===
def key_rows_ok(pool: pd.DataFrame, key_cols: Sequence[str]) -> np.ndarray:
    """Rows with every key column present (benchmark_matching.R drops NA age bands)."""
    ok = np.ones(len(pool), dtype=bool)
    for c in key_cols:
        ok &= pool[c].notna().to_numpy()
    return ok


def greedy_match_k(pool: pd.DataFrame, y: np.ndarray, key_cols: Sequence[str], k: int,
                   rng: np.random.Generator) -> Tuple[np.ndarray, np.ndarray, Dict]:
    """Greedy nearest-Mahalanobis(PC1-10) 1:k matching without replacement inside exact
    key strata, in k ROUNDS (every case gets its first control before any case gets a
    second, MatchIt's ratio semantics, so the 1:k case set is the 1:1 case set); cases
    visited in random order.  Returns (rows, set_id per row, summary)."""
    pcs = [c for c in H.PCS10 if c in pool.columns]
    P = pool[pcs].to_numpy(dtype=np.float64)
    cov = np.cov(P.T) + 1e-12 * np.eye(len(pcs))
    Li = np.linalg.inv(np.linalg.cholesky(cov))
    Z = P @ Li.T
    strata = pool[list(key_cols)].astype(str).agg("|".join, axis=1).to_numpy()
    ok = key_rows_ok(pool, key_cols)
    rows: List[int] = []
    sets: List[int] = []
    sid = 0
    n_ctrl_hist = np.zeros(k + 1, dtype=np.int64)
    for _, idx in pd.Series(np.flatnonzero(ok)).groupby(strata[ok]).indices.items():
        idx = np.flatnonzero(ok)[np.asarray(idx)]
        idx = np.asarray(idx)
        ci = idx[y[idx] == 1]
        ki = idx[y[idx] == 0]
        if ci.size == 0 or ki.size == 0:
            continue
        ci = ci[rng.permutation(ci.size)]
        avail = np.ones(ki.size, dtype=bool)
        Zk = Z[ki]
        got: List[List[int]] = [[] for _ in range(ci.size)]
        for _round in range(k):                      # round r: every case takes ONE more control
            for ii, c in enumerate(ci):
                if not avail.any():
                    break
                d = np.sum((Zk - Z[c]) ** 2, axis=1)
                d[~avail] = np.inf
                j = int(np.argmin(d))
                if not np.isfinite(d[j]):
                    break
                avail[j] = False
                got[ii].append(int(ki[j]))
        for ii, c in enumerate(ci):
            if not got[ii]:
                continue
            rows.append(int(c))
            sets.append(sid)
            for r_ in got[ii]:
                rows.append(r_)
                sets.append(sid)
            n_ctrl_hist[len(got[ii])] += 1
            sid += 1
    rows_a = np.asarray(rows, dtype=np.int64)
    sets_a = np.asarray(sets, dtype=np.int64)
    summ = dict(n_sets=int(sid), n_rows=int(rows_a.size), n_cases=int(sid), n_controls=int(rows_a.size - sid),
                controls_per_set_hist={int(i): int(v) for i, v in enumerate(n_ctrl_hist) if i > 0})
    return rows_a, sets_a, summ


# ================================================================== arms =====
def base_frame(pool: pd.DataFrame, y: np.ndarray) -> pd.DataFrame:
    base = pool.copy()
    base["Y"] = y.astype(np.int64)
    base["Sex"] = pd.to_numeric(base["Sex"], errors="coerce").fillna(0).astype(np.int64)
    base["Age"] = pd.to_numeric(base["Age"], errors="coerce")
    base["Age"] = base["Age"].fillna(base["Age"].mean()).astype(np.float64)
    base["lpa_cluster"] = base["lpa_cluster"].astype(str)
    base["center"] = base["center"].astype(str)
    base["Study_strat"] = base["Study"].astype(str)
    return base


def mixed_rows(labels: np.ndarray, y: np.ndarray) -> np.ndarray:
    d = pd.DataFrame({"l": pd.Series(labels).astype(str).to_numpy(), "y": y})
    ct = d.groupby("l").y.agg(["sum", "size"])
    ok = ct.index[(ct["sum"] > 0) & (ct["sum"] < ct["size"])]
    return np.flatnonzero(d.l.isin(ok).to_numpy())


def build_frames(pool: pd.DataFrame, y: np.ndarray, matched: Dict[str, np.ndarray], py_sets: Dict[str, Tuple[np.ndarray, np.ndarray]],
                 arms: Sequence[str]) -> Tuple[Dict, Dict[str, CLR]]:
    """Unconditional frames (df, rows, covs, n_pairs) and CLR objects per arm."""
    frames: Dict = {}
    h_arms = [a for a in arms if a in H_ARMS]
    if h_arms:
        frames.update(H.arm_frames(pool, y, matched, h_arms))
    base = base_frame(pool, y)
    pcs20 = [c for c in H.PCS20 if c in base.columns]
    pcs10 = [c for c in H.PCS10 if c in base.columns]
    pair33 = ["Age", "Sex", "center", "lpa_cluster", "cluster_prob", "entropy_norm"] + pcs10
    if "saturated_cell_pcs" in arms:
        # the oracle under a POLYGENIC alternative: cell dummies AND PC1-20.  A polygenic score built from real SNPs is correlated
        # with ancestry inside the cells (causal allele frequencies differ by ancestry), so the cell-only model (N0) leaves a
        # within-cell ancestry channel open (its non-causal SNPs ran at 38x FPR@1e-5 under the alternative, 1.3x under the null).
        rows = mixed_rows(pool["cell_A"].astype(str).to_numpy(), y)
        d = base.iloc[rows].reset_index(drop=True)
        d["cell_A"] = d["cell_A"].astype(str)
        frames["saturated_cell_pcs"] = (d, rows, ["Age", "Sex"] + pcs20 + ["cell_A"], None)
    if "strat_project" in arms:
        rows = mixed_rows(pool["cell_SC"].astype(str).to_numpy() if "cell_SC" in pool.columns
                          else (pool["Study"].astype(str) + "|" + pool["center"].astype(str)).to_numpy(), y)
        d = base.iloc[rows].reset_index(drop=True)
        d["cell_SC_strat"] = (d["Study"].astype(str) + "|" + d["center"].astype(str))
        frames["strat_project"] = (d, rows, ["Age", "Sex"] + pcs20 + ["cell_SC_strat"], None)
    for arm, k in PY_MATCH_K.items():
        if arm in arms and arm in py_sets:
            rows, _ = py_sets[arm]
            d = base.iloc[rows].reset_index(drop=True)
            d["Study_cov"] = H.collapse_levels(d["Study"], d["Y"].to_numpy())
            frames[arm] = (d, rows, pair33 + ["Study_cov"], int((y[rows] == 1).sum()))
    clr: Dict[str, CLR] = {}
    all_rows = np.arange(len(pool))
    if "strat_cell_clr" in arms:
        clr["strat_cell_clr"] = CLR(all_rows, pool["cell_A"].astype(str).to_numpy(), y)
    if "strat_project_clr" in arms:
        lab = (pool["Study"].astype(str) + "|" + pool["center"].astype(str)).to_numpy()
        clr["strat_project_clr"] = CLR(all_rows, lab, y)
    if "strat_keyP_clr" in arms:
        lab = pool[H.MATCH_KEYS["KEY_P"]].astype(str).agg("|".join, axis=1).to_numpy()
        okr = np.flatnonzero(key_rows_ok(pool, H.MATCH_KEYS["KEY_P"]))
        clr["strat_keyP_clr"] = CLR(okr, lab[okr], y)
    if "riga_projectkey_clr" in arms and matched.get("KEY_P") is not None:
        rows = np.asarray(matched["KEY_P"], dtype=np.int64)
        h = rows.size // 2
        clr["riga_projectkey_clr"] = CLR(rows, np.concatenate([np.arange(h), np.arange(h)]).astype(str), y)
    for arm, k in PY_MATCH_K.items():
        carm = arm + "_clr"
        if carm in arms and arm in py_sets:
            rows, sets = py_sets[arm]
            clr[carm] = CLR(rows, sets.astype(str), y)
    return frames, clr


# ================================================================ metrics ====
def chi2_from_p(p: np.ndarray) -> np.ndarray:
    p64 = p.astype(np.float64)
    with np.errstate(divide="ignore", invalid="ignore"):
        q = chi2.isf(np.clip(p64, 1e-300, 1.0), 1)
    q = np.where(np.isfinite(p64), q, np.nan)
    q = np.where(p64 <= 0, CHI2_CAP, q)
    return np.minimum(q, CHI2_CAP)


def power_block(p: np.ndarray) -> Dict:
    ok = np.isfinite(p) & (p > 0) & (p <= 1)
    pp = p[ok]
    out: Dict = dict(n=int(pp.size), n_nan=int((~ok).sum()))
    q = chi2_from_p(pp) if pp.size else np.zeros(0)
    out["mean_chi2"] = float(q.mean()) if pp.size else np.nan
    out["median_chi2"] = float(np.median(q)) if pp.size else np.nan
    for a, tag in POWER_ALPHAS:
        out[f"n_lt_{tag}"] = int((pp < a).sum())
        out[f"power_{tag}"] = float((pp < a).mean()) if pp.size else np.nan
    return out


# ================================================================ run task ===
def run_task(args) -> int:
    replicate = int(args.replicate)
    null_mode = bool(args.ystar_from)
    spike = args.spike
    tag = args.tag or ("nullarms" if null_mode else ("burden" if spike == "burden" else "power"))
    out_dir = Path(args.out_root) / (args.out_sub or ("nullarms" if null_mode else ("burden" if spike == "burden" else "native")))
    skip_scan = bool((not null_mode) and spike == "burden")
    npz_dir = out_dir / "npz"
    npz_dir.mkdir(parents=True, exist_ok=True)
    manifest = Path(args.manifest)
    if not null_mode:
        write_power_manifest(manifest)
    arms = list(ARM_ORDER) if not args.arms else [a for a in args.arms.split(",") if a]
    if null_mode and not args.arms:
        arms = list(NEW_ARMS)
    bad = [a for a in arms if a not in ARM_ORDER]
    if bad:
        raise SystemExit(f"unknown arm(s) {bad}")

    pool = pd.read_csv(args.pool_csv, low_memory=False)
    N = len(pool)
    bfile = Path(args.bfile)
    bed = H.Bed(bfile)
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
    M_all = bed.n_snp
    bim = bed.bim.reset_index(drop=True)
    ann = H.snp_annotation(bim, Path(args.persnp) if args.persnp else None)
    maf_pool = pool_maf_from_acount(Path(args.acount), bim)
    artifact = np.zeros(M_all, dtype=bool)
    if Path(args.artifact_list).exists():
        art = pd.read_csv(args.artifact_list, usecols=["ID"], dtype={"ID": str}).ID.astype(str)
        artifact = bim.id.astype(str).isin(set(art)).to_numpy()
    burden = H.load_burden_masks(Path(args.burden_dir) if args.burden_dir else None)
    burden = {m: b for m, b in burden.items() if b["B"].shape[1] == N}
    tables = unit_spans(load_gene_tables(Path(args.burden_dir), burden)) if burden else {}
    log(f"pool N={N}; bed {bfile.name}: {bed.n_fam} samples, {M_all} SNPs; annotation excluded={int(ann.excluded.sum())} "
        f"tgated={int(ann.tgated.sum())} clean={int(ann.clean.sum())} artifact={int(artifact.sum())}; burden masks {len(burden)}")

    # ---------------------------------------------------------------- Y ----
    ad = pool["AD"].to_numpy().astype(np.int8)
    causal = pd.DataFrame({"snp_index": pd.Series(dtype=np.int64), "snp": pd.Series(dtype=str), "chr": pd.Series(dtype=np.int64),
                           "pos": pd.Series(dtype=np.int64), "causal_class": pd.Series(dtype=str), "maf_bin": pd.Series(dtype=str),
                           "maf_pool": pd.Series(dtype=float), "level": pd.Series(dtype=np.int64), "c": pd.Series(dtype=float),
                           "sign": pd.Series(dtype=np.int64), "beta": pd.Series(dtype=float), "ncp_expected": pd.Series(dtype=float)})
    causal_b = pd.DataFrame()
    linked: Dict[str, np.ndarray] = {}
    near = np.zeros(M_all, dtype=bool)
    eta = np.zeros(N, dtype=np.float64)
    p_true = np.full(N, np.nan)
    cells_df = pd.DataFrame()
    vjson: Dict = {}
    fails: List[str] = []
    matched_from_npz: Dict[str, np.ndarray] = {}
    keep = np.ones(N, dtype=bool)
    if null_mode:
        z = np.load(args.ystar_from, allow_pickle=False)
        y = z["ystar"].astype(np.int8)
        seed = int(z["seed"]) if "seed" in z.files else int(args.seed)
        if y.size != N:
            raise SystemExit(f"{args.ystar_from}: ystar {y.size} != pool {N}")
        for key, arm in (("KEY0", "riga_pair33"), ("KEY3", "riga_cohortkey"), ("KEY_P", "riga_projectkey")):
            k = f"rows__{arm}"
            if k in z.files:
                matched_from_npz[key] = z[k].astype(np.int64)
        vjson = dict(mode="nullarms", ystar_from=str(args.ystar_from), n_case_y=int(y.sum()),
                     matched_keys_from_npz=sorted(matched_from_npz))
        log(f"NULL-ARM mode: Y* from {Path(args.ystar_from).name} (seed {seed}, {int(y.sum())} cases); "
            f"MatchIt rows for {sorted(matched_from_npz)}")
    else:
        seed = read_power_manifest(manifest, replicate) if args.seed is None else int(args.seed)
        ss = np.random.SeedSequence(seed)
        ch = ss.spawn(4)
        rng_c, rng_y, rng_b, rng_m = (np.random.default_rng(c) for c in ch)
        if spike in ("snp", "both"):
            causal = draw_causal_snps(ann, maf_pool, artifact, rng_c, classes=[c for c in args.snp_classes.split(",") if c])
            near = near_causal_mask(ann, causal)
        # genotypes of the causal SNPs (same decoder + imputation as the scan)
        betas = np.zeros(len(causal))
        sds = np.zeros(len(causal))
        Gc = np.zeros((len(causal), N), dtype=np.float64)
        for i, j in enumerate(causal.snp_index.to_numpy() if len(causal) else []):
            g8 = bed.read_chunk(int(j), 1)[:, pool_idx]
            g, _, _ = H.impute_chunk(g8)
            Gc[i] = g[0]
            sds[i] = float(g[0].std())
            betas[i] = float(causal.sign.iloc[i]) * float(causal.c.iloc[i]) / max(sds[i], 1e-6)
            eta += betas[i] * Gc[i]
        if len(causal):
            causal["sd_geno"] = sds
            causal["beta"] = betas
            causal["maf_scan"] = np.minimum(Gc.mean(1) / 2.0, 1 - Gc.mean(1) / 2.0)
        eta_snp_var = float(eta.var())
        eta_b = np.zeros(N, dtype=np.float64)
        if burden and spike in ("burden", "both"):
            causal_b = draw_causal_burden(tables, N, rng_c, n_per=int(args.burden_per_stratum))
            for _, r in causal_b.iterrows():
                eta_b += float(r.gamma) * burden[r["mask"]]["B"][int(r.unit_index)].astype(np.float64)
            eta += eta_b
            linked = burden_linked_flags(tables, causal_b, window=int(args.burden_link_window))
        log(f"causal design: {len(causal)} SNPs ({dict(causal.causal_class.value_counts())}), "
            f"{len(causal_b)} burden units; var(eta) SNP {eta_snp_var:.3f} total {eta.var():.3f} (sd {eta.std():.3f})")
        # ---- optional ascertainment emulation in the single-label cells (--mono-ascertain RATE) ----
        # Cells with only cases or only controls cannot be given labels by the within-cell model, so by default their
        # subjects keep their labels and carry NO signal (conservative for the pooled designs, which use them).  With
        # --mono-ascertain r > 0 every subject of such a cell gets a LATENT status ~ Bernoulli(expit(logit(r) + eta)) and
        # only the subjects whose latent status equals the project's label are kept (an all-control project keeps its
        # latent controls, an all-case project its latent cases): the retained subjects carry the case / control
        # enrichment of a recruited case-control sample, while the batch confounding stays.  Matched, stratified and
        # oracle designs never use these cells, so they are unchanged; the pooled designs gain their real information.
        keep = np.ones(N, dtype=bool)
        asc_rows = []
        if args.mono_ascertain > 0:
            r_lat = float(args.mono_ascertain)
            # eta is NOT centred (uncentred dosages x random-sign effects give a replicate-specific offset, -3.5 in the
            # first smoke): solve the latent intercept so that the MEAN latent case rate over the whole pool equals r_lat
            # (inside the mixed cells alpha_cell absorbs the offset; here it must be removed explicitly).
            a_lat = brentq(lambda a: float(expit(a + eta).mean() - r_lat), -80.0, 80.0, xtol=1e-10, maxiter=500)
            log(f"ascertainment: latent intercept {a_lat:.3f} (mean eta {eta.mean():.3f}) -> mean latent rate {expit(a_lat + eta).mean():.3f}")
            cells_all = pool["cell_A"].astype(str).to_numpy()
            for c, idx in pd.Series(np.arange(N)).groupby(cells_all).indices.items():
                idx = np.asarray(idx)
                k = int(ad[idx].sum())
                if 0 < k < idx.size:
                    continue
                p_lat = expit(a_lat + eta[idx])
                y_lat = (rng_y.random(idx.size) < p_lat).astype(np.int8)
                kp = y_lat == ad[idx]
                keep[idx] = kp
                asc_rows.append(dict(cell=c, label=int(ad[idx][0]), n=int(idx.size), n_kept=int(kp.sum()),
                                     mean_eta_all=float(eta[idx].mean()), mean_eta_kept=float(eta[idx][kp].mean()) if kp.any() else np.nan))
            n_drop = int((~keep).sum())
            log(f"ascertainment (latent rate {r_lat:.3f}): {len(asc_rows)} single-label cells, {n_drop} of "
                f"{sum(r['n'] for r in asc_rows)} subjects dropped -> pool N {int(keep.sum())}")
            pool = pool.iloc[keep].reset_index(drop=True)
            ad = ad[keep]
            eta = eta[keep]
            pool_idx = pool_idx[keep]
            pool_iid = pool_iid[keep]
            if len(causal):
                Gc = Gc[:, keep]
            for m_ in burden:
                burden[m_]["B"] = np.ascontiguousarray(burden[m_]["B"][:, keep])
            N = int(keep.sum())
        y, p_true, cells_df = generate_y(pool, eta, rng_y)
        vjson, fails = validate_power(pool, y, p_true, cells_df, eta)
        vjson.update(mono_ascertain_rate=float(args.mono_ascertain), n_pool_kept=int(N), n_pool_dropped=int((~keep).sum()),
                     ascertainment_cells=asc_rows)
        vjson.update(mode="power", spike=spike, eta_var_snp=eta_snp_var, eta_var_burden=float(eta_b.var()), eta_var_total=float(eta.var()), eta_sd=float(eta.std()),
                     n_causal_snps=int(len(causal)), n_causal_burden=int(len(causal_b)),
                     n_near_causal_snps=int(near.sum()),
                     alpha_range=[float(cells_df.alpha.min()), float(cells_df.alpha.max())] if len(cells_df) else None)
        # expected information on the causal SNPs (cell-stratified, at the truth)
        cells = pool["cell_A"].astype(str).to_numpy()
        w = p_true * (1 - p_true)
        codes, _ = pd.factorize(cells)
        S = sps.csr_matrix((np.ones(N), (np.arange(N), codes)))
        ST = S.T.tocsr()
        sw = np.asarray(ST @ w).ravel()
        vjson["mean_w_truth"] = float(w.mean())
        if len(causal):
            Gw = Gc * w[None, :]
            gbar = np.asarray(ST @ Gw.T).T / np.maximum(sw, 1e-12)[None, :]
            ss = np.asarray(ST @ (Gw * Gc).T).T - sw[None, :] * gbar * gbar
            V_true = ss.sum(1)
            causal["info_cellstrat"] = V_true
            causal["ncp_expected"] = betas ** 2 * V_true
        if len(causal_b):
            # the same for the causal burden units (their B row, at the truth)
            vb = []
            for _, r in causal_b.iterrows():
                b = burden[r["mask"]]["B"][int(r.unit_index)].astype(np.float64)
                bw = b * w
                bbar = np.asarray(ST @ bw).ravel() / np.maximum(sw, 1e-12)
                vb.append(float((np.asarray(ST @ (bw * b)).ravel() - sw * bbar * bbar).sum()))
            causal_b["info_cellstrat"] = vb
            causal_b["ncp_expected"] = causal_b.gamma.to_numpy() ** 2 * np.asarray(vb)
        ncp_msg = (f"expected NCP mean {causal.ncp_expected.mean():.1f} [{causal.ncp_expected.min():.1f}, {causal.ncp_expected.max():.1f}]"
                   if len(causal) else "no causal SNPs")
        if len(causal_b):
            ncp_msg += (f"; burden expected NCP mean {causal_b.ncp_expected.mean():.1f} "
                        f"[{causal_b.ncp_expected.min():.1f}, {causal_b.ncp_expected.max():.1f}]")
        log(f"VP1 cells={vjson['VP1_n_cells']} mono={vjson['VP1_n_mono_cells']} max|z|={vjson['VP1_max_abs_z_count']:.2f} | "
            f"VP2 chi2 {vjson['VP2_chi2_center_y']:.1f} (real {vjson['VP2_chi2_center_y_real']:.1f}) "
            f"LRT centre|PCs {vjson['VP2_lrt_center_pcs_y']:.1f} (real {vjson['VP2_lrt_center_pcs_real']:.1f}) | "
            f"VP4 marginal max|z|={vjson['VP4_max_abs_z']:.2f}" + (f" | VP4b given eta {vjson['VP4b_max_abs_z']:.2f}" if 'VP4b_max_abs_z' in vjson else "")
            + f" | VP5 leakage={vjson['VP5_identity_leakage']:.4f} | mean w {w.mean():.3f} | "
            f"cases {int(y.sum())} (real {int(ad.sum())}); {ncp_msg}")
        del Gc
    vjson.update(replicate=replicate, seed=seed, tag=tag, bfile=str(bfile), pool_csv=str(args.pool_csv), driver=DRIVER_PATH.name)
    npz_path = npz_dir / f"rep{replicate:03d}_seed{seed}.npz"
    json_path = npz_dir / f"rep{replicate:03d}_seed{seed}.json"
    failed_marker = npz_dir / f"FAILED_rep{replicate:03d}_seed{seed}.json"
    if npz_path.exists() and json_path.exists() and not args.force:
        log(f"replicate {replicate} seed {seed} already done ({npz_path.name}) -> skip")
        return 0
    json_path.write_text(json.dumps(vjson, indent=1, default=float))
    if fails:
        failed_marker.write_text(json.dumps(dict(replicate=replicate, seed=seed, fails=fails), indent=1))
        log("VALIDATION FAILED: " + " | ".join(fails))
        return 3
    if null_mode:
        rng_m = np.random.default_rng(np.random.SeedSequence(seed).spawn(5)[4])
        rng_b = np.random.default_rng(np.random.SeedSequence(seed).spawn(5)[3])

    # ------------------------------------------------------------ matching --
    matched: Dict[str, np.ndarray] = dict(matched_from_npz)
    need_keys = [k for a, k in H.MATCHED_ARMS.items() if a in arms]
    if "riga_projectkey_clr" in arms:
        need_keys.append("KEY_P")
    need_keys = [k for k in dict.fromkeys(need_keys) if k not in matched]
    matcher = "npz" if matched else "none"
    if need_keys:
        mroot = Path(args.match_dir) if args.match_dir else (out_dir / "match")
        mdir = mroot / f"rep{replicate:03d}_seed{seed}"
        mdir.mkdir(parents=True, exist_ok=True)
        mtag = f"rep{replicate:03d}"
        pheno_cols = [c for c in ("IID", "row_idx", "fam_pos", "Study", "Cohort", "center", "AD") if c in pool.columns]
        pheno_csv = mdir / f"y_pool_{mtag}.csv"
        pp = pool[pheno_cols].copy()
        pp["Y"] = y.astype(np.int64)
        if "row_idx" in pp.columns:
            pp["row_idx"] = np.arange(N)          # benchmark_matching.R carries row_idx through; it must be the position in THIS csv
        pp.to_csv(pheno_csv, index=False)
        script = Path(args.match_script)
        if args.match_fallback == "python" or not script.exists():
            fb = H.python_fallback_match(pool, y, seed=0)
            matched.update({k: fb[k] for k in need_keys if k in fb})
            matcher = "python_fallback"
        else:
            matched.update(H.run_match_r(script, pheno_csv, pheno_csv, mdir, mtag, need_keys, args.match_cmd, seed=seed))
            matcher = "matchit"
    for k, v in matched.items():
        log(f"  matched {k}: {len(v) // 2} pairs ({len(v)} rows)")
    py_sets: Dict[str, Tuple[np.ndarray, np.ndarray]] = {}
    py_summ: Dict[str, Dict] = {}
    for arm, k in PY_MATCH_K.items():
        if arm in arms or (arm + "_clr") in arms:
            rows, sets, summ = greedy_match_k(pool, y, H.MATCH_KEYS["KEY_P"], k, rng_m)
            py_sets[arm] = (rows, sets)
            py_summ[arm] = summ
            log(f"  greedy 1:{k} on KEY_P: {summ['n_sets']} sets, {summ['n_controls']} controls, hist {summ['controls_per_set_hist']}")

    # ----------------------------------------------------------- null fits --
    frames, clr = build_frames(pool, y, matched, py_sets, arms)
    nfs: Dict[str, Dict] = {}
    arm_info: Dict[str, Dict] = {}
    for arm in arms:
        if arm in clr:
            c = clr[arm]
            arm_info[arm] = dict(error="", family="clr", n_subjects=c.n_rows, n_cases_arm=c.n_cases, n_strata=c.n_strata,
                                 n_effective=c.n_effective, info_weight=c.info_weight, n_pairs=np.nan,
                                 converged=True, inestimable_design=False,
                                 matcher=("greedy" if arm.startswith("riga_projectkey_py") else
                                          (matcher if arm == "riga_projectkey_clr" else "none")))
            log(f"  [{arm:26s}] CLR rows={c.n_rows:6d} cases={c.n_cases:5d} strata={c.n_strata:5d} n_eff={c.n_effective:.0f}")
            continue
        if arm not in frames or frames[arm][0] is None:
            arm_info[arm] = dict(error="rows unavailable", n_subjects=0)
            log(f"  [{arm}] SKIPPED (no rows)")
            continue
        d, rows, covs, n_pairs = frames[arm]
        t1 = time.time()
        try:
            nf = H.null_fit_local(d, covs)
        except Exception as exc:  # noqa: BLE001
            arm_info[arm] = dict(error=f"{type(exc).__name__}: {exc}", n_subjects=int(len(d)))
            log(f"  [{arm}] null fit FAILED: {arm_info[arm]['error']}")
            continue
        nfs[arm] = nf
        y_arm = d["Y"].to_numpy()
        n1, n0 = int(y_arm.sum()), int((1 - y_arm).sum())
        arm_info[arm] = dict(error="", family="uncond", n_subjects=nf["n_subjects"], n_cases_arm=nf["n_cases_arm"],
                             n_pairs=(n_pairs if n_pairs is not None else np.nan), rank=nf["rank"], converged=nf["converged"],
                             n_params=nf["n_params"], max_abs_coef_design=nf["max_abs_coef_design"],
                             max_abs_coef_any=nf["max_abs_coef_any"], inestimable_design=nf["inestimable_design"],
                             n_effective=float(4.0 * n1 * n0 / max(n1 + n0, 1)), n_eff_w=nf["n_effective"],
                             warn_classes=nf["warn_classes"], perfect_separation_warning=nf["perfect_separation_warning"],
                             null_fit_s=float(time.time() - t1),
                             matcher=("greedy" if arm in PY_MATCH_K else (matcher if arm in H.MATCHED_ARMS else "none")),
                             covariates=";".join(c for c in covs if not c.upper().startswith("PC"))
                             + (";PC1..PC20" if any(c in H.PCS20[10:] for c in covs) else (";PC1..PC10" if any(c in H.PCS10 for c in covs) else "")))
        if arm in py_summ:
            arm_info[arm].update(n_sets=py_summ[arm]["n_sets"], n_controls=py_summ[arm]["n_controls"])
        log(f"  [{arm:26s}] N={nf['n_subjects']:6d} cases={nf['n_cases_arm']:5d} rank={nf['rank']:4d} conv={nf['converged']} "
            f"max|design coef|={nf['max_abs_coef_design']:.2f} inest={nf['inestimable_design']} sep={nf['perfect_separation_warning']} ({time.time() - t1:.1f}s)")
    live_u = [a for a in arms if a in nfs]
    live_c = [a for a in arms if a in clr]
    live = [a for a in arms if a in nfs or a in clr]
    if not live:
        raise SystemExit("no arm is live")

    # ---------------------------------------------------------- decoder ----
    scratch = Path(args.scratch) / f"rep{replicate:03d}"
    dec = None
    if not args.no_verify_decoder and not skip_scan:
        dec = H.verify_decoder_plink2(bed, pool_idx, pool_iid, scratch, plink2=args.plink2)
        if dec["bad"]:
            failed_marker.write_text(json.dumps(dict(replicate=replicate, seed=seed, fails=["decoder mismatch"], decoder=dec)))
            return 5

    # ---------------------------------------------------------- bed pass ---
    M = 0 if skip_scan else (M_all if not args.max_chunks else min(M_all, args.max_chunks * args.chunk))
    P = {arm: {"chi2": np.full(M, np.nan, dtype=np.float32), "spa": np.full(M, np.nan, dtype=np.float32),
               "flags": np.zeros(M, dtype=np.int8)} for arm in live_u}
    P.update({arm: {"chi2": np.full(M, np.nan, dtype=np.float32)} for arm in live_c})
    maf = np.full(M, np.nan, dtype=np.float32)
    missr = np.full(M, np.nan, dtype=np.float32)
    t_pass = time.time()
    t_arm = {arm: 0.0 for arm in live}
    n_chunks = (M + args.chunk - 1) // args.chunk
    for start in range(0, M, args.chunk):
        m = min(args.chunk, M - start)
        G8 = bed.read_chunk(start, m)[:, pool_idx]
        G, maf_c, miss_c = H.impute_chunk(G8)
        del G8
        maf[start:start + m] = maf_c
        missr[start:start + m] = miss_c
        for arm in live_u:
            ta = time.time()
            rows = frames[arm][1]
            Ga = G if (rows.size == N and np.array_equal(rows, np.arange(N))) else G[:, rows]
            p_spa, p_chi2, fl = H.score_from_nullfit(Ga, nfs[arm], args.fastspa_cutoff)
            P[arm]["chi2"][start:start + m] = p_chi2
            P[arm]["spa"][start:start + m] = p_spa
            P[arm]["flags"][start:start + m] = fl
            t_arm[arm] += time.time() - ta
        for arm in live_c:
            ta = time.time()
            p, _, _ = clr[arm].score(G[:, clr[arm].rows])
            P[arm]["chi2"][start:start + m] = p
            t_arm[arm] += time.time() - ta
        if (start // args.chunk) % 10 == 0:
            log(f"  chunk {start // args.chunk + 1}/{n_chunks} ({time.time() - t_pass:.0f}s elapsed)")
    log(f"bed pass done in {time.time() - t_pass:.0f}s; per-arm s: " + ", ".join(f"{ARM_LABEL.get(a, a)}={t_arm[a]:.0f}" for a in live))

    # ------------------------------------------------------------ metrics --
    ann_m = ann.iloc[:M].reset_index(drop=True)
    near_m = near[:M]
    excluded = ann_m.excluded.to_numpy()
    null_mask = (~excluded) & (~near_m)
    if len(causal):
        causal_in = causal[causal.snp_index < M].reset_index(drop=True)
    else:
        causal_in = causal
    common = dict(design=tag, scheme=("nullarms" if null_mode else "power"), scenario="adsp_pool", seed=seed,
                  replicate=replicate, n_pool=int(N), mono_ascertain=float(args.mono_ascertain if not null_mode else 0.0), n_snps_bed=int(M), n_causal=int(len(causal_in)),
                  n_near_causal=int(near_m.sum()), bfile=bfile.name, tool="native", is_external=False,
                  driver=DRIVER_PATH.name, eta_sd=(float(eta.std()) if not null_mode else 0.0))
    unit_rows: List[Dict] = []
    power_rows: List[Dict] = []
    causal_rows: List[Dict] = []
    for arm in live:
        info = arm_info[arm]
        tests = ("chi2", "spa") if arm in live_u else ("chi2",)
        base = dict(common, arm=arm, arm_label=ARM_LABEL.get(arm, arm), family=info.get("family", ""),
                    is_matched=bool(arm in H.MATCHED_ARMS or arm.startswith("riga_projectkey_py")),
                    primary_test=("spa" if arm in live_u else "chi2"),
                    **{k: v for k, v in info.items() if k not in ("family",)}, covset=arm)
        for test in tests:
            p = P[arm][test].astype(np.float64)
            m_null = H.metrics(p[null_mask])
            m_null_clean = H.metrics(p[null_mask & ann_m.clean.to_numpy()])
            m_null_tg = H.metrics(p[null_mask & ann_m.tgated.to_numpy()])
            row = dict(base, test=test, **{f"null_{k}": v for k, v in m_null.items()},
                       **{f"nullclean_{k}": v for k, v in m_null_clean.items()},
                       **{f"nulltgated_{k}": v for k, v in m_null_tg.items()})
            if len(causal_in):
                pc = p[causal_in.snp_index.to_numpy()]
                row.update({f"causal_{k}": v for k, v in power_block(pc).items()})
                for (cls, mb, lev), sub in causal_in.groupby(["causal_class", "maf_bin", "level"], sort=True):
                    pcs_ = p[sub.snp_index.to_numpy()]
                    power_rows.append(dict(base, test=test, causal_class=cls, maf_bin=mb, level=int(lev),
                                           c=float(SNP_LEVELS[int(lev)]), ncp_expected_mean=float(sub.ncp_expected.mean()),
                                           **power_block(pcs_)))
                for cls, sub in causal_in.groupby("causal_class", sort=True):
                    power_rows.append(dict(base, test=test, causal_class=cls, maf_bin="all", level=-1, c=np.nan,
                                           ncp_expected_mean=float(sub.ncp_expected.mean()),
                                           **power_block(p[sub.snp_index.to_numpy()])))
                q = chi2_from_p(pc)
                for i, r in causal_in.iterrows():
                    causal_rows.append(dict(seed=seed, replicate=replicate, arm=arm, arm_label=ARM_LABEL.get(arm, arm), test=test,
                                            snp=r.snp, snp_index=int(r.snp_index), chr=int(r.chr), pos=int(r.pos),
                                            causal_class=r.causal_class, maf_bin=r.maf_bin, level=int(r.level), c=float(r.c),
                                            sign=int(r.sign), beta=float(r.beta), maf_pool=float(r.maf_pool),
                                            ncp_expected=float(r.ncp_expected), p=float(pc[i]), chi2=float(q[i])))
            unit_rows.append(row)
            msg = (f"  [{ARM_LABEL.get(arm, arm):3s} {arm:26s} {test}] null lam={m_null['lambda_gc']:.3f} "
                   f"1e-3={m_null['fpr_1e3_x']:.2f}x 1e-5={m_null['fpr_1e5_x']:.1f}x")
            if len(causal_in):
                msg += (f" | causal mean chi2={row['causal_mean_chi2']:.1f} power@1e-5={row['causal_power_1e5']:.2f} "
                        f"@5e-8={row['causal_power_5e8']:.2f}")
            log(msg)

    # ------------------------------------------------------------- burden --
    burden_rows: List[Dict] = []
    burden_out: Dict[str, np.ndarray] = {}
    for mask, bd in burden.items():
        B = bd["B"]
        cmaf = bd["cmac"] / (2.0 * N)
        lk = linked.get(mask, np.zeros(B.shape[0], dtype=bool))
        clev = np.full(B.shape[0], -1, dtype=np.int64)
        cgam = np.zeros(B.shape[0])
        cncp = np.zeros(B.shape[0])
        if len(causal_b):
            sub = causal_b[causal_b["mask"] == mask]
            clev[sub.unit_index.to_numpy()] = sub.level.to_numpy()
            cgam[sub.unit_index.to_numpy()] = sub.gamma.to_numpy()
            cncp[sub.unit_index.to_numpy()] = sub.ncp_expected.to_numpy()
        burden_out[f"causal_level__{mask}"] = clev
        burden_out[f"gamma__{mask}"] = cgam
        burden_out[f"ncp_expected__{mask}"] = cncp
        burden_out[f"linked__{mask}"] = lk
        for arm in live:
            if arm in live_u:
                rows = frames[arm][1]
                p_spa, p_chi2, fl = H.score_from_nullfit(B[:, rows], nfs[arm], args.fastspa_cutoff)
                tests = (("chi2", p_chi2), ("spa", p_spa))
                burden_out[f"chi2__{mask}__{arm}"] = p_chi2.astype(np.float32)
                burden_out[f"spa__{mask}__{arm}"] = p_spa.astype(np.float32)
            else:
                p_c, _, _ = clr[arm].score(B[:, clr[arm].rows].astype(np.float64))
                tests = (("chi2", p_c),)
                burden_out[f"chi2__{mask}__{arm}"] = p_c.astype(np.float32)
            for test, p in tests:
                base = dict(common, arm=arm, arm_label=ARM_LABEL.get(arm, arm), test=test, mask=mask,
                            n_subjects=arm_info[arm]["n_subjects"], primary_test=("spa" if arm in live_u else "chi2"))
                nul = (~lk) & (clev < 0)
                burden_rows.append(dict(base, group="null", n_units=int(nul.sum()), **H.tail_metrics(p[nul])))
                for lo, hi in H.CMAF_BINS:
                    mk = nul & (cmaf >= lo) & (cmaf < hi)
                    burden_rows.append(dict(base, group=f"null_cmaf[{lo:g},{hi:g})", n_units=int(mk.sum()), **H.tail_metrics(p[mk])))
                if len(causal_b):
                    for li in range(len(BURDEN_LEVELS)):
                        mk = clev == li
                        if mk.any():
                            burden_rows.append(dict(base, group=f"causal_L{li}", level=li, a=float(BURDEN_LEVELS[li]),
                                                    n_units=int(mk.sum()), gamma_mean=float(cgam[mk].mean()),
                                                    cmac_mean=float(bd["cmac"][mk].mean()), **power_block(p[mk])))
        log(f"  burden {mask}: {B.shape[0]} units ({int((clev >= 0).sum())} causal, {int(lk.sum())} linked) scored on {len(live)} arms")

    # ------------------------------------------------------ VP6 + json ----
    if len(causal_in) and "saturated_cell" in live_u:
        q = chi2_from_p(P["saturated_cell"]["chi2"][causal_in.snp_index.to_numpy()].astype(np.float64))
        exp_ncp = causal_in.ncp_expected.to_numpy()
        vjson["VP6_saturated_cell_mean_chi2_minus1"] = float(np.nanmean(q) - 1.0)
        vjson["VP6_expected_ncp_mean"] = float(exp_ncp.mean())
        vjson["VP6_ratio_realized_over_expected"] = float((np.nanmean(q) - 1.0) / max(exp_ncp.mean(), 1e-9))
        log(f"VP6 realized (mean chi2 - 1) / expected NCP on the causal SNPs, saturated_cell: {vjson['VP6_ratio_realized_over_expected']:.3f}")
    if len(causal_b) and "saturated_cell" in live_u:
        qs, es = [], []
        for mask in burden:
            clev = burden_out[f"causal_level__{mask}"]
            if (clev >= 0).any():
                qs.append(chi2_from_p(burden_out[f"chi2__{mask}__saturated_cell"][clev >= 0].astype(np.float64)))
                es.append(burden_out[f"ncp_expected__{mask}"][clev >= 0])
        if qs:
            qs = np.concatenate(qs)
            es = np.concatenate(es)
            vjson["VP6b_burden_saturated_cell_mean_chi2_minus1"] = float(np.nanmean(qs) - 1.0)
            vjson["VP6b_burden_expected_ncp_mean"] = float(es.mean())
            vjson["VP6b_burden_ratio_realized_over_expected"] = float((np.nanmean(qs) - 1.0) / max(es.mean(), 1e-9))
            log(f"VP6b realized / expected NCP on the causal burden units, saturated_cell: {vjson['VP6b_burden_ratio_realized_over_expected']:.3f}")
    vjson.update(arm_info=arm_info, decoder_verification=dec, t_bed_pass_s=float(time.time() - t_pass), t_arm_s=t_arm,
                 n_snps=int(M), burden_masks=list(burden), matcher=matcher, py_match=py_summ,
                 n_pairs={k: int(len(v) // 2) for k, v in matched.items()}, arms_live=live, t_total_s=float(time.time() - T0))
    json_path.write_text(json.dumps(vjson, indent=1, default=float))

    # ---------------------------------------------------------- outputs ----
    # one CSV per replicate (10 concurrent array tasks appending to one file would race)
    rep_sfx = f"_rep{replicate:03d}"
    units_csv = out_dir / f"units_plasmode_{tag}{rep_sfx}.csv"
    power_csv = out_dir / f"power_plasmode_{tag}{rep_sfx}.csv"
    causal_csv = out_dir / f"causal_plasmode_{tag}{rep_sfx}.csv"
    burden_csv = out_dir / f"burden_plasmode_{tag}{rep_sfx}.csv"
    pd.DataFrame(unit_rows).to_csv(units_csv, index=False)
    if power_rows:
        pd.DataFrame(power_rows).to_csv(power_csv, index=False)
    if causal_rows:
        pd.DataFrame(causal_rows).to_csv(causal_csv, index=False)
    if burden_rows:
        pd.DataFrame(burden_rows).to_csv(burden_csv, index=False)
    npz = dict(y=y.astype(np.int8), ystar=y.astype(np.int8), ad=ad, eta=eta.astype(np.float32), p_true=p_true.astype(np.float32),
               keep=(keep if not null_mode else np.ones(N, dtype=bool)),
               maf=maf, miss_rate=missr, seed=np.uint64(seed), replicate=np.int64(replicate),
               chr=ann_m.chr.to_numpy().astype(np.int16), pos=ann_m.pos.to_numpy().astype(np.int64),
               excluded=excluded, near_causal=near_m, arms=np.array(live), mode=np.array([vjson["mode"]]))
    if len(causal_in):
        npz.update(causal_idx=causal_in.snp_index.to_numpy().astype(np.int64),
                   causal_class=causal_in.causal_class.to_numpy().astype(str),
                   causal_maf_bin=causal_in.maf_bin.to_numpy().astype(str),
                   causal_level=causal_in.level.to_numpy().astype(np.int8), causal_c=causal_in.c.to_numpy().astype(np.float32),
                   causal_beta=causal_in.beta.to_numpy().astype(np.float32), causal_sign=causal_in.sign.to_numpy().astype(np.int8),
                   causal_ncp_expected=causal_in.ncp_expected.to_numpy().astype(np.float32))
    else:
        npz["causal_idx"] = np.zeros(0, dtype=np.int64)
    for arm in live_u:
        npz[f"p_chi2__{arm}"] = P[arm]["chi2"]
        npz[f"p_spa__{arm}"] = P[arm]["spa"]
        npz[f"flags__{arm}"] = P[arm]["flags"]
        npz[f"rows__{arm}"] = frames[arm][1].astype(np.int32)
    for arm in live_c:
        npz[f"p_chi2__{arm}"] = P[arm]["chi2"]
        npz[f"rows__{arm}"] = clr[arm].rows.astype(np.int32)
        npz[f"strata__{arm}"] = clr[arm].codes.astype(np.int32)
    for k, v in burden_out.items():
        npz[f"burden_{k}"] = v
    np.savez_compressed(npz_path, **npz)
    if len(causal_b):
        causal_b.assign(seed=seed, replicate=replicate).to_csv(npz_dir / f"rep{replicate:03d}_seed{seed}_causal_burden.csv", index=False)
    log(f"wrote {units_csv.name} (+{len(unit_rows)}), {power_csv.name} (+{len(power_rows)}), {causal_csv.name} (+{len(causal_rows)}), "
        f"{burden_csv.name} (+{len(burden_rows)}), {npz_path.name}; total {time.time() - T0:.0f}s")
    return 0


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(description=__doc__.split("Arms")[0], formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--replicate", type=int, required=True)
    ap.add_argument("--spike", choices=["snp", "burden", "both"], default="snp",
                    help="snp: causal single variants + the 281k scan (out subdir native); burden: causal chr22 burden units only, "
                         "no scan (out subdir burden); both: everything in one replicate")
    ap.add_argument("--seed", type=int, default=None, help="override the manifest seed")
    ap.add_argument("--out-sub", default="", help="output subdirectory under --out-root (default native | burden | nullarms)")
    ap.add_argument("--snp-classes", default="clean,tgated",
                    help="causal SNP classes to draw (clean = unflagged, missingness-clean, not centre-batch t-gated; tgated = the "
                         "centre/project-differentiated probe class).  Technically driven causal dosages make Y track the sequencing "
                         "metrics within cells and inflate every technically driven null SNP (2026-09-21: 5 of 20 replicates with the probe "
                         "class had the oracle's internal null at 70-270x at 1e-5), so 'clean' alone gives a usable internal null.")
    ap.add_argument("--mono-ascertain", type=float, default=0.0,
                    help="latent case rate for the single-label cells: subjects whose latent status (from the same model) differs "
                         "from the project's label are dropped (ascertainment emulation); 0 = keep their labels (no signal)")
    ap.add_argument("--burden-per-stratum", type=int, default=N_PER_BURDEN_STRATUM, help="causal burden units per (mask x level)")
    ap.add_argument("--burden-link-window", type=int, default=0,
                    help="bp: null burden units within this distance of a causal gene are 'linked' (excluded from the internal null)")
    ap.add_argument("--tag", default="", help="output tag (default power | nullarms)")
    ap.add_argument("--ystar-from", default="", help="null-arm mode: a deliverable-B npz whose ystar / MatchIt rows are reused")
    ap.add_argument("--pool-csv", default=str(H.OUT_ROOT / "plasmode_pool.csv"))
    ap.add_argument("--manifest", default=str(OUT_ROOT / "seed_manifest_power.csv"))
    ap.add_argument("--bfile", default=str(H.BFILE_POOL))
    ap.add_argument("--acount", default=str(ACOUNT))
    ap.add_argument("--artifact-list", default=str(ARTIFACT_LIST))
    ap.add_argument("--persnp", default=str(H.PERSNP))
    ap.add_argument("--burden-dir", default=str(H.OUT_ROOT / "burden_chr22"))
    ap.add_argument("--out-root", default=str(OUT_ROOT))
    ap.add_argument("--scratch", default=os.environ.get("SIM_SCRATCH", "${SCRATCH_DIR}/"
                    "a773d3e1-08d5-4d13-8c1c-7004109bf225/scratchpad/build/plasmode_power"))
    ap.add_argument("--match-script", default=str(H.MATCH_SCRIPT))
    ap.add_argument("--match-cmd", default="Rscript {script} --pool {pool} --out {outdir} --tag {tag} --ystar {ystar} --keys {keys}")
    ap.add_argument("--match-dir", default="")
    ap.add_argument("--match-fallback", choices=["none", "python"], default="none")
    ap.add_argument("--arms", default="", help="comma list (default: all 19; null-arm mode: the 11 new arms)")
    ap.add_argument("--chunk", type=int, default=2000)
    ap.add_argument("--max-chunks", type=int, default=0)
    ap.add_argument("--fastspa-cutoff", type=float, default=0.05)
    ap.add_argument("--no-verify-decoder", action="store_true")
    ap.add_argument("--plink2", default=H.PLINK2)
    ap.add_argument("--force", action="store_true")
    return ap


def main() -> int:
    args = build_parser().parse_args()
    if args.persnp == "":
        args.persnp = None
    return run_task(args)


if __name__ == "__main__":
    sys.exit(main())
