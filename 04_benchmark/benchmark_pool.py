#!/usr/bin/env python
"""
benchmark_pool.py -- DELIVERABLE B data layer, step 1 (run ONCE via plasmode_pool.sbatch).

Builds the master table of the 27,306-subject unrelated phenotyped pool that every
plasmode task reads, the cell table of the Y* schemes, the per-task seed manifest and
(unless --skip-plink) the pool-restricted PLINK1 fileset for REGENIE.

Pool construction is a COPY of adsp_center_batch_common.py L70-88 (that module is not
importable: argparse + file reads at module level).  The pool is
    pheno/adsp_combined.csv  x  processed/output.pruned.unrelated.unrels
merged with the QC file (SeqCenter, Study, technical metrics) and the LPA table, ordered by
position in processed/chrall.common.pruned.fam (fam_pos), and cross-checked against
revision/tables/adsp_center_batch/rare_chr22_per_sample.csv (same 27,306 IIDs).

Outputs (outputs/plasmode/):
  plasmode_pool.csv    one row per pool subject, in .fam order (row i <-> fam_pos[i]);
                       columns documented in COLUMN_DOC below
  plasmode_cells.csv   cells of scheme A (Cohort|Study|center) and A2 (lpa|Cohort|Study|center):
                       scheme, cell, n, n_case, rate, mono
  seed_manifest.csv    SeedSequence(20260915).spawn(60): task, task_type, spawn_key, seed
                       (seed = first uint64 word of child.generate_state(1) -> use
                       numpy.random.SeedSequence(seed).spawn(k) inside the task)
  pool_keep.txt        FID IID of the pool for plink2 --keep
  pool_281k.{bed,bim,fam}   plink2 --keep pool --make-bed of processed/chrall.common.pruned
                       (sample order == plasmode_pool.csv row order, asserted)
  plasmode_pool_build.json  counts and checks
Nothing existing is modified.
"""
import argparse, os, sys, json, time, subprocess
import numpy as np
import pandas as pd

ROOT = "${ADSP_ROOT}"
PHENO = f"{ROOT}/pheno/adsp_combined.csv"
QCF = f"{ROOT}/pheno/gcad.r5.wgs.58507.2024.11.03.quality.metrics.ALL.csv"
LPA = f"{ROOT}/adsp_lpa_clusters_20PC.tsv"
UNREL = f"{ROOT}/processed/output.pruned.unrelated.unrels"
BFILE = f"{ROOT}/processed/chrall.common.pruned"
PERSAMPLE = f"{ROOT}/revision/tables/adsp_center_batch/rare_chr22_per_sample.csv"
PLINK2 = "${HOME_DIR}/bin/plink2"
OUTDIR_DEFAULT = f"{ROOT}/outputs/plasmode"
ROOT_SEED = 20260915
N_TASKS = 60

TECH_NUM = ["depth_of_coverage_30x", "average.genome.coverage", "Q30_score_pcent",
            "variantEval_novelsites_perc", "variantEval_totalvariants_TiTv",
            "variantEval_novelvariants_TiTv", "variantEval_total_Variants", "Freemix", "Chipmix"]
TECH_FAC = ["Sequencer", "PCR.free", "Read.length..bp."]
MAC_BINS = ["MAC1", "MAC2_5", "MAC6_10", "MAC11_50", "MAC51_500", "MACgt500"]
PCS = [f"PC{i}" for i in range(1, 21)]

COLUMN_DOC = {
    "IID": "sample id (matches .fam IID, GDS sample.id, adsp_combined.csv IID)",
    "row_idx": "0-based row index in this file == 0-based sample index in pool_281k.bed",
    "fam_pos": "0-based sample index in processed/chrall.common.pruned.fam (bed decoder subset index)",
    "SUBJID": "subject id from adsp_combined.csv",
    "SeqCenter": "raw QC-file SeqCenter label (17 raw labels)",
    "center": "11-level merge: riga_design_harness.R L76-80 rule (MEDGENOME, Unknown, Genentech/GENENTECH, "
              "HudsonAlpha, Macrogen/MACROGEN -> Other) + KCCG normalization (Kinghorn... -> KCCG)",
    "center_cb": "adsp_center_batch_common.py rule (raw centers with >= 100 pool subjects kept, else Other; "
                 "MEDGENOME and DNALINK are own levels here) -- the level set of revision/tables/adsp_center_batch/",
    "Cohort": "adsp_combined.csv Cohort (50 levels in pool)",
    "Study": "GCAD project label = QC-file Study (65 levels in pool)",
    "AD": "observed AD label 0/1 (real Y)",
    "Sex": "adsp_combined.csv Sex (1/2 coding as in the file)",
    "Age": "numeric age, '+' stripped",
    "age_band": "R cut(Age, seq(30, 100, 3), right = TRUE) label, e.g. '(75,78]'; NA outside (30, 100]",
    "Race": "adsp_combined.csv Race code",
    "lpa_cluster": "LPA cluster id (1..13) from adsp_lpa_clusters_20PC.tsv",
    "cluster_prob": "raw posterior of the assigned cluster",
    "entropy": "posterior entropy",
    "entropy_norm": "entropy / log(K), K = 13 post_k columns",
    "PC1..PC20": "adsp_combined.csv PCs (the deployed in-sample PCs)",
    "MAC1..MACgt500, n_missing, missing_rate": "chr22 QC7 per-sample minor-allele counts by pool-MAC bin "
              "from rare_chr22_per_sample.csv (built on QC_label7 PASS over the same pool)",
    "<metric>": "raw technical metric from the QC file (9 numeric: " + ", ".join(TECH_NUM) + ")",
    "<metric>_z": "standardized over the pool non-missing values; missing -> 0 (mean)",
    "<metric>_na": "1 if the raw metric is missing (Freemix 29, Chipmix 3,569 in the pool)",
    "Sequencer, PCR.free, Read.length..bp.": "raw QC-file factor labels",
    "Sequencer_norm": "upper-cased Sequencer (merges 'Novaseq'/'NovaSeq'): 6 levels",
    "PCR_free": "PCR.free relabelled (PCRFree/PCRAmplified/Unknown)",
    "Read_length": "Read.length..bp. as string (150/100/101/90)",
    "cell_A": "Cohort|Study|center  (scheme A permutation cell)",
    "cell_A2": "lpa_cluster|Cohort|Study|center  (scheme A2 cell)",
    "cell_SC": "Study|center",
    "SeqCenter_norm": "raw SeqCenter with Kinghorn... -> KCCG only (14 levels)",
    "cell_A_raw": "Cohort|Study|SeqCenter_norm (129 cells, 44 mono, 7,194 subjects: the designer's count)",
    "cell_A_cb": "Cohort|Study|center_cb (128 cells)",
}


def norm_center(s):
    # COPY of adsp_center_batch_common.py norm_center (L77-81)
    u = s.strip().upper()
    if u.startswith("KINGHORN"): return "KCCG"
    if u.startswith("HUDSONALPHA"): return "HUDSONALPHA"
    return u


RIGA_OTHER = {"MEDGENOME", "Unknown", "Genentech", "GENENTECH",
              "HudsonAlpha Institute for Biotechnology", "Macrogen", "MACROGEN"}


def riga_center(s):
    # riga_design_harness.R L75-80: NA -> Unknown; listed labels -> Other; else raw label.
    # Plus the KCCG normalization of the long Kinghorn label (label only; the same subjects).
    if pd.isna(s): s = "Unknown"
    if s in RIGA_OTHER: return "Other"
    if s.startswith("Kinghorn"): return "KCCG"
    return s


def age_band(a):
    # R: cut(Age, breaks = seq(30, 100, by = 3), right = TRUE)  ->  (30,33], (33,36], ..., (97,100]
    if pd.isna(a) or a <= 30 or a > 100: return np.nan
    k = int(np.ceil((a - 30) / 3.0)) - 1
    lo = 30 + 3 * k
    return f"({lo},{lo + 3}]"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default=OUTDIR_DEFAULT)
    ap.add_argument("--skip-plink", action="store_true", help="table only (login-node dry run)")
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--plink-mem-mb", type=int, default=6000)
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    T0 = time.time()
    def log(*a): print(f"[{time.time()-T0:7.1f}s]", *a, flush=True)
    checks = {}

    # ------------------------------------------------------------ pool (COPY of L70-88)
    fam = pd.read_csv(f"{BFILE}.fam", sep=r"\s+", header=None,
                      names=["FID", "IID", "PAT", "MAT", "SEX", "PHE"], dtype=str)
    n_fam = len(fam)
    ph = pd.read_csv(PHENO)
    unrel = set(pd.read_csv(UNREL, sep=r"\s+", header=None, dtype=str)[1])
    qc_all = pd.read_csv(QCF, dtype=str)
    qc = qc_all[["Sample_name", "SeqCenter", "Study"]]
    lpa = pd.read_csv(LPA, sep="\t")
    K_LPA = len([c for c in lpa.columns if c.startswith("post_k")])
    pool = ph[ph.IID.isin(unrel)].copy()
    pool = pool.merge(qc, left_on="IID", right_on="Sample_name", how="left")
    pool = pool.merge(lpa[["IID", "lpa_cluster", "cluster_prob", "entropy"]], on="IID", how="left")
    assert pool.SeqCenter.notna().all() and pool.lpa_cluster.notna().all()
    pool["center_raw"] = pool.SeqCenter.map(norm_center)
    fam_pos = pd.Series(np.arange(n_fam), index=fam.IID.values)
    pool["fam_pos"] = fam_pos.reindex(pool.IID.values).values
    assert pool.fam_pos.notna().all()
    pool = pool.sort_values("fam_pos").reset_index(drop=True)
    n = len(pool)
    cc = pool.center_raw.value_counts()
    keep_centers = [c for c in cc.index if cc[c] >= 100]
    pool["center_cb"] = np.where(pool.center_raw.isin(keep_centers), pool.center_raw, "Other")
    # ------------------------------------------------------------ end of copied block
    pool["fam_pos"] = pool.fam_pos.astype(int)
    pool["row_idx"] = np.arange(n)
    assert not pool.IID.duplicated().any()
    log(f"pool N={n} (fam {n_fam}); center_cb levels {pool.center_cb.nunique()}")
    checks["n_pool"] = int(n)

    # cross-check against the adsp_center_batch per-sample table (same pool, GDS order there)
    ps = pd.read_csv(PERSAMPLE)
    assert set(ps.IID) == set(pool.IID), "pool differs from rare_chr22_per_sample.csv"
    ps = ps.set_index("IID").reindex(pool.IID.values)
    assert (ps.center.values == pool.center_cb.values).all(), "center_cb differs from per-sample table"
    for c in MAC_BINS + ["n_missing", "missing_rate"]:
        pool[c] = ps[c].values
    checks["persample_crosscheck"] = "PASS (IID set and center_cb identical)"

    # centre (11 levels, riga rule + KCCG normalization)
    pool["center"] = pool.SeqCenter.map(riga_center)
    ctab = pool.center.value_counts()
    log("center (riga 11-level):", ctab.to_dict())
    assert pool.center.nunique() == 11, ctab
    checks["center_levels"] = ctab.to_dict()

    # phenotype columns
    pool["Age"] = pd.to_numeric(pool.Age.astype(str).str.replace("+", "", regex=False), errors="coerce")
    assert pool.Age.notna().all()
    pool["age_band"] = pool.Age.map(age_band)
    checks["n_age_band_na"] = int(pool.age_band.isna().sum())
    pool["AD"] = pool.AD.astype(int)
    pool["entropy_norm"] = pool.entropy / np.log(K_LPA)
    pool["lpa_cluster"] = pool.lpa_cluster.astype(int)
    checks["K_lpa"] = int(K_LPA)

    # technical metrics
    qt = qc_all.set_index("Sample_name").reindex(pool.IID.values)
    for c in TECH_NUM:
        v = pd.to_numeric(qt[c], errors="coerce")
        pool[c] = v.values
        pool[c + "_na"] = v.isna().astype(int).values
        z = (v - v.mean()) / v.std(ddof=1)
        pool[c + "_z"] = z.fillna(0.0).values
    for c in TECH_FAC:
        pool[c] = qt[c].values
    pool["Sequencer_norm"] = pool.Sequencer.astype(str).str.upper()
    pool["PCR_free"] = pool["PCR.free"].astype(str)
    pool["Read_length"] = pool["Read.length..bp."].astype(str)
    checks["tech_na"] = {c: int(pool[c + "_na"].sum()) for c in TECH_NUM}
    checks["tech_factor_levels"] = {"Sequencer_norm": pool.Sequencer_norm.value_counts().to_dict(),
                                    "PCR_free": pool.PCR_free.value_counts().to_dict(),
                                    "Read_length": pool.Read_length.value_counts().to_dict()}

    # cells
    pool["cell_A"] = pool.Cohort.astype(str) + "|" + pool.Study.astype(str) + "|" + pool.center
    pool["cell_A2"] = pool.lpa_cluster.astype(str) + "|" + pool.cell_A
    pool["cell_SC"] = pool.Study.astype(str) + "|" + pool.center
    # raw-centre variants (the designer's 129-cell / 44-mono / 7,194 count is reproduced ONLY by the
    # raw SeqCenter label with the Kinghorn -> KCCG normalization, 14 levels; recorded, not primary)
    pool["SeqCenter_norm"] = pool.SeqCenter.astype(str).str.replace(r"^Kinghorn.*", "KCCG", regex=True)
    pool["cell_A_raw"] = pool.Cohort.astype(str) + "|" + pool.Study.astype(str) + "|" + pool.SeqCenter_norm
    pool["cell_A_cb"] = pool.Cohort.astype(str) + "|" + pool.Study.astype(str) + "|" + pool.center_cb
    cells = []
    for scheme, col in [("A", "cell_A"), ("A2", "cell_A2"), ("SC", "cell_SC"),
                        ("A_raw", "cell_A_raw"), ("A_cb", "cell_A_cb")]:
        g = pool.groupby(col).AD.agg(["size", "sum"]).reset_index()
        g.columns = ["cell", "n", "n_case"]
        g["scheme"] = scheme
        g["rate"] = g.n_case / g.n
        g["mono"] = ((g.n_case == 0) | (g.n_case == g.n)).astype(int)
        cells.append(g[["scheme", "cell", "n", "n_case", "rate", "mono"]])
        yb = pool.AD.values
        r2 = 1 - ((yb - pool.groupby(col).AD.transform("mean").values) ** 2).sum() / ((yb - yb.mean()) ** 2).sum()
        checks[f"cells_{scheme}"] = {"n_cells": int(len(g)), "n_mono": int(g.mono.sum()),
                                     "n_subjects_mono": int(g.n[g.mono == 1].sum()),
                                     "n_cells_lt20": int((g.n < 20).sum()),
                                     "between_cell_R2_AD": float(r2)}
        log(f"scheme {scheme}: {checks[f'cells_{scheme}']}")
    cells = pd.concat(cells, ignore_index=True)
    cells.to_csv(f"{args.outdir}/plasmode_cells.csv", index=False)

    # centre x AD dependence on the real Y (V2 anchor numbers)
    from scipy import stats
    for ccol in ["center", "center_cb", "SeqCenter_norm"]:
        ct = pd.crosstab(pool[ccol], pool.AD)
        chi2, p, dof, _ = stats.chi2_contingency(ct, correction=False)
        checks[f"{ccol}_x_AD_realY"] = {"levels": int(pool[ccol].nunique()), "chi2": float(chi2), "df": int(dof),
                                        "cramers_v": float(np.sqrt(chi2 / n))}
        log(f"{ccol} x AD (real Y):", checks[f"{ccol}_x_AD_realY"])

    # write pool
    cols = (["IID", "row_idx", "fam_pos", "SUBJID", "SeqCenter", "center", "center_cb", "Cohort", "Study",
             "AD", "Sex", "Age", "age_band", "Race", "lpa_cluster", "cluster_prob", "entropy", "entropy_norm"]
            + PCS + MAC_BINS + ["n_missing", "missing_rate"]
            + TECH_NUM + [c + "_z" for c in TECH_NUM] + [c + "_na" for c in TECH_NUM]
            + TECH_FAC + ["Sequencer_norm", "PCR_free", "Read_length", "SeqCenter_norm",
                          "cell_A", "cell_A2", "cell_SC", "cell_A_raw", "cell_A_cb"])
    pool[cols].to_csv(f"{args.outdir}/plasmode_pool.csv", index=False)
    log(f"wrote plasmode_pool.csv ({n} x {len(cols)})")

    # seed manifest
    ss = np.random.SeedSequence(ROOT_SEED)
    children = ss.spawn(N_TASKS)
    man = pd.DataFrame({
        "task": np.arange(1, N_TASKS + 1),
        "task_type": ["A_native"] * 20 + ["A2_native"] * 20 + ["spare"] * 20,
        "spawn_key": [c.spawn_key[0] for c in children],
        "seed": [int(c.generate_state(1, dtype=np.uint64)[0]) for c in children],
    })
    man.to_csv(f"{args.outdir}/seed_manifest.csv", index=False)
    assert man.seed.nunique() == N_TASKS

    # keep file
    keep = pd.DataFrame({"FID": fam.FID.values[pool.fam_pos.values], "IID": pool.IID.values})
    assert (keep.IID.values == fam.IID.values[pool.fam_pos.values]).all()
    keep.to_csv(f"{args.outdir}/pool_keep.txt", sep="\t", header=False, index=False)

    with open(f"{args.outdir}/COLUMNS.json", "w") as f:
        json.dump(COLUMN_DOC, f, indent=1)

    if not args.skip_plink:
        out = f"{args.outdir}/pool_281k"
        cmd = [PLINK2, "--bfile", BFILE, "--keep", f"{args.outdir}/pool_keep.txt", "--make-bed",
               "--out", out, "--threads", str(args.threads), "--memory", str(args.plink_mem_mb)]
        log("plink2:", " ".join(cmd))
        r = subprocess.run(cmd, capture_output=True, text=True)
        print(r.stdout[-3000:]); print(r.stderr[-2000:])
        assert r.returncode == 0, "plink2 failed"
        f2 = pd.read_csv(f"{out}.fam", sep=r"\s+", header=None, dtype=str)
        assert len(f2) == n and (f2[1].values == pool.IID.values).all(), \
            "pool_281k.fam order != plasmode_pool.csv order"
        b2 = pd.read_csv(f"{out}.bim", sep=r"\s+", header=None, dtype=str)
        b1 = pd.read_csv(f"{BFILE}.bim", sep=r"\s+", header=None, dtype=str)
        assert len(b2) == len(b1) and (b2[1].values == b1[1].values).all(), "bim changed"
        checks["pool_281k"] = {"n_samples": int(len(f2)), "n_snps": int(len(b2)),
                               "bed_bytes": int(os.path.getsize(out + ".bed")),
                               "expected_bytes": int(3 + len(b2) * ((n + 3) // 4)),
                               "fam_order_equals_pool_csv": True}
        assert checks["pool_281k"]["bed_bytes"] == checks["pool_281k"]["expected_bytes"]
        log("pool_281k:", checks["pool_281k"])
    checks["elapsed_s"] = time.time() - T0
    with open(f"{args.outdir}/plasmode_pool_build.json", "w") as f:
        json.dump(checks, f, indent=1, default=str)
    log("DONE")


if __name__ == "__main__":
    main()
