"""
synthetic_comparator_designs.py — external GWAS-tool comparators for the RIGA simulation
=======================================================================

Adds REGENIE, GMMAT and GENESIS to the six pure-Python covariate-set
variants already implemented in ``synthetic_null_simulation.py``.

Everything here is opt-in: ``simulation_power.run_all_methods_raw`` only
calls into this module when the environment variable ``COMPARATORS`` is
set (e.g. ``COMPARATORS=regenie,gmmat,genesis``).  With the variable
unset the six-method path is bit-for-bit unchanged.

Why a separate dosage matrix
----------------------------
The simulation carries **standardized** genotypes ``X0`` (per-SNP
z-scores, continuous, roughly -3.5 .. +4.5 after LD-Cholesky and
batch-effect injection).  REGENIE reads PLINK .bed (hard calls only),
GMMAT's fast C score test reads .bed, and GENESIS's
``MatrixGenotypeReader`` wants integer 0/1/2.  So we derive a hard-call
dosage matrix ``D`` from ``X0`` by a **monotone per-SNP quantile map**
onto HWE genotype proportions at a per-SNP MAF drawn from U(0.05, 0.5)
— the same marginal MAF distribution the underlying arrays were built
from (``genotype_simulation.simulate_genotypes``).

Because the map is monotone within each SNP and applied jointly across
all N subjects, every between-subject ordering that carries signal
(population mean shift, sequencing-center batch shift, causal effect)
survives; only the marginal is re-discretised.

Anchor method
-------------
Discretisation is a real (if mild) transformation, so we also report
``combined_pcs_dosage``: the *same* vectorised logistic score test used
by ``combined_pcs``, run on ``D`` instead of ``X0``.  That is the
apples-to-apples anchor for REGENIE / GMMAT / GENESIS — same data,
different method — while ``combined_pcs`` vs ``combined_pcs_dosage``
quantifies the cost of discretisation itself.

Relatedness / GRM
-----------------
The simulation contains **no cryptic relatedness by construction**:
subjects are drawn i.i.d. within population.  Handing GMMAT/GENESIS an
identity (or no) kinship would collapse both to a re-parameterised GLM
and make the comparison uninformative.  We therefore supply an
**empirical GCTA-style GRM** computed from the simulated genotypes,
``G = Z Z' / M`` with ``Z_ij = (D_ij - 2 p_j) / sqrt(2 p_j (1 - p_j))``.
That matrix is not noise: it encodes the two-population ancestry split
and the center batch structure, which is exactly the confounding a
mixed model is supposed to absorb.  This is stated in the run README.
"""
from __future__ import annotations

import os
import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Dict, List, Optional, Sequence

import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

VALID_COMPARATORS = ("regenie", "gmmat", "genesis")

# Extra method keys produced when comparators are enabled.  ``combined_pcs_dosage``
# is always produced (it is the same-data anchor and costs ~1 s).
DOSAGE_ANCHOR = "combined_pcs_dosage"

REGENIE_BIN = os.environ.get("REGENIE_BIN", "${HOME_DIR}/bin/regenie")
RSCRIPT_BIN = os.environ.get("RSCRIPT_BIN", "Rscript")
R_MODULE = os.environ.get("R_MODULE", "r/4.4.0")
COMPARATORS_R = str(Path(__file__).resolve().parent / "synthetic_comparator_designs.R")

# LD block size — kept in sync with simulation_power.LD_BLOCK_SIZE.  Blocks are
# mapped onto pseudo-chromosomes 1..n_blocks so that REGENIE's LOCO makes sense
# (a single-chromosome dataset would leave step-2 with no leave-one-out
# prediction) and so that no causal-block signal leaks into the null-block
# background through the step-1 ridge.
LD_BLOCK_SIZE = 2000

MAF_LO, MAF_HI = 0.05, 0.50


def comparators_from_env(value: Optional[str] = None) -> List[str]:
    """Parse the COMPARATORS env var into a validated list (may be empty)."""
    raw = os.environ.get("COMPARATORS", "") if value is None else value
    out: List[str] = []
    for tok in raw.split(","):
        tok = tok.strip().lower()
        if not tok:
            continue
        if tok not in VALID_COMPARATORS:
            raise ValueError(
                f"Unknown comparator {tok!r}; valid: {VALID_COMPARATORS}"
            )
        if tok not in out:
            out.append(tok)
    return out


def scratch_root() -> Path:
    """Node-local scratch for the per-config genotype files."""
    for key in ("SIM_SCRATCH", "TMPDIR", "SLURM_TMPDIR"):
        v = os.environ.get(key)
        if v and os.path.isdir(v):
            return Path(v)
    return Path("/tmp")


# ---------------------------------------------------------------------------
# 1. Dosage path
# ---------------------------------------------------------------------------

def dosage_from_standardized(
    X0: np.ndarray,
    seed: int,
    maf_lo: float = MAF_LO,
    maf_hi: float = MAF_HI,
    row_chunk: int = 2000,
) -> np.ndarray:
    """
    Monotone per-SNP quantile map from standardized values to 0/1/2 hard calls.

    For SNP j a MAF ``p_j ~ U(maf_lo, maf_hi)`` is drawn, giving HWE counts
    ``n0 = round(N (1-p)^2)`` and ``n1 = round(N 2p(1-p))``.  The n0 smallest
    values of row j become 0, the next n1 become 1, the rest become 2.

    Assignment is by RANK, not by comparison against the n0-th and n1-th order
    statistics.  The thresholding form is wrong whenever a row has ties: X0 is
    not tie-free — the first SNP of each LD-Cholesky block takes only about
    three distinct values — and when the k0-th order statistic equals the row
    maximum, a strict ``>`` comparison collapses the whole row to a constant.
    A monomorphic column makes REGENIE step 1 abort outright
    ("SNP snpXXXX has low variance"), which sets the ENTIRE regenie array for
    that unit to NaN and silently deletes all 14,000 of its tests.  In the
    first launch this hit 5 of 7 array tasks, every one of them in the
    ``2admixed`` scenarios, which would have left REGENIE's per-scenario
    numbers computed over a non-randomly thinned subset of replicates.

    Ranking makes the counts exact by construction: ``k0`` and ``k1`` are
    clipped so every row gets at least one of each genotype, so no row can be
    monomorphic and the realised MAF matches the target exactly.

    Returns (M, N) int8.
    """
    X0 = np.asarray(X0)
    M, N = X0.shape
    rng = np.random.default_rng(int(seed))
    maf = rng.uniform(maf_lo, maf_hi, size=M)

    n0 = np.rint(N * (1.0 - maf) ** 2).astype(np.int64)
    n1 = np.rint(N * 2.0 * maf * (1.0 - maf)).astype(np.int64)
    k0 = np.clip(n0, 1, N - 2)                 # rank of the first "1"
    k1 = np.clip(n0 + n1, k0 + 1, N - 1)       # rank of the first "2"

    # Tie-break jitter: X0 is standardized, so genuine differences are O(1e-2)
    # and 1e-9 orders exact duplicates arbitrarily without perturbing any real
    # ordering. Seeded separately so it cannot alias the MAF draw.
    jrng = np.random.default_rng(int(seed) + 982451653)

    D = np.empty((M, N), dtype=np.int8)
    for start in range(0, M, row_chunk):
        stop = min(start + row_chunk, M)
        chunk = np.asarray(X0[start:stop], dtype=np.float64)
        chunk = chunk + 1e-9 * jrng.standard_normal(chunk.shape)
        rank = np.argsort(np.argsort(chunk, axis=1), axis=1)
        D[start:stop] = (
            (rank >= k0[start:stop, None]).astype(np.int8)
            + (rank >= k1[start:stop, None]).astype(np.int8)
        )
    return D


def empirical_grm(D: np.ndarray, block: int = 4000) -> np.ndarray:
    """
    GCTA-style empirical GRM: G = Z Z' / M, Z = (D - 2p) / sqrt(2p(1-p)).

    Computed in row blocks to bound memory.  Returns (N, N) float64.
    """
    M, N = D.shape
    G = np.zeros((N, N), dtype=np.float64)
    used = 0
    for start in range(0, M, block):
        stop = min(start + block, M)
        sub = np.asarray(D[start:stop], dtype=np.float32)
        p = sub.mean(axis=1, keepdims=True) / 2.0
        sd = np.sqrt(2.0 * p * (1.0 - p))
        keep = (sd[:, 0] > 1e-8)
        if not keep.any():
            continue
        Z = (sub[keep] - 2.0 * p[keep]) / sd[keep]
        G += (Z.T @ Z).astype(np.float64)
        used += int(keep.sum())
    if used == 0:
        raise ValueError("empirical_grm: no polymorphic SNPs")
    G /= used
    # Symmetrise (guards float32 accumulation asymmetry) — PSD by construction.
    G = 0.5 * (G + G.T)
    return G


# ---------------------------------------------------------------------------
# 2. PLINK1 binary writer
# ---------------------------------------------------------------------------

# PLINK1 .bed 2-bit codes, A1 = allele in bim column 5 = the counted allele.
#   00 -> hom A1 (dosage 2)   01 -> missing   10 -> het (1)   11 -> hom A2 (0)
_BED_LUT = np.array([0b11, 0b10, 0b00], dtype=np.uint8)  # index by dosage 0,1,2


def write_plink(
    D: np.ndarray,
    df: pd.DataFrame,
    prefix: str,
    block_size: int = LD_BLOCK_SIZE,
) -> None:
    """Write .bed/.bim/.fam for the (M, N) int8 dosage matrix."""
    M, N = D.shape
    pfx = Path(prefix)

    idx = np.arange(M)
    chrom = idx // block_size + 1
    pos = idx % block_size + 1
    bim = pd.DataFrame(
        {
            "chr": chrom,
            "id": [f"snp{i}" for i in idx],
            "cm": 0,
            "bp": pos,
            "a1": "A",
            "a2": "G",
        }
    )
    bim.to_csv(pfx.with_suffix(".bim"), sep="\t", header=False, index=False)

    iid = _sample_ids(N)
    fam = pd.DataFrame(
        {
            "fid": iid,
            "iid": iid,
            "pid": 0,
            "mid": 0,
            "sex": np.where(df["gender"].to_numpy() == 1, 1, 2),
            "phe": -9,
        }
    )
    fam.to_csv(pfx.with_suffix(".fam"), sep=" ", header=False, index=False)

    npad = (-N) % 4
    with open(pfx.with_suffix(".bed"), "wb") as fh:
        fh.write(bytes([0x6C, 0x1B, 0x01]))     # magic + SNP-major
        for start in range(0, M, 2000):
            stop = min(start + 2000, M)
            codes = _BED_LUT[np.asarray(D[start:stop], dtype=np.int8)]
            if npad:
                codes = np.concatenate(
                    [codes, np.zeros((stop - start, npad), dtype=np.uint8)], axis=1
                )
            c = codes.reshape(stop - start, -1, 4)
            packed = (
                c[:, :, 0]
                | (c[:, :, 1] << 2)
                | (c[:, :, 2] << 4)
                | (c[:, :, 3] << 6)
            ).astype(np.uint8)
            fh.write(packed.tobytes(order="C"))


def _sample_ids(n: int) -> List[str]:
    return [f"I{i:06d}" for i in range(n)]


def _pc_cols(df: pd.DataFrame) -> List[str]:
    pcs = [c for c in df.columns if c.upper().startswith("PC")]
    return sorted(pcs, key=lambda x: int("".join(ch for ch in x if ch.isdigit()) or 0))


def write_pheno_covar(df: pd.DataFrame, workdir: Path) -> Dict[str, Path]:
    """Write REGENIE pheno/covar files and a combined CSV for the R tools."""
    iid = _sample_ids(len(df))
    pcs = _pc_cols(df)

    pheno = pd.DataFrame({"FID": iid, "IID": iid, "Y": df["Y"].to_numpy().astype(int)})
    ppath = workdir / "pheno.txt"
    pheno.to_csv(ppath, sep=" ", index=False)

    cov = pd.DataFrame({"FID": iid, "IID": iid})
    cov["center"] = df["center"].to_numpy().astype(int)
    cov["age"] = df["age"].to_numpy()
    cov["gender"] = df["gender"].to_numpy().astype(int)
    for c in pcs:
        cov[c] = df[c].to_numpy()
    cpath = workdir / "covar.txt"
    cov.to_csv(cpath, sep=" ", index=False)

    rdat = cov.copy()
    rdat.insert(0, "id", iid)
    rdat["Y"] = df["Y"].to_numpy().astype(int)
    rpath = workdir / "rpheno.csv"
    rdat.drop(columns=["FID", "IID"]).to_csv(rpath, index=False)

    return {"pheno": ppath, "covar": cpath, "rpheno": rpath, "pcs": pcs}


# ---------------------------------------------------------------------------
# 3. REGENIE
# ---------------------------------------------------------------------------

def run_regenie(
    workdir: Path,
    prefix: Path,
    files: Dict,
    M: int,
    threads: int = 1,
    step1_bsize: int = 1000,
    step2_bsize: int = 1000,
) -> np.ndarray:
    """
    REGENIE step 1 (--bt --loocv) then step 2 (--bt --firth --approx).

    Returns an M-length p-value array aligned to variant order (NaN where
    REGENIE dropped the variant).
    """
    out1 = workdir / "fit1"
    cmd1 = [
        REGENIE_BIN, "--step", "1",
        "--bed", str(prefix),
        "--phenoFile", str(files["pheno"]),
        "--covarFile", str(files["covar"]),
        "--catCovarList", "center",
        "--bt", "--loocv",
        "--bsize", str(step1_bsize),
        "--lowmem", "--lowmem-prefix", str(workdir / "tmp_rg"),
        "--threads", str(threads),
        "--out", str(out1),
    ]
    _run(cmd1, workdir / "regenie_step1.log")

    pred = Path(str(out1) + "_pred.list")
    if not pred.exists():
        raise RuntimeError("regenie step 1 produced no _pred.list")

    out2 = workdir / "assoc"
    cmd2 = [
        REGENIE_BIN, "--step", "2",
        "--bed", str(prefix),
        "--phenoFile", str(files["pheno"]),
        "--covarFile", str(files["covar"]),
        "--catCovarList", "center",
        "--bt", "--firth", "--approx", "--pThresh", "0.05",
        "--pred", str(pred),
        "--bsize", str(step2_bsize),
        "--threads", str(threads),
        "--out", str(out2),
    ]
    _run(cmd2, workdir / "regenie_step2.log")

    res_path = Path(str(out2) + "_Y.regenie")
    if not res_path.exists():
        cands = sorted(workdir.glob("assoc_*.regenie"))
        if not cands:
            raise RuntimeError("regenie step 2 produced no .regenie output")
        res_path = cands[0]

    res = pd.read_csv(res_path, sep=r"\s+", comment="#")
    p = np.full(M, np.nan, dtype=np.float64)
    ids = res["ID"].astype(str).str.slice(3).astype(int).to_numpy()
    logp = res["LOG10P"].to_numpy(dtype=np.float64)
    ok = np.isfinite(logp) & (ids >= 0) & (ids < M)
    p[ids[ok]] = np.power(10.0, -logp[ok])
    return p


def _run(cmd: Sequence[str], logfile: Path) -> None:
    with open(logfile, "w") as fh:
        proc = subprocess.run(
            list(cmd), stdout=fh, stderr=subprocess.STDOUT, check=False
        )
    if proc.returncode != 0:
        tail = ""
        try:
            tail = "".join(open(logfile).readlines()[-30:])
        except OSError:
            pass
        raise RuntimeError(
            f"command failed (rc={proc.returncode}): {' '.join(map(str, cmd))}\n{tail}"
        )


# ---------------------------------------------------------------------------
# 4. GMMAT + GENESIS (single Rscript call — they share the GRM)
# ---------------------------------------------------------------------------

def run_r_comparators(
    workdir: Path,
    prefix: Path,
    files: Dict,
    D: np.ndarray,
    G: np.ndarray,
    which: Sequence[str],
    threads: int = 1,
) -> Dict[str, np.ndarray]:
    """Run GMMAT::glmm.score and/or GENESIS::assocTestSingle. Returns dict."""
    M, N = D.shape

    geno_path = workdir / "geno.i8"
    with open(geno_path, "wb") as fh:
        fh.write(np.ascontiguousarray(D, dtype=np.int8).tobytes(order="C"))
    grm_path = workdir / "grm.f64"
    with open(grm_path, "wb") as fh:
        fh.write(np.ascontiguousarray(G, dtype=np.float64).tobytes(order="C"))

    cmd = [
        RSCRIPT_BIN, "--vanilla", COMPARATORS_R,
        "--workdir", str(workdir),
        "--bfile", str(prefix),
        "--geno", str(geno_path),
        "--grm", str(grm_path),
        "--pheno", str(files["rpheno"]),
        "--m", str(M),
        "--n", str(N),
        "--block-size", str(LD_BLOCK_SIZE),
        "--which", ",".join(which),
    ]
    env = os.environ.copy()
    env.setdefault("OMP_NUM_THREADS", str(threads))
    env.setdefault("OPENBLAS_NUM_THREADS", str(threads))
    logfile = workdir / "R_comparators.log"
    with open(logfile, "w") as fh:
        proc = subprocess.run(
            cmd, stdout=fh, stderr=subprocess.STDOUT, check=False, env=env
        )
    if proc.returncode != 0:
        tail = "".join(open(logfile).readlines()[-40:])
        raise RuntimeError(f"Rscript comparators failed (rc={proc.returncode})\n{tail}")

    out: Dict[str, np.ndarray] = {}
    for name in which:
        f = workdir / f"{name}_pvals.txt"
        if not f.exists():
            raise RuntimeError(f"{name}: expected output {f} missing")
        arr = np.loadtxt(f, dtype=np.float64)
        if arr.shape[0] != M:
            raise RuntimeError(f"{name}: got {arr.shape[0]} p-values, expected {M}")
        out[name] = arr
    return out


# ---------------------------------------------------------------------------
# 5. Top-level driver
# ---------------------------------------------------------------------------

def run_comparators(
    df: pd.DataFrame,
    X0: np.ndarray,
    which: Sequence[str],
    dosage_seed: int,
    tag: str = "cfg",
    threads: Optional[int] = None,
    keep_workdir: bool = False,
    return_dosage: bool = False,
    grm_snp_idx: Optional[np.ndarray] = None,
):
    """
    Build the dosage view of ``X0`` and run the requested external tools.

    Returns ``dict method -> (M,) p-value array``.  Always includes
    ``combined_pcs_dosage`` (the Python score test on the same hard calls).
    On tool failure the offending method's array is all-NaN and a warning is
    printed; the rest of the replicate still completes.

    ``grm_snp_idx`` restricts the SNPs used to build the empirical GRM.  Pass the
    null-block index (SNPs outside the causal LD block).  Without it the GRM is
    built from every SNP including the causal block, which puts GMMAT and
    GENESIS at a systematic disadvantage against REGENIE: REGENIE's step-2 test
    uses a leave-one-chromosome-out predictor, so the causal block never enters
    the predictor for a variant tested on another pseudo-chromosome, whereas a
    whole-genome GRM carries the phenotype-correlated causal signal into every
    test.  The result is over-correction — in the smoke run the two mixed models
    returned null-block lambda_GC = 0.749 and roughly half the power of the
    other methods, which is an artefact of the asymmetry rather than a property
    of either method.
    """
    try:
        from synthetic_null_simulation import assoc_Y_on_SNP_plus_covars
    except ImportError:                                       # imported as a package
        from simulation.simulation_2pop_qqplot import assoc_Y_on_SNP_plus_covars

    which = [w for w in which if w in VALID_COMPARATORS]
    M, N = X0.shape
    if threads is None:
        threads = int(os.environ.get("SLURM_CPUS_PER_TASK", "1"))

    D = dosage_from_standardized(X0, seed=dosage_seed)

    out: Dict[str, np.ndarray] = {}
    pcs = _pc_cols(df)
    out[DOSAGE_ANCHOR] = assoc_Y_on_SNP_plus_covars(
        D.astype(np.float64), df, covariate_cols=["center", "age", "gender"] + pcs
    ).values

    if not which:
        return (out, D) if return_dosage else out

    root = scratch_root()
    workdir = Path(tempfile.mkdtemp(prefix=f"riga_{tag}_", dir=str(root)))
    try:
        prefix = workdir / "geno"
        write_plink(D, df, str(prefix))
        files = write_pheno_covar(df, workdir)

        if "regenie" in which:
            try:
                out["regenie"] = run_regenie(workdir, prefix, files, M, threads=threads)
            except Exception as exc:                      # noqa: BLE001
                print(f"  [comparators] regenie FAILED: {exc}")
                out["regenie"] = np.full(M, np.nan)

        r_which = [w for w in which if w in ("gmmat", "genesis")]
        if r_which:
            try:
                if grm_snp_idx is None:
                    print("  [comparators] WARNING: GRM built from ALL SNPs "
                          "(causal block included) - mixed models will be "
                          "over-corrected relative to REGENIE's LOCO")
                    G = empirical_grm(D)
                else:
                    gi = np.asarray(grm_snp_idx, dtype=np.int64)
                    print(f"  [comparators] GRM from {gi.size} null-block SNPs "
                          f"of {D.shape[0]} (causal block excluded)")
                    G = empirical_grm(D[gi])
                res = run_r_comparators(
                    workdir, prefix, files, D, G, r_which, threads=threads
                )
                out.update(res)
            except Exception as exc:                      # noqa: BLE001
                print(f"  [comparators] R tools FAILED: {exc}")
                for w in r_which:
                    out.setdefault(w, np.full(M, np.nan))
    finally:
        if not keep_workdir:
            shutil.rmtree(workdir, ignore_errors=True)
        else:
            print(f"  [comparators] kept workdir {workdir}")

    return (out, D) if return_dosage else out
