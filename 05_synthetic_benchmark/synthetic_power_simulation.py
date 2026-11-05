"""
synthetic_power_simulation.py — Power analysis simulation for multi-ancestry GWAS methods
================================================================================

Overview
--------
This module extends the null simulation (``synthetic_null_simulation.py``) by
introducing causal variants with pre-specified effect sizes.  The same seven
ancestry-composition scenarios and six GWAS methods are evaluated, now for
statistical *power* (fraction of causal SNPs detected) rather than type-I
error calibration.

All dataset-building infrastructure (phenotype simulation, cohort assembly,
batch-effect injection, PCA, GMM labelling, association tests) is imported
directly from ``synthetic_null_simulation.py``.  Only the phenotype generation
is replaced by ``simulate_population_with_preset_genetics()``, which accepts
a pre-computed per-SNP effect vector instead of sampling one internally.

Cohort design (inherited from simulation_2pop_qqplot)
-----------------------------------------------------
Subjects are assembled into N_COHORTS=3 multi-ancestry sequencing cohorts
(study labels 1–3) via ``sample_subpopulations()``.  Each cohort contains
subjects from *both* Pop A and Pop B, so study membership is not a proxy
for ancestry.  See the ``sample_subpopulations`` docstring in
``synthetic_null_simulation.py`` for full design rationale.

Effect-size parameterisation
-----------------------------
Per-variant log-odds effect sizes are derived from liability-scale
heritability h² and causal variant count K via:

    β_causal = sqrt( h² · π² / (3 · K · (1 − h²)) )

where π²/3 is the logistic residual variance.  This gives an approximate
liability-scale h².  The default sweep is h² ∈ {0.3, 0.5, 0.7} with K = 10,
corresponding to per-variant odds ratios of ≈ 1.45, 1.78, and 2.40.

Causal variant placement
------------------------
Causal variants are restricted to a single randomly chosen LD block
(contiguous window of LD_BLOCK_SIZE=2 000 SNPs).  FPR and lambda_GC_null
are then computed exclusively on SNPs *outside* that block, so calibration
metrics reflect the test's behaviour under the null rather than LD tagging
of the causal signal.

Genetic architecture modes
--------------------------
  shared        Identical effect vectors in both populations (β_A = β_B).
  pop_specific  A fraction ``pop_specific_frac`` of causal variants are
                private to one population (zero effect in the other).
  heterogeneous All variants shared, but per-population effects are drawn
                from a bivariate normal with cross-population genetic
                correlation rG (default 0.8).
  mixed         Population-specific causal sets AND heterogeneous effects
                for the shared portion.

Six scenarios
---------------
  2pure_homo        AFR vs EUR, homogeneous age/gender
  2pure_hetero      AFR vs EUR, heterogeneous age/gender
  1pure_1admixed_homo   AFR vs CEU50%-YRI50%, homogeneous
  1pure_1admixed_hetero AFR vs CEU50%-YRI50%, heterogeneous
  2admixed_homo         CEU50%-YRI50% vs CEU50%-CHS50%, homogeneous
  2admixed_hetero       CEU50%-YRI50% vs CEU50%-CHS50%, heterogeneous

  Note on the 2-admixed scenarios (``2admixed_homo``, ``2admixed_hetero``)
  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  CEU5-YRI5 (50 % CEU + 50 % YRI) and CEU5-CHS5 (50 % CEU + 50 % CHS)
  share half their ancestry and are simulated from chromosome 22 only.
  Chr22 lacks sufficient ancestry-informative markers to separate the YRI
  vs CHS components: the two populations overlap extensively in PC space,
  GMM accuracy approaches 50 % (chance), and all six methods perform
  similarly regardless of effect size.  This is an expected consequence of
  the genetic data available — not a software defect.  The 2-admixed
  scenarios serve as a stress test for the methods, not as a scenario where
  meaningful method differences are expected.

Power metrics (per method, scenario, h², K)
-------------------------------------------
  power_{5e8,1e5,1e3}      Fraction of causal SNPs with p < threshold
  fpr_{5e8,1e5,1e3}        Fraction of null SNPs (outside causal block)
                            with p < threshold
  lambda_gc                 Genomic inflation factor (full set)
  lambda_gc_null            Lambda_GC restricted to null SNPs
  power_gc_{5e8,1e5,1e3}   Power after genomic-control correction
                            (test statistics divided by lambda_gc_null)
  pop_accuracy              GMM population assignment accuracy vs true label
  n_null_idx                Number of SNPs outside the causal LD block

Output files (per replicate)
-----------------------------
  power_run_{id:03d}_metrics.csv   One row per scenario × config × method
  power_run_{id:03d}_pvals.npz     Raw M-dimensional p-value arrays

Usage (command line)
--------------------
  python synthetic_power_simulation.py --replicate-id 1 --seed 42 \\
      --h2-values 0.3,0.5,0.7 --num-causal-values 10 \\
      --genetic-arch shared --out-dir outputs/power_runs

SLURM array jobs: see run_power_array.sbatch.
"""

import os
import sys
import numpy as np
import pandas as pd
from typing import Dict, List, Tuple, Optional
from pathlib import Path
import argparse
from scipy.stats import chi2 as chi2_dist
from sklearn.decomposition import PCA

_src_dir = Path(__file__).resolve().parent.parent
if str(_src_dir) not in sys.path:
    sys.path.insert(0, str(_src_dir))
from simulation.simulation_2pop_qqplot import (
    _r_truncnorm,
    sample_subpopulations,
    assign_population_labels,
    _assign_centers,
    inject_batch_effects,
    maybe_inject_study_batch,
    verify_batch_heterogeneity,
    pvals_combined_study,
    pvals_combined_pcs,
    pvals_combined_population,
    pvals_combined_study_pcs,
    pvals_combined_population_pcs,
    pvals_matched_population,
    pvals_matched_study,
    pvals_matched_study_center,
    pvals_matched_population_pcs,
    pvals_matched_study_pcs,
    match_cases_controls_with_center,
    match_cases_controls_with_study,
    match_cases_controls_with_study_center,
    INTERACTION_MODE,
    _load_default_arrays,
    N_PCS,
    RNG_SEED,
    POP_A_PREV,
    POP_B_PREV,
)
from plot_simulation_main_figure import compute_lambda_gc

def h2_to_beta_causal(h2: float, num_causal: int) -> float:
    """Per-variant effect size on standardised genotypes for a target h²."""
    return np.sqrt(h2 * np.pi**2 / (3.0 * num_causal * (1.0 - h2)))


# Default configs use strong per-SNP effects (OR ≈ 1.45–2.4 for h2=0.3–0.7, K=10).
# That makes the causal LD block light up strongly (many tagging SNPs get small p-values).
# For less LD-tagging dominance / more realistic effect sizes, use lower h² or larger K
# (e.g. h2=0.1 K=10 → OR≈1.2; h2=0.3 K=100 → OR≈1.04).
POWER_CONFIGS = [
    {"h2": 0.3, "num_causal": 10, "label": "h2=0.3, K=10"},
    {"h2": 0.5, "num_causal": 10, "label": "h2=0.5, K=10"},
    {"h2": 0.7, "num_causal": 10, "label": "h2=0.7, K=10"},
]
# Optional: more moderate effect sizes (less LD tagging)
# POWER_CONFIGS_MODERATE = [
#     {"h2": 0.1, "num_causal": 10, "label": "h2=0.1, K=10"},   # OR ~1.2
#     {"h2": 0.2, "num_causal": 10, "label": "h2=0.2, K=10"},   # OR ~1.3
#     {"h2": 0.3, "num_causal": 100, "label": "h2=0.3, K=100"},  # OR ~1.04
# ]

METHOD_ORDER = [
    "combined_study",
    "combined_pcs",
    "combined_population",
    "combined_study_pcs",
    "combined_population_pcs",
    "matched_population",
    "matched_study",
    "matched_study_center",
    "matched_population_pcs",
    "matched_study_pcs",
]
# NOTE: METHOD_ORDER and METHOD_ORDER_FULL are documentation only.  Nothing in
# the repo imports either of them; the run is driven entirely by the dict that
# run_all_methods_raw returns.  Adding a name here wires nothing.

# Arena A comparators (REGENIE / GMMAT / GENESIS).  Opt-in via
#   COMPARATORS=regenie,gmmat,genesis
# When the variable is unset, _COMPARATORS is [] and the six-method path above
# is completely unchanged.
_sim_dir = Path(__file__).resolve().parent
if str(_sim_dir) not in sys.path:
    sys.path.insert(0, str(_sim_dir))
try:
    from synthetic_comparator_designs import comparators_from_env, DOSAGE_ANCHOR
    _COMPARATORS = comparators_from_env()
except Exception as _exc:                                     # noqa: BLE001
    print(f"[simulation_power] comparators unavailable: {_exc}")
    _COMPARATORS, DOSAGE_ANCHOR = [], "combined_pcs_dosage"

METHOD_ORDER_FULL = METHOD_ORDER + ([DOSAGE_ANCHOR] + _COMPARATORS if _COMPARATORS else [])

# Number of LD blocks reserved as GRM markers for the Arena A comparators and
# excluded from evaluation.  2 of 7 blocks leaves 4 evaluable null blocks
# (~8,000 SNPs/unit) and gives the GRM ~4,000 markers, ample for the
# two-population split and batch structure it has to encode.
N_GRM_BLOCKS = 2

SCENARIO_ORDER = [
    "2pure_homo",
    "2pure_hetero",
    "1pure_1admixed_homo",
    "1pure_1admixed_hetero",
    "2admixed_homo",
    "2admixed_hetero",
]


# ---------------------------------------------------------------------------
# Genetic architecture modes
# ---------------------------------------------------------------------------

GENETIC_ARCHS = ["shared", "pop_specific", "heterogeneous", "mixed"]

# LD block size (SNPs per block) for restricting causal variants to one block.
# FPR and lambda_GC null are computed on SNPs outside this block only.
LD_BLOCK_SIZE = 2000


def get_ld_block_assignments(M: int, block_size: int = LD_BLOCK_SIZE) -> np.ndarray:
    """Return block index per SNP (0, 1, 2, ...) for M SNPs."""
    return np.arange(M, dtype=np.intp) // block_size


def generate_shared_causal_effects(
    M: int,
    beta_causal: float,
    num_causal: int,
    rng: np.random.Generator,
    allowed_indices: Optional[np.ndarray] = None,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Select causal SNP indices and assign per-variant log-odds effects.

    Each causal variant receives a random sign (+/-) times beta_causal,
    so that effects are not all in the same direction.  Non-causal
    variants have beta = 0.

    Parameters
    ----------
    M : total number of SNPs
    beta_causal : absolute per-variant effect size on *standardised*
                  genotypes (log-odds scale).  OR ≈ exp(beta_causal).
    num_causal : number of causal variants
    rng : numpy Generator for reproducibility
    allowed_indices : if provided, causal SNPs are chosen only from these
                      indices (e.g. one LD block).

    Returns
    -------
    causal_idx : array of causal SNP indices
    beta_snp   : M-length array of per-SNP effects
    """
    if allowed_indices is not None:
        if len(allowed_indices) < num_causal:
            raise ValueError(
                f"allowed_indices has {len(allowed_indices)} SNPs but num_causal={num_causal}; "
                "use a larger LD block or fewer causal variants."
            )
        causal_idx = rng.choice(allowed_indices, size=num_causal, replace=False)
    else:
        causal_idx = rng.choice(M, size=num_causal, replace=False)
    signs = rng.choice([-1, 1], size=num_causal)
    beta_snp = np.zeros(M)
    beta_snp[causal_idx] = signs * beta_causal
    return causal_idx, beta_snp


def generate_ancestry_aware_causal_effects(
    M: int,
    beta_causal: float,
    num_causal: int,
    rng: np.random.Generator,
    genetic_arch: str = "shared",
    genetic_corr: float = 0.8,
    pop_specific_frac: float = 0.5,
    allowed_indices: Optional[np.ndarray] = None,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Generate per-population causal effect vectors under various genetic
    architecture models.

    Parameters
    ----------
    M : total number of SNPs
    beta_causal : per-variant effect magnitude (from h2_to_beta_causal)
    num_causal : total causal variants across both populations
    rng : numpy Generator
    genetic_arch : one of 'shared', 'pop_specific', 'heterogeneous', 'mixed'
    genetic_corr : cross-population genetic correlation rG
        (used by 'heterogeneous' and 'mixed')
    pop_specific_frac : fraction of causal variants that are
        population-specific (used by 'pop_specific' and 'mixed')
    allowed_indices : if provided, causal SNPs are chosen only from these
                      indices (e.g. one LD block).

    Returns
    -------
    causal_idx : union of all causal SNP indices
    beta_snp_A : M-length effect vector for population A
    beta_snp_B : M-length effect vector for population B
    """
    if genetic_arch not in GENETIC_ARCHS:
        raise ValueError(
            f"Unknown genetic_arch '{genetic_arch}'; "
            f"choose from {GENETIC_ARCHS}"
        )

    if allowed_indices is not None:
        if len(allowed_indices) < num_causal:
            raise ValueError(
                f"allowed_indices has {len(allowed_indices)} SNPs but num_causal={num_causal}; "
                "use a larger LD block or fewer causal variants."
            )
        causal_idx = rng.choice(allowed_indices, size=num_causal, replace=False)
    else:
        causal_idx = rng.choice(M, size=num_causal, replace=False)
    beta_snp_A = np.zeros(M)
    beta_snp_B = np.zeros(M)

    if genetic_arch == "shared":
        signs = rng.choice([-1, 1], size=num_causal)
        beta_snp_A[causal_idx] = signs * beta_causal
        beta_snp_B[causal_idx] = signs * beta_causal
        return causal_idx, beta_snp_A, beta_snp_B

    # -- split causal set into shared / popA-only / popB-only --
    use_pop_specific = genetic_arch in ("pop_specific", "mixed")
    n_shared = (
        max(1, round(num_causal * (1.0 - pop_specific_frac)))
        if use_pop_specific else num_causal
    )
    n_pop_specific = num_causal - n_shared
    n_A_only = n_pop_specific // 2
    n_B_only = n_pop_specific - n_A_only

    rng.shuffle(causal_idx)
    shared_idx = causal_idx[:n_shared]
    A_only_idx = causal_idx[n_shared : n_shared + n_A_only]
    B_only_idx = causal_idx[n_shared + n_A_only :]

    # -- shared variants --
    use_hetero = genetic_arch in ("heterogeneous", "mixed")
    if n_shared > 0:
        if use_hetero:
            rG = float(np.clip(genetic_corr, -1.0, 1.0))
            cov = beta_causal ** 2 * np.array([[1.0, rG], [rG, 1.0]])
            effects = rng.multivariate_normal([0.0, 0.0], cov, size=n_shared)
            beta_snp_A[shared_idx] = effects[:, 0]
            beta_snp_B[shared_idx] = effects[:, 1]
        else:
            signs = rng.choice([-1, 1], size=n_shared)
            beta_snp_A[shared_idx] = signs * beta_causal
            beta_snp_B[shared_idx] = signs * beta_causal

    # -- population-specific variants (zero in the other population) --
    if n_A_only > 0:
        signs_A = rng.choice([-1, 1], size=n_A_only)
        beta_snp_A[A_only_idx] = signs_A * beta_causal

    if n_B_only > 0:
        signs_B = rng.choice([-1, 1], size=n_B_only)
        beta_snp_B[B_only_idx] = signs_B * beta_causal

    return causal_idx, beta_snp_A, beta_snp_B


# ---------------------------------------------------------------------------
# Phenotype simulation with directly-specified genetic effects
# ---------------------------------------------------------------------------

def simulate_population_with_preset_genetics(
    X_std: np.ndarray,
    beta_snp: np.ndarray,
    target_prev: float = 0.109,
    beta_age: float = 0.03,
    beta_gender: float = -0.2,
    mode: str = "homogeneous",
    study_label: int = 0,
    population_label: Optional[str] = None,
    unmeasured_env_effect: float = 0.0,
    rng: Optional[np.random.Generator] = None,
) -> pd.DataFrame:
    """
    Simulate one population with a pre-computed beta_snp vector.

    The genetic linear predictor g = X_std^T @ beta_snp is added directly
    to the log-odds (no rescaling).  This gives exact control over
    per-variant effect sizes as specified.

    unmeasured_env_effect
        Log-odds shift from an unmeasured environmental factor that
        differs between populations (e.g. diet, SES, environmental
        exposure).  NOT included in the GWAS regression covariates,
        so it creates residual confounding in combined analyses.
        Matching on population balances it because it is constant
        within a population.
    """
    rng = rng or np.random.default_rng()
    _M, N = X_std.shape

    geffc = X_std.T @ beta_snp

    if mode == "homogeneous":
        age = rng.uniform(30.0, 80.0, size=N)
        gender = rng.binomial(1, 0.5, size=N)
    elif mode == "heterogeneous_A":
        age = _r_truncnorm(N, 30, 60, 40, 10, rng)
        gender = rng.binomial(1, 0.7, size=N)
    elif mode == "heterogeneous_B":
        age = _r_truncnorm(N, 55, 80, 70, 10, rng)
        gender = rng.binomial(1, 0.3, size=N)
    else:
        raise ValueError(f"Unknown mode: {mode}")

    logit_target = np.log(target_prev / (1 - target_prev))
    mean_cov_effect = beta_age * np.mean(age) + beta_gender * np.mean(gender)
    conf_eff = logit_target - mean_cov_effect - np.mean(geffc)

    log_odds = (
        geffc + conf_eff + beta_age * age + beta_gender * gender
        + unmeasured_env_effect
    )
    prob = 1.0 / (1.0 + np.exp(-log_odds))
    Y = rng.binomial(1, prob)

    df = pd.DataFrame({"age": age, "gender": gender, "Y": Y, "study": study_label})
    if population_label is not None:
        df["population"] = population_label
    return df


# ---------------------------------------------------------------------------
# Dataset builder
# ---------------------------------------------------------------------------

def _build_power_setting_dataset(
    X_std_A: np.ndarray,
    X_std_B: np.ndarray,
    heterogeneity: bool,
    n_components: int,
    seed: int,
    beta_causal: float,
    num_causal: int,
    batch_effect: bool = True,
    batch_frac: float = 0.05,
    batch_size: float = 0.3,
    genetic_arch: str = "shared",
    genetic_corr: float = 0.8,
    pop_specific_frac: float = 0.5,
    scenario_key: str = None,
    use_path_c: bool = False,
) -> Tuple[pd.DataFrame, np.ndarray, np.ndarray, np.ndarray, float, np.ndarray]:
    """
    Build one analysis-ready dataset with ancestry-aware causal variants.

    Causal variants are restricted to a single LD block so that FPR can be
    reported for SNPs not in LD with causal (null_idx = SNPs outside that block).

    Mirrors _build_setting_dataset from the null simulation, adding
    batch effects, center assignment, PCA, and GMM population labels.

    Returns (df, Xstd_s, causal_idx, beta_snp_A, pop_accuracy, null_idx).
    null_idx : indices of SNPs not in the causal LD block (for FPR / lambda_GC null).
    """
    rng = np.random.default_rng(seed)

    # ---- Path C: synthesise per-pop genotypes from real per-pop MAFs/LD
    # In Path C every replicate gets fresh genotypes built from the panel,
    # carrying real per-SNP allele-frequency divergence between Pop A and
    # Pop B (the signal needed for classical stratification).
    if use_path_c:
        from path_c import build_path_c_geno, joint_standardise
        # Pool size: needs to cover the per-pop case-fraction demand from
        # COHORT_*_SIZES × COHORT_CASE_FRACTIONS_{A,B}.  Pop B demands
        # ~609 cases at prev=0.08 → need N ≥ 7600.  Use 8000 with slack.
        N_POP_POOL = 8000
        X_std_A, X_std_B = build_path_c_geno(scenario_key, n_per_pop=N_POP_POOL, seed=seed)

    M = X_std_A.shape[0]

    block_assignments = get_ld_block_assignments(M)
    n_blocks = int(block_assignments.max()) + 1
    causal_block_id = rng.integers(0, n_blocks)
    in_causal_block = block_assignments == causal_block_id
    allowed_indices = np.where(in_causal_block)[0]

    causal_idx, beta_snp_A, beta_snp_B = generate_ancestry_aware_causal_effects(
        M, beta_causal, num_causal, rng,
        genetic_arch=genetic_arch,
        genetic_corr=genetic_corr,
        pop_specific_frac=pop_specific_frac,
        allowed_indices=allowed_indices,
    )

    null_idx = np.where(~in_causal_block)[0]

    mode_A = "heterogeneous_A" if heterogeneity else "homogeneous"
    mode_B = "heterogeneous_B" if heterogeneity else "homogeneous"

    df_A = simulate_population_with_preset_genetics(
        X_std_A, beta_snp_A, mode=mode_A,
        target_prev=POP_A_PREV,
        study_label=0, population_label="PopA",
        rng=rng,
    )
    df_B = simulate_population_with_preset_genetics(
        X_std_B, beta_snp_B, mode=mode_B,
        target_prev=POP_B_PREV,
        study_label=1, population_label="PopB",
        rng=rng,
    )

    df_s, idx_a, idx_b = sample_subpopulations(df_A, df_B, rng=rng)
    Xstd_s = np.concatenate((X_std_A[:, idx_a], X_std_B[:, idx_b]), axis=1)
    # Path C: standardise jointly (across the combined cohort) so the
    # per-SNP Pop signal is preserved for the methods to interact with.
    if use_path_c:
        from path_c import joint_standardise
        Xstd_s = joint_standardise(Xstd_s)

    if batch_effect:
        df_s = _assign_centers(df_s)
        center_labels = df_s["center"].values
        Xstd_s = inject_batch_effects(
            Xstd_s, center_labels,
            frac_affected=batch_frac,
            effect_size=batch_size,
            clip_dosage=False,
            rng=rng,
        )
        # Mechanism B (INTERACTION_MODE 3/4 only).  No-op in every other mode
        # and it never draws from `rng`, so the existing code path is unchanged.
        Xstd_s, _mechb = maybe_inject_study_batch(Xstd_s, df_s, seed)
        verify_batch_heterogeneity(Xstd_s, center_labels, rng=rng)

    pca = PCA(n_components=N_PCS, svd_solver="randomized", random_state=seed)
    pcs = pca.fit_transform(Xstd_s.T)
    pcs_df = pd.DataFrame(pcs, columns=[f"PC{i+1}" for i in range(pcs.shape[1])])
    df_s = pd.merge(df_s, pcs_df, left_index=True, right_index=True)
    df_s = assign_population_labels(df_s, n_components=n_components, random_state=seed)

    pop_accuracy = float((df_s["population"] == df_s["assigned_pop"]).mean())

    return df_s, Xstd_s, causal_idx, beta_snp_A, pop_accuracy, null_idx


# ---------------------------------------------------------------------------
# GWAS runner returning raw (M-dimensional) p-values
# ---------------------------------------------------------------------------

def run_all_methods_raw(
    df: pd.DataFrame, X0: np.ndarray, age_tolerance: int = 1,
    comparators: Optional[List[str]] = None,
    dosage_seed: Optional[int] = None,
    tag: str = "cfg",
    grm_snp_idx: Optional[np.ndarray] = None,
) -> Dict[str, np.ndarray]:
    """
    Return method -> raw M-dimensional p-value array (NaN for failed fits).

    The six pure-Python covariate-set variants are always run on the
    standardized matrix ``X0`` and are unaffected by ``comparators``.

    When ``comparators`` is a non-empty list (from ``COMPARATORS=...``), a
    hard-call dosage view of ``X0`` is derived (see ``synthetic_comparator_designs.py``) and
    the requested external tools are run on it, together with
    ``combined_pcs_dosage`` — the same Python score test on the same hard
    calls, which is the same-data anchor for the external tools.
    """
    out: Dict[str, np.ndarray] = {}
    M = X0.shape[0]
    out["combined_study"] = pvals_combined_study(df, X0).values
    out["combined_pcs"] = pvals_combined_pcs(df, X0).values
    out["combined_population"] = pvals_combined_population(df, X0).values
    out["combined_study_pcs"] = pvals_combined_study_pcs(df, X0).values
    out["combined_population_pcs"] = pvals_combined_population_pcs(df, X0).values
    try:
        out["matched_population"] = pvals_matched_population(
            df, X0, age_tolerance=age_tolerance
        ).values
    except ValueError:
        out["matched_population"] = np.full(M, np.nan)

    # Study-resolution matched arms.  These MUST fill with NaN of length M on
    # failure, not with an empty array: run_single_power_job indexes every
    # method array by causal_idx / eval_null_idx, so a short array raises
    # IndexError and kills the whole replicate.  Exception (not ValueError) is
    # caught deliberately and asymmetrically: study cells are half the size of
    # centre cells, so statsmodels' PerfectSeparationError and singular-QR
    # paths are live risks, and nothing should be allowed to kill a 38-minute
    # replicate at unit 17 of 18.  Failures are not hidden -- an arm that fails
    # often shows up as units_all_nan in S10d_arenaA_tool_status.csv.
    for _name, _fn in (("matched_study", pvals_matched_study),
                       ("matched_study_center", pvals_matched_study_center),
                       ("matched_population_pcs", pvals_matched_population_pcs),
                       ("matched_study_pcs", pvals_matched_study_pcs)):
        try:
            out[_name] = _fn(df, X0, age_tolerance=age_tolerance).values
        except Exception as exc:                                  # noqa: BLE001
            print(f"    [{_name}] unit failed: {type(exc).__name__}: {exc}")
            out[_name] = np.full(M, np.nan)

    if comparators:
        from synthetic_comparator_designs import run_comparators
        out.update(
            run_comparators(
                df, X0, which=comparators,
                dosage_seed=int(dosage_seed if dosage_seed is not None else 0),
                tag=tag,
                grm_snp_idx=grm_snp_idx,
            )
        )
    return out


# ---------------------------------------------------------------------------
# Power / FPR metrics
# ---------------------------------------------------------------------------

# P-value thresholds used for power/FPR reporting (single tuple, one pass per method).
# 5e-7 is added so that pooled-run empirical FPR can be computed at this level;
# ~10,000 null runs provide ~70 expected FPs at 5e-7 (CV ≈ 12%) per scenario.
# 5e-8 cannot be reliably estimated by direct simulation (would require ~80K runs);
# use the pooled-run 95% CI upper bound or QQ extrapolation for that level.
POWER_THRESHOLDS = (5e-8, 5e-7, 1e-5, 1e-3)


def compute_power(raw_pvals: np.ndarray, causal_idx: np.ndarray, threshold: float) -> float:
    """Fraction of causal variants with p < threshold."""
    causal_p = raw_pvals[causal_idx]
    valid = causal_p[np.isfinite(causal_p) & (causal_p > 0) & (causal_p <= 1)]
    if len(valid) == 0:
        return np.nan
    return float(np.mean(valid < threshold))


def compute_fpr(
    raw_pvals: np.ndarray,
    causal_idx: np.ndarray,
    threshold: float,
    null_idx: Optional[np.ndarray] = None,
) -> float:
    """
    Fraction of null variants with p < threshold (false positive rate).

    If null_idx is provided (SNPs not in LD with causal), FPR is computed
    on that set only. Otherwise FPR uses all non-causal SNPs.
    """
    if null_idx is not None:
        nc_p = raw_pvals[null_idx]
    else:
        M = len(raw_pvals)
        non_causal = np.setdiff1d(np.arange(M), causal_idx)
        nc_p = raw_pvals[non_causal]
    valid = nc_p[np.isfinite(nc_p) & (nc_p > 0) & (nc_p <= 1)]
    if len(valid) == 0:
        return np.nan
    return float(np.mean(valid < threshold))


def _lambda_gc_subset(raw_pvals: np.ndarray, indices: np.ndarray) -> float:
    """Lambda_GC computed on a subset of SNP indices."""
    p = raw_pvals[indices]
    valid = p[np.isfinite(p) & (p > 0) & (p <= 1)]
    return compute_lambda_gc(valid)


def compute_power_gc_corrected(
    raw_pvals: np.ndarray,
    causal_idx: np.ndarray,
    non_causal_idx: np.ndarray,
    threshold: float,
    lgc: Optional[float] = None,
) -> float:
    """
    Power after genomic-control correction.

    Standard GWAS practice: estimate lambda_GC from non-causal (null)
    variants, then divide every test statistic by lambda_GC before
    applying the significance threshold.  Well-calibrated methods
    (lambda_GC ≈ 1) keep their raw power; inflated methods pay a
    penalty proportional to their inflation.

    If lgc is provided, it is used instead of estimating from non_causal_idx
    (avoids recomputing when computing at multiple thresholds).
    """
    if lgc is None:
        nc_p = raw_pvals[non_causal_idx]
        nc_valid = nc_p[np.isfinite(nc_p) & (nc_p > 0) & (nc_p <= 1)]
        if nc_valid.size == 0:
            return np.nan
        lgc = compute_lambda_gc(nc_valid)
    if not np.isfinite(lgc) or lgc <= 0:
        return np.nan

    causal_p = raw_pvals[causal_idx]
    ok = np.isfinite(causal_p) & (causal_p > 0) & (causal_p <= 1)
    if ok.sum() == 0:
        return np.nan

    causal_chi2 = chi2_dist.isf(causal_p[ok], df=1)
    corrected_chi2 = causal_chi2 / lgc
    corrected_p = chi2_dist.sf(corrected_chi2, df=1)
    return float(np.mean(corrected_p < threshold))


def _metrics_for_method(
    raw_pvals: np.ndarray,
    causal_idx: np.ndarray,
    null_idx: np.ndarray,
    thresholds: Tuple[float, ...] = POWER_THRESHOLDS,
    fst: np.ndarray = None,
) -> Dict[str, float]:
    """
    Compute all power/FPR/lambda metrics for one method in a single pass.

    Extracts causal and null p-values once, computes valid masks once,
    lambda_GC (full and null) once, then power/FPR/GC-power at all thresholds.
    Returns a dict with keys matching the row columns used in run_single_power_job.
    """
    causal_p = raw_pvals[causal_idx]
    null_p = raw_pvals[null_idx]

    valid_mask_full = np.isfinite(raw_pvals) & (raw_pvals > 0) & (raw_pvals <= 1)
    causal_valid_mask = np.isfinite(causal_p) & (causal_p > 0) & (causal_p <= 1)
    null_valid_mask = np.isfinite(null_p) & (null_p > 0) & (null_p <= 1)

    valid_pvals = raw_pvals[valid_mask_full]
    causal_valid = causal_p[causal_valid_mask]
    null_valid = null_p[null_valid_mask]

    n_valid = int(np.sum(valid_mask_full))
    n_causal_valid = int(np.sum(causal_valid_mask))

    lambda_gc = compute_lambda_gc(valid_pvals) if valid_pvals.size else np.nan
    lambda_gc_null = compute_lambda_gc(null_valid) if null_valid.size else np.nan

    # Power at all thresholds (one pass)
    if causal_valid.size == 0:
        powers = [np.nan] * len(thresholds)
    else:
        powers = [float(np.mean(causal_valid < t)) for t in thresholds]

    # FPR at all thresholds (one pass)
    if null_valid.size == 0:
        fprs = [np.nan] * len(thresholds)
    else:
        fprs = [float(np.mean(null_valid < t)) for t in thresholds]

    # GC-corrected power at all thresholds (reuse lambda_gc_null)
    if (
        null_valid.size == 0
        or not np.isfinite(lambda_gc_null)
        or lambda_gc_null <= 0
        or causal_valid_mask.sum() == 0
    ):
        power_gcs = [np.nan] * len(thresholds)
    else:
        causal_chi2 = chi2_dist.isf(causal_p[causal_valid_mask], df=1)
        corrected_chi2 = causal_chi2 / lambda_gc_null
        corrected_p = chi2_dist.sf(corrected_chi2, df=1)
        power_gcs = [float(np.mean(corrected_p < t)) for t in thresholds]

    # Map to the exact keys used in the CSV (5e8, 5e7, 1e5, 1e3)
    thresh_keys = ["5e8", "5e7", "1e5", "1e3"]
    out = {
        "n_valid_pvals": n_valid,
        "n_causal_valid": n_causal_valid,
        "n_null_idx": len(null_idx),
        "lambda_gc": lambda_gc,
        "lambda_gc_null": lambda_gc_null,
        "median_pval_causal": float(np.nanmedian(causal_p)),
        "median_pval_null": float(np.nanmedian(null_p)),
    }
    for i, k in enumerate(thresh_keys):
        if i < len(thresholds):
            out[f"power_{k}"] = powers[i]
            out[f"fpr_{k}"] = fprs[i]
            out[f"power_gc_{k}"] = power_gcs[i]

    # FST-stratified null calibration (Path C only — fst=None elsewhere).
    # Lambda_GC and FPR are computed separately on (a) low-FST null SNPs
    # ("clean") and (b) high-FST null SNPs ("strat") so we can attribute
    # any inflation to ancestry-driven stratification rather than
    # numerical artefact.
    if fst is not None and null_valid.size > 0:
        null_fst = fst[null_idx][null_valid_mask]
        # Per-SNP MAF-difference threshold of 0.10 ⇒ FST ≈ 0.04 for typical
        # MAFs.  Calibrated empirically to split null SNPs roughly 70/30.
        clean_mask = null_fst < 0.01
        strat_mask = null_fst > 0.10
        for tag, m in [("clean", clean_mask), ("strat", strat_mask)]:
            sub = null_valid[m]
            out[f"n_null_{tag}"] = int(m.sum())
            out[f"lambda_gc_null_{tag}"] = (
                compute_lambda_gc(sub) if sub.size else np.nan
            )
            for i, k in enumerate(thresh_keys):
                if i < len(thresholds):
                    out[f"fpr_{k}_{tag}"] = (
                        float(np.mean(sub < thresholds[i])) if sub.size else np.nan
                    )

    return out


# ---------------------------------------------------------------------------
# Scenario definitions (same eight as the main figure)
# ---------------------------------------------------------------------------

def _scenario_configs_power(arrays: Dict[str, np.ndarray]) -> List[Dict]:
    """
    Six scenarios matching the null simulation in ``_scenario_configs()``,
    with batch effects enabled by default.

    Cohort structure
    ----------------
    Within each scenario, ``_build_power_setting_dataset`` calls
    ``sample_subpopulations`` from ``simulation_2pop_qqplot``, which
    assembles three multi-ancestry cohorts (study 1–3), each containing
    subjects from both Pop A and Pop B.  See that function's docstring for
    the full cohort-design rationale.

    2-admixed scenarios
    -------------------
    The ``2admixed_homo`` and ``2admixed_hetero`` scenarios pair
    CEU5-YRI5 (50 % CEU + 50 % YRI) with CEU5-CHS5 (50 % CEU + 50 % CHS).
    Both populations share half their ancestry and are simulated from chr22
    only, which lacks sufficient ancestry-informative markers to separate
    their non-CEU components in PCA.  Consequently:
      • GMM population assignment approaches 50 % accuracy (chance level).
      • All six methods produce similar power and FPR — there is no
        detectable structure for any method to exploit.
    This is an expected data limitation, not a software defect.  The
    2-admixed scenarios serve as stress tests, not as scenarios where method
    differences are anticipated.
    """
    return [
        {
            "key": "2pure_homo",
            "label": "2 pure populations, homogeneous covariates",
            "X_std_A": arrays["Xstd_afr"],
            "X_std_B": arrays["Xstd_eur"],
            "heterogeneity": False, "n_components": 2,
            "batch_effect": True, "batch_frac": 0.05, "batch_size": 0.3,
        },
        {
            "key": "2pure_hetero",
            "label": "2 pure populations, heterogeneous covariates",
            "X_std_A": arrays["Xstd_afr"],
            "X_std_B": arrays["Xstd_eur"],
            "heterogeneity": True, "n_components": 2,
            "batch_effect": True, "batch_frac": 0.05, "batch_size": 0.3,
        },
        {
            "key": "1pure_1admixed_homo",
            "label": "1 pure + 1 admixed, homogeneous covariates",
            "X_std_A": arrays["Xstd_afr"],
            "X_std_B": arrays["Xstd_ceu5yri5"],
            "heterogeneity": False, "n_components": 2,
            "batch_effect": True, "batch_frac": 0.05, "batch_size": 0.3,
        },
        {
            "key": "1pure_1admixed_hetero",
            "label": "1 pure + 1 admixed, heterogeneous covariates",
            "X_std_A": arrays["Xstd_afr"],
            "X_std_B": arrays["Xstd_ceu5yri5"],
            "heterogeneity": True, "n_components": 2,
            "batch_effect": True, "batch_frac": 0.05, "batch_size": 0.3,
        },
        # ── 2-admixed scenarios ────────────────────────────────────────────
        # CEU5-YRI5 and CEU5-CHS5 share 50 % CEU ancestry and are simulated
        # from chr22 only.  PCA cannot separate their non-CEU components;
        # GMM accuracy ≈ 50 %; all six methods behave similarly.
        # See module docstring for full explanation.
        # ──────────────────────────────────────────────────────────────────
        {
            "key": "2admixed_homo",
            "label": "2 admixed populations, homogeneous covariates",
            "X_std_A": arrays["Xstd_ceu5yri5"],
            "X_std_B": arrays["Xstd_ceu5chs5"],
            "heterogeneity": False, "n_components": 2,
            "batch_effect": True, "batch_frac": 0.05, "batch_size": 0.3,
        },
        {
            "key": "2admixed_hetero",
            "label": "2 admixed populations, heterogeneous covariates",
            "X_std_A": arrays["Xstd_ceu5yri5"],
            "X_std_B": arrays["Xstd_ceu5chs5"],
            "heterogeneity": True, "n_components": 2,
            "batch_effect": True, "batch_frac": 0.05, "batch_size": 0.3,
        },
    ]


# ---------------------------------------------------------------------------
# Main driver
# ---------------------------------------------------------------------------

def _config_key(
    h2: float,
    num_causal: int,
    genetic_arch: str = "shared",
    genetic_corr: float = 0.8,
    pop_specific_frac: float = 0.5,
    causal_frac: Optional[float] = None,
) -> str:
    if causal_frac is not None:
        base = f"h2_{h2:.1f}_pct{causal_frac:.4f}"
    else:
        base = f"h2_{h2:.1f}_K{num_causal}"
    if genetic_arch == "shared":
        return base
    if genetic_arch == "pop_specific":
        return f"{base}_popspec{pop_specific_frac:.1f}"
    if genetic_arch == "heterogeneous":
        return f"{base}_hetero_rg{genetic_corr:.1f}"
    if genetic_arch == "mixed":
        return f"{base}_mixed_ps{pop_specific_frac:.1f}_rg{genetic_corr:.1f}"
    return base


def _resolve_num_causal(pc: Dict, M: int) -> int:
    """
    Resolve number of causal variants from config.
    If config has causal_frac (e.g. 0.001 for 0.1%), K = causal_frac * M, capped by LD block size.
    Otherwise use num_causal.
    """
    if pc.get("causal_frac") is not None:
        pct = float(pc["causal_frac"])
        K = round(pct * M)
        K = max(1, min(K, LD_BLOCK_SIZE))
        return K
    return int(pc["num_causal"])


def run_single_power_job(
    arrays: Dict[str, np.ndarray],
    replicate_id: int,
    seed: int,
    age_tolerance: int,
    out_dir: str,
    power_configs: Optional[List[Dict]] = None,
):
    """
    Run one power-analysis replicate and save per-run output files.
    """
    if power_configs is None:
        power_configs = POWER_CONFIGS

    out_path = Path(out_dir)
    out_path.mkdir(parents=True, exist_ok=True)

    scenarios = _scenario_configs_power(arrays)
    metrics_rows: List[Dict] = []
    pvals_archive: Dict[str, np.ndarray] = {}

    n_configs = len(power_configs)

    for s_idx, cfg in enumerate(scenarios):
        M_scenario = cfg["X_std_A"].shape[0]
        for pc_idx, pc in enumerate(power_configs):
            h2 = pc["h2"]
            num_causal = _resolve_num_causal(pc, M_scenario)
            genetic_arch = pc.get("genetic_arch", "shared")
            genetic_corr = pc.get("genetic_corr", 0.8)
            pop_specific_frac = pc.get("pop_specific_frac", 0.5)
            causal_frac = pc.get("causal_frac")

            beta_causal = h2_to_beta_causal(h2, num_causal)
            ck = _config_key(
                h2, num_causal, genetic_arch,
                genetic_corr, pop_specific_frac,
                causal_frac=causal_frac,
            )

            scenario_seed = seed * 10000 + s_idx * n_configs + pc_idx

            print(f"  Scenario: {cfg['key']}, config: {ck}")

            use_path_c = os.environ.get("USE_PATH_C", "0") == "1"
            df_s, Xstd_s, causal_idx, beta_snp, pop_accuracy, null_idx = _build_power_setting_dataset(
                X_std_A=cfg["X_std_A"], X_std_B=cfg["X_std_B"],
                heterogeneity=cfg["heterogeneity"],
                n_components=cfg["n_components"],
                seed=scenario_seed,
                beta_causal=beta_causal, num_causal=num_causal,
                batch_effect=cfg.get("batch_effect", True),
                batch_frac=cfg.get("batch_frac", 0.05),
                batch_size=cfg.get("batch_size", 0.3),
                genetic_arch=genetic_arch,
                genetic_corr=genetic_corr,
                pop_specific_frac=pop_specific_frac,
                scenario_key=cfg["key"],
                use_path_c=use_path_c,
            )

            # ---- GRM marker reservation (Arena A comparators only) ----------
            # GMMAT and GENESIS must not be handed a GRM built from the variants
            # they are about to test.  Measured on this pipeline's own data, a
            # variant whose own 2,000-SNP LD block is in the GRM is tested at
            # lambda_GC 0.788; remove that block from the GRM and the same
            # variants calibrate at 1.006.  (Removing the CAUSAL block is a
            # separate, additive effect worth ~2x in causal chi-square.)  Both
            # harms come from the same source: 92% of a SNP's total r^2 lies
            # inside its own block, so a whole-genome GRM here is ~14% "the
            # variant being tested" -- an assumption real GRMs satisfy and this
            # simulation does not, because one of 7 pseudo-chromosomes is a
            # seventh of the genome.
            #
            # Rather than refit the null model 7 times per unit (leave-one-block
            # -out, ~2.6x the wall clock), reserve whole blocks as GRM markers
            # and never evaluate on them.  No tested-and-evaluated variant then
            # contributes to the random effect at all, which is a stronger
            # guarantee than LOCO and mirrors real practice, where the GRM comes
            # from a pruned marker set.  REGENIE keeps its own LOCO on the same
            # .bed and is unaffected.
            n_blocks = int(np.ceil(M_scenario / LD_BLOCK_SIZE))
            causal_block = int(causal_idx[0]) // LD_BLOCK_SIZE
            avail = [b for b in range(n_blocks) if b != causal_block]
            grm_blocks = sorted(avail)[-N_GRM_BLOCKS:] if len(avail) > N_GRM_BLOCKS else avail
            block_of = np.arange(M_scenario) // LD_BLOCK_SIZE
            grm_snp_idx = np.where(np.isin(block_of, grm_blocks))[0]
            # Evaluate only where neither the causal signal nor the GRM reaches.
            eval_null_idx = np.setdiff1d(null_idx, grm_snp_idx, assume_unique=False)
            print(f"    blocks: {n_blocks} total | causal {causal_block} | "
                  f"GRM {grm_blocks} ({grm_snp_idx.size} SNPs) | "
                  f"evaluated null {eval_null_idx.size}")

            raw_pvals = run_all_methods_raw(
                df_s, Xstd_s, age_tolerance=age_tolerance,
                comparators=_COMPARATORS,
                dosage_seed=scenario_seed,
                tag=f"r{replicate_id:03d}s{s_idx}c{pc_idx}",
                grm_snp_idx=grm_snp_idx,
            )

            # Per-SNP FST for FST-stratified null reporting (Path C only).
            scenario_fst = None
            if use_path_c:
                from path_c import fst_for_scenario
                scenario_fst = fst_for_scenario(cfg["key"])

            # Matched-sample sizes, recorded per arm so that no calibration
            # number from a matched arm can ever be read without the N it was
            # computed on.  The finer matching keys retain materially fewer
            # pairs than matched_population, so "better calibrated" and
            # "smaller and less powered" are confounded by construction and the
            # reader has to be able to see both.  ~0.5 s/unit to recompute.
            n_pairs_by_method: Dict[str, float] = {}
            for _mname, _mfn in (
                ("matched_population", match_cases_controls_with_center),
                ("matched_study", match_cases_controls_with_study),
                ("matched_study_center", match_cases_controls_with_study_center),
            ):
                try:
                    _cs, _ct = _mfn(df_s, age_tolerance=age_tolerance)
                    n_pairs_by_method[_mname] = float(min(len(_cs), len(_ct)))
                except Exception as exc:                          # noqa: BLE001
                    print(f"    [{_mname}] pair count failed: "
                          f"{type(exc).__name__}: {exc}")
                    n_pairs_by_method[_mname] = np.nan
            print("    matched pairs: " + " | ".join(
                f"{k}={v:.0f}" for k, v in n_pairs_by_method.items()
                if np.isfinite(v)
            ) + f"  (N={len(df_s)})")

            for method, pvals in raw_pvals.items():
                pvals_archive[f"{cfg['key']}__{ck}__{method}"] = pvals

                m = _metrics_for_method(pvals, causal_idx, eval_null_idx, fst=scenario_fst)
                row = {
                    "replicate": replicate_id,
                    "seed": seed,
                    "scenario_key": cfg["key"],
                    "scenario_label": cfg["label"],
                    "h2": h2,
                    "beta_causal": beta_causal,
                    "odds_ratio": np.exp(beta_causal),
                    "num_causal": num_causal,
                    "genetic_arch": genetic_arch,
                    "genetic_corr": genetic_corr,
                    "pop_specific_frac": pop_specific_frac,
                    "config_key": ck,
                    "method": method,
                    "pop_accuracy": pop_accuracy,
                    # NaN for every non-matched arm, which analyse on all N.
                    "n_matched_pairs": n_pairs_by_method.get(method, np.nan),
                    "n_subjects_total": len(df_s),
                    "interaction_mode": INTERACTION_MODE,
                    **m,
                }
                if causal_frac is not None:
                    row["causal_frac"] = causal_frac
                metrics_rows.append(row)

            pvals_archive[f"{cfg['key']}__{ck}__causal_idx"] = causal_idx
            # Persist the reserved GRM markers so the aggregator excludes the
            # same variants it would otherwise silently score.
            pvals_archive[f"{cfg['key']}__{ck}__grm_idx"] = grm_snp_idx

    metrics_df = pd.DataFrame(metrics_rows)
    metrics_path = out_path / f"power_run_{replicate_id:03d}_metrics.csv"
    pvals_path = out_path / f"power_run_{replicate_id:03d}_pvals.npz"
    metrics_df.to_csv(metrics_path, index=False)
    np.savez_compressed(pvals_path, **pvals_archive)
    print(f"Saved: {metrics_path}, {pvals_path}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run one power-analysis replicate for multi-ancestry GWAS methods."
    )
    parser.add_argument("--replicate-id", type=int, default=1)
    parser.add_argument("--seed", type=int, default=RNG_SEED)
    parser.add_argument("--age-tolerance", type=int, default=1)
    parser.add_argument("--out-dir", type=str, default="outputs/power_runs")
    parser.add_argument(
        "--h2-values", type=str, default="0.3,0.5,0.7",
        help="Comma-separated heritability values (default: 0.3,0.5,0.7)",
    )
    parser.add_argument(
        "--num-causal-values", type=str, default="10",
        help="Comma-separated numbers of causal variants (default: 10). Ignored if --causal-frac is set.",
    )
    parser.add_argument(
        "--causal-frac",
        type=str,
        default=None,
        help=(
            "Comma-separated fraction of SNPs that are causal (e.g. 0.001,0.01 for 0.1%% and 1%%). "
            "K = causal_frac * M per scenario, capped by LD block size. Overrides --num-causal-values."
        ),
    )
    parser.add_argument(
        "--genetic-arch", type=str, default="shared",
        help=(
            "Comma-separated genetic architecture modes "
            "(shared, pop_specific, heterogeneous, mixed). "
            "Default: shared"
        ),
    )
    parser.add_argument(
        "--pop-specific-frac", type=float, default=0.5,
        help=(
            "Fraction of causal variants that are population-specific "
            "(for pop_specific/mixed). Default: 0.5"
        ),
    )
    parser.add_argument(
        "--effect-corr", type=float, default=0.8,
        help=(
            "Cross-population genetic correlation rG "
            "(for heterogeneous/mixed). Default: 0.8"
        ),
    )
    return parser.parse_args()


def _build_power_configs(
    h2_str: str,
    num_causal_str: str,
    genetic_arch_str: str = "shared",
    pop_specific_frac: float = 0.5,
    effect_corr: float = 0.8,
    causal_frac_str: Optional[str] = None,
) -> List[Dict]:
    h2_values = [float(x.strip()) for x in h2_str.split(",")]
    genetic_archs = [x.strip() for x in genetic_arch_str.split(",")]
    configs = []
    if causal_frac_str is not None:
        causal_frac_values = [float(x.strip()) for x in causal_frac_str.split(",")]
        for h2 in h2_values:
            for pct in causal_frac_values:
                for ga in genetic_archs:
                    configs.append({
                        "h2": h2,
                        "causal_frac": pct,
                        "genetic_arch": ga,
                        "genetic_corr": effect_corr,
                        "pop_specific_frac": pop_specific_frac,
                        "label": f"h2={h2:.1f}, pct={pct:.4f}, {ga}",
                    })
    else:
        num_causal_values = [int(x.strip()) for x in num_causal_str.split(",")]
        for h2 in h2_values:
            for nc in num_causal_values:
                for ga in genetic_archs:
                    configs.append({
                        "h2": h2,
                        "num_causal": nc,
                        "genetic_arch": ga,
                        "genetic_corr": effect_corr,
                        "pop_specific_frac": pop_specific_frac,
                        "label": f"h2={h2:.1f}, K={nc}, {ga}",
                    })
    return configs


if __name__ == "__main__":
    args = _parse_args()
    power_configs = _build_power_configs(
        args.h2_values,
        args.num_causal_values,
        args.genetic_arch,
        args.pop_specific_frac,
        args.effect_corr,
        causal_frac_str=args.causal_frac,
    )
    arrays = _load_default_arrays()
    run_single_power_job(
        arrays=arrays,
        replicate_id=args.replicate_id,
        seed=args.seed,
        age_tolerance=args.age_tolerance,
        out_dir=args.out_dir,
        power_configs=power_configs,
    )
