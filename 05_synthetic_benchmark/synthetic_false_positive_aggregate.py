"""
synthetic_false_positive_aggregate.py — Pool null p-values across all simulation runs to compute
empirical FPR at alpha = 5e-5, 5e-7, and 5e-8.

Uses a streaming approach: for each .npz file, counts false positives at each
alpha threshold without storing all p-values in memory.  Memory usage is O(1)
per (scenario, config, method) accumulator, not O(N_total).

Usage
-----
Null simulation runs (all SNPs are null):

    python synthetic_false_positive_aggregate.py --run-dir outputs/null_runs \\
        --run-type null --out outputs/fpr_null.csv

Power simulation runs (null SNPs are those OUTSIDE the causal LD block):

    python synthetic_false_positive_aggregate.py --run-dir outputs/power_runs \\
        --run-type power --out outputs/fpr_power.csv

Design rationale and sample-size requirements
---------------------------------------------
Direct empirical FPR at alpha requires N_null × alpha >= ~50 expected false
positives per scenario x method to keep Poisson CV <= 14%.

With M_null = 14 000 null SNPs per null-simulation run (or ~12 000 for power
runs), the required number of replicates is:

    alpha   null-sim runs needed   power-sim runs needed
    -----   --------------------   ---------------------
    5e-5    71   (100 already done)  83   (100 already done)
    5e-7    7143 (need ~10 000)      8333 (need ~10 000)
    5e-8    71429 -- infeasible; use 95% CI upper bound instead

At 5e-8 with 10 000 power runs (120 M null tests per scenario), if zero FPs
are observed the 95% upper bound is -ln(0.05) / 120e6 ~ 2.5e-8 < 5e-8, which
establishes calibration without requiring a direct observation.

Output columns
--------------
  scenario_key, method, alpha, n_null_total,
  n_fp,               # observed false positives
  fpr_empirical,      # n_fp / n_null_total
  ci_lo_95, ci_hi_95, # 95% Wilson CI
  expected_fp         # n_null_total * alpha (expected under perfect calibration)
"""

import argparse
import sys
import warnings
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

# Alpha levels to evaluate
ALPHAS = [5e-5, 5e-7, 5e-8]

# LD block size (must match synthetic_power_simulation.py)
LD_BLOCK_SIZE = 2000


# ---------------------------------------------------------------------------
# Wilson 95% confidence interval for a binomial proportion
# ---------------------------------------------------------------------------

def wilson_ci(n_success: int, n_total: int, z: float = 1.96) -> Tuple[float, float]:
    """Wilson score interval for a binomial proportion."""
    if n_total == 0:
        return (0.0, 1.0)
    p_hat = n_success / n_total
    denom = 1.0 + z**2 / n_total
    centre = (p_hat + z**2 / (2 * n_total)) / denom
    half_width = (z / denom) * np.sqrt(p_hat * (1 - p_hat) / n_total + z**2 / (4 * n_total**2))
    return (max(0.0, centre - half_width), min(1.0, centre + half_width))


def one_sided_upper_bound(n_total: int, confidence: float = 0.95) -> float:
    """
    Exact Clopper-Pearson 95% one-sided upper bound when n_success = 0.
    P(X >= 1 | n, p) <= 1 - confidence  =>  p <= 1 - confidence^(1/n).
    """
    if n_total == 0:
        return 1.0
    return 1.0 - confidence ** (1.0 / n_total)


# ---------------------------------------------------------------------------
# Streaming counters: count FPs at each alpha without storing all p-values
# ---------------------------------------------------------------------------

class FPCounter:
    """Accumulates FP counts and total null tests across files."""
    __slots__ = ("n_total", "fp_counts")

    def __init__(self, alphas: List[float]):
        self.n_total = 0
        self.fp_counts = {a: 0 for a in alphas}

    def update(self, pvals: np.ndarray) -> None:
        """Count valid null p-values below each alpha threshold."""
        valid = pvals[np.isfinite(pvals) & (pvals > 0) & (pvals <= 1)]
        self.n_total += len(valid)
        for alpha, count in self.fp_counts.items():
            self.fp_counts[alpha] = count + int(np.sum(valid < alpha))


def _get_null_idx_from_causal(causal_idx: np.ndarray, M: int,
                               block_size: int = LD_BLOCK_SIZE) -> np.ndarray:
    """
    Given causal variant indices and total SNP count M, return indices of SNPs
    that are NOT in the same LD block as any causal variant.
    """
    causal_blocks = set(int(i) // block_size for i in causal_idx)
    block_ids = np.arange(M) // block_size
    null_mask = np.ones(M, dtype=bool)
    for b in causal_blocks:
        null_mask[block_ids == b] = False
    return np.where(null_mask)[0]


# ---------------------------------------------------------------------------
# Streaming aggregation: null simulation
# ---------------------------------------------------------------------------

def _stream_null_runs(
    run_dir: Path, alphas: List[float],
) -> Dict[Tuple[str, str], FPCounter]:
    """Stream through run_NNN_pvals.npz files, counting FPs per (scenario, method)."""
    pval_files = sorted(run_dir.glob("run_*_pvals.npz"))
    if not pval_files:
        raise FileNotFoundError(f"No run_*_pvals.npz files found in {run_dir}")

    counters: Dict[Tuple[str, str], FPCounter] = defaultdict(lambda: FPCounter(alphas))
    n_files = len(pval_files)

    for i, fpath in enumerate(pval_files):
        if (i + 1) % 100 == 0 or i == 0:
            print(f"  Processing file {i+1}/{n_files}: {fpath.name}", flush=True)
        archive = np.load(fpath, allow_pickle=False)
        for key in archive.files:
            parts = key.split("__")
            if len(parts) != 2:
                continue
            scenario_key, method = parts
            counters[(scenario_key, method)].update(archive[key])
        archive.close()

    return dict(counters)


# ---------------------------------------------------------------------------
# Streaming aggregation: power simulation
# ---------------------------------------------------------------------------

def _stream_power_runs(
    run_dir: Path, alphas: List[float],
) -> Dict[Tuple[str, str, str], FPCounter]:
    """Stream through power_run_NNN_pvals.npz files, counting FPs for null SNPs only."""
    pval_files = sorted(run_dir.glob("power_run_*_pvals.npz"))
    if not pval_files:
        raise FileNotFoundError(f"No power_run_*_pvals.npz files found in {run_dir}")

    counters: Dict[Tuple[str, str, str], FPCounter] = defaultdict(lambda: FPCounter(alphas))
    n_files = len(pval_files)

    for i, fpath in enumerate(pval_files):
        if (i + 1) % 100 == 0 or i == 0:
            print(f"  Processing file {i+1}/{n_files}: {fpath.name}", flush=True)

        archive = np.load(fpath, allow_pickle=False)

        # Collect causal_idx arrays keyed by (scenario, config)
        causal_map: Dict[Tuple[str, str], np.ndarray] = {}
        for key in archive.files:
            if key.endswith("__causal_idx"):
                stem = key[: -len("__causal_idx")]
                parts = stem.split("__")
                if len(parts) == 2:
                    causal_map[tuple(parts)] = archive[key]

        # Process pval arrays: extract null SNPs, count FPs
        for key in archive.files:
            if key.endswith("__causal_idx"):
                continue
            parts = key.split("__")
            if len(parts) != 3:
                continue
            scenario_key, config_key, method = parts
            pvals = archive[key]
            M = len(pvals)

            causal_key = (scenario_key, config_key)
            if causal_key not in causal_map:
                warnings.warn(
                    f"No causal_idx found for {causal_key} in {fpath.name}; skipping."
                )
                continue

            null_idx = _get_null_idx_from_causal(causal_map[causal_key], M)
            counters[(scenario_key, config_key, method)].update(pvals[null_idx])

        archive.close()

    return dict(counters)


# ---------------------------------------------------------------------------
# Convert counters to result rows
# ---------------------------------------------------------------------------

def counter_to_rows(
    counter: FPCounter,
    scenario_key: str = "",
    method: str = "",
    extra_cols: Optional[Dict] = None,
) -> List[Dict]:
    """Convert an FPCounter into one dict per alpha level."""
    rows = []
    n_total = counter.n_total
    for alpha, n_fp in counter.fp_counts.items():
        fpr = n_fp / n_total if n_total > 0 else np.nan
        ci_lo, ci_hi = wilson_ci(n_fp, n_total)

        if n_fp == 0 and n_total > 0:
            ub_1sided = one_sided_upper_bound(n_total)
        else:
            ub_1sided = np.nan

        row = {
            "scenario_key": scenario_key,
            "method": method,
            "alpha": alpha,
            "n_null_total": n_total,
            "n_fp": n_fp,
            "fpr_empirical": fpr,
            "ci_lo_95": ci_lo,
            "ci_hi_95": ci_hi,
            "ci_hi_95_onesided": ub_1sided,
            "expected_fp": n_total * alpha,
            "fpr_ratio_vs_nominal": fpr / alpha if fpr > 0 and np.isfinite(fpr) else np.nan,
        }
        if extra_cols:
            row.update(extra_cols)
        rows.append(row)
    return rows


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Pool null p-values across simulation runs and compute empirical FPR."
    )
    parser.add_argument(
        "--run-dir", required=True,
        help="Directory containing run_*_pvals.npz (null sim) or power_run_*_pvals.npz (power sim).",
    )
    parser.add_argument(
        "--run-type", choices=["null", "power"], default="null",
        help="'null' to pool all SNPs (null sim); 'power' to pool only null SNPs outside LD block.",
    )
    parser.add_argument(
        "--out", default=None,
        help="Output CSV path (default: <run-dir>/fpr_pooled.csv).",
    )
    parser.add_argument(
        "--alphas", default="5e-5,5e-7,5e-8",
        help="Comma-separated alpha levels (default: 5e-5,5e-7,5e-8).",
    )
    return parser.parse_args()


def main() -> None:
    args = _parse_args()
    run_dir = Path(args.run_dir)
    alphas = [float(x.strip()) for x in args.alphas.split(",")]
    out_path = Path(args.out) if args.out else run_dir / "fpr_pooled.csv"

    print(f"Run directory : {run_dir}")
    print(f"Run type      : {args.run_type}")
    print(f"Alpha levels  : {alphas}")

    M_per_run = 14_000 if args.run_type == "null" else 12_000
    for alpha in alphas:
        n_runs_needed = 50 / (M_per_run * alpha)
        print(
            f"  alpha={alpha:.0e}: need ~{n_runs_needed:.0f} runs for 50 expected FPs "
            f"(CV ~ 14%)"
        )

    rows: List[Dict] = []

    if args.run_type == "null":
        print("\nStreaming null-simulation p-values (all M SNPs are null)...")
        counters = _stream_null_runs(run_dir, alphas)
        for (scenario_key, method) in sorted(counters):
            c = counters[(scenario_key, method)]
            print(f"  {scenario_key} x {method}: {c.n_total:,} null tests")
            rows.extend(counter_to_rows(c, scenario_key=scenario_key, method=method))

    else:  # power
        print("\nStreaming power-simulation null p-values (SNPs outside causal LD block)...")
        counters = _stream_power_runs(run_dir, alphas)
        for key in sorted(counters):
            scenario_key, config_key, method = key
            c = counters[key]
            print(f"  {scenario_key} x {config_key} x {method}: {c.n_total:,} null tests")
            rows.extend(
                counter_to_rows(
                    c, scenario_key=scenario_key, method=method,
                    extra_cols={"config_key": config_key},
                )
            )

    df = pd.DataFrame(rows)
    df.to_csv(out_path, index=False, float_format="%.6g")
    print(f"\nSaved: {out_path}")
    print(f"Rows: {len(df)}")

    # Quick summary
    for alpha in alphas:
        sub = df[df["alpha"] == alpha]
        n_runs_approx = sub["n_null_total"].median() / M_per_run if not sub.empty else 0
        well_cal = sub[sub["fpr_ratio_vs_nominal"].notna() & (sub["fpr_ratio_vs_nominal"] <= 2.0)]
        no_ratio = sub[sub["fpr_ratio_vs_nominal"].isna()]  # n_fp == 0
        print(
            f"\nalpha={alpha:.0e}  (~{n_runs_approx:.0f} runs pooled):"
            f"  {len(well_cal)}/{len(sub)} combos with FPR/alpha <= 2"
            f"  ({len(no_ratio)} with 0 FPs)"
        )
        if not sub.empty and sub["n_fp"].sum() == 0:
            sample_n = int(sub["n_null_total"].iloc[0])
            ub = one_sided_upper_bound(sample_n)
            print(
                f"  No FPs observed anywhere; 95% one-sided UB ~ {ub:.2e}"
                f" ({'< alpha' if ub < alpha else '> alpha'})"
            )


if __name__ == "__main__":
    main()
