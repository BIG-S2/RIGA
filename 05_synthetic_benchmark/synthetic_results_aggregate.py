"""
synthetic_results_aggregate.py — pool the Arena A benchmark and report it honestly.

Reads outputs/power_runs_bench/power_run_*_{metrics.csv,pvals.npz} produced by
revision/scripts/arenaA_bench.sbatch (COMPARATORS=regenie,gmmat,genesis) and
emits the false-positive-rate, power, inflation and tool-status tables for the
ten methods:

  combined_study, combined_pcs, combined_population, combined_study_pcs,
  combined_population_pcs, matched_population      (the six native variants;
                                                    matched_population is RIGA)
  combined_pcs_dosage                              (same test, hard-call view —
                                                    the discretisation anchor)
  regenie, gmmat, genesis                          (the external comparators)


Why this script exists rather than src/simulation/synthetic_false_positive_aggregate.py
------------------------------------------------------------------
1. synthetic_false_positive_aggregate.py's one_sided_upper_bound() is wrong.  It returns
   1 - confidence^(1/n) = 1 - 0.95^(1/n) ~ 0.0513/n.  The 95% one-sided
   Clopper-Pearson upper bound for zero events is the rule of three,
   1 - 0.05^(1/n) ~ 3/n — about 58x larger.  run_power_array.sbatch's own
   design comment quotes the correct value ("120M tests, 0 FPs => UB ~ 2.5e-8",
   i.e. 3/120e6), so the code contradicts the design it implements.  Every
   bound this script prints comes from scipy.stats.beta and is unit-tested
   against the rule of three at the bottom of the file.

2. The revision requires the bound at alpha = 5e-8 *whether or not* zero events
   were seen.  synthetic_false_positive_aggregate.py only computes a bound when n_fp == 0.  The
   general exact upper limit is  Beta^{-1}(1 - a; x + 1, n - x).

3. Binomial intervals assume independent tests.  Null SNPs inside one replicate
   share 3,000 subjects and sit in LD blocks, so the effective number of
   independent units is nearer the replicate count than the SNP count.  A
   nominal binomial CI is therefore too narrow.  This script also reports a
   cluster (replicate-level) bootstrap interval, and the paper should quote the
   wider of the two.

Usage
-----
  python synthetic_results_aggregate.py --run-dir outputs/power_runs_bench \\
      --out-dir revision/tables --fig revision/figures
"""
from __future__ import annotations

import argparse
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd
from scipy.stats import beta as beta_dist

LD_BLOCK_SIZE = 2000
ALPHAS = [1e-3, 1e-5, 5e-7, 5e-8]

# Below this many expected false positives a point estimate is not reportable:
# the Poisson coefficient of variation exceeds 20%, so the estimate carries
# less information than the interval does.  At and above it the point estimate
# is quoted with its exact interval.
MIN_EXPECTED_FP_FOR_POINT = 25.0

METHOD_LABEL = {
    "combined_study": "Pooled + study",
    "combined_pcs": "Pooled + PCs",
    "combined_population": "Pooled + population",
    "combined_study_pcs": "Pooled + study + PCs",
    "combined_population_pcs": "Pooled + population + PCs",
    "matched_population": "RIGA (matched)",
    "matched_study": "RIGA (matched, study key)",
    "matched_study_center": "RIGA (matched, study+centre key)",
    "matched_population_pcs": "RIGA (matched + PCs)",
    "matched_study_pcs": "RIGA (matched, study key + PCs)",
    "combined_pcs_dosage": "Pooled + PCs (hard calls)",
    "regenie": "REGENIE",
    "gmmat": "GMMAT",
    "genesis": "GENESIS",
}
EXTERNAL = ("regenie", "gmmat", "genesis")


# ---------------------------------------------------------------------------
# Exact binomial (Clopper-Pearson) intervals
# ---------------------------------------------------------------------------

def cp_interval(x: int, n: int, conf: float = 0.95) -> Tuple[float, float]:
    """Exact two-sided Clopper-Pearson interval for x successes in n trials."""
    if n == 0:
        return (0.0, 1.0)
    a = 1.0 - conf
    lo = 0.0 if x == 0 else float(beta_dist.ppf(a / 2.0, x, n - x + 1))
    hi = 1.0 if x == n else float(beta_dist.ppf(1.0 - a / 2.0, x + 1, n - x))
    return (lo, hi)


def cp_upper_one_sided(x: int, n: int, conf: float = 0.95) -> float:
    """
    Exact one-sided Clopper-Pearson upper confidence limit.

    Solves  P(X <= x | n, p) = 1 - conf  for p, i.e. the largest p under which
    the observed count would not be surprisingly small.  Equals the rule of
    three, 1 - (1-conf)^(1/n), when x = 0.
    """
    if n == 0:
        return 1.0
    if x >= n:
        return 1.0
    return float(beta_dist.ppf(conf, x + 1, n - x))


def cp_lower_one_sided(x: int, n: int, conf: float = 0.95) -> float:
    """
    Exact one-sided Clopper-Pearson LOWER confidence limit.

    This is the statistic that decides inflation.  "Upper bound <= alpha" is the
    wrong test and is not even attainable here: a perfectly calibrated method
    pooled over 108M null tests produces ~5.4 false positives at alpha = 5e-8,
    whose 95% upper limit is ~2x alpha.  Demanding upper <= alpha would require
    observing FEWER false positives than calibration predicts.  A method is
    declared inflated only when the lower limit exceeds alpha.
    """
    if n == 0 or x == 0:
        return 0.0
    return float(beta_dist.ppf(1.0 - conf, x, n - x + 1))


def _self_test() -> None:
    """The rule of three, and a textbook value, as executable assertions."""
    n = 120_000_000
    # 1 - 0.05^(1/n) loses most of its significant digits to cancellation at this
    # n, so the reference is computed with expm1 and compared relatively.
    ref = -np.expm1(np.log(0.05) / n)
    got = cp_upper_one_sided(0, n)
    assert abs(got - ref) / ref < 1e-6, f"rule of three: {got} vs {ref}"
    assert abs(got * n - 2.9957) < 1e-3, "rule of three ~ 3/n"
    lo, hi = cp_interval(0, 100)
    assert lo == 0.0 and abs(hi - 0.03621669) < 1e-6
    lo, hi = cp_interval(5, 100)
    assert abs(lo - 0.01643) < 1e-4 and abs(hi - 0.11284) < 1e-4


# ---------------------------------------------------------------------------
# Streaming accumulation
# ---------------------------------------------------------------------------

class Counter:
    """
    False-positive counts, with per-unit detail retained.

    A "unit" is one (replicate, scenario, config) cell.  Units are kept
    individually for two reasons: the cluster bootstrap must resample
    REPLICATES (the independent sampling unit) rather than units, and the
    common-unit denominator needs to know which units each method actually
    produced.
    """
    __slots__ = ("n_total", "fp", "units")

    def __init__(self) -> None:
        self.n_total = 0
        self.fp = {a: 0 for a in ALPHAS}
        # unit_key -> (n_valid, {alpha: count})
        self.units: Dict[Tuple[int, str, str], Tuple[int, Dict[float, int]]] = {}

    def update(self, p: np.ndarray, unit_key: Tuple[int, str, str]) -> None:
        finite = np.isfinite(p)
        # p == 0 is an underflowed p-value, not an invalid one.  Dropping it
        # would delete the single most extreme false positive from both the
        # numerator and the denominator — precisely the observations the tail
        # thresholds exist to count.  Keep them and count them at every alpha.
        v = p[finite & (p >= 0) & (p <= 1)]
        n_zero = int(np.sum(v == 0))
        self.n_total += v.size
        counts = {}
        for a in ALPHAS:
            c = int(np.sum(v < a))     # p == 0 < a for every alpha considered
            self.fp[a] += c
            counts[a] = c
        counts["_n_zero"] = n_zero
        self.units[unit_key] = (v.size, counts)

    def merge(self, other: "Counter") -> None:
        self.n_total += other.n_total
        for a in ALPHAS:
            self.fp[a] += other.fp[a]
        self.units.update(other.units)

    def restricted(self, keep: set) -> "Counter":
        """A copy holding only the units in `keep`."""
        c = Counter()
        for k, (n, counts) in self.units.items():
            if k not in keep:
                continue
            c.n_total += n
            for a in ALPHAS:
                c.fp[a] += counts[a]
            c.units[k] = (n, counts)
        return c


def null_index(causal_idx: np.ndarray, M: int) -> np.ndarray:
    """SNPs not in any LD block containing a causal variant."""
    blocks = set(int(i) // LD_BLOCK_SIZE for i in causal_idx)
    bid = np.arange(M) // LD_BLOCK_SIZE
    mask = np.ones(M, dtype=bool)
    for b in blocks:
        mask[bid == b] = False
    return np.where(mask)[0]


def frozen_manifest(run_dir: Path) -> List[int]:
    """
    Replicate ids having BOTH a pvals archive and a metrics CSV, resolved once.

    The benchmark array writes continuously, and an earlier version globbed the
    pvals files at the start of streaming and the metrics files ~40 s later,
    after ~1.2 GB of I/O. A replicate that landed in between appeared in the
    lambda/power tables but not the false-positive table, so the two described
    different replicate sets while both were labelled the same N. The effect was
    small (lambda shifted <= 3.2e-4) but it is a race, not a transcription slip,
    and it would silently recur on the final 500-replicate aggregation.
    """
    pv = {int(f.name.split("_")[2]) for f in run_dir.glob("power_run_*_pvals.npz")}
    mt = {int(f.name.split("_")[2]) for f in run_dir.glob("power_run_*_metrics.csv")}
    both = sorted(pv & mt)
    if pv ^ mt:
        print(f"  note: {len(pv ^ mt)} replicate(s) had only one of the two files "
              f"and are excluded from BOTH tables: {sorted(pv ^ mt)[:8]}")
    return both


def stream(run_dir: Path, keep_ids: List[int]):
    """Accumulate counters keyed by (scenario, config, method), and tool status."""
    ids = set(keep_ids)
    files = sorted(f for f in run_dir.glob("power_run_*_pvals.npz")
                   if int(f.name.split("_")[2]) in ids)
    if not files:
        raise FileNotFoundError(f"no power_run_*_pvals.npz in {run_dir}")

    counters: Dict[Tuple[str, str, str], Counter] = defaultdict(Counter)
    # status[method] = [n_units_total, n_units_all_nan, n_units_partial_nan]
    status: Dict[str, List[int]] = defaultdict(lambda: [0, 0, 0])
    # GMMAT and GENESIS are two implementations of one estimator (logistic mixed
    # model + score test on the same GRM). Their agreement is a wiring check, not
    # corroboration, so it is measured rather than asserted: a reviewer who
    # correlates the two comparator columns must find the number already stated.
    agree: List[Dict] = []

    for i, f in enumerate(files):
        if i % 25 == 0:
            print(f"  [{i+1}/{len(files)}] {f.name}", flush=True)
        try:
            rep_id = int(f.name.split("_")[2])
        except (IndexError, ValueError):
            rep_id = i
        try:
            ar = np.load(f, allow_pickle=False)
        except Exception as exc:                                  # noqa: BLE001
            print(f"    SKIP unreadable {f.name}: {exc}")
            continue

        pair: Dict[Tuple[str, str, str], np.ndarray] = {}
        causal, grm = {}, {}
        for k in ar.files:
            for suffix, tgt in (("__causal_idx", causal), ("__grm_idx", grm)):
                if k.endswith(suffix):
                    parts = k[: -len(suffix)].split("__")
                    if len(parts) == 2:
                        tgt[tuple(parts)] = ar[k]

        for k in ar.files:
            if k.endswith("__causal_idx") or k.endswith("__grm_idx"):
                continue
            parts = k.split("__")
            if len(parts) != 3:
                continue
            scen, cfg, meth = parts
            p = ar[k]
            M = len(p)
            if (scen, cfg) not in causal:
                continue

            fin = np.isfinite(p)
            status[meth][0] += 1
            if not fin.any():
                status[meth][1] += 1
                # An all-NaN unit is a tool failure. It must NOT simply vanish
                # from that tool's denominator: silently dropping the units a
                # tool failed on flatters the tool that failed. It is recorded
                # in the status table and excluded from the common-unit
                # stratum below.
                continue
            if not fin.all():
                status[meth][2] += 1

            # Evaluate only where neither the causal signal nor the GRM
            # reaches. The reserved GRM markers are in the GRM handed to
            # GMMAT/GENESIS, so scoring them would reintroduce exactly the
            # proximal contamination the reservation exists to remove.
            nidx = null_index(causal[(scen, cfg)], M)
            gidx = grm.get((scen, cfg))
            if gidx is not None and gidx.size:
                nidx = np.setdiff1d(nidx, gidx, assume_unique=False)
            counters[(scen, cfg, meth)].update(p[nidx], (rep_id, scen, cfg))
            if meth in ("gmmat", "genesis"):
                pair[(scen, cfg, meth)] = p[nidx]

        for (scen, cfg, _m) in list(pair):
            a = pair.get((scen, cfg, "gmmat"))
            b = pair.get((scen, cfg, "genesis"))
            if a is None or b is None:
                continue
            ok = np.isfinite(a) & np.isfinite(b) & (a > 0) & (b > 0)
            if ok.sum() < 100:
                continue
            la, lb = -np.log10(a[ok]), -np.log10(b[ok])
            agree.append({
                "replicate": rep_id, "scenario_key": scen, "config_key": cfg,
                "n": int(ok.sum()),
                "pearson_r_neglog10p": float(np.corrcoef(la, lb)[0, 1]),
                "max_abs_diff_neglog10p": float(np.max(np.abs(la - lb))),
                "median_abs_diff_neglog10p": float(np.median(np.abs(la - lb))),
            })
            pair.pop((scen, cfg, "gmmat"), None)
            pair.pop((scen, cfg, "genesis"), None)
        pair.clear()
        ar.close()

    return counters, status, len(files), agree


# ---------------------------------------------------------------------------
# Reporting
# ---------------------------------------------------------------------------

def cluster_bootstrap(units: Dict, alpha: float, n_boot: int = 2000,
                      seed: int = 12345, conf: float = 0.95) -> Tuple[float, float]:
    """
    Replicate-level cluster bootstrap of the FPR.

    Null SNPs are not independent: within one replicate they share 3,000
    subjects and sit in LD blocks, and the SAME replicate contributes a unit to
    every (scenario, config) cell.  The independent sampling unit is therefore
    the REPLICATE, not the SNP and not the (replicate, scenario, config) unit.
    Resampling anything finer understates the variance — which is the opposite
    of what this interval exists to do.

    Returns (lo, hi) percentile limits, or (nan, nan) with fewer than 2
    replicates.
    """
    by_rep: Dict[int, List[float]] = defaultdict(lambda: [0.0, 0.0])
    for (rep_id, _scen, _cfg), (n, counts) in units.items():
        by_rep[rep_id][0] += n
        by_rep[rep_id][1] += counts[alpha]
    if len(by_rep) < 2:
        return (np.nan, np.nan)
    sizes = np.array([v[0] for v in by_rep.values()], dtype=np.float64)
    hits = np.array([v[1] for v in by_rep.values()], dtype=np.float64)
    rng = np.random.default_rng(seed)
    n_rep = len(sizes)
    idx = rng.integers(0, n_rep, size=(n_boot, n_rep))
    num = hits[idx].sum(axis=1)
    den = sizes[idx].sum(axis=1)
    with np.errstate(invalid="ignore", divide="ignore"):
        r = np.where(den > 0, num / den, np.nan)
    a = (1.0 - conf) / 2.0
    return (float(np.nanpercentile(r, 100 * a)), float(np.nanpercentile(r, 100 * (1 - a))))


def rows_for(counter: Counter, **keys) -> List[Dict]:
    out = []
    n = counter.n_total
    for a in ALPHAS:
        x = counter.fp[a]
        lo, hi = cp_interval(x, n)
        ub = cp_upper_one_sided(x, n)
        lb = cp_lower_one_sided(x, n)
        exp_fp = n * a
        blo, bhi = cluster_bootstrap(counter.units, a)
        # How much wider the clustered interval is than the binomial one. This
        # is the design effect the non-independence buys, and it is the number
        # that decides whether the nominal Clopper-Pearson bound may be quoted
        # on its own.
        deff = ((bhi - blo) / (hi - lo)) if (np.isfinite(blo) and hi > lo) else np.nan
        n_zero = sum(c.get("_n_zero", 0) for _n, c in counter.units.values())
        # Inflation is decided by the LOWER limit, not the upper one — see
        # cp_lower_one_sided().  The upper limit answers the different and
        # equally necessary question "how much inflation can we rule out?".
        if n == 0:
            verdict = "no data"
        elif lb > a:
            verdict = "inflated"
        elif hi < a:
            verdict = "conservative"
        else:
            verdict = "consistent with nominal"
        out.append({
            **keys,
            "alpha": a,
            "n_null_tests": n,
            "n_fp": x,
            "expected_fp": exp_fp,
            # The point estimate is emitted always but flagged: at alpha = 5e-8
            # the design does not support it and the paper must quote the bound.
            "fpr_point": (x / n) if n else np.nan,
            "point_estimate_reportable": bool(exp_fp >= MIN_EXPECTED_FP_FOR_POINT),
            "cp_lo_95": lo,
            "cp_hi_95": hi,
            "cp_lower_95_onesided": lb,
            "cp_upper_95_onesided": ub,
            "cp_upper_over_alpha": ub / a,
            "verdict": verdict,
            "boot_lo_95": blo,
            "boot_hi_95": bhi,
            "design_effect_width_ratio": deff,
            "n_units": len(counter.units),
            "n_pvals_underflow_to_zero": n_zero,
            "fpr_ratio_point": (x / n / a) if (n and x) else np.nan,
        })
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--run-dir", default="outputs/power_runs_bench")
    ap.add_argument("--out-dir", default="revision/tables")
    ap.add_argument("--fig", default="revision/figures")
    ap.add_argument("--min-replicates", type=int, default=20)
    args = ap.parse_args()

    _self_test()
    print("Clopper-Pearson self-test passed (rule of three reproduced).\n")

    run_dir = Path(args.run_dir)
    out_dir = Path(args.out_dir); out_dir.mkdir(parents=True, exist_ok=True)
    fig_dir = Path(args.fig); fig_dir.mkdir(parents=True, exist_ok=True)

    print(f"streaming {run_dir}")
    manifest = frozen_manifest(run_dir)
    print(f"frozen manifest: {len(manifest)} replicates with both files")
    counters, status, n_files, agree = stream(run_dir, manifest)
    if n_files < args.min_replicates:
        print(f"WARNING: only {n_files} replicates available "
              f"(< --min-replicates {args.min_replicates}); "
              f"tail thresholds will be uninformative.", file=sys.stderr)

    methods = sorted({k[2] for k in counters})
    scenarios = sorted({k[0] for k in counters})
    configs = sorted({k[1] for k in counters})
    print(f"\n{n_files} replicates | {len(scenarios)} scenarios | "
          f"{len(configs)} configs | {len(methods)} methods")
    missing = [m for m in EXTERNAL if m not in methods]
    if missing:
        print(f"WARNING: external comparator(s) absent from the archive: {missing}",
              file=sys.stderr)

    # ---- pooled over the whole design (the stratum that supports alpha=5e-8)
    pooled: Dict[str, Counter] = {}
    for (scen, cfg, meth), c in counters.items():
        pooled.setdefault(meth, Counter()).merge(c)

    rows = []
    for meth, c in pooled.items():
        rows.extend(rows_for(c, stratum="pooled", scenario="ALL", config="ALL",
                             method=meth, method_label=METHOD_LABEL.get(meth, meth)))

    # ---- pooled, restricted to units every method completed
    # Without this, a method that failed on some units is scored only on the
    # units it survived, while its competitors are scored on all of them. That
    # is a comparison between different designs, and it favours whichever tool
    # failed most.
    common_units = None
    for meth, c in pooled.items():
        u = set(c.units)
        common_units = u if common_units is None else (common_units & u)
    common_units = common_units or set()
    n_all = len(set().union(*[set(c.units) for c in pooled.values()])) if pooled else 0
    print(f"units completed by every method: {len(common_units)} of {n_all}")
    if common_units and len(common_units) < n_all:
        for meth, c in pooled.items():
            rows.extend(rows_for(c.restricted(common_units), stratum="pooled_common_units",
                                 scenario="ALL", config="ALL", method=meth,
                                 method_label=METHOD_LABEL.get(meth, meth)))

    # ---- pooled over the scenarios where ancestry is actually resolvable
    # 2admixed_homo and 2admixed_hetero pair CEU5-YRI5 with CEU5-CHS5 from chr22
    # only. synthetic_power_simulation.py documents these as stress tests: the two
    # populations share half their ancestry, PCA cannot separate the non-CEU
    # components, and GMM assignment sits at chance. No method can exploit
    # structure that is not recoverable, so pooling them in dilutes every
    # method's inflation toward the same value and understates the contrast the
    # benchmark exists to measure. Reported both ways rather than chosen.
    structured: Dict[str, Counter] = {}
    for (scen, cfg, meth), c in counters.items():
        if scen.startswith("2admixed"):
            continue
        structured.setdefault(meth, Counter()).merge(c)
    for meth, c in structured.items():
        rows.extend(rows_for(c, stratum="pooled_resolvable_ancestry", scenario="4 of 6",
                             config="ALL", method=meth,
                             method_label=METHOD_LABEL.get(meth, meth)))

    # ---- per scenario
    per_scen: Dict[Tuple[str, str], Counter] = {}
    for (scen, cfg, meth), c in counters.items():
        per_scen.setdefault((scen, meth), Counter()).merge(c)
    for (scen, meth), c in per_scen.items():
        rows.extend(rows_for(c, stratum="scenario", scenario=scen, config="ALL",
                             method=meth, method_label=METHOD_LABEL.get(meth, meth)))

    # ---- fully stratified
    for (scen, cfg, meth), c in counters.items():
        rows.extend(rows_for(c, stratum="scenario_x_config", scenario=scen, config=cfg,
                             method=meth, method_label=METHOD_LABEL.get(meth, meth)))

    fpr = pd.DataFrame(rows)
    fpr_path = out_dir / "S10_arenaA_fpr.csv"
    fpr.to_csv(fpr_path, index=False, float_format="%.6g")
    print(f"\nwrote {fpr_path}  ({len(fpr)} rows)")

    # ---- tool status
    st = pd.DataFrame(
        [{"method": m, "method_label": METHOD_LABEL.get(m, m),
          "units_attempted": v[0], "units_all_nan": v[1], "units_partial_nan": v[2],
          "unit_failure_rate": v[1] / v[0] if v[0] else np.nan}
         for m, v in sorted(status.items())])
    st_path = out_dir / "S10d_arenaA_tool_status.csv"
    st.to_csv(st_path, index=False, float_format="%.6g")
    print(f"wrote {st_path}")

    # ---- GMMAT vs GENESIS: one estimator, two implementations
    if agree:
        ag = pd.DataFrame(agree)
        ag_path = out_dir / "S10e_arenaA_gmmat_genesis_agreement.csv"
        ag.to_csv(ag_path, index=False, float_format="%.6g")
        print(f"wrote {ag_path}")
        print("\nGMMAT vs GENESIS across "
              f"{len(ag)} units: median r(-log10 p) = {ag.pearson_r_neglog10p.median():.6f}, "
              f"median |diff| = {ag.median_abs_diff_neglog10p.median():.2e}, "
              f"worst |diff| = {ag.max_abs_diff_neglog10p.max():.2e} -log10 units.")
        print("These are two implementations of ONE estimator (logistic mixed "
              "model + score test, same GRM), so this is a wiring check, not "
              "independent corroboration. Report them as one comparator family.")

    # ---- power and lambda, from the per-replicate metrics CSVs
    mfiles = sorted(f for f in run_dir.glob("power_run_*_metrics.csv")
                    if int(f.name.split("_")[2]) in set(manifest))
    if mfiles:
        met = pd.concat([pd.read_csv(f) for f in mfiles], ignore_index=True)
        gcols = ["scenario_key", "config_key", "h2", "method"]
        agg = {c: ["mean", "std"] for c in
               [c for c in ("power_5e8", "power_5e7", "power_1e5", "power_1e3",
                            "lambda_gc", "lambda_gc_null", "pop_accuracy")
                if c in met.columns]}
        pw = met.groupby(gcols, dropna=False).agg(agg)
        pw.columns = ["_".join(c) for c in pw.columns]
        pw = pw.reset_index()
        pw["n_replicates"] = met.groupby(gcols, dropna=False).size().values
        pw["method_label"] = pw["method"].map(lambda m: METHOD_LABEL.get(m, m))
        pw_path = out_dir / "S10b_arenaA_power.csv"
        pw.to_csv(pw_path, index=False, float_format="%.6g")
        print(f"wrote {pw_path}  ({len(mfiles)} metrics files)")

        lam = (met.groupby(["method"], dropna=False)["lambda_gc_null"]
                  .agg(["mean", "std", "median", "min", "max", "count"]).reset_index())
        lam["method_label"] = lam["method"].map(lambda m: METHOD_LABEL.get(m, m))
        lam_path = out_dir / "S10c_arenaA_lambda.csv"
        lam.to_csv(lam_path, index=False, float_format="%.6g")
        print(f"wrote {lam_path}")
    else:
        pw = lam = None
        print("no metrics CSVs found; power/lambda tables skipped", file=sys.stderr)

    # ---- console summary, in the form the paper needs
    print("\n" + "=" * 96)
    print("ARENA A — pooled false-positive rate, all scenarios and configs")
    print("=" * 96)
    pool = fpr[fpr.stratum == "pooled"]
    for a in ALPHAS:
        sub = pool[pool.alpha == a].sort_values("method")
        n = int(sub.n_null_tests.iloc[0]) if len(sub) else 0
        print(f"\nalpha = {a:.0e}   ({n:,} null tests per method, "
              f"{n * a:.1f} expected false positives)")
        hdr = (f"  {'method':<26} {'obs':>6} {'FPR':>14} {'95% CP upper':>13} "
               f"{'UB/alpha':>9}  verdict")
        print(hdr); print("  " + "-" * (len(hdr) - 2))
        for _, r in sub.iterrows():
            pt = f"{r.fpr_point:.3e}" if r.point_estimate_reportable else "not reportable"
            print(f"  {r.method_label:<26} {int(r.n_fp):>6} {pt:>14} "
                  f"{r.cp_upper_95_onesided:>13.3e} {r.cp_upper_over_alpha:>9.2f}  {r.verdict}")

    if lam is not None:
        print("\n" + "=" * 96)
        print("ARENA A — null-block lambda_GC by method (mean over replicates)")
        print("=" * 96)
        for _, r in lam.sort_values("mean").iterrows():
            print(f"  {r.method_label:<26} {r['mean']:.4f} +/- {r['std']:.4f} "
                  f"(median {r['median']:.4f}, n = {int(r['count'])})")

    print("\n" + "=" * 96)
    print("ARENA A — external tool status")
    print("=" * 96)
    for _, r in st.iterrows():
        print(f"  {r.method_label:<26} attempted {r.units_attempted:>6}  "
              f"all-NaN {r.units_all_nan:>4}  partial-NaN {r.units_partial_nan:>4}")

    print("\nNOTE for the manuscript: at alpha = 5e-8 the point estimate is "
          f"{'reportable' if (pool[(pool.alpha == 5e-8)].expected_fp.iloc[0] >= MIN_EXPECTED_FP_FOR_POINT) else 'NOT reportable'}"
          " at this replicate budget; quote the Clopper-Pearson upper bound.")
    print("Binomial intervals assume independent tests; the replicate-level "
          "bootstrap columns (boot_lo_95/boot_hi_95) are the conservative "
          "alternative and should be quoted where they are wider.")
    print("\nDONE")


if __name__ == "__main__":
    main()
