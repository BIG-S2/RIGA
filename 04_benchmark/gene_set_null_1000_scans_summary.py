#!/usr/bin/env python
"""gene_set_null_1000_scans_summary.py (2026-10-19, package geneset_null_1000_20261019; revised the same day after review) -- copy of
revision/scripts/gene_set_null_native_tests.py (2026-10-13) with the output path, the result root and the design list as
arguments, a --max-rep filter, a per-design bootstrap generator, a meta.json requirement and a short TXT summary. The per-scan
computation is unchanged: the primary SB (Cauchy combination of the SKAT(1,25), burden(1,25), SKAT(1,1) and burden(1,1) p values
= STAAR-O without ACAT-V), STAAR-O, burden and SKAT, with and without the cMAC >= 10 floor; per test the observed/expected ratio
at 1e-3, 1e-4, 1e-5, 2.5e-6 and 1e-6 with a percentile bootstrap over replicates, and the share of null scans with any BH FDR-5%
discovery (Clopper-Pearson 95%). Read-only on inputs.

Differences from the frozen summariser, and why
  * --designs defaults to strat_project ONLY. The rows of the other designs (pooled_project 100 scans, saturated_cell 50) are
    copied verbatim from the frozen 100-scan table fpr_rv_native_tests.csv (Supplementary Table 15b) into the output, marked in
    a 'source' column (--no-append-frozen omits them). Reason: the frozen script draws every bootstrap from ONE generator
    consumed in design order (saturated_cell, strat_project, pooled_project), so recomputing pooled_project after 1,000
    stratified files would change its boot_lo/boot_hi although nothing about that design changed.
  * The bootstrap generator is seeded per design, np.random.default_rng([seed, design_index]) with the index in the canonical
    order (saturated_cell 0, strat_project 1, pooled_project 2), so a design's rows no longer depend on which other designs are
    run in the same call. Consequence: on the same input files the deterministic columns (n_rep, n_tests_mean, observed,
    expected, ratio, share_any_fdr, share_lo, share_hi) reproduce the frozen table exactly, while boot_lo/boot_hi are a
    different random draw of the same bootstrap (verified 19 Oct 2026 on the frozen 100 stratified files; see README.md).
  * --max-rep (default 1000) and a meta.json requirement: a genes file without its meta.json is an unfinished replicate and is
    skipped (listed in the TXT).
Default output revision/tables/plasmode_power/fpr_rv_native_tests_1000.csv (+ fpr_rv_native_tests_1000_summary.txt); the 100-scan
table is never written (asserted); an existing output is moved aside with a timestamp, never overwritten.
"""
import argparse, glob, os, re, sys, time
import numpy as np, pandas as pd
from scipy.stats import beta

ROOT_DEFAULT = "${ADSP_ROOT}/outputs/plasmode/rv_tools/genome/tools"
TAB = "${ADSP_ROOT}/revision/tables/plasmode_power"
OUT_DEFAULT = f"{TAB}/fpr_rv_native_tests_1000.csv"
TXT_DEFAULT = f"{TAB}/fpr_rv_native_tests_1000_summary.txt"
FROZEN_100 = f"{TAB}/fpr_rv_native_tests.csv"
ALPHAS = [1e-3, 1e-4, 1e-5, 2.5e-6, 1e-6]
DESIGN_INDEX = {"saturated_cell": 0, "strat_project": 1, "pooled_project": 2}
COLS = ["design", "test", "cmac_floor", "n_rep", "n_tests_mean", "alpha", "observed", "expected", "ratio", "boot_lo", "boot_hi", "share_any_fdr", "share_lo", "share_hi"]


def cauchy(P):
    P = np.clip(np.asarray(P, float), 1e-300, 1 - 1e-16)
    t = np.where(P < 1e-15, 1.0 / (P * np.pi), np.tan((0.5 - P) * np.pi)).mean(axis=1)
    return np.where(t > 1e15, 1.0 / (t * np.pi), 0.5 - np.arctan(t) / np.pi)


def bh_n(p, q=0.05):
    p = np.sort(p[np.isfinite(p)]); m = len(p); ok = np.nonzero(p <= q * np.arange(1, m + 1) / m)[0]
    return 0 if len(ok) == 0 else ok[-1] + 1


def clopper_pearson(k, n):
    lo = beta.ppf(0.025, k, n - k + 1) if k > 0 else 0.0
    hi = beta.ppf(0.975, k + 1, n - k) if k < n else 1.0
    return lo, hi


def safe_out(path):
    if os.path.exists(path):
        base, ext = os.path.splitext(path)
        moved = f"{base}.{time.strftime('%Y%m%dT%H%M%S')}{ext}"
        os.replace(path, moved)
        print(f"existing {path} moved to {moved}")
    return path


def rep_of(f):
    return int(re.search(r"_rep(\d+)_genes", os.path.basename(f)).group(1))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", default=ROOT_DEFAULT)
    ap.add_argument("--out", default=OUT_DEFAULT)
    ap.add_argument("--txt", default=TXT_DEFAULT)
    ap.add_argument("--designs", default="strat_project", help="designs to recompute (comma list from saturated_cell,strat_project,pooled_project)")
    ap.add_argument("--max-rep", type=int, default=1000, help="use replicates <= this (default 1000; 100 reproduces the frozen deterministic columns)")
    ap.add_argument("--seed", type=int, default=20261013, help="bootstrap seed; the generator is default_rng([seed, design index])")
    ap.add_argument("--no-append-frozen", action="store_true", help="do not copy the frozen table's rows of the designs that are not recomputed")
    ap.add_argument("--frozen", default=FROZEN_100, help="the frozen 100-scan table whose other-design rows are copied")
    ap.add_argument("--limit", type=int, default=0, help="debug: use at most this many files per design")
    a = ap.parse_args()
    assert os.path.abspath(a.out) != os.path.abspath(FROZEN_100), "refusing to write the frozen 100-scan table"
    designs = [d for d in a.designs.split(",") if d]
    assert all(d in DESIGN_INDEX for d in designs), f"unknown design in {designs}"
    rows, found, shares, skipped = [], {}, {}, {}
    for design in designs:
        rng = np.random.default_rng([a.seed, DESIGN_INDEX[design]])      # per-design stream (see the module docstring)
        files = sorted(f for f in glob.glob(f"{a.root}/{design}_rep*_genes.tsv.gz") if "_syn_" not in f)
        files = [f for f in files if rep_of(f) <= a.max_rep]
        skipped[design] = sorted(rep_of(f) for f in files if not os.path.exists(f.replace("_genes.tsv.gz", "_meta.json")))
        files = [f for f in files if rep_of(f) not in skipped[design]]
        if a.limit:
            files = files[:a.limit]
        found[design] = sorted(rep_of(f) for f in files)
        acc = {}
        for f in files:
            d = pd.read_csv(f, sep="\t"); d = d[np.isfinite(d.p_staar_O)]
            tests = {"SB": cauchy(d[["p_staar_S_1_25", "p_staar_B_1_25", "p_staar_S_1_1", "p_staar_B_1_1"]].to_numpy()), "STAAR-O": d.p_staar_O.to_numpy(),
                     "burden": cauchy(d[["p_staar_B_1_25", "p_staar_B_1_1"]].to_numpy()), "SKAT": cauchy(d[["p_staar_S_1_25", "p_staar_S_1_1"]].to_numpy())}
            for floor in (0, 10):
                keep = (d.cMAC >= floor).to_numpy()
                for t, p in tests.items():
                    pk = p[keep]; ac = acc.setdefault((t, floor), dict(n=[], any=[], k={al: [] for al in ALPHAS}))
                    ac["n"].append(pk.size); ac["any"].append(int(bh_n(pk) > 0))
                    for al in ALPHAS: ac["k"][al].append(int(np.sum(pk < al)))
        for (t, floor), ac in acc.items():
            n = np.asarray(ac["n"], float); nrep = n.size; nany = int(np.sum(ac["any"]))
            lo_s, hi_s = clopper_pearson(nany, nrep)
            shares[(design, t, floor)] = (nany, nrep, nany / nrep, lo_s, hi_s)
            for al in ALPHAS:
                k = np.asarray(ac["k"][al], float); idx = rng.integers(0, nrep, size=(2000, nrep))
                boot = k[idx].sum(axis=1) / (n[idx].sum(axis=1) * al)
                rows.append(dict(design=design, test=t, cmac_floor=floor, n_rep=nrep, n_tests_mean=n.mean(), alpha=al, observed=int(k.sum()),
                                 expected=float((n * al).sum()), ratio=float(k.sum() / (n * al).sum()), boot_lo=np.percentile(boot, 2.5), boot_hi=np.percentile(boot, 97.5),
                                 share_any_fdr=nany / nrep, share_lo=lo_s, share_hi=hi_s))
    res = pd.DataFrame(rows, columns=COLS)
    res["source"] = f"recomputed {time.strftime('%Y-%m-%d')} by geneset_null_1000_20261019/gene_set_null_1000_scans_summary.py (per-design bootstrap seed [{a.seed}, design index]) from {a.root}"
    res.to_csv(safe_out(a.out), index=False)
    copied, n_copied = [], 0
    if not a.no_append_frozen and os.path.exists(a.frozen):
        # the frozen rows are appended as TEXT (no float parsing/re-printing: pandas' fast parser is not round-trip exact in the
        # last digit), so they are byte-identical to the frozen table's lines; the recomputed rows above come first
        src = f"copied verbatim from the frozen 100-scan table {os.path.basename(a.frozen)} (2026-10-13; single-stream bootstrap, seed {a.seed})"
        with open(a.frozen) as fh:
            hdr = fh.readline().rstrip("\n").split(",")
            assert hdr == COLS, f"frozen table columns {hdr} != {COLS}"
            keep = [ln.rstrip("\n") for ln in fh if ln.strip() and ln.split(",")[0] not in designs]
        copied = sorted({ln.split(",")[0] for ln in keep})
        n_copied = len(keep)
        with open(a.out, "a") as fh:
            for ln in keep:
                fh.write(f'{ln},"{src}"\n')
        res = pd.concat([res, pd.read_csv(a.out).iloc[len(res):]], ignore_index=True)
    print(res[(res.alpha == 2.5e-6)][["design", "test", "cmac_floor", "n_rep", "observed", "ratio", "boot_lo", "boot_hi", "share_any_fdr"]].round(3).to_string(index=False))
    # ---- short TXT summary
    L = [f"Permuted gene-set null (within-cell label permutation), STAAR-family tests -- {time.strftime('%Y-%m-%d %H:%M')}",
         f"Result root: {a.root}", f"Table: {a.out}", f"Recomputed designs: {', '.join(designs)} (replicates <= {a.max_rep})",
         f"Rows copied unchanged (as text) from the frozen 100-scan table ({a.frozen}): {', '.join(copied) if copied else 'none'} ({n_copied} rows)", ""]
    for design in designs:
        reps = found.get(design, [])
        if design == "strat_project":
            missing = sorted(set(range(1, a.max_rep + 1)) - set(reps))
            L.append(f"[{design}] replicates found: {len(reps)} of {a.max_rep} requested; missing {len(missing)}: {missing if len(missing) <= 60 else str(missing[:60]) + ' ...'}")
        else:
            L.append(f"[{design}] replicates found: {len(reps)}")
        if skipped.get(design):
            L.append(f"  skipped (genes file without meta.json = unfinished): {skipped[design]}")
        L.append("  share of null scans with any BH 5% discovery = realised FDR under the complete null (Clopper-Pearson 95%)")
        L.append(f"  {'test':10s} {'cMAC floor':10s} {'scans with hit':>15s} {'share':>8s} {'CP 95% interval':>20s}")
        for floor in (10, 0):
            for t in ("SB", "STAAR-O", "burden", "SKAT"):
                if (design, t, floor) in shares:
                    k, n, s, lo, hi = shares[(design, t, floor)]
                    L.append(f"  {t:10s} {('>= ' + str(floor)) if floor else 'none':10s} {f'{k} / {n}':>15s} {s:8.3f} {f'[{lo:.3f}, {hi:.3f}]':>20s}")
        sub = res[(res.design == design) & (res.cmac_floor == 10) & (res.alpha == 2.5e-6)]
        for _, r in sub.iterrows():
            L.append(f"  tail at 2.5e-6, cMAC >= 10, {r.test:8s}: observed {int(r.observed)} expected {r.expected:.1f} ratio {r.ratio:.2f} [{r.boot_lo:.2f}, {r.boot_hi:.2f}]")
        L.append("")
    L.append("Reading: under a complete null every discovery is false, so a scan's false-discovery proportion is 1 if it has any BH")
    L.append("discovery and 0 otherwise; the share of scans with any discovery is therefore the realised FDR of the procedure, to be")
    L.append("compared with the nominal 5%. With 1,000 scans the Clopper-Pearson interval at 8% is about 6.4-9.9% (+/- 1.7 points).")
    L.append("Bootstrap intervals (boot_lo/boot_hi) use one generator per design, default_rng([seed, design index]); the frozen 100-scan")
    L.append("table used a single stream in design order, so its boot columns are a different draw of the same bootstrap. The")
    L.append("deterministic columns of a design are identical for identical input files.")
    txt = "\n".join(L) + "\n"
    with open(safe_out(a.txt), "w") as fh:
        fh.write(txt)
    print(txt)
    return 0


if __name__ == "__main__":
    sys.exit(main())
