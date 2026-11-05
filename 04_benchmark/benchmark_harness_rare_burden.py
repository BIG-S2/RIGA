#!/usr/bin/env python
"""
benchmark_harness_rare_burden.py -- benchmark_harness.py with the MAF < 0.1 % burden masks (2026-09-22 user rule:
rare-variant tests include only variants with MAF < 0.1 %, i.e. pool MAC <= 54).

The harness file is untouched; this wrapper patches it in memory:
  * BURDEN_MASKS (the masks that receive causal burden units) -> coding_maf001, gene_maf001, coding_mac6_54, gene_mac6_54
    (every mask found in --burden-dir is still scored as null units: syn_maf001 and win100k_maf001 included);
  * --burden-dir defaults to outputs/plasmode/burden_chr22_maf001 (plasmode_burden_chr22_maf001.R);
  * --all-arms expands to --arms <every arm of the harness> (null-arm mode scores only the 11 new arms by default);
  * --strip-snp (null-arm mode) deletes the single-variant keys (p_chi2__*, p_spa__*, flags__*) from the npz this run
    wrote, so that plasmode_power_aggregate.py --nullarms-npz picks up ONLY the burden null results of this run; the
    single-variant null pools already hold these Y* (outputs/plasmode/native, power/nullarms, power/nullarms2).
    Use with --max-chunks 1 to make the single-variant scan a 2,000-SNP formality.
Every other option is passed through unchanged.
"""
from __future__ import annotations

import glob
import re
import sys
from pathlib import Path

import numpy as np

SCRIPTS = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPTS))
import benchmark_harness as P  # noqa: E402

P.BURDEN_MASKS = ("coding_maf001", "gene_maf001", "coding_mac6_54", "gene_mac6_54")
MAF_DIR = str(P.H.OUT_ROOT / "burden_chr22_maf001")


def opt(argv, name, default=None):
    if name in argv:
        i = argv.index(name)
        return argv[i + 1] if i + 1 < len(argv) else default
    return default


def main() -> int:
    argv = list(sys.argv[1:])
    strip = "--strip-snp" in argv
    argv = [a for a in argv if a != "--strip-snp"]
    if "--all-arms" in argv:
        argv = [a for a in argv if a != "--all-arms"] + ["--arms", ",".join(P.ARM_ORDER)]
    if "--burden-dir" not in argv:
        argv += ["--burden-dir", MAF_DIR]
    sys.argv = [str(SCRIPTS / "benchmark_harness.py")] + argv
    print(f"[maf001 wrapper] BURDEN_MASKS={P.BURDEN_MASKS} burden_dir={opt(argv, '--burden-dir')} strip_snp={strip}", flush=True)
    rc = P.main()
    if rc not in (0, None):
        return int(rc)
    if strip:
        rep = int(opt(argv, "--replicate", "0"))
        out_root = Path(opt(argv, "--out-root", str(P.OUT_ROOT)))
        null_mode = bool(opt(argv, "--ystar-from", ""))
        sub = opt(argv, "--out-sub", "") or ("nullarms" if null_mode else "native")
        hits = sorted(glob.glob(str(out_root / sub / "npz" / f"rep{rep:03d}_seed*.npz")))
        if not hits:
            print(f"[maf001 wrapper] WARNING: no npz to strip under {out_root / sub / 'npz'} for replicate {rep}", flush=True)
            return 0
        for f in hits:
            z = dict(np.load(f, allow_pickle=False))
            drop = [k for k in z if re.match(r"^(p_chi2__|p_spa__|flags__)", k)]
            for k in drop:
                del z[k]
            np.savez_compressed(f, **z)
            print(f"[maf001 wrapper] stripped {len(drop)} single-variant keys from {Path(f).name}; kept {len(z)} keys "
                  f"({sum(k.startswith('burden_') for k in z)} burden keys)", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
