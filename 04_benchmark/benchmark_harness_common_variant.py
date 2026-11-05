#!/usr/bin/env python
"""
benchmark_harness_common_variant.py -- benchmark_harness.py for the LARGE common-variant simulation (2026-09-23):
  * SNP effect levels c = 0.09 / 0.12 / 0.15 log-odds per SD of genotype (per-allele OR 1.15 / 1.20 / 1.26 at MAF 0.30),
    one step up from the 0.06 / 0.09 / 0.12 of the first run, so that the 5e-8 read-out sits in the informative range;
  * its own seed manifest (SeedSequence(20260924).spawn(40) -> outputs/plasmode/power/seed_manifest_bigpower.csv), so the
    replicates are independent of tags main / clean / asc_clean;
  * no burden masks (--burden-dir /nonexistent unless given): single variants only.
The harness file is untouched; everything else is passed through.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

SCRIPTS = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPTS))
import benchmark_harness as P  # noqa: E402

P.SNP_LEVELS = (0.09, 0.12, 0.15)
MANIFEST = P.OUT_ROOT / "seed_manifest_bigpower.csv"


def write_manifest(path: Path, n: int = 40) -> None:
    if path.exists():
        return
    children = np.random.SeedSequence(20260924).spawn(n)   # own root: the null manifest uses 20260923
    rows = [dict(task=i + 1, task_type="bigpower_native", spawn_key=i,
                 seed=int(c.generate_state(1, dtype=np.uint64)[0])) for i, c in enumerate(children)]
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(".tmp.csv")
    pd.DataFrame(rows).to_csv(tmp, index=False)
    tmp.rename(path)


def main() -> int:
    argv = list(sys.argv[1:])
    if "--manifest" not in argv:
        argv += ["--manifest", str(MANIFEST)]
    if "--burden-dir" not in argv:
        argv += ["--burden-dir", "/nonexistent"]
    mi = argv.index("--manifest")
    write_manifest(Path(argv[mi + 1]))
    sys.argv = [str(SCRIPTS / "benchmark_harness.py")] + argv
    print(f"[big wrapper] SNP_LEVELS={P.SNP_LEVELS} manifest={argv[mi + 1]}", flush=True)
    return P.main()


if __name__ == "__main__":
    sys.exit(main())
