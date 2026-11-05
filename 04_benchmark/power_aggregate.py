#!/usr/bin/env python
"""
power_aggregate.py -- plasmode_power_aggregate.py for the LARGE common-variant simulation (2026-09-23):
tag big = power replicates outputs/plasmode/power/bigpower (six native designs, effect levels 0.09 / 0.12 / 0.15),
permuted-null pools outputs/plasmode/power/bignull/npz (up to 1,000 Y*), external tools (REGENIE Firth + LOCO, GMMAT,
GENESIS; all with centre + collapsed project + PC1-20) from outputs/plasmode/power/tools_{bignull,bigpower}/npz.
The aggregator file is untouched: its REGENIE set is widened in memory to the covset names gmmat_project and
genesis_project (arm names regenie_gmmat_project / regenie_genesis_project; the test label "firth" is bookkeeping:
GMMAT is a plain score test and GENESIS a saddlepoint-corrected score test), and no burden tables are requested.
Usage: python power_aggregate.py [extra aggregator options]
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import plasmode_power_aggregate as A  # noqa: E402

ROOT = Path("${ADSP_ROOT}")
A.REG_SETS = ["pooled_project", "gmmat_project", "genesis_project"]
A.ARM_LABEL.update({"regenie_gmmat_project": "G5", "regenie_genesis_project": "E5"})
A.ARM_DESC.update({"regenie_gmmat_project": "GMMAT glmmkin + score test, sparse GRM, + centre + project",
                   "regenie_genesis_project": "GENESIS fitNullModel + Score.SPA, sparse GRM, + centre + project"})
sys.argv = [sys.argv[0], "--tag", "big",
            "--power-dir", str(ROOT / "outputs/plasmode/power/bigpower"),
            "--burden-dir", "/nonexistent", "--burden-dense-dir", "",
            "--null-npz", "/nonexistent",
            "--nullarms-npz", str(ROOT / "outputs/plasmode/power/bignull/npz"),
            "--reg-null-npz", str(ROOT / "outputs/plasmode/power/tools_bignull/npz"),
            "--reg-power-npz", str(ROOT / "outputs/plasmode/power/tools_bigpower/npz"),
            "--max-rep", "1000",
            "--generator-note", "LARGE common-variant run: permuted null = observed labels permuted inside the 128 cohort x project x "
            "centre cells (new seeds, up to 1,000 replicates; external tools on the first 200); power = within-cell generator, clean "
            "causal class, effect levels c = 0.09 / 0.12 / 0.15 per SD (per-allele OR 1.15 / 1.20 / 1.26 at MAF 0.30), 40 replicates; "
            "designs: oracle N0, pooled + PCs N1, pooled + centre + project N5, RIGA 1:1 centre key M1 and project key M3, stratified S1, "
            "REGENIE Firth + LOCO (R5), GMMAT (G5) and GENESIS (E5) with centre + project + PC1-20 and the pool's sparse GRM"] + sys.argv[1:]
sys.exit(A.main())
