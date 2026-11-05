#!/usr/bin/env python
# ukb_build_replication_tables.py (2026-10-19): copy of revision/scripts/ukb_build_S11_S13_reconstructed.py reading the
# root-checked lookups (tables/review/ukb_rootcheck_genesets.csv, ukb_rootcheck_sv.csv) and writing
# supp_v3/S11_ukb_geneset_replication_rootcheck.csv and S13_ukb_sv_direct_rootcheck.csv. The only changes are those
# input/output names and the dropped backup-copy line; nothing existing is overwritten.
"""ukb_build_S11_S13_reconstructed.py (2026-10-13): Supplementary Tables 11 (UK Biobank replication of the 14 primary gene sets)
and 13 (direct test of the single variants) from the corrected-cohort saddlepoint lookups (ukb_03_reconstructed_spa_lookups.R).
Writes revision/tables/supp_v3/S11_ukb_geneset_replication.csv and S13_ukb_sv_direct.csv (old files kept with a date suffix)."""
import shutil
from pathlib import Path
import numpy as np, pandas as pd
REV = Path("${ADSP_ROOT}/revision"); R = REV / "tables/review"; S = REV / "tables/supp_v3"
g = pd.read_csv(R / "ukb_rootcheck_genesets.csv")
fam = {"coding": "Coding", "noncoding": "Noncoding", "ncRNA": "ncRNA"}
s11 = pd.DataFrame({"Family": g.family.map(fam), "Gene": g.gene, "Mask": g.cat, "Chr": g.chr, "#SNV": g.nSNV, "cMAC": g.cMAC, "ADSP SB P": g.adsp_SB,
                    "UKB burden P (SPA)": g.SPA_burden, "Verdict": np.where(g.SPA_burden < 0.05 / 14, "replicates (P < 0.05/14)", np.where(g.SPA_burden < 0.05, "nominal P < 0.05", "not replicated"))})
s11 = s11.sort_values(["Family", "ADSP SB P"]); s11.to_csv(S / "S11_ukb_geneset_replication_rootcheck.csv", index=False)
print(s11.to_string(index=False))
sv = pd.read_csv(R / "ukb_rootcheck_sv.csv"); sv = sv[sv.gene != "APOE"]; s2 = pd.read_csv(S / "S2_fdr_single_variants.csv")
s2 = s2.rename(columns={"CHR": "chr", "POS": "pos", "REF": "ref", "ALT": "alt"})[["chr", "pos", "ref", "alt", "MAF", "OR", "P"]].rename(columns={"MAF": "ADSP MAF", "OR": "ADSP OR", "P": "ADSP P"})
m = sv.merge(s2, on=["chr", "pos", "ref", "alt"], how="left")
m["locus"] = m.chr.astype(str) + ":" + (m.pos // 1_000_000).astype(str)                       # 2026-10-14: 1-Mb loci; the lead (smallest ADSP P) carries the replication test
m["is_lead"] = m.groupby("locus")["ADSP P"].transform("min") == m["ADSP P"]; nloci = m.locus.nunique()
verdict = np.where(m.pvalue.isna(), "absent", np.where(~m.is_lead, "secondary variant of the locus", np.where(m.pvalue < 0.05 / nloci, f"replicates (P < 0.05/{nloci})", np.where(m.pvalue < 0.05, "nominal P < 0.05", "not replicated"))))
s13 = pd.DataFrame({"Gene": m.gene, "Variant": m.chr.astype(str) + ":" + m.pos.astype(str) + ":" + m.ref + ":" + m.alt, "Locus lead": np.where(m.is_lead, "lead", ""), "ADSP MAF": m["ADSP MAF"], "ADSP OR": m["ADSP OR"], "ADSP P": m["ADSP P"],
                    "UKB MAF": m.MAF, "UKB N": m.N, "UKB P (SPA)": m.pvalue, "Verdict": verdict})
s13 = s13.sort_values(["ADSP P"])
# 2026-10-19 rootcheck: the backup-copy line of the original is dropped; nothing existing is written
s13.to_csv(S / "S13_ukb_sv_direct_rootcheck.csv", index=False); print(s13.to_string(index=False))
