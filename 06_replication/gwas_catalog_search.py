#!/usr/bin/env python3
"""Search the full GWAS Catalog association file for every manuscript gene, using BOTH
MAPPED_GENE and 'REPORTED GENE(S)' columns, and flag AD/dementia-relevant traits.
"""
import csv, json, re, sys
from collections import defaultdict

CAT = "${SCRATCH_DIR}/0ecc26ea-f4e5-4737-9951-5392110419ec/scratchpad/gwas-catalog-download-associations-alt-full.tsv"
OUT = "${ADSP_ROOT}/revision/logs/catalog_hits.json"

GENES = {
 # printed symbol -> set of symbols to match (current + deprecated + the printed one)
 "CACNA1D": {"CACNA1D"}, "PRKAR1B": {"PRKAR1B"}, "CKMT1A": {"CKMT1A", "CKMT1"},
 "CATSPER2": {"CATSPER2"}, "ZNF532": {"ZNF532"}, "DCP1A": {"DCP1A"},
 "PIK3R6": {"PIK3R6"}, "SEMA4A": {"SEMA4A"}, "ADAMTS14": {"ADAMTS14"},
 "MYO1E": {"MYO1E"}, "APOE": {"APOE"}, "KIF18B": {"KIF18B"}, "CASP4": {"CASP4"},
 "USH2A": {"USH2A"}, "NEB": {"NEB"}, "PCDHB14": {"PCDHB14"}, "RUNX1": {"RUNX1"},
 "APMAP": {"APMAP", "C20orf3"}, "GSC": {"GSC"}, "HIST1H4F": {"HIST1H4F", "H4C6", "HIST1H4F"},
 "ZC3H12A": {"ZC3H12A"}, "TNFRSF1B": {"TNFRSF1B"}, "EBF2": {"EBF2"},
 "SNCAIP": {"SNCAIP"}, "PRKCZ": {"PRKCZ"}, "SIRT2": {"SIRT2"}, "CNTNAP1": {"CNTNAP1"},
 "FGG": {"FGG"}, "CST1": {"CST1"}, "ZNF337": {"ZNF337"},
 "C19orf71": {"C19orf71", "TEKTIP1"}, "QPRT": {"QPRT"}, "ANP32B": {"ANP32B"},
 "SLC6A15": {"SLC6A15"},
 "PSEN1": {"PSEN1"}, "TOMM40": {"TOMM40"}, "ELMSAN1": {"ELMSAN1", "MIDEAS", "MIER5", "C14orf117"},
 "MIR4505": {"MIR4505"},
 # confusables, for contrast
 "PIK3R5": {"PIK3R5"}, "CKMT1B": {"CKMT1B"}, "SNCA": {"SNCA"}, "CNTNAP2": {"CNTNAP2"},
 "FGB": {"FGB"}, "CST3": {"CST3"}, "ACE": {"ACE"}, "H4C5": {"H4C5"},
}

SYM2KEY = defaultdict(set)
for k, syns in GENES.items():
    for s in syns:
        SYM2KEY[s.upper()].add(k)

AD_RE = re.compile(
    r"alzheimer|dementia|cognitive declin|cognitive impair|"
    r"amyloid|neurofibrillary|\btau\b|tauopath|braak|cerad|"
    r"lewy body|neurodegenerat|mild cognitive|"
    r"neuropatholog|apolipoprotein e|apoe|hippocampal atroph|"
    r"cognitive function|memory performance|episodic memory|"
    r"age.related cognitive", re.I)

# tighter core-AD regex (for the strict verdict)
CORE_RE = re.compile(r"alzheimer|dementia|mild cognitive impair|"
                     r"neurofibrillary|braak|cerad|lewy body|amyloid", re.I)

SPLIT = re.compile(r"[,;]| - |/")


def syms_in(cell):
    if not cell:
        return set()
    return {t.strip().upper() for t in SPLIT.split(cell) if t.strip()}


res = defaultdict(lambda: {"n_assoc": 0, "traits": defaultdict(int),
                           "ad_hits": [], "core_ad_hits": []})

with open(CAT, newline="", encoding="utf-8", errors="replace") as f:
    rd = csv.DictReader(f, delimiter="\t")
    for i, row in enumerate(rd):
        mapped = row.get("MAPPED_GENE") or ""
        reported = row.get("REPORTED GENE(S)") or ""
        hit = syms_in(mapped) | syms_in(reported)
        keys = set()
        for s in hit:
            if s in SYM2KEY:
                keys |= SYM2KEY[s]
        if not keys:
            continue
        trait = row.get("DISEASE/TRAIT", "")
        rec = {"trait": trait, "pmid": row.get("PUBMEDID"),
               "author": row.get("FIRST AUTHOR"), "date": row.get("DATE"),
               "snp": row.get("SNPS"), "p": row.get("P-VALUE"),
               "ptext": row.get("P-VALUE (TEXT)", ""),
               "mapped": mapped, "reported": reported,
               "n": (row.get("INITIAL SAMPLE SIZE") or "")[:110]}
        isad = bool(AD_RE.search(trait))
        iscore = bool(CORE_RE.search(trait))
        for k in keys:
            r = res[k]
            r["n_assoc"] += 1
            r["traits"][trait] += 1
            if isad:
                r["ad_hits"].append(rec)
            if iscore:
                r["core_ad_hits"].append(rec)

final = {}
for k, r in res.items():
    tr = sorted(r["traits"].items(), key=lambda x: -x[1])
    # dedupe AD hits by (pmid, snp, trait)
    def dedupe(lst):
        seen, o = set(), []
        for x in lst:
            key = (x["pmid"], x["snp"], x["trait"])
            if key in seen:
                continue
            seen.add(key)
            o.append(x)
        return o
    final[k] = {"n_assoc": r["n_assoc"], "n_traits": len(r["traits"]),
                "top_traits": tr[:8],
                "n_ad": len(dedupe(r["ad_hits"])),
                "n_core_ad": len(dedupe(r["core_ad_hits"])),
                "ad_hits": dedupe(r["ad_hits"])[:40],
                "core_ad_hits": dedupe(r["core_ad_hits"])[:40]}

json.dump(final, open(OUT, "w"), indent=1)

for k in sorted(final, key=lambda x: -final[x]["n_core_ad"]):
    v = final[k]
    print(f"{k:10s} n_assoc={v['n_assoc']:5d} n_traits={v['n_traits']:4d} "
          f"n_AD_broad={v['n_ad']:3d} n_AD_core={v['n_core_ad']:3d}")
print("WROTE", OUT, file=sys.stderr)
