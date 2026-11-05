#!/usr/bin/env python3
"""Build novelty_audit.tsv and novelty_audit.md from the verified resolution/catalog/lit JSONs."""
import json, io

BASE = "${ADSP_ROOT}/revision"
RES = json.load(open(f"{BASE}/logs/gene_resolution.json"))
NEI = json.load(open(f"{BASE}/logs/neighbors.json"))
CAT = json.load(open(f"{BASE}/logs/catalog_hits.json"))

# ---------------------------------------------------------------- curated fields
# prior_gwas : verified GWAS Catalog AD/dementia rows (rsID, p, PMID) -- "" if none
# prior_lit  : verified primary rare-variant/functional AD literature (PMID) -- "" if none
# verdict    : novel | not previously reported in ADSP | previously reported - reframe as replication | drop
C = {
"CACNA1D": dict(cls="single-variant",
  prior_gwas="none AD/dementia (329 assoc., 0 AD; region dominated by blood-pressure signals)",
  prior_lit="No direct AD genetic report. Cav1.3 calcium-signalling literature is Parkinson-centred (PMID 24781149; PMID 36592226).",
  verdict="not previously reported in ADSP",
  note="OVERLAPS DCP1A by 18,720 bp -- the draft counts DCP1A and CACNA1D as two separate 3p21.1 loci. They are one locus."),
"DCP1A": dict(cls="single-variant",
  prior_gwas="none AD/dementia (24 assoc., 0 AD)",
  prior_lit="none AD-specific (4 PubMed AD-context papers, none genetic)",
  verdict="not previously reported in ADSP",
  note="OVERLAPS CACNA1D by 18,720 bp. SYMBOL HAZARD: DCP1 is the previous HGNC symbol of ACE (HGNC:2707), a top-tier AD gene with 14 AD Catalog rows; DCP1A is HGNC:18714. Any pathway/text-mining tool fed 'DCP1' will return ACE."),
"PRKAR1B": dict(cls="single-variant",
  prior_gwas="none AD/dementia (72 assoc., 0 AD)",
  prior_lit="ESTABLISHED dementia gene: PRKAR1B p.Leu50Arg causes a neurodegenerative dementia with unique pathology (PMID 24722252, Brain 2014); screened in early-onset dementia cohorts (PMID 25108559); mechanism confirmed across species (PMID 38743596, Brain 2024); NLPD-PKA allostery (PMID 40244081).",
  verdict="previously reported - reframe as replication",
  note="Catalog-clean but NOT novel. The largest single novelty correction in the draft."),
"CKMT1A": dict(cls="single-variant",
  prior_gwas="none AD/dementia (1 assoc. total: diastolic BP, PMID 27841878, and that row is MAPPED to the pseudogenes STRCP1-CATSPER2P1)",
  prior_lit="none AD-specific",
  verdict="not previously reported in ADSP",
  note="24,768 bp from CATSPER2 -- same 15q15.3 locus, not two. 87,985 bp from its segmental duplicate CKMT1B. 15q15.3 is a recurrent segmental-duplication/CNV region (STRC/STRCP1, CATSPER2/CATSPER2P1) with documented short-read mapping failure (PMID 35022556; PMID 40341240). MAPPABILITY CAVEAT REQUIRED."),
"CATSPER2": dict(cls="single-variant",
  prior_gwas="rs146163343, p=9e-6, proxy AD/ADRD in African-ancestry (PMID 40708016) -- suggestive, NOT genome-wide",
  prior_lit="none AD-specific (0 PubMed AD papers)",
  verdict="not previously reported in ADSP",
  note="Same 15q15.3 locus as CKMT1A (24,768 bp). Same segmental-duplication mappability caveat."),
"ZNF532": dict(cls="single-variant",
  prior_gwas="none AD/dementia (13 assoc., 0 AD)",
  prior_lit="none (0 PubMed AD papers)",
  verdict="novel",
  note="chr18 -- NOT among the chromosomes the draft lists for its novel signals (1,3,7,10,13,15,16,22)."),
"PIK3R6": dict(cls="single-variant",
  prior_gwas="rs60256690, p=1e-8, 'Alzheimer's disease or family history of Alzheimer's disease' (PMID 39998322, Willett 2025) -- GENOME-WIDE SIGNIFICANT",
  prior_lit="none additional (0 PubMed AD papers)",
  verdict="previously reported - reframe as replication",
  note="11,215 bp from PIK3R5. chr17 -- NOT among the chromosomes the draft lists."),
"SEMA4A": dict(cls="single-variant",
  prior_gwas="rs2842879, p=7e-7, CSF p-tau at 24 months (PMID 35086473), MAPPED intergenic SEMA4A-SLC25A44 -- suggestive",
  prior_lit="none AD-specific genetic",
  verdict="not previously reported in ADSP",
  note="Overlaps LMNA neighbourhood (7,281 bp)."),
"ADAMTS14": dict(cls="single-variant",
  prior_gwas="none AD/dementia (39 assoc., 0 AD)",
  prior_lit="none AD-specific",
  verdict="novel", note="3,613 bp from PALD1."),
"MYO1E": dict(cls="single-variant",
  prior_gwas="rs7180978, p=6e-6, proxy AD/ADRD (PMID 40708016) -- suggestive, NOT genome-wide",
  prior_lit="none AD-specific",
  verdict="not previously reported in ADSP",
  note="Overlaps LDHAL6B/FAM81A."),
"APOE": dict(cls="single-variant (positive control)",
  prior_gwas="4,763 assoc., 244 core-AD rows -- the canonical AD locus",
  prior_lit="vast", verdict="previously reported - reframe as replication",
  note="Top variant chr19:44,909,967 lies 571 bp 3' of APOE, in the APOE-APOC1 interval. TOMM40 ends 92 bp before APOE starts."),

"KIF18B": dict(cls="coding (live draft)",
  prior_gwas="none AD/dementia (22 assoc., 0 AD)",
  prior_lit="1 PubMed AD-context paper; the cited stout_kif18b_2011 (PMID 21820313 family) describes a MITOTIC astral-microtubule kinesin -- a weak mechanistic fit for post-mitotic neurons.",
  verdict="drop",
  note="NULL in the corrected analysis (STAAR-O 0.316-0.935). Prior project work attributes it to ultra-rare cMAC + ACAT-V double counting."),
"CASP4": dict(cls="coding (live draft)",
  prior_gwas="none AD/dementia (5 assoc., 0 AD)",
  prior_lit="STRONG prior AD functional evidence: human-specific CASP4 drives Alzheimer-related synaptic/behavioural deficits (PMID 27516385, Hum Mol Genet 2016); cell-type-specific pyroptosis in AD (PMID 36481964, Acta Neuropathol 2023); caspase-4 non-canonical inflammasome drives AD pathology (PMID 40497980, 2025).",
  verdict="drop",
  note="NULL in the corrected analysis (STAAR-O 0.543-0.991). Also never novel biologically -- the non-canonical inflammasome link is prior work, not this study's discovery."),
"USH2A": dict(cls="coding (live draft)",
  prior_gwas="rs111033333, p=3e-8, 'Alzheimer's disease (age of onset) in APOE e4 non-carriers', MAPPED to USH2A (PMID 33637690, He 2021, Transl Psychiatry) -- GENOME-WIDE SIGNIFICANT. NB: that paper's headline genes are ERN1 and SPPL2C; the USH2A row is a secondary result within it, so cite it as a Catalog entry, not as 'USH2A was reported as an AD gene'.",
  prior_lit="PMID 33794991 (age-informed AD association analysis)",
  verdict="drop",
  note="Not significant in the corrected analysis; and it was already reported, so it could never have been claimed as novel."),
"NEB": dict(cls="coding (live draft)",
  prior_gwas="none AD/dementia (70 assoc., 0 AD)",
  prior_lit="none AD-specific (neuromuscular literature only, PMID 22468856)",
  verdict="drop",
  note="Not significant in the corrected analysis. NEB spans 249 kb -- a very large mutational target and a known rare-variant-burden false-positive magnet."),
"PCDHB14": dict(cls="coding (live draft)",
  prior_gwas="none AD/dementia (4 assoc., 0 AD)",
  prior_lit="none AD-specific",
  verdict="drop",
  note="Not significant in the corrected analysis. Sits in the tandem protocadherin-beta array (PCDHB13 4,364 bp; PCDHB12 10,772 bp; PCDHB15 17,636 bp) -- high paralog mis-assignment risk."),
"RUNX1": dict(cls="coding (live draft)",
  prior_gwas="rs7281993, p=8e-6, 'Dementia' (PMID 39046104); rs2834655, p=3e-6, age-related cognitive decline (PMID 30954325) -- both suggestive",
  prior_lit="29 PubMed AD-context papers, none rare-variant",
  verdict="drop",
  note="Not significant in the corrected analysis. The cited bertran_identification_2022 is about RUNX1 in NON-ALCOHOLIC FATTY LIVER DISEASE via systems-biology modelling -- it does not support an interleukin/neuroimmune claim in AD."),

"APMAP": dict(cls="noncoding (live draft)",
  prior_gwas="none AD/dementia (28 assoc., 0 AD); several rows are CST7 pQTLs",
  prior_lit="ESTABLISHED AD modifier: APMAP is an endogenous suppressor of Abeta production (PMID 25180020, Hum Mol Genet 2015 -- the draft's own mosser_adipocyte_2015); APMAP interactome alters APP processing in AD brain (PMID 30704515).",
  verdict="previously reported - reframe as replication",
  note="Nearest gene is CST7 at 1,532 bp -- cystatin F, a canonical disease-associated-microglia marker and itself a pQTL in this region. An enhancer/DHS signal here cannot be assigned to APMAP rather than CST7 by proximity alone."),
"GSC": dict(cls="noncoding (live draft)",
  prior_gwas="none AD/dementia (19 assoc., 0 AD)",
  prior_lit="none AD-specific",
  verdict="novel",
  note="3.6 kb gene; 143 bp from GSC-DT. Symbol hazard: GSC2 (previous symbol GSCL) is a different gene on chr22."),
"HIST1H4F": dict(cls="noncoding (live draft)",
  prior_gwas="none AD/dementia. 23 Catalog rows name HIST1H4F, but every one is either a bulk histone-cluster/HLA REPORTED-gene list (PMID 30664634, 28540026, 32066663) or an intergenic MAPPED assignment (H1-3 - H4C6, H4C6 - H4C7). No locus-specific, no AD.",
  prior_lit="none AD-specific for H4C6 itself",
  verdict="novel",
  note="DEPRECATED SYMBOL: HIST1H4F -> H4C6 (HGNC:4783). NOT H4C5 (HGNC:4790, 35,372 bp away). 400 bp gene inside the histone cluster (H4C7 at 5,818 bp) AND inside the extended MHC (chr6:25.7-33.4 Mb) -- long-range LD and paralog ambiguity caveat required. The draft's H4K16 ACETYLATION claim is supported by nativio_dysregulation_2018 but NOT by fu_histone_2024, which is about histone METHYLATION."),
"ZC3H12A": dict(cls="noncoding (live draft)",
  prior_gwas="none AD/dementia (9 assoc., 0 AD)",
  prior_lit="none AD-specific genetic (4 AD-context papers)",
  verdict="novel", note="Regnase-1 immune biology is well established (PMID 19701189 / matsushita_zc3h12a_2009) but not AD-specific."),
"TNFRSF1B": dict(cls="noncoding (live draft)",
  prior_gwas="none AD/dementia (32 assoc., 0 AD)",
  prior_lit="DIRECT prior AD variant-level report: 'TNFRSF1B Gene Variants and Related Soluble TNFR2 Levels Impact Resilience in Alzheimer's Disease' (PMID 33716716, Front Aging Neurosci 2021). TNFR2 is an established AD therapeutic target (PMID 30778285 = the draft's own orti-casan_targeting_2019; PMID 40527437).",
  verdict="previously reported - reframe as replication",
  note="Catalog-clean but variant-level prior evidence exists outside the Catalog."),
"EBF2": dict(cls="noncoding (live draft)",
  prior_gwas="none core AD; 2 broad rows are PHF-tau SNPxSNP interaction terms MAPPED intergenic CDCA2-EBF2 (PMID 32450446)",
  prior_lit="none AD-specific",
  verdict="novel", note="683 Catalog associations overall -- a highly pleiotropic region; 246 kb to the nearest neighbour."),

"SNCAIP": dict(cls="SKATMC (PROVISIONALLY INVALID)",
  prior_gwas="ZERO AD/dementia rows (23 assoc., 0 AD)",
  prior_lit="No SNCAIP variant-dementia association could be verified. Existing SNCAIP genetics is Parkinson/Lewy-body and mostly expression or methylation (PMID 18366718, 17467279, 32568194).",
  verdict="drop",
  note="CITATIONS DO NOT SUPPORT THE CLAIM. The draft asserts 'recent literature establishes a robust correlation between SNCAIP variants and broader cognitive impairment and dementia' citing (a) khandelwal_relationship_2010 -- a review of PARKIN and protein aggregation that reports no SNCAIP variant; (b) dashtipour_hypermethylation_2017 -- SNCAIP DNA METHYLATION in Parkinson cortex, not variants, not AD; (c) fan_alpha-synuclein_2021 -- alpha-synuclein and cognition in PARKINSON disease, not SNCAIP. Different chromosome from SNCA (chr5 vs chr4)."),
"PRKCZ": dict(cls="SKATMC (PROVISIONALLY INVALID)",
  prior_gwas="none AD/dementia (51 assoc., 0 AD)", prior_lit="31 AD-context papers, none rare-variant/genetic",
  verdict="drop", note="Result provisionally invalid pending SKATMC indexing verification."),
"SIRT2": dict(cls="SKATMC (PROVISIONALLY INVALID)",
  prior_gwas="none AD/dementia (3 assoc., 0 AD -- all SIRT2 protein-level pQTLs)",
  prior_lit="ESTABLISHED AD target: PMID 37488853 (the draft's own sola-sevilla_sirt2_2024); microglial SIRT2 deficiency aggravates amyloid pathology (PMID 40499846); dual HDAC6/SIRT2 inhibition in AD (PMID 40381221).",
  verdict="drop", note="Result provisionally invalid; and the gene was never novel for AD."),
"CNTNAP1": dict(cls="SKATMC (PROVISIONALLY INVALID)",
  prior_gwas="none AD/dementia (10 assoc., 0 AD)",
  prior_lit="none AD-specific; CNTNAP1 literature is congenital hypomyelinating neuropathy (PMID 29882456)",
  verdict="drop",
  note="SYMBOL HAZARD: CNTNAP2 (chr7, different chromosome) HAS 10 core-AD Catalog rows including rs114360492 p=2e-9 (PMID 30617256, Jansen). Must be explicitly disambiguated. The cited altuna_dna_2019 is about hippocampal DNA methylation and neurogenesis -- it does not support the node-of-Ranvier claim."),
"FGG": dict(cls="SKATMC (PROVISIONALLY INVALID)",
  prior_gwas="none AD/dementia (91 assoc., 0 AD)",
  prior_lit="Fibrinogen-AD biology established (PMID 31056353 / merlini_fibrinogen_2019) but not FGG-variant-specific",
  verdict="drop",
  note="Tandem fibrinogen cluster: FGA 13,307 bp, FGB 31,251 bp. FGA is NEARER than FGB. A rare-variant set at this locus cannot be attributed to FGG specifically."),
"CST1": dict(cls="SKATMC (PROVISIONALLY INVALID)",
  prior_gwas="rs147525344, p=5e-9, 'Alzheimer disease and age of onset' (PMID 26830138, Herold 2016, Mol Psychiatry), MAPPED intergenic CST1-CSTP2 -- GENOME-WIDE SIGNIFICANT. NB: that paper's headline genes are OSBPL6, PTPRG and PDCL3; the CST1-CSTP2 row is a secondary result and the assignment is intergenic, so it is evidence that the region is not virgin territory rather than that CST1 is an established AD gene.",
  prior_lit="none additional",
  verdict="previously reported - reframe as replication",
  note="Cystatin cluster: CST4 58,522 bp, CST2 72,501 bp, CST3 109,004 bp. CST3 (cystatin C) is itself AD-associated (PMID 40708016). Cluster is unresolvable by proximity."),
"ZNF337": dict(cls="SKATMC (PROVISIONALLY INVALID)",
  prior_gwas="rs6050790, p=1e-6, 'Psychosis and affective disturbance in Alzheimer's disease' (PMID 38784964) -- suggestive",
  prior_lit="none (0 PubMed AD papers)",
  verdict="drop",
  note="Cannot be described as a 'novel' gene: an AD-phenotype Catalog entry already exists."),
"C19orf71": dict(cls="SKATMC (PROVISIONALLY INVALID)",
  prior_gwas="none AD/dementia (5 assoc., 0 AD)", prior_lit="none (0 PubMed AD papers)",
  verdict="drop", note="DEPRECATED SYMBOL: C19orf71 -> TEKTIP1 (HGNC:34496)."),
"QPRT": dict(cls="SKATMC (PROVISIONALLY INVALID)",
  prior_gwas="none AD/dementia (8 assoc., 0 AD)", prior_lit="16 AD-context papers (kynurenine pathway), none genetic",
  verdict="drop",
  note="QPRT lies inside the recurrent 16p11.2 BP4-BP5 CNV interval (with SPN, KIF22, MAZ, PRRT2, KCTD13, TAOK2, MAPK3). Overlaps SPN by 7,920 bp. Structural-variation confounding caveat required."),
"ANP32B": dict(cls="SKATMC (PROVISIONALLY INVALID)",
  prior_gwas="none AD/dementia (111 assoc., 0 AD)", prior_lit="none AD-specific genetic",
  verdict="drop", note="Result provisionally invalid."),
"SLC6A15": dict(cls="SKATMC (PROVISIONALLY INVALID)",
  prior_gwas="AD-related: AD polygenic-risk-score extremes rs6539809 p=4e-10 (PMID 35589863); cognitive-decline rate in late MCI rs189010983 p=5e-7 (PMID 26252872); plasma Abeta42/40 rs55851056 p=5e-6 (PMID 35623891). DEPRESSION PLEIOTROPY CONFIRMED: major depressive disorder rs2125716 p=3e-8, REPORTED GENE(S)=SLC6A15 (PMID 27479909).",
  prior_lit="3 AD-context papers",
  verdict="previously reported - reframe as replication",
  note="DEPRESSION-PLEIOTROPY CAVEAT REQUIRED. Depression is both an AD risk factor and an AD prodrome, so an SLC6A15 signal against a graded cognitive outcome is confounded by mood phenotype. All AD-related Catalog rows are MAPPED to the intergenic RPL6P25-SLC6A15 interval, not to SLC6A15 itself."),

"FRMD4A": dict(cls="single-variant (current draft)",
  prior_gwas="ESTABLISHED AD LOCUS: rs7081208 / rs2446581 / rs17314229, p=1e-10, 'Alzheimer's disease', MAPPED to FRMD4A (PMID 22430674, Lambert 2012, Mol Psychiatry). 142 associations overall.",
  prior_lit="Replicated AD susceptibility gene since 2012.",
  verdict="previously reported - reframe as replication",
  note="POSITIVE FINDING, MISLABELLED. FRMD4A currently appears inside a list under a heading that says the scan 'identifies novel loci'. It is not novel -- it is a replicated AD locus, and recovering it is evidence the matched scan retains real signal. Promote it, do not bury it."),
"TUBA3C": dict(cls="single-variant (current draft)",
  prior_gwas="rs1306405136, p=4e-8, 'Alzheimer's disease or family history of Alzheimer's disease', MAPPED TUBA3C - PSPC1P1 (PMID 39998322, Willett 2025) -- GENOME-WIDE SIGNIFICANT",
  prior_lit="none AD-specific",
  verdict="previously reported - reframe as replication",
  note="Reported by the SAME recent multi-ancestry AD meta-analysis (Willett 2025) that reported PIK3R6. Two of the draft's named single-variant genes come from that one paper."),
"TAFA5": dict(cls="single-variant (current draft, MISLABELLED as Z84468.1)",
  prior_gwas="SEVEN AD-relevant rows: late-onset AD rs1034435 p=3e-6 (PMID 27770636); neurofibrillary tangles rs75814742 p=6e-6 (PMID 31497858); age at dementia onset in PSEN1 E280A carriers rs11705431 p=3e-6 (PMID 36951251); memory performance in cognitive impairment rs3810640/rs3827414 p=7e-6 (PMID 37985223); total PHF-tau interaction (PMID 32450446).",
  prior_lit="Brain-expressed TAFA/FAM19A chemokine-like family.",
  verdict="previously reported - reframe as replication",
  note="RELABEL REQUIRED. The chr22:48,539,023 lead variant is printed as 'Z84468.1', a clone-based Ensembl identifier with no HGNC record. The variant lies INSIDE TAFA5 (protein-coding). Relabelling converts an uninterpretable clone ID into a replication of prior AD and neurofibrillary-tangle work."),
"FAM234B": dict(cls="single-variant (current draft)",
  prior_gwas="none AD/dementia (90 assoc., 0 AD)", prior_lit="none AD-specific",
  verdict="novel", note="chr12:13,044,185-13,142,521."),
"GRIA2": dict(cls="single-variant (current draft)",
  prior_gwas="none AD/dementia (6 assoc., 0 AD)", prior_lit="AMPA-receptor biology well established but no prior AD genetic association",
  verdict="novel", note="chr4:157,204,182-157,387,146."),
"MAPK4": dict(cls="single-variant (current draft)",
  prior_gwas="none AD/dementia (97 assoc., 0 AD; blood-pressure dominated)", prior_lit="none AD-specific",
  verdict="novel", note="chr18:50,559,451-50,731,830."),
"MARCH8": dict(cls="single-variant (current draft)",
  prior_gwas="none AD/dementia (432 assoc., 0 AD)", prior_lit="none AD-specific",
  verdict="novel",
  note="DEPRECATED SYMBOL: MARCH8 -> MARCHF8 (HGNC:23356). The chr10:45,594,673 lead variant lies inside MARCHF8. All MARCH-family symbols were renamed to MARCHF in 2020 to stop spreadsheet software coercing them to dates."),
"VAMP1": dict(cls="single-variant (current draft)",
  prior_gwas="none AD/dementia (9 assoc., 0 AD)", prior_lit="none AD-specific",
  verdict="novel", note="chr12:6,462,237-6,470,948."),
"TSC22D1-AS1": dict(cls="single-variant (current draft)",
  prior_gwas="none AD/dementia (22 assoc., 0 AD)", prior_lit="none",
  verdict="novel", note="lncRNA, chr13:44,575,766-44,755,728."),

"PSEN1": dict(cls="coding (CORRECTED result)",
  prior_gwas="ZERO Catalog AD rows (26 assoc., 0 AD) -- expected: PSEN1 is a Mendelian gene, invisible to common-variant GWAS",
  prior_lit="CANONICAL early-onset familial AD gene (PMID 7596406, Sherrington, Nature 1995; PMID 7651536). 4,240 AD papers, 2,443 with rare-variant/mutation terms.",
  verdict="previously reported - reframe as replication",
  note="POSITIVE CONTROL, and a strong one: the absence of PSEN1 from the common-variant Catalog while rare-variant aggregation recovers it at p=1.7e-8 is exactly the behaviour a well-calibrated rare-variant framework should show. Do not present as novel."),
"TOMM40": dict(cls="noncoding (CORRECTED result)",
  prior_gwas="1,787 assoc., 158 core-AD rows -- one of the most-reported AD loci",
  prior_lit="TOMM40's AD signal is long argued to reflect LD with APOE: PMID 22869030 ('TOMM40 association with Alzheimer disease: tales of APOE and linkage disequilibrium'); PMID 34625307 ('Definitive roles of TOMM40-APOE-APOC1 variants').",
  verdict="previously reported - reframe as replication",
  note="TOMM40 ends 92 bp before APOE begins. The chr19 enhancer/DHS result is NOT independent of APOE and must be reported conditional on APOE genotype or explicitly labelled as the APOE locus."),
"ELMSAN1": dict(cls="noncoding (CORRECTED result)",
  prior_gwas="none AD/dementia (146 assoc., 0 AD)",
  prior_lit="none AD-specific (1 PubMed paper)",
  verdict="novel",
  note="DEPRECATED SYMBOL: ELMSAN1 -> MIDEAS (HGNC:19853). Two caveats: (i) MIR4505 lies ENTIRELY INSIDE MIDEAS, so the ncRNA and enhancer/DHS hits are ONE locus, not two; (ii) MIDEAS starts 491,431 bp from the end of PSEN1 on chr14 -- independence from the PSEN1 signal must be demonstrated, not assumed."),
"MIR4505": dict(cls="ncRNA (CORRECTED result)",
  prior_gwas="no Catalog rows",
  prior_lit="1 paper, Lewy-body/RBD serum miRNA panel (PMID 35962561) -- not AD-specific",
  verdict="novel",
  note="NOT AN INDEPENDENT LOCUS: MIR4505 (chr14:73,758,747-73,758,819) is contained within MIDEAS (73,715,122-73,790,285), 43,625 bp inside its 5' end. Reporting MIDEAS and MIR4505 as separate findings double-counts one signal."),
}

ORDER = ["CACNA1D","DCP1A","PRKAR1B","CKMT1A","CATSPER2","ZNF532","PIK3R6","SEMA4A",
         "ADAMTS14","MYO1E","APOE",
         "FRMD4A","TUBA3C","TAFA5","FAM234B","GRIA2","MAPK4","MARCH8","VAMP1","TSC22D1-AS1",
         "KIF18B","CASP4","USH2A","NEB","PCDHB14","RUNX1",
         "APMAP","GSC","HIST1H4F","ZC3H12A","TNFRSF1B","EBF2",
         "SNCAIP","PRKCZ","SIRT2","CNTNAP1","FGG","CST1","ZNF337","C19orf71","QPRT",
         "ANP32B","SLC6A15",
         "PSEN1","TOMM40","ELMSAN1","MIR4505"]


def _one(c):
    """Distances use the difference convention: gap = next_start - prev_end;
    overlap = min(end) - max(start).  Matches the numbers already established."""
    d = c["dist"]
    if d < 0:
        return f"{c['sym']} (OVERLAPS by {abs(d) - 1:,} bp)"
    return f"{c['sym']} ({d:,} bp)"


def fmt_neighbour(g):
    n = NEI.get(g)
    if not n:
        return "n/a"
    picks = n["confusable"][:2] if n["confusable"] else n["nearest_any"][:1]
    # always surface a same-family paralog if one exists within the window
    fam = [c for c in n["confusable"] if c.get("same_stem")]
    for c in fam:
        if c not in picks:
            picks = picks[:1] + [c]
            break
    return "; ".join(_one(c) for c in picks)


rows = []
for g in ORDER:
    r, c = RES[g], C[g]
    cur = r.get("hgnc_symbol")
    coord = f"chr{r['chr']}:{r['start']:,}-{r['end']:,}" if r.get("chr") else "n/a"
    rows.append(dict(
        printed=g, current=cur, deprecated=("YES" if cur != g else ""),
        hgnc=r.get("hgnc_id", ""), ens=r.get("ens_id", ""),
        coord=coord, cls=c["cls"], neighbour=fmt_neighbour(g),
        prior_gwas=c["prior_gwas"], prior_lit=c["prior_lit"],
        verdict=c["verdict"], note=c["note"]))

# --------------------------------------------------------------------- TSV
cols = ["printed_symbol", "current_hgnc_symbol", "symbol_deprecated", "hgnc_id",
        "ensembl_id", "grch38_coordinates", "result_class",
        "nearest_confusable_neighbour_and_distance",
        "prior_AD_dementia_gwas_catalog", "prior_rare_variant_or_functional_AD_literature",
        "verdict", "notes"]
with open(f"{BASE}/tables/novelty_audit.tsv", "w") as f:
    f.write("\t".join(cols) + "\n")
    for r in rows:
        f.write("\t".join(str(r[k]).replace("\t", " ").replace("\n", " ") for k in
                ["printed", "current", "deprecated", "hgnc", "ens", "coord", "cls",
                 "neighbour", "prior_gwas", "prior_lit", "verdict", "note"]) + "\n")

print("wrote novelty_audit.tsv:", len(rows), "rows")
json.dump(rows, open(f"{BASE}/logs/novelty_rows.json", "w"), indent=1)
