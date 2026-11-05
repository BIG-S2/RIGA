#!/usr/bin/env python3
"""supplementary_tables_specifications.py (2026-10-04) -- render the CSV tables in revision/tables/supp_v3/
as LaTeX longtable fragments for revision/manuscript/supplementary_information_v3.tex.
Usage: python3 supplementary_tables_specifications.py [--csvdir DIR] [--outdir DIR]
Every number is read from the CSV; nothing is typed in here."""
import argparse, csv, math, pathlib, re

def esc(s):
    s = "" if s is None else str(s)
    if s in ("nan", "NA", "None"): return ""
    for a, b in [("\\", r"\textbackslash{}"), ("&", r"\&"), ("%", r"\%"), ("#", r"\#"), ("_", r"\_"), ("$", r"\$"), ("~", r"\textasciitilde{}"), ("^", r"\^{}")]:
        s = s.replace(a, b)
    return s

def is_na(v):
    return v is None or str(v).strip() in ("", "nan", "NA", "None")

def f_int(v):
    if is_na(v): return ""
    x = float(v); return f"{int(round(x)):,}"

def f_num(d):
    def f(v):
        if is_na(v): return ""
        return f"{float(v):,.{d}f}"
    return f

def f_sci(v, sig=2):
    if is_na(v): return ""
    x = float(v)
    if x == 0: return "0"
    if 0.001 <= x < 1: return f"{x:.{max(sig,2)}g}" if x >= 0.01 else f"{x:.{sig}g}"
    e = int(math.floor(math.log10(abs(x)))); m = x / 10 ** e
    m = round(m, sig - 1)
    if m >= 10: m /= 10; e += 1
    return f"${m:.{sig-1}f}\\times10^{{{e}}}$"

def f_str(v): return esc(v)
def f_test(v):
    """Expand the display label without changing the source test identifier."""
    return "SKAT--burden omnibus test" if v == "SB (SKAT + burden omnibus)" else esc(v)
def f_variant(v):
    """chr:pos:ref:alt with break points after each colon and inside long alleles"""
    if is_na(v): return ""
    parts = str(v).split(":")
    out = []
    for i, x in enumerate(parts):
        if i >= 2 and len(x) > 6:
            out.append("\\allowbreak{}".join(esc(x[j:j+6]) for j in range(0, len(x), 6)))
        else:
            out.append(esc(x))
    return ":\\allowbreak{}".join(out)
def f_gene(v): return "" if is_na(v) else r"\textit{" + esc(v) + "}"
def f_locus(v):   # locus label: italicize gene symbols (all-caps tokens) and leave the descriptive words upright
    return "" if is_na(v) else " ".join(r"\textit{" + esc(w) + "}" if re.fullmatch(r"[A-Z][A-Z0-9-]+", w) else esc(w) for w in str(v).split())
def f_family(v): return esc(str(v).replace("Gene-centric ", ""))
def f_source(v): return "FDR 5\%" if str(v).startswith("FDR") else "text (n.s.)"
def f_tier(v): return "FDR" if str(v).startswith("FDR") else "suggestive"
def f_verdict(v):
    v = str(v)
    return {"APOE region": "APOE", "variant absent from FinnGen": "absent"}.get(v, "replicated" if v.startswith("replicated") else ("nominal" if v.startswith("nominal") else esc(v)))
def f_novelty(v):
    v = str(v)
    return "APOE" if v == "APOE region" else ("known AD locus" if v.startswith("AD association") else ("none" if v.startswith("no AD") else esc(v)))
def f_bool(v): return "" if is_na(v) else ("yes" if str(v).lower() in ("true", "1", "yes") else "no")

def read_csv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh))

def longtable(rows, cols, caption, label, font=r"\footnotesize", colspec=None, landscape=False, note=None):
    """cols: list of (header, formatter) where formatter is a column name + fmt function, or a callable on the row."""
    heads = [c[0] for c in cols]
    if colspec is None: colspec = "l" + "r" * (len(cols) - 1)
    out = []
    if landscape: out.append(r"\begin{landscape}")
    out.append(font)
    out.append(r"\begin{longtable}{" + colspec + "}")
    out.append(r"\multicolumn{" + str(len(cols)) + r"}{@{}p{\linewidth}@{}}{" + caption + r"}\label{" + label + r"}\\[4pt]")
    out.append(r"\toprule " + " & ".join(r"\textbf{" + h + "}" for h in heads) + r"\\ \midrule")
    out.append(r"\endfirsthead")
    out.append(r"\multicolumn{" + str(len(cols)) + r"}{l}{\textit{(continued)}}\\ \toprule " + " & ".join(r"\textbf{" + h + "}" for h in heads) + r"\\ \midrule")
    out.append(r"\endhead")
    out.append(r"\midrule \multicolumn{" + str(len(cols)) + r"}{r}{\textit{continued on next page}}\\")
    out.append(r"\endfoot")
    out.append(r"\bottomrule" + (r" \multicolumn{" + str(len(cols)) + r"}{p{0.98\linewidth}}{\footnotesize " + note + r"}\\" if note else ""))
    out.append(r"\endlastfoot")
    for r in rows:
        cells = []
        for h, f in cols:
            cells.append(f(r) if callable(f) and not isinstance(f, tuple) else f[1](r.get(f[0])))
        out.append(" & ".join(cells) + r"\\")
    out.append(r"\end{longtable}")
    out.append(r"\normalsize")
    if landscape: out.append(r"\end{landscape}")
    return "\n".join(out) + "\n"

def col(name, fmt, header=None):
    return (header if header is not None else esc(name), (name, fmt))

def or_ci(r):
    if is_na(r.get("OR")): return ""
    return f"{float(r['OR']):.2f} ({float(r['OR 95% CI lower']):.2f}--{float(r['OR 95% CI upper']):.2f})"

def nearest(r):
    g = r.get("Nearest protein-coding gene")
    if is_na(g): return ""
    d = r.get("Distance to nearest gene (bp)")
    return f_gene(g) + ("" if is_na(d) or float(d) == 0 else f" ({f_int(d)} bp)")

def ukb_verdict(r):
    disc = {("TREM2", "missense"), ("NECTIN2", "missense")}
    if (r.get("Gene name"), r.get("Category")) not in disc or r.get("arm") != "noRIGA": return ""
    return "replicated ($P < 0.05/2$)" if float(r["STAAR-B"]) < 0.025 else "not replicated"

def carriers(which):
    def f(r):
        return f"{f_int(r[f'Carriers ({which})'])} / {f_int(r[f'N {which}'])}"
    return f

SPECS = {
 "S1a": dict(csv="S1a_analytic_dataset_flow.csv", colspec="P{0.68\\linewidth}rr", font=r"\small",
    cols=[col("Step", f_str), col("N", f_int, "N"), col("N removed at this step", f_int, "Removed at this step")],
    caption=r"\textbf{Supplementary Table 1a. From the ADSP Release 5 whole-genome sequencing release to the analytic set of the stratified design.} The counts removed at each step; no latent-class confidence filter and no matching step is applied. Source: \texttt{S1a\_analytic\_dataset\_flow.csv}.", label="tab:s1a"),
 "S1b": dict(csv="S1b_analytic_race.csv", font=r"\small", colspec="P{0.08\\linewidth}P{0.36\\linewidth}rrrr",
    cols=[col("Race code", f_str, "Code"), col("Self-reported group", f_str), col("N", f_int), col("Cases", f_int), col("Controls", f_int), col("% of analytic set", f_num(2), r"\%")],
    caption=r"\textbf{Supplementary Table 1b. Self-reported race/ethnicity of the analytic set} (30{,}894 participants).", label="tab:s1b"),
 "S1c": dict(csv="S1c_analytic_seqcenter.csv", font=r"\small", colspec="P{0.44\\linewidth}rrrrr",
    cols=[col("Sequencing center", f_str), col("N", f_int), col("Cases", f_int), col("Controls", f_int), col("Cells", f_int, "Cells"), col("% of analytic set", f_num(2), r"\%")],
    caption=r"\textbf{Supplementary Table 1c. Sequencing centers of the analytic set.} Two center labels that differ only in capitalization are merged; centers with fewer than 100 analytic participants are pooled as Other.", label="tab:s1c"),
 "S1d": dict(csv="S1d_analytic_projects.csv", font=r"\scriptsize", colspec="P{0.26\\linewidth}rrrrrr",
    cols=[col("Sequencing project", f_str), col("N in pool", f_int, "Pool N"), col("Cases in pool", f_int, "Pool cases"), col("N retained", f_int, "Kept N"), col("Cases retained", f_int, "Kept cases"), col("Controls retained", f_int, "Kept controls"), col("Retention", f_num(3))],
    caption=r"\textbf{Supplementary Table 1d. Sequencing projects: the phenotyped pool and the analytic set.} Projects with a retention of 0 contain only cases or only controls in every cell and carry no within-cell contrast (Proposition 1b).", label="tab:s1d"),
 "S1e": dict(csv="S1e_analytic_cohorts.csv", font=r"\scriptsize", colspec="P{0.34\\linewidth}rrrrrr",
    cols=[col("Cohort", f_str), col("N in pool", f_int, "Pool N"), col("Cases in pool", f_int, "Pool cases"), col("N retained", f_int, "Kept N"), col("Cases retained", f_int, "Kept cases"), col("Controls retained", f_int, "Kept controls"), col("Retention", f_num(3))],
    caption=r"\textbf{Supplementary Table 1e. Contributing cohorts: the phenotyped pool and the analytic set.}", label="tab:s1e"),
 "S1f": dict(csv="S1f_project_center_cells.csv", font=r"\scriptsize", colspec="P{0.26\\linewidth}P{0.22\\linewidth}rrrl",
    cols=[col("Sequencing project", f_str), col("Sequencing center", f_str), col("N", f_int), col("Cases", f_int), col("Controls", f_int), col("Retained (both labels present)", f_bool, "Retained")],
    caption=r"\textbf{Supplementary Table 1f. The 79 sequencing project $\times$ center cells of the phenotyped pool.} The 64 cells that contain both cases and controls form the strata of the design; the 15 single-label cells are set aside.", label="tab:s1f"),
 "S1g": dict(csv="S1g_analytic_latent_class.csv", font=r"\small",
    cols=[col("Latent class", f_str), col("N", f_int), col("Cases", f_int), col("Controls", f_int)],
    caption=r"\textbf{Supplementary Table 1g. Latent genetic classes of the analytic set} (13-class Gaussian mixture on PC1--PC20; modal class).", label="tab:s1g"),
 "S2": dict(csv="S2_fdr_single_variants.csv", font=r"\scriptsize", colspec="P{0.17\\linewidth}P{0.13\\linewidth}rP{0.12\\linewidth}rrP{0.13\\linewidth}P{0.18\\linewidth}", landscape=True,
    cols=[col("Variant (chr:pos:ref:alt)", f_variant, "Variant (GRCh38)"), col("Locus label (Figure 4a)", f_str, "Locus label"), col("MAF", f_num(4)),
          ("OR (95\% CI)", or_ci), col("P", f_sci, "$P$"), col("BH q value", f_sci, "$q$"), ("Nearest protein-coding gene", nearest), col("Variant-to-gene relation (FAVOR)", f_str, "Relation (FAVOR)")],
    caption=r"\textbf{Supplementary Table 2. The 296 single variants admitted at FDR 5\% in the stratified design.} Score test with the normal approximation on the 37{,}900{,}595 variants with minor allele count $\geq 50$ that carry \texttt{QC\_label7}, scanned on that label directly; odds ratios from the score-based estimate; $q$ = Benjamini--Hochberg adjusted $P$; locus labels as in Figure 4a; nearest protein-coding gene (distance when the variant is not inside it) and the FAVOR variant-to-gene relation. Full columns in \texttt{S2\_fdr\_single\_variants.csv}.", label="tab:s2"),
 "S3": dict(csv="S3_fdr_gene_mask_sets.csv", font=r"\tiny", landscape=True,
    cols=[col("Analysis", f_family, "Family"), col("Gene name", f_gene, "Gene"), col("Chr", f_str), col("Category", f_str, "Mask"), col("#SNV", f_int, r"\#SNV"), col("cMAC", f_int),
          col("Burden(1,25)", f_sci), col("SKAT(1,25)", f_sci), col("ACAT-V(1,25)", f_sci), col("STAAR-B(1,25)", f_sci, "STAAR-B"), col("STAAR-S(1,25)", f_sci, "STAAR-S"), col("STAAR-A(1,25)", f_sci, "STAAR-A"), col("ACAT-O", f_sci), col("STAAR-O", f_sci), col("SB", f_sci, r"\shortstack{SKAT--burden\\omnibus test\\(primary)}"), col("BH q value", f_sci, r"\shortstack{$q$ (SKAT--burden\\omnibus test)}")],
    caption=r"\textbf{Supplementary Table 3. Every gene $\times$ mask set admitted at FDR 5\% within its family} under the primary SKAT--burden omnibus test (coding: 69{,}555 sets with cMAC $\geq10$; noncoding: 101{,}641; ncRNA: 7{,}704). Component $P$ values of STAAR with Beta(1,25) weights, the annotation-weighted STAAR-B/S/A components, ACAT-O, the STAAR-O omnibus and the SKAT--burden omnibus test, the Cauchy combination of all SKAT and burden $P$ values (STAAR-O without ACAT-V); $q$ = Benjamini--Hochberg adjusted SKAT--burden omnibus test $P$ values within the family. The \textit{NECTIN2} missense set is carried by a single allele near \textit{APOE} (Results). All components in \texttt{S3\_fdr\_gene\_mask\_sets.csv}.", label="tab:s3"),
 "S4": dict(csv="S4_carrier_counts.csv", font=r"\tiny", landscape=True,
    cols=[col("Gene name", f_gene, "Gene"), col("Chr", f_str), col("Category", f_str, "Mask"), col("Source", f_source, "Set"), col("#SNV in set", f_int, r"\#SNV"), col("cMAC", f_int),
          ("Carriers / N, cases", carriers("cases")), ("Carriers / N, controls", carriers("controls")),
          col("Carrier freq. (cases)", f_num(4), "Freq. cases"), col("Carrier freq. (controls)", f_num(4), "Freq. controls"), col("Carrier OR (unadjusted)", f_num(2), "Carrier OR"), col("Fisher P (carrier, unadjusted)", f_sci, "Fisher $P$")],
    caption=r"\textbf{Supplementary Table 4. Rare-variant carrier counts in cases and controls of the analytic set} for every gene $\times$ mask set of Supplementary Table 3 and for the five coding masks of \textit{PSEN1}. Variants selected exactly as in the deployed scan (\texttt{QC\_label7}, SNVs, MAF $<1\%$ in the analytic set); a carrier holds at least one minor allele of the set. The carrier odds ratio and Fisher $P$ are unadjusted descriptive quantities, not the test statistics of the scan. Set: admitted at FDR 5\% (Supplementary Table~3) or discussed in the text without reaching significance. Minor allele counts by case status and per-variant counts are in \texttt{S4\_carrier\_counts.csv} and \texttt{S4b\_variant\_level\_carriers.csv}.", label="tab:s4"),
 "S8": dict(csv="S8_replication_finngen.csv", font=r"\tiny", landscape=True, colspec="lrrP{0.13\\linewidth}rrrrrlrP{0.12\\linewidth}rP{0.15\\linewidth}",
    cols=[col("Tier", f_tier), col("CHR", f_str), col("POS", f_int), col("Locus label", f_str), col("MAF", f_num(4)), col("ADSP P", f_sci, "ADSP $P$"), col("ADSP OR (ALT)", f_num(2), "ADSP OR"),
          col("FinnGen P (exact variant)", f_sci, "FinnGen $P$"), col("FinnGen beta (ALT, aligned)", f_num(3), r"FinnGen $\beta$"), col("Direction concordant", f_bool, "Same dir."),
          col("FinnGen best P in window", f_sci, "Best $P$ $\\pm$250 kb"), col("Verdict", f_verdict), col("Catalog AD associations in window", f_int, "Catalog AD"), col("Novelty", f_novelty)],
    caption=r"\textbf{Supplementary Table 8. FinnGen replication and GWAS Catalog audit of the single-variant leads of the stratified scan.} The 213 leads with $P < 10^{-5}$ pruned at 500 kb (FDR-5\% tier and suggestive tier), the exact-variant FinnGen release-12 G6\_ALZHEIMER result with the effect aligned to the ADSP alternative allele, the best FinnGen $P$ within 250 kb, the pre-specified verdict (replicated: same direction and $P < 0.05/11$ over the lead variants of the 11 non-\textit{APOE} discovery loci (the suggestive leads carry no verdict); nominal: same direction and $P < 0.05$; absent: variant not in the FinnGen panel) and the number of genome-wide-significant Alzheimer's disease or dementia associations within 250 kb in the GWAS Catalog, with the novelty call. Density-matched control loci carry such an association 5.5\% of the time. Summary statistics in \texttt{S8b\_replication\_summary.csv}.", label="tab:s8"),
 "S9": dict(csv="S9_catalog_genes.csv", font=r"\scriptsize", colspec="lP{0.17\\linewidth}rrrP{0.3\\linewidth}",
    cols=[col("Gene", f_gene), col("Masks admitted at FDR 5%", f_str, r"Masks at FDR 5\%"), col("Catalog associations (any trait)", f_int, "Catalog assoc. (any)"), col("Strict AD/dementia associations", f_int, "AD/dementia, genome-wide"),
          col("Best AD P", f_sci, "Best AD $P$"), col("First AD report", f_str, "Strongest AD report")],
    caption=r"\textbf{Supplementary Table 9. GWAS Catalog audit of the genes with a gene $\times$ mask set admitted at FDR 5\%, and \textit{PSEN1}.} Associations whose mapped or reported gene is the symbol, with the strict Alzheimer's disease or dementia filter of Supplementary Table 8 (genome-wide significant, main-effect disease traits). Genes without an entry (\textit{MMP17}, \textit{PGAM5}, \textit{MIPOL1}, \textit{RN7SL293P}) have no reported AD association and are presented as candidates requiring replication; \textit{PSEN1} is a Mendelian gene without a GWAS Catalog AD entry.", label="tab:s9"),
 "S10": dict(csv="S10_ukb_coding_lookup_rootcheck_snx22.csv", font=r"\scriptsize",
    cols=[col("Gene name", f_gene, "Gene"), col("Category", f_str, "Mask"), col("#SNV", f_int, r"\#SNV"), col("cMAC", f_int), col("STAAR-B", f_sci, "Burden $P$ (saddlepoint)"), ("Discovery set", ukb_verdict)],
    caption=r"\textbf{Supplementary Table 10. UK Biobank whole-genome coding results for the coding genes named in this study.} The 200k whole-genome release with the corrected phenotype (182{,}238 participants, 1{,}753 Alzheimer's disease cases, phecode 290.11); saddlepoint burden test on the five coding masks with the annotations of the ADSP scan, from the genome-wide UK Biobank coding scan that also provides the synonymous calibration reference of the Methods. The three coding sets admitted at FDR 5\% in ADSP are judged in Supplementary Table~11; the genes here are shown for reference. Saddlepoint $P$ values come from the root-verified solver described in the Methods. The final row gives the \textit{SNX22} synonymous set, included because its $P$ value before root verification was $2.8\times10^{-16}$; with a verified root it is 0.052.", label="tab:s10"),
 "S11": dict(csv="S11_ukb_geneset_replication_rootcheck.csv", font=r"\scriptsize",
    cols=[col("Family", f_str), col("Gene", f_gene), col("Mask", f_str), col("Chr", f_str), col("#SNV", f_int, r"\#SNV"), col("cMAC", f_int), col("ADSP SB P", f_sci, r"\shortstack{ADSP SKAT--burden\\omnibus test $P$}"), col("UKB burden P (SPA)", f_sci, "UKB burden $P$"), col("Verdict", f_str)],
    caption=r"\textbf{Supplementary Table 11. UK Biobank whole-genome replication of the 14 gene sets admitted at FDR 5\% in ADSP.} The 200k whole-genome release with the corrected phenotype (182{,}238 participants, 1{,}753 Alzheimer's disease cases, phecode 290.11), the same masks, annotations and MAF cutoff as the ADSP scan. Because the cohort is 1:103 case--control, the replication statistic is the saddlepoint burden test; STAAR provides no saddlepoint SKAT, so sets that ADSP admits through their SKAT component face a weaker test here. Saddlepoint $P$ values come from the root-verified solver described in the Methods. Replication requires $P < 0.05/14$ over the 14 sets.", label="tab:s11"),
 "S12": dict(csv="S12_lambda_by_ancestry.csv", font=r"\small",
    cols=[col("Group", f_str), col("N", f_int), col("Bin", f_str), col("Tests", f_int), col("lambda_GC", f_num(3), r"$\lambda_{GC}$"), col("n < 5e-8 outside APOE", f_int, r"$n < 5\times10^{-8}$ outside \textit{APOE}"), col("n < 1e-5 outside APOE", f_int, r"$n < 10^{-5}$ outside \textit{APOE}")],
    caption=r"\textbf{Supplementary Table 12. Calibration of the stratified single-variant scan within each self-reported ancestry group.} The design refitted inside each group (project $\times$ center cells formed within the group, mixed cells only; age, sex, cell and PC1--PC10 with the sparse GRM; Methods), the same scan and the same \texttt{QC\_label7} filter with minor allele count $\geq 50$ within the group; genomic-control $\lambda_{GC}$ over all tests and by minor-allele-frequency bin, with the full stratified scan for reference. The Asian group has no variant with MAC $\geq 50$ below 1\% frequency.", label="tab:s12"),
 "S12b": dict(csv="S12b_race_discovery_variants.csv", font=r"\scriptsize", landscape=True,
    cols=[col("Locus", f_locus), col("key", f_variant, "Variant (chr:pos:ref:alt)"), col("Group", f_str), col("N", f_int), col("MAF", f_num(4)), col("Est", f_num(3), r"$\hat\beta$"), col("SE", f_num(3)), col("P", f_sci, "$P$"), col("Est_full", f_num(3), r"Full-scan $\hat\beta$"), col("P_full", f_sci, "Full-scan $P$")],
    caption=r"\textbf{Supplementary Table 12b. The 21 non-\textit{APOE} single variants (11 loci) of the full stratified scan inside each self-reported ancestry group.} The \textit{APOE} lead and the 21 non-\textit{APOE} variants admitted at FDR 5\%, tested within each group by the per-ancestry refit of Supplementary Table 12 (minor allele count $\geq 50$ within the group; blank where the variant does not reach that count). Five rare variants (near \textit{MMP17}, \textit{KLHL1}, \textit{CRTAC1} and \textit{BAIAP2L1}) reach the count in no group, \textit{TUBB8} only in the White group and \textit{KCNV1} only in the Black or African American group. Inverse-variance fixed-effect meta-analysis of the group estimates: of the 14 variants tested in two or more groups, 12 have the same direction in every group tested, with Cochran's $Q$ $P \geq 0.27$ except chr2:127{,}124{,}439 ($Q$ $P = 0.043$), and meta-analytic estimates close to the full-scan estimates (\textit{CELF4} 0.243 against 0.249, \textit{MS4A4A} $-0.116$ against $-0.108$). \textit{SHARPIN} has an opposite, non-significant estimate in the Black or African American group ($Q$ $P = 0.11$). The \textit{BIN1} lead chr2:127{,}133{,}851 is heterogeneous ($Q$ $P = 1.0\times10^{-4}$): $\hat\beta = 0.127$, 0.078 and 0.183 in the White, Hispanic/Latino and Asian groups and $-0.108$ ($P = 0.023$) in the Black or African American group, whereas chr2:127{,}135{,}234, 1.4 kb away, is concordant in all four ($Q$ $P = 0.49$). \textit{APOE}: 0.859 ($Q$ $P = 3.1\times10^{-6}$, the Hispanic/Latino estimate being smaller).", label="tab:s12b"),
 "S12c": dict(csv="S12c_race_gw_hits.csv", font=r"\scriptsize", landscape=True,
    cols=[col("Group", f_str), col("key", f_variant, "Variant (chr:pos:ref:alt)"), col("N", f_int), col("MAF", f_num(4)), col("Est", f_num(3), r"$\hat\beta$"), col("SE", f_num(3)), col("P", f_sci, "$P$"), col("MAF_full", f_num(4), "Full-scan MAF"), col("Est_full", f_num(3), r"Full-scan $\hat\beta$"), col("P_full", f_sci, "Full-scan $P$")],
    caption=r"\textbf{Supplementary Table 12c. Variants below $5\times10^{-8}$ outside \textit{APOE} in the per-ancestry refits.} Seven variants, each significant in one group only (one White, one Hispanic/Latino, three Black or African American, two Asian) and none below the FDR threshold of the full scan (full-sample $P \geq 6.5\times10^{-7}$). Two are indels in a homopolymer or GC-rich context, two lie in pericentromeric sequence and two, 2.6 kb apart at 11p13, have near-identical estimates; one, in the Asian group, is common (within-group minor allele frequency 0.41). These are the classes that the call-artifact audit of this study flags; they are reported for completeness and carry no claim.", label="tab:s12c"),
 "S13": dict(csv="S13_ukb_sv_direct_rootcheck.csv", font=r"\scriptsize", landscape=True,
    cols=[col("Gene", f_gene), col("Variant", f_variant, "Variant (chr:pos:ref:alt)"), col("Locus lead", f_str, "Lead"), col("ADSP MAF", f_num(4)), col("ADSP OR", f_num(2)), col("ADSP P", f_sci, "ADSP $P$"), col("UKB MAF", f_sci, "UKB MAF"), col("UKB N", f_int), col("UKB P (SPA)", f_sci, "UKB $P$"), col("Verdict", f_str)],
    caption=r"\textbf{Supplementary Table 13. Direct test of the non-\textit{APOE} single variants of the ADSP replication set in the UK Biobank whole genomes.} The 200k whole-genome release with the corrected phenotype (182{,}238 participants, 1{,}753 Alzheimer's disease cases, phecode 290.11) on the UK Biobank null model used for the gene-set lookups; variants matched by position and alleles with no minor-allele-count floor; $P$ from the saddlepoint score test with the root-verified solver described in the Methods. The 21 non-\textit{APOE} variants fall in 11 one-megabase loci; the lead variant of each locus (smallest ADSP $P$) carries the replication test at $P < 0.05/11$, and the other variants are listed for completeness.", label="tab:s13"),
 "S14a": dict(csv="S14a_split_summary.csv", font=r"\small",
    cols=[col("half", f_str, "Half"), col("N", f_int, "Participants"), col("cases", f_int, "Cases"), col("controls", f_int, "Controls"), col("cells", f_int, "Cells"), col("clusters", f_int, "Relatedness clusters"), col("subjects_in_clusters_gt1", f_int, "In clusters of size $> 1$"), col("tested", f_int, "Variants tested"), col("lambda", f_num(3), r"$\lambda_{GC}$"), col("fdr_variants", f_int, "FDR 5\% variants"), col("fdr_non_apoe", f_int, r"outside \textit{APOE}")],
    caption=r"\textbf{Supplementary Table 14a. The two held-out halves of the stratified design.} Relatedness clusters (pairs related at third degree or closer in the sparse GRM, or members of the same declared family) were allocated whole to one half, balancing cases and controls within every project $\times$ center cell; each half was analyzed with the design of the main scan (own null model, own sparse GRM, \texttt{QC\_label7}, minor allele count $\geq 50$ within the half).", label="tab:s14a"),
 "S14b": dict(csv="S14b_split_crosstest.csv", font=r"\scriptsize", landscape=True,
    cols=[col("discovery_half", f_str, "Discovery half"), col("key", f_variant, "Variant (chr:pos:ref:alt)"), col("MAF", f_num(4)), col("Est", f_num(3), r"$\hat\beta$"), col("pvalue", f_sci, "$P$"), col("q", f_sci, "BH $q$"), col("Est_rep", f_num(3), r"$\hat\beta$ (other half)"), col("pvalue_rep", f_sci, "$P$ (other half)"), col("verdict", f_str, "Verdict")],
    caption=r"\textbf{Supplementary Table 14b. Held-out replication across the two halves.} Leads admitted at FDR 5\% in one half (pruned at 500 kb, the \textit{APOE} region as one locus) tested in the other half; replication requires the same direction and $P < 0.05$ divided by the number of non-\textit{APOE} leads of that discovery half.", label="tab:s14b"),
 "S14c": dict(csv="S14c_split_discovery_variants.csv", font=r"\scriptsize", landscape=True,
    cols=[col("key", f_variant, "Variant (chr:pos:ref:alt)"), col("Locus label (Figure 3a)", f_locus, "Locus"), col("Beta (log-odds)", f_num(3), r"Full-sample $\hat\beta$"), col("P", f_sci, "Full-sample $P$"), col("Est_A", f_num(3), r"$\hat\beta$ half A"), col("P_A", f_sci, "$P$ half A"), col("Est_B", f_num(3), r"$\hat\beta$ half B"), col("P_B", f_sci, "$P$ half B")],
    caption=r"\textbf{Supplementary Table 14c. The 21 non-\textit{APOE} single variants (11 loci) of the full stratified scan in each held-out half.} The \textit{APOE} lead and the 21 non-\textit{APOE} variants admitted at FDR 5\% in the full scan, with their estimates in each half (minor allele count $\geq 50$ within the half; blank where the variant does not reach that count; descriptive, each half having about half the sample). Of the 21, 15 can be tested in both halves; all 15 are concordant in direction, with $P < 0.01$ in each half.", label="tab:s14c"),



 "S15a": dict(csv="S15a_benchmark_sv_replicate_ci.csv", font=r"\scriptsize", landscape=True,
    cols=[col("Panel filter", f_str), col("Design", f_str), col("Test", f_str), col("alpha", f_sci, r"$\alpha$"), col("Null tests", f_int), col("Observed", f_int), col("Expected", f_num(1)),
          col("Ratio", f_num(3)), col("Bootstrap 95% lower", f_num(3), "Bootstrap lower"), col("Bootstrap 95% upper", f_num(3), "Bootstrap upper")],
    caption=r"\textbf{Supplementary Table 15a. Single-variant false-positive counts on the permuted null of the benchmark, with intervals from the variation between replicates.} For each design, the observed number of null tests below $\alpha$ over 1{,}000 permuted-null replicates, its expectation and their ratio, with a percentile bootstrap interval over replicates (2{,}000 resamples), for the normal-approximation score test used in the ADSP scans and for the saddlepoint test. The benchmark panel is shown with and without the removal of heterozygote-excess sites (pooled heterozygote excess at $P<10^{-4}$, the direction and threshold of \texttt{QC\_label7}; 76{,}734 of 281{,}156 sites). Without the filter, all-heterozygous mapping artifacts produce repeated false positives under the normal approximation (174 of 188 stratified events at $5\times10^{-8}$ fall on 19 such sites); with it, the two tests give the same counts. The \textit{APOE} and known-locus windows are excluded throughout.", label="tab:s15a"),
 "S15b": dict(csv="S15b_geneset_null_native_tests.csv", font=r"\scriptsize", landscape=True,
    cols=[col("Floor", f_str), col("Design", f_str), col("Test", f_test), col("Replicates", f_int), col("Gene tests per replicate", f_int, "Tests per replicate"), col("Observed at 2.5e-6", f_int, "Observed"), col("Ratio at 2.5e-6", f_num(2), r"Ratio at $2.5\times10^{-6}$"),
          col("Bootstrap 95% lower", f_num(2), "Lower"), col("Bootstrap 95% upper", f_num(2), "Upper"), col("Share of scans with an FDR-5% discovery", f_num(2), r"Scans with an FDR discovery"), col("Clopper-Pearson lower", f_num(3), "Lower"), col("Clopper-Pearson upper", f_num(3), "Upper")],
    caption=r"\textbf{Supplementary Table 15b. The gene-set statistics of this study on the permuted null, with and without the cMAC floor.} For the oracle, stratified and pooled designs, the false-positive ratio at the gene-scan threshold with a percentile bootstrap over replicates, and the share of null scans in which the Benjamini--Hochberg procedure at FDR 5\% over all coding genes admits at least one gene, which under a complete null equals the realized false discovery rate, with a Clopper--Pearson interval. The SKAT--burden omnibus test is the Cauchy combination of the SKAT and burden $p$ values (STAAR-O without ACAT-V); burden and SKAT are the corresponding Cauchy combinations of the two weight families. Genes with cMAC $<10$ have discrete, skewed statistics; removing them brings the SKAT--burden omnibus test and the burden test to nominal and leaves STAAR-O above it.", label="tab:s15b"),
 "S16": dict(csv="S1h_diversity_accounting.csv", font=r"\scriptsize", landscape=True,
    cols=[col("dimension", f_str, "Dimension"), col("group", f_str, "Group"), col("n_pool", f_int, "Eligible pool"), col("n_analytic", f_int, "Analytic set"), col("pct_retained", f_num(1), r"\% retained"),
          col("pool_share", f_num(2), r"\% of pool"), col("analytic_share", f_num(2), r"\% of analytic set"), col("info_share", f_num(2), r"\% of information")],
    caption=r"\textbf{Supplementary Table 16. Who the stratified design retains, and how much information each group carries.} For each age band, contributing cohort, genetic-ancestry super-population and self-reported group: the eligible pool (38{,}129 participants with recorded status, sex, age and project), the analytic set (30{,}894), the percentage retained, each group's share of the pool and of the analytic set, and its share of the design's information, where a participant's information weight is $\hat\pi_c(1-\hat\pi_c)$ for the case fraction $\hat\pi_c$ of its project $\times$ sequencing center cell. Groups whose cells contain one outcome only, such as the all-control LASI-DAD project, contribute no within-cell contrast and are not retained.", label="tab:s16"),
 "S17": dict(csv="S17_conditional_apoe_trem2.csv", font=r"\scriptsize", landscape=True,
    cols=[col("Gene", f_gene), col("Mask", f_str), col("Conditioned on", f_str), col("Variants conditioned on", f_int, "Known variants"), col("SB P", f_sci, r"\shortstack{SKAT--burden\\omnibus test $P$}"), col("SB P conditional", f_sci, r"\shortstack{SKAT--burden\\omnibus test $P$,\\conditional}"),
          col("STAAR-O P", f_sci, "STAAR-O $P$"), col("STAAR-O P conditional", f_sci, "STAAR-O $P$, conditional")],
    caption=r"\textbf{Supplementary Table 17. Conditional analysis of the \textit{APOE}-region and \textit{TREM2} gene sets.} Each set admitted at FDR 5\% in the \textit{APOE} region is re-tested conditional on \textit{APOE} $\varepsilon$4 (rs429358) and $\varepsilon$2 (rs7412), and additionally on the lead indel of the single-variant scan (chr19:44{,}909{,}967 TGG$>$T); the \textit{TREM2} missense and promoter sets are conditioned on R47H and R62H. STAARpipeline's conditional functions with the same masks, annotations and MAF cutoff as the discovery scan; the SKAT--burden omnibus test is the primary statistic. Every \textit{APOE}-region set loses significance, whereas the \textit{TREM2} sets are attenuated but retain association.", label="tab:s17"),
 "S18": dict(csv="S18_robustness_specifications.csv", font=r"\tiny", landscape=True,
    cols=[col("Signal", f_str), col("Family", f_str), col("FDR boundary", f_sci, "FDR boundary"), col("full", f_sci, "Full model"), col("nocls", f_sci, "No latent-class terms"), col("pc20", f_sci, "PC1--20"),
          col("noAPOEsel", f_sci, r"No \textit{APOE}-extreme projects"), col("noSelected", f_sci, r"No \textit{APOE}-extreme or early-onset projects"), col("LOO n", f_int, "LOO fits"), col("LOO min p", f_sci, "LOO min $P$"), col("LOO max p", f_sci, "LOO max $P$"),
          col("LOO n below boundary", f_int, "LOO fits below boundary")],
    caption=r"\textbf{Supplementary Table 18. Robustness of every reported signal to the null-model specification and to leaving out any one project.} Each signal of the stratified design is re-tested under refitted null models: the full model; without the latent-class terms; with 20 instead of 10 principal components; without the two \textit{APOE}-extreme projects; without the \textit{APOE}-extreme and early-onset projects; and leaving out each project in turn (LOO). Gene sets use the primary SKAT--burden omnibus test; single variants the score test. The FDR boundary is the critical value of the full scan in that family, so a LOO $P$ below it means the signal would still be admitted with that project removed, which is a conservative yardstick because every exclusion also reduces the sample.", label="tab:s18"),
 "S19": dict(csv="S19_sv_two_tests.csv", font=r"\tiny", landscape=True,
    cols=[col("CHR", f_str, "Chr"), col("POS", f_int, "Position"), col("REF", f_str), col("ALT", f_str), col("MAC", f_int), col("MAF", f_num(4)), col("P normal approximation", f_sci, "$P$, normal approximation"), col("P saddlepoint", f_sci, "$P$, saddlepoint"),
          col("Ratio", f_num(2), "Ratio"), col("Below full-scan boundary (saddlepoint)", f_str, "Below boundary (saddlepoint)")],
    caption=r"\textbf{Supplementary Table 19. The 296 single variants admitted at FDR 5\% re-tested with STAARpipeline's saddlepoint score test.} Both tests use the deployed null model of the stratified design. The saddlepoint $P$ values are larger by a near-constant factor across the allele-frequency range, including common variants where the two approximations should coincide, so the difference reflects the implementations' treatment of the null model rather than a tail correction; the benchmark, which implements both approximations on the same score, shows the normal approximation calibrated (Supplementary Table~15a). Under the saddlepoint variant, 274 of the 296 variants remain below the full-scan boundary ($3.9\times10^{-7}$), including the \textit{BIN1} and \textit{SHARPIN} leads and the rare variants near \textit{MMP17}, \textit{KLHL1}, \textit{CRTAC1} and \textit{TUBB8}; the \textit{MS4A4A}, \textit{CELF4} and \textit{DGKB} variants and the rare \textit{KCNV1} and \textit{BAIAP2L1} variants do not.", label="tab:s19"),
 "S20": dict(csv="S20_adsp_gvc_catalog.csv", font=r"\scriptsize", landscape=True,
    cols=[col("Gene", f_gene), col("Family in this study", f_str, "Family"), col("ADSP GVC status", f_str, "Status in the GVC catalog"), col("Entries as reported locus", f_int, "Entries"), col("Best tier", f_str, "Best tier"), col("Evidence type", f_str), col("Adjudicated studies", f_str), col("Entries within 500 kb", f_int, "Entries within 500 kb")],
    caption=r"\textbf{Supplementary Table 20. The reported genes in the ADSP Gene Verification Committee catalog.} The Top Genes Repository (topgenes.niagads.org, accessed 13 October 2026) lists 423 locus entries from 29 adjudicated studies, with committee-assigned evidence tiers 1 to 7 (195 entries at tier 1), separate single-variant and variant-set (gene-based) evidence, and the population of each study; 357 entries come from European-ancestry analyses. A gene is a reported locus when it appears in an entry's reported or grouped loci; otherwise it may lie within 500 kb of a reported locus. \textit{NECTIN2} has one entry, a tier-2 gene-based report in European ancestry (Nazarian 2019), and otherwise appears only as a neighbor of the \textit{APOE} locus, consistent with the conditional analysis of this study.", label="tab:s20"),
 "S6": dict(csv="S6_variant_counts_by_chromosome.csv", font=r"\scriptsize", landscape=True,
    cols=[col("Chromosome", f_str), col("Variant records in the annotated GDS (bi-allelic after decomposition)", f_int, "Records in aGDS"), col("Pass QC_label7 (AN > 3,636; RUTH signed log10 p > -4; F_MISSING < 0.05)", f_int, "QC pass"), col("QC_label7 pass fraction", f_num(4), "Fraction"),
          col("Single-variant tests (MAC >= 50, QC_label7)", f_int, "SV tests"), col("Gene-centric coding gene x mask tests", f_int, "Coding tests"), col("Gene-centric noncoding gene x mask tests", f_int, "Noncoding tests"), col("Gene-centric ncRNA gene x mask tests", f_int, "ncRNA tests")],
    caption=r"\textbf{Supplementary Table 6. Per-chromosome counts of variant records, of records carrying the variant QC label and of tests performed.} Records are bi-allelic after multi-allelic decomposition; the QC label requires AN $> 3{,}636$, no excess of heterozygotes at $p < 10^{-4}$ under \textsc{Ruth} and genotype missingness $< 5\%$; single-variant tests require minor allele count $\geq 50$ in the analytic set; gene $\times$ mask tests are those with at least two rare variants (MAF $< 1\%$) in the set.", label="tab:s6"),
 "S7": dict(csv="S7_threshold_tiers.csv", font=r"\scriptsize", landscape=True,
    cols=[col("Family", f_str), col("Tests (m)", f_int, "$m$"), col("Masks", f_int), col("lambda_GC", f_num(3), r"$\lambda_{GC}$"), col("Fixed threshold", f_sci, "Fixed"), col("Admitted at fixed", f_int, "$n$"),
          col("Bonferroni 0.05/m", f_sci, "$0.05/m$"), col("Admitted at Bonferroni", f_int, "$n$"), col("Admitted at FDR 5% (k)", f_int, r"FDR $k$"), col("FDR critical value 0.05 k/m", f_sci, "$0.05\\,k/m$"), col("Largest admitted p", f_sci, "Largest $P$")],
    caption=r"\textbf{Supplementary Table 7. What each threshold admits in every scan of the stratified design.} The fixed thresholds are the conventions ($5\times10^{-8}$ single variants; $5\times10^{-7}$ coding, $3.57\times10^{-7}$ noncoding and $2.5\times10^{-6}$ ncRNA, i.e.\ Bonferroni over 20{,}000 genes times the number of masks); Bonferroni is over the tests actually performed; the FDR rule is Benjamini--Hochberg at 5\% within the family, the rule used throughout the paper. $\lambda_{GC}$ is the bulk genomic-control factor of the family.", label="tab:s7"),
}

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--csvdir", default="${ADSP_ROOT}/revision/tables/supp_v3")
    ap.add_argument("--outdir", default="${ADSP_ROOT}/revision/manuscript/supp_tables")
    a = ap.parse_args()
    outdir = pathlib.Path(a.outdir); outdir.mkdir(parents=True, exist_ok=True)
    for key, sp in SPECS.items():
        path = pathlib.Path(a.csvdir) / sp["csv"]
        if not path.exists():
            print(f"{key}: MISSING {path.name}"); continue
        rows = read_csv(path)
        tex = longtable(rows, sp["cols"], sp["caption"], sp["label"], font=sp.get("font", r"\footnotesize"),
                        colspec=sp.get("colspec"), landscape=sp.get("landscape", False), note=sp.get("note"))
        (outdir / f"{key}.tex").write_text(tex)
        print(f"{key}: {len(rows)} rows -> {key}.tex")

if __name__ == "__main__":
    main()
