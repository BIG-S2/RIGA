#!/usr/bin/env Rscript
## supplementary_table2_single_variants.R (2026-10-14): Supplementary Table 2 (single variants at FDR 5%) from the QC_label7 rescan of the
## stratified design, with the single-variant rows of Tables 6 and 7 updated in place. The S2 section is taken verbatim from
## supplementary_table1_sample_description.R (2026-10-04) with the input paths changed to individual_analysis_qc7/.
suppressPackageStartupMessages({ library(data.table); library(SeqArray); library(gdsfmt); library(STAARpipelineSummary) })
OUT  <- "${ADSP_ROOT}/revision/tables/supp_v3"; ROOT <- "${ADSP_ROOT}/output/matched_lpa_pc_pair33_projcell"
log_lines <- character(0); LOG <- function(...) { s <- sprintf(...); message(s); log_lines <<- c(log_lines, s) }
file.copy(file.path(OUT, "S2_fdr_single_variants.csv"), file.path(OUT, "S2_fdr_single_variants.qc2scan_2026-10-04.csv"), overwrite = FALSE)
## ---------------- S2: FDR-significant single variants ----------------
LOG(""); LOG("=== S2: single variants at FDR 5%% ===")
g <- as.data.table(get(load(file.path(ROOT, "individual_analysis_qc7", "results_individual_analysis_genome.Rdata"))))   # QC_label7 rescan (2026-10-14)
LOG("single-variant results: %d rows; columns: %s", nrow(g), paste(names(g), collapse = ", "))
g <- g[is.finite(pvalue) & pvalue > 0]
m_sv <- nrow(g); g[, q := p.adjust(pvalue, "BH")]
lambda_sv <- qchisq(median(g$pvalue), df = 1, lower.tail = FALSE) / qchisq(0.5, df = 1)
sig <- g[q < 0.05][order(pvalue)]; k_sv <- nrow(sig)
LOG("tests %d; lambda_GC %.4f; BH-5%% discoveries %d; critical value %.3e; largest admitted p %.3e",
    m_sv, lambda_sv, k_sv, 0.05 * k_sv / m_sv, max(sig$pvalue))
sv_by_chr <- g[, .(`Single-variant tests (MAC >= 50, QC_label7)` = .N), by = .(Chromosome = CHR)][order(Chromosome)]
tiers_sv <- data.table(Family = "Single variant", `Tests (m)` = m_sv, Masks = 1L, lambda_GC = round(lambda_sv, 4),
  `Fixed threshold` = 5e-8, `Admitted at fixed` = sum(g$pvalue < 5e-8),
  `Bonferroni 0.05/m` = 0.05 / m_sv, `Admitted at Bonferroni` = sum(g$pvalue < 0.05 / m_sv),
  `Admitted at FDR 5% (k)` = k_sv, `FDR critical value 0.05 k/m` = 0.05 * k_sv / m_sv,
  `Largest admitted p` = max(sig$pvalue))
rm(g); invisible(gc())

## annotated-Manhattan labels for the same variants
lab <- fread(file.path(ROOT, "individual_analysis_qc7", "manhattan_MAC_20_sig_annotated_v4.csv"))
sig <- merge(sig, lab[, .(CHR, POS, REF, ALT, locus_label = label)],
             by = c("CHR", "POS", "REF", "ALT"), all.x = TRUE)
LOG("labels matched for %d of %d variants", sum(!is.na(sig$locus_label)), nrow(sig))

## nearest protein-coding gene and FAVOR variant-to-gene relation (as in build_supp_tables.R)
agds_dir  <- get(load("${ADSP_ROOT}/output/agds_dir.Rdata"))
genes_tbl <- get("genes", envir = asNamespace("STAARpipelineSummary"))
genes_tbl$chromosome_name <- as.integer(genes_tbl$chromosome_name)
genes_tbl <- genes_tbl[!is.na(genes_tbl$chromosome_name) & genes_tbl$chromosome_name %in% 1:22, ]
nearest_hgnc <- function(chr, pos) {
  g <- genes_tbl[genes_tbl$chromosome_name == chr, , drop = FALSE]
  if (nrow(g) == 0) return(c(NA_character_, NA_real_))
  d <- pmax(pmax(g$start_position - pos, 0), pmax(pos - g$end_position, 0))
  d[pos >= g$start_position & pos <= g$end_position] <- 0
  i <- which.min(d); c(g$hgnc_symbol[i], d[i])
}
first_gene <- function(s) {
  if (is.na(s) || s == "") return(NA_character_)
  tok <- trimws(sub("\\(.*$", "", strsplit(s, "[,;]")[[1]][1]))
  if (tok == "" || tok == "NONE") NA_character_ else tok
}
genehancer_gene <- function(s) {
  if (is.na(s) || s == "") return(NA_character_)
  m <- regmatches(s, regexec("connected_gene=([A-Za-z0-9._-]+)", s))[[1]]
  if (length(m) >= 2) m[2] else NA_character_
}
fmt_exonic <- function(s) {
  s <- sub(" SNV$", "", s)
  s <- sub("^stopgain$", "Stop-gain", s); s <- sub("^stoploss$", "Stop-loss", s)
  s <- sub("^nonsynonymous$", "Nonsynonymous", s); sub("^synonymous$", "Synonymous", s)
}
make_relation_label <- function(gc_cat, gc_info, ex_cat, gh, cage, dhs) {
  gene_main <- first_gene(gc_info)
  nz <- function(x) if (is.na(x)) "" else x
  ex_cat <- nz(ex_cat); gc_cat <- nz(gc_cat); gh <- nz(gh); cage <- nz(cage); dhs <- nz(dhs)
  if (ex_cat != "") {
    if (!is.na(gene_main)) return(sprintf("%s of %s", fmt_exonic(ex_cat), gene_main))
    return(fmt_exonic(ex_cat))
  }
  lc <- tolower(gc_cat)
  if (lc %in% c("intronic", "ncrna_intronic") && !is.na(gene_main)) return(sprintf("Intronic to %s", gene_main))
  if (lc %in% c("utr3", "utr5", "ncrna_utr3", "ncrna_utr5") && !is.na(gene_main))
    return(sprintf("%s of %s", toupper(sub("ncrna_", "", lc)), gene_main))
  if (lc %in% c("upstream", "downstream", "upstream;downstream") && !is.na(gene_main))
    return(sprintf("%s of %s", paste0(toupper(substr(gc_cat, 1, 1)), substr(gc_cat, 2, nchar(gc_cat))), gene_main))
  if (lc %in% c("splicing", "ncrna_splicing") && !is.na(gene_main)) return(sprintf("Splicing of %s", gene_main))
  if (lc %in% c("ncrna_exonic", "ncrna") && !is.na(gene_main)) return(sprintf("ncRNA of %s", gene_main))
  if (cage != "") return(if (!is.na(gene_main)) sprintf("CAGE promoter near %s", gene_main) else "CAGE promoter")
  if (gh != "") { eg <- genehancer_gene(gh); if (is.na(eg)) eg <- gene_main
    return(if (!is.na(eg)) sprintf("Predicted enhancer of %s", eg) else "Predicted enhancer") }
  if (dhs != "") return(if (!is.na(gene_main)) sprintf("DHS near %s", gene_main) else "DHS")
  if (!is.na(gene_main)) return(sprintf("Intergenic near %s", gene_main))
  "Intergenic"
}
sig[, varinfo := paste(CHR, POS, REF, ALT, sep = "-")]
sig[, `:=`(relation_label = NA_character_, gencode_category = NA_character_, gencode_exonic_category = NA_character_)]
for (chrN in sort(unique(sig$CHR))) {
  rows <- which(sig$CHR == chrN); gp <- agds_dir[chrN]
  if (is.na(gp) || !file.exists(gp)) { LOG("  chr%d: aGDS missing, skipped", chrN); next }
  gg <- seqOpen(gp)
  varinfo_all <- seqGetData(gg, "annotation/info/FunctionalAnnotation/VarInfo")
  idx <- match(sig$varinfo[rows], varinfo_all); found <- !is.na(idx)
  LOG("  chr%d: %d FDR variants, %d matched in aGDS annotation", chrN, length(rows), sum(found))
  if (any(found)) {
    seqSetFilter(gg, variant.sel = idx[found], verbose = FALSE)
    A <- "annotation/info/FunctionalAnnotation/"
    gc_cat <- seqGetData(gg, paste0(A, "genecode_comprehensive_category"))
    gc_info <- seqGetData(gg, paste0(A, "genecode_comprehensive_info"))
    ex_cat <- seqGetData(gg, paste0(A, "genecode_comprehensive_exonic_category"))
    gh <- seqGetData(gg, paste0(A, "genehancer")); cage <- seqGetData(gg, paste0(A, "cage_tc"))
    dhs <- seqGetData(gg, paste0(A, "rdhs"))
    labs <- vapply(seq_along(gc_cat), function(i)
      make_relation_label(gc_cat[i], gc_info[i], ex_cat[i], gh[i], cage[i], dhs[i]), character(1))
    set(sig, i = rows[found], j = "relation_label", value = labs)
    set(sig, i = rows[found], j = "gencode_category", value = gc_cat)
    set(sig, i = rows[found], j = "gencode_exonic_category", value = ex_cat)
  }
  seqClose(gg); rm(varinfo_all); invisible(gc())
}
nn <- t(mapply(nearest_hgnc, sig$CHR, sig$POS))
sig[, nearest_gene := nn[, 1]]; sig[, nearest_gene_distance_bp := as.numeric(nn[, 2])]
has_est <- all(c("Est", "Est_se") %in% names(sig))
if (has_est) {
  sig[, `:=`(OR = exp(Est), OR_L95 = exp(Est - 1.959964 * Est_se), OR_U95 = exp(Est + 1.959964 * Est_se))]
} else LOG("no Est/Est_se columns in the saddlepoint single-variant scan: effect-size columns omitted")
setorder(sig, pvalue)
s2 <- data.table(CHR = sig$CHR, POS = sig$POS, REF = sig$REF, ALT = sig$ALT,
                 `Variant (chr:pos:ref:alt)` = sprintf("%d:%d:%s:%s", sig$CHR, sig$POS, sig$REF, sig$ALT))
for (cc in c("ALT_AF", "MAF", "N", "MAC", "Score", "Score_se", "Est", "Est_se", "OR", "OR_L95", "OR_U95", "pvalue", "pvalue_log10", "q"))
  if (cc %in% names(sig)) s2[, (cc) := sig[[cc]]]
setnames(s2, old = intersect(c("Est", "Est_se", "OR_L95", "OR_U95", "pvalue", "pvalue_log10", "q", "Score_se"), names(s2)),
         new = c("Beta (log-odds)", "SE", "OR 95% CI lower", "OR 95% CI upper", "P", "-log10 P", "BH q value", "Score SE")[
           match(intersect(c("Est", "Est_se", "OR_L95", "OR_U95", "pvalue", "pvalue_log10", "q", "Score_se"), names(s2)),
                 c("Est", "Est_se", "OR_L95", "OR_U95", "pvalue", "pvalue_log10", "q", "Score_se"))])
s2[, `Locus label (Figure 3a)` := sig$locus_label]
s2[, `Nearest protein-coding gene` := sig$nearest_gene]
s2[, `Distance to nearest gene (bp)` := sig$nearest_gene_distance_bp]
s2[, `Variant-to-gene relation (FAVOR)` := sig$relation_label]
s2[, `GENCODE category` := sig$gencode_category]
s2[, `GENCODE exonic category` := sig$gencode_exonic_category]
fwrite(s2, file.path(OUT, "S2_fdr_single_variants.csv"))
LOG("S2 written: %d variants (%d in the APOE region chr19:44.4-46.4 Mb)", nrow(s2),
    sum(s2$CHR == 19 & s2$POS >= 44.4e6 & s2$POS <= 46.4e6))
## ---- Table 7: replace the single-variant row; Table 6: replace the single-variant test counts ----
s7 <- fread(file.path(OUT, "S7_threshold_tiers.csv")); s7 <- rbind(tiers_sv, s7[Family != "Single variant"], fill = TRUE); fwrite(s7, file.path(OUT, "S7_threshold_tiers.csv"))
s6 <- fread(file.path(OUT, "S6_variant_counts_by_chromosome.csv")); cc <- grep("^Single-variant tests", names(s6), value = TRUE)[1]
s6 <- merge(s6[, setdiff(names(s6), cc), with = FALSE], sv_by_chr, by = "Chromosome", all.x = TRUE); setnames(s6, "Single-variant tests (MAC >= 50, QC_label7)", cc)
fwrite(s6, file.path(OUT, "S6_variant_counts_by_chromosome.csv")); LOG("S6 and S7 single-variant rows updated")
writeLines(log_lines, file.path(OUT, "build_log_S2_qc7.txt")); cat("DONE\n")
