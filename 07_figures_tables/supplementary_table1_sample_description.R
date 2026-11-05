###########################################################
# supplementary_table1_sample_description.R (2026-10-04)
# Supplementary Tables for the STRATIFIED design D_S
# (output/matched_lpa_pc_pair33_projcell; QC_label7; FDR 5% rule).
# Derived from revision/scripts/build_supp_tables.R (matched design,
# 2026-08-06). Writes revision/tables/supp_v3/. Every column comes
# from a result object on disk; anything absent is logged and omitted.
###########################################################
suppressPackageStartupMessages({
  library(data.table); library(SeqArray); library(gdsfmt)
  library(STAARpipelineSummary)
})
OUT   <- "${ADSP_ROOT}/revision/tables/supp_v3"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
ROOT  <- "${ADSP_ROOT}/output/matched_lpa_pc_pair33_projcell"
MATCH <- "${ADSP_ROOT}/output/matched/matched_samples_lpa_pc.pair33_projcell.csv"
TAB   <- "${ADSP_ROOT}/revision/tables"
SRC   <- "${ADSP_ROOT}/src"
log_lines <- character(0)
LOG <- function(...) { s <- sprintf(...); message(s); log_lines <<- c(log_lines, s) }

## ---------------- S1: from the release to the analytic set ----------------
LOG("=== S1: analytic dataset flow and composition ===")
ph_all <- fread("${ADSP_ROOT}/pheno/adsp_combined.csv")
n_release <- uniqueN(fread("${ADSP_ROOT}/pheno/gcad.r5.wgs.58507.2024.11.03.quality.metrics.ALL.csv",
                           select = "Sample_name", colClasses = "character")$Sample_name)
LOG("R5 release quality-metrics manifest: %d sequenced individuals; adsp_combined.csv rows: %d", n_release, nrow(ph_all))
n_pheno <- nrow(ph_all[!is.na(AD) & !is.na(Sex) & !is.na(Age)])
LOG("non-missing AD/Sex/Age: %d", n_pheno)
## the pre-design pool exactly as riga_study_units_null_model.R built it
src <- readLines(file.path(SRC, "prepare_sample_and_one_to_one_matching.R"))
cut <- grep("^m\\.out <- matchit", src)[1]
eval(parse(text = paste(src[1:(cut - 1)], collapse = "\n")))
dfx <- as.data.table(dfx)
if (!"Cohort" %in% names(dfx))
  dfx <- merge(dfx, ph_all[, .(IID, Cohort)], by = "IID", all.x = TRUE)
qc <- fread("${ADSP_ROOT}/pheno/gcad.r5.wgs.58507.2024.11.03.quality.metrics.ALL.csv",
            select = c("Sample_name", "Study"), colClasses = "character")
setnames(qc, "Sample_name", "IID"); qc <- unique(qc, by = "IID")
dfx <- merge(dfx, qc, by = "IID", all.x = TRUE)
dfx[is.na(Study) | Study == "", Study := "Unknown"]
dfx[, cell_SC := paste(as.character(Study), as.character(SeqCenter), sep = "|")]
n_pool <- nrow(dfx)
LOG("pre-design pool (phenotype + latent class + center + project): %d (cases %d)",
    n_pool, sum(dfx$AD == 1))
cells <- dfx[, .(N = .N, Cases = sum(AD == 1), Controls = sum(AD == 0)),
             by = .(`Project x center cell` = cell_SC, `Sequencing project` = Study,
                    `Sequencing center` = as.character(SeqCenter))]
cells[, `Retained (both labels present)` := Cases > 0 & Controls > 0]
setorder(cells, -`Retained (both labels present)`, -N)
LOG("project x center cells: %d; informative %d; subjects in single-label cells %d",
    nrow(cells), sum(cells$`Retained (both labels present)`),
    sum(cells[`Retained (both labels present)` == FALSE]$N))
matched <- fread(MATCH)
n_an <- nrow(matched); n_case <- sum(matched$AD == 1); n_ctrl <- sum(matched$AD == 0)
LOG("analytic set: %d (cases %d / controls %d)", n_an, n_case, n_ctrl)
stopifnot(n_an == sum(cells[`Retained (both labels present)` == TRUE]$N))

s1 <- data.table(
  Step = c("1. ADSP R5 joint-called WGS release (GRCh38; sequenced individuals in the release quality-metrics manifest)",
           "2. Harmonized phenotype table with non-missing AD status, sex and age",
           "3. Latent-class assignment, sequencing center and sequencing-project label available (pre-design pool)",
           "4. Member of a project x center cell that contains both cases and controls (stratified design; Proposition 1b)",
           "4a. Analytic set -- cases", "4b. Analytic set -- controls"),
  N = c(n_release, n_pheno, n_pool, n_an, n_case, n_ctrl),
  `N removed at this step` = c(NA_integer_, n_release - n_pheno, n_pheno - n_pool,
                               n_pool - n_an, NA_integer_, NA_integer_))
fwrite(s1, file.path(OUT, "S1a_analytic_dataset_flow.csv"))

## composition of the analytic set
old_race <- fread(file.path(TAB, "S1b_matched_ancestry.csv"))
race_lab <- setNames(old_race$`Self-reported group`, as.character(old_race$`Race code`))
rt <- matched[, .(N = .N, Cases = sum(AD == 1), Controls = sum(AD == 0)), by = .(`Race code` = Race)][order(-N)]
rt[, `Self-reported group` := race_lab[as.character(`Race code`)]]
rt[, `% of analytic set` := round(100 * N / n_an, 2)]
setcolorder(rt, c("Race code", "Self-reported group", "N", "Cases", "Controls", "% of analytic set"))
fwrite(rt, file.path(OUT, "S1b_analytic_race.csv"))
ct <- matched[, .(N = .N, Cases = sum(AD == 1), Controls = sum(AD == 0), `Cells` = uniqueN(cell_SC)),
              by = .(`Sequencing center` = SeqCenter)][order(-N)]
ct[, `% of analytic set` := round(100 * N / n_an, 2)]
fwrite(ct, file.path(OUT, "S1c_analytic_seqcenter.csv"))
proj <- merge(dfx[, .(`N in pool` = .N, `Cases in pool` = sum(AD == 1)), by = .(`Sequencing project` = Study)],
              matched[, .(`N retained` = .N, `Cases retained` = sum(AD == 1), `Controls retained` = sum(AD == 0)),
                      by = .(`Sequencing project` = Study)], all.x = TRUE)
for (cc in c("N retained", "Cases retained", "Controls retained")) proj[is.na(get(cc)), (cc) := 0L]
proj[, Retention := round(`N retained` / `N in pool`, 3)]
setorder(proj, -`N in pool`)
fwrite(proj, file.path(OUT, "S1d_analytic_projects.csv"))
coh <- merge(dfx[, .(`N in pool` = .N, `Cases in pool` = sum(AD == 1)), by = .(Cohort)],
             matched[, .(`N retained` = .N, `Cases retained` = sum(AD == 1), `Controls retained` = sum(AD == 0)),
                     by = .(Cohort)], all.x = TRUE)
for (cc in c("N retained", "Cases retained", "Controls retained")) coh[is.na(get(cc)), (cc) := 0L]
coh[, Retention := round(`N retained` / `N in pool`, 3)]
setorder(coh, -`N in pool`)
fwrite(coh, file.path(OUT, "S1e_analytic_cohorts.csv"))
fwrite(cells, file.path(OUT, "S1f_project_center_cells.csv"))
lc <- matched[, .(N = .N, Cases = sum(AD == 1), Controls = sum(AD == 0)), by = .(`Latent class` = lpa_cluster)][order(`Latent class`)]
fwrite(lc, file.path(OUT, "S1g_analytic_latent_class.csv"))
LOG("S1 written: race %d rows, center %d, projects %d (retained %d), cohorts %d (retained %d), cells %d, classes %d",
    nrow(rt), nrow(ct), nrow(proj), sum(proj$`N retained` > 0), nrow(coh), sum(coh$`N retained` > 0), nrow(cells), nrow(lc))
rm(dfx); invisible(gc())

## ---------------- S2: FDR-significant single variants ----------------
LOG(""); LOG("=== S2: single variants at FDR 5%% ===")
g <- as.data.table(get(load(file.path(ROOT, "individual_analysis",
                                      "results_individual_analysis_genome.QC7.Rdata"))))
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
lab <- fread(file.path(ROOT, "individual_analysis", "manhattan_MAC_20_sig_annotated.csv"))
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

## ---------------- S3: gene x mask sets at FDR 5% (per family) ----------------
LOG(""); LOG("=== S3: gene-centric sets at FDR 5%% ===")
to_dt <- function(obj, cat) {
  ## the merged STAARpipelineSummary objects are matrices of mode list (is.list TRUE): keep them whole
  if (is.list(obj) && !is.data.frame(obj) && !is.matrix(obj))
    obj <- do.call(rbind, Filter(function(z) is.matrix(z) || is.data.frame(z), obj))
  x <- as.data.table(as.data.frame(obj, stringsAsFactors = FALSE, check.names = FALSE))
  for (cc in names(x)) if (is.list(x[[cc]])) set(x, j = cc, value = unlist(x[[cc]], use.names = FALSE))
  for (cc in setdiff(names(x), c("Gene name", "Category"))) {
    v <- suppressWarnings(as.numeric(x[[cc]]))
    if (all(is.na(v) == is.na(x[[cc]]))) set(x, j = cc, value = v)
  }
  if (!"Category" %in% names(x)) x[, Category := cat]
  x
}
load_cat <- function(dir, cats, files = paste0(cats, ".Rdata")) rbindlist(lapply(seq_along(cats), function(i) {
  x <- to_dt(get(load(file.path(dir, files[i]))), cats[i]); LOG("  %s/%s: %d rows", basename(dir), cats[i], nrow(x)); x }),
  fill = TRUE, use.names = TRUE)
fam <- list(
  coding    = load_cat(file.path(ROOT, "coding"), c("plof", "plof_ds", "missense", "disruptive_missense", "synonymous")),
  noncoding = load_cat(file.path(ROOT, "noncoding"), c("UTR", "upstream", "downstream", "promoter_CAGE", "promoter_DHS", "enhancer_CAGE", "enhancer_DHS")),
  ncRNA     = load_cat(file.path(ROOT, "ncRNA"), "ncRNA", "results_ncRNA_genome.Rdata"))
fixed_alpha <- c(coding = 5e-7, noncoding = 3.57e-7, ncRNA = 2.5e-6)
fam_label <- c(coding = "Gene-centric coding", noncoding = "Gene-centric noncoding", ncRNA = "Gene-centric ncRNA")
s3 <- list(); tiers <- list(tiers_sv); gm_by_chr <- list()
for (f in names(fam)) {
  d <- fam[[f]]
  LOG("%s columns: %s", f, paste(names(d), collapse = ", "))
  d[, p := as.numeric(`STAAR-O`)]
  d <- d[is.finite(p) & p > 0]
  m <- nrow(d); d[, q := p.adjust(p, "BH")]
  lam <- qchisq(median(d$p), df = 1, lower.tail = FALSE) / qchisq(0.5, df = 1)
  hit <- d[q < 0.05][order(p)]; k <- nrow(hit)
  LOG("%s: %d tests (%d masks), lambda %.3f, fixed %.2e -> %d, Bonferroni %.2e -> %d, FDR 5%% -> %d (critical %.3e, largest p %.3e)",
      f, m, uniqueN(d$Category), lam, fixed_alpha[f], sum(d$p < fixed_alpha[f]), 0.05 / m, sum(d$p < 0.05 / m), k,
      0.05 * k / m, if (k) max(hit$p) else NA_real_)
  tiers[[f]] <- data.table(Family = fam_label[f], `Tests (m)` = m, Masks = uniqueN(d$Category), lambda_GC = round(lam, 4),
    `Fixed threshold` = fixed_alpha[f], `Admitted at fixed` = sum(d$p < fixed_alpha[f]),
    `Bonferroni 0.05/m` = 0.05 / m, `Admitted at Bonferroni` = sum(d$p < 0.05 / m),
    `Admitted at FDR 5% (k)` = k, `FDR critical value 0.05 k/m` = 0.05 * k / m,
    `Largest admitted p` = if (k) max(hit$p) else NA_real_)
  if ("Chr" %in% names(d)) gm_by_chr[[f]] <- d[, .N, by = .(Chromosome = as.integer(Chr))][order(Chromosome)]
  if (k) { hit[, Analysis := fam_label[f]]; hit[, `BH q value` := q]; hit[, p := NULL]; hit[, q := NULL]; s3[[f]] <- hit }
}
s3 <- rbindlist(s3, fill = TRUE, use.names = TRUE)
keep <- c("Analysis", "Gene name", "Chr", "Category", "#SNV", "cMAC",
          "SKAT(1,25)", "SKAT(1,1)", "Burden(1,25)", "Burden(1,1)", "ACAT-V(1,25)", "ACAT-V(1,1)",
          "STAAR-S(1,25)", "STAAR-S(1,1)", "STAAR-B(1,25)", "STAAR-B(1,1)", "STAAR-A(1,25)", "STAAR-A(1,1)",
          "ACAT-O", "STAAR-O", "BH q value")
keep <- keep[keep %in% names(s3)]
s3 <- s3[, ..keep]
setorder(s3, Analysis, `STAAR-O`)
fwrite(s3, file.path(OUT, "S3_fdr_gene_mask_sets.csv"))
LOG("S3 written: %d rows (%s)", nrow(s3), paste(sprintf("%s/%s %.2g", s3$`Gene name`, s3$Category, s3$`STAAR-O`), collapse = "; "))

## ---------------- S7: threshold tiers per family ----------------
s7 <- rbindlist(tiers, use.names = TRUE)
fwrite(s7, file.path(OUT, "S7_threshold_tiers.csv"))
LOG("S7 written: %d families", nrow(s7))

## ---------------- S6: per-chromosome variant and test counts ----------------
LOG(""); LOG("=== S6: per-chromosome counts ===")
pf <- fread(file.path(TAB, "qclabel_pass_fractions.csv"))[node == "QC_label7", .(Chromosome = chr,
        `Variant records in the annotated GDS (bi-allelic after decomposition)` = n_variants,
        `Pass QC_label7 (AN > 3,636; RUTH signed log10 p > -4; F_MISSING < 0.05)` = n_pass,
        `QC_label7 pass fraction` = pass_frac)]
s6 <- merge(pf, sv_by_chr, by = "Chromosome", all.x = TRUE)
for (f in names(gm_by_chr)) {
  x <- copy(gm_by_chr[[f]]); setnames(x, "N", sprintf("%s gene x mask tests", fam_label[f]))
  s6 <- merge(s6, x, by = "Chromosome", all.x = TRUE)
}
setorder(s6, Chromosome)
tot <- s6[, lapply(.SD, function(z) if (is.numeric(z)) sum(z, na.rm = TRUE) else NA), .SDcols = names(s6)]
tot[, Chromosome := NA_integer_]; tot[, `QC_label7 pass fraction` := round(tot[[3]] / tot[[2]], 5)]
s6 <- rbind(s6, tot); s6[, Chromosome := ifelse(is.na(Chromosome), "Total", as.character(Chromosome))]
fwrite(s6, file.path(OUT, "S6_variant_counts_by_chromosome.csv"))
LOG("S6 written: totals %s", paste(names(tot)[-1], unlist(tot[, -1]), sep = " = ", collapse = "; "))

writeLines(log_lines, file.path(OUT, "build_log_S1_S2_S3_S6_S7.txt"))
LOG("done.")
