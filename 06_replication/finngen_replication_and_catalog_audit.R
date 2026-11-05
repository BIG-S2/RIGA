#!/usr/bin/env Rscript
## finngen_replication_and_catalog_audit.R (2026-10-14; reads the QC_label7 rescan and its v4 labels) -- from replication_projcell_fdr.R (2026-10-06) -- external support for the stratified-design (D_S)
## single-variant findings and GWAS Catalog audit of every named gene (revision plan Phases 1 & 4).
##   (1) FinnGen R12 G6_ALZHEIMER lookup of the FDR-5% single variants (exact variant, allele-aligned
##       direction) and of the suggestive leads (p < 1e-5, pruned at 500 kb), with a genome-wide
##       direction sign test and the density-matched control-locus enrichment of qc7_fwd_replication_5arm.R.
##   (2) GWAS Catalog (strict AD/dementia, genome-wide significant) within 250 kb of each lead, and a
##       gene-symbol search for the gene x mask sets of Supplementary Table 3 plus PSEN1.
## Outputs: revision/tables/supp_v3/S8_replication_finngen.csv, S8b_replication_summary.csv,
##          S9_catalog_genes.csv, build_log_replication.txt
suppressPackageStartupMessages(library(data.table)); setDTthreads(4)
ROOT <- "${ADSP_ROOT}"
OUTD <- "${ADSP_ROOT}/output/matched_lpa_pc_pair33_projcell"
TAB  <- file.path(ROOT, "revision/tables/supp_v3"); WORK <- file.path(ROOT, "revision/work/replication")
CATF <- "${WORK_DIR}/gwas/gwas-catalog-associations_ontology-annotated.tsv"
WINDOW <- 250000L; APOE <- c(19L, 44400000L, 46400000L); SUG <- 1e-5
set.seed(42)
log_lines <- character(0); LOG <- function(...) { s <- sprintf(...); message(s); log_lines <<- c(log_lines, s) }

## ---- D_S scan ----
G <- as.data.table(get(load(file.path(OUTD, "individual_analysis_qc7", "results_individual_analysis_genome.Rdata"))))
G <- G[is.finite(pvalue) & pvalue > 0]; G[, q := p.adjust(pvalue, "BH")]
G[, `:=`(CHR = as.integer(CHR), POS = as.integer(POS))]
G[, apoe := CHR == APOE[1] & POS >= APOE[2] & POS <= APOE[3]]
LOG("D_S scan: %d tests; FDR 5%%: %d variants (%d in APOE region); p<1e-5: %d", nrow(G), sum(G$q < 0.05), sum(G$q < 0.05 & G$apoe), sum(G$pvalue < SUG))
prune <- function(d, dist = 500000L) {
  d <- d[order(pvalue)]; keep <- logical(nrow(d)); kc <- integer(0); kp <- integer(0)
  for (i in seq_len(nrow(d))) if (!length(kc) || !any(kc == d$CHR[i] & abs(kp - d$POS[i]) < dist)) {
    keep[i] <- TRUE; kc <- c(kc, d$CHR[i]); kp <- c(kp, d$POS[i]) }
  d[keep] }
L <- prune(G[pvalue < SUG]); L[, tier := ifelse(q < 0.05, "FDR 5%", "suggestive (p < 1e-5)")]
L[, id := .I]
LOG("pruned leads: %d (FDR %d, suggestive %d; APOE-region leads %d)", nrow(L), sum(L$tier == "FDR 5%"), sum(L$tier != "FDR 5%"), sum(L$apoe))
lab <- fread(file.path(OUTD, "individual_analysis_qc7", "manhattan_MAC_20_sig_annotated_v4.csv"))
L <- merge(L, lab[, .(CHR, POS, REF, ALT, label)], by = c("CHR", "POS", "REF", "ALT"), all.x = TRUE)
rm(G); invisible(gc())

## ---- FinnGen: exact variant (allele-aligned) and window ----
FG <- readRDS(file.path(WORK, "fg_G6_ALZHEIMER.rds"))[chrom <= 22L]
LOG("FinnGen G6_ALZHEIMER autosomal variants: %d; columns: %s", nrow(FG), paste(names(FG), collapse = ","))
FGk <- FG[, .(CHR = as.integer(chrom), POS = as.integer(pos), ref, alt, fg_p = pval, fg_beta = beta, fg_af = af_alt)]
ex <- merge(L[, .(id, CHR, POS, REF, ALT)], FGk, by = c("CHR", "POS"), allow.cartesian = TRUE)
ex[, aligned := fifelse(REF == ref & ALT == alt, 1, fifelse(REF == alt & ALT == ref, -1, NA_real_))]
ex <- ex[!is.na(aligned)][order(id, fg_p)][!duplicated(id)]
ex[, fg_beta_alt := aligned * fg_beta]
L <- merge(L, ex[, .(id, fg_exact_p = fg_p, fg_beta_alt, fg_af_alt = fifelse(aligned == 1, fg_af, 1 - fg_af))], by = "id", all.x = TRUE)
FGP <- FG[, .(CHR = as.integer(chrom), start = as.integer(pos), end = as.integer(pos), pval)]; setkey(FGP, CHR, start, end)
win_stats <- function(pts) {
  h <- pts[, .(CHR, start = POS - WINDOW, end = POS + WINDOW, id)]; setkey(h, CHR, start, end)
  ov <- foverlaps(FGP, h, type = "within", nomatch = NULL, by.x = c("CHR", "start", "end"), by.y = c("CHR", "start", "end"))
  ov[, .(n_fg = .N, fg_window_best_p = min(pval)), by = id] }
L <- merge(L, win_stats(L[, .(CHR, POS, id)]), by = "id", all.x = TRUE); L[is.na(n_fg), n_fg := 0L]
L[, concordant := sign(Est) == sign(fg_beta_alt)]
## replication set = the discovery signals (FDR 5%, non-APOE); correction over that set only
n_disc <- sum(L$tier == "FDR 5%" & !L$apoe)
LOG("discovery replication set: %d non-APOE FDR signals; Bonferroni threshold 0.05/%d = %.3g", n_disc, n_disc, 0.05 / n_disc)
L[, verdict := fifelse(apoe, "APOE region", fifelse(tier != "FDR 5%", "suggestive lead (direction reported, no claim)",
             fifelse(is.na(fg_exact_p), "variant absent from FinnGen",
             fifelse(concordant & fg_exact_p < 0.05 / n_disc, "replicated (same direction, p < 0.05 over the discovery set)",
             fifelse(concordant & fg_exact_p < 0.05, "nominal support (same direction, p < 0.05)", "not supported")))))]
## density-matched control loci (as in qc7_fwd_replication_5arm.R)
CHRLEN <- c(248956422,242193529,198295559,190214555,181538259,170805979,159345973,145138636,138394717,133797422,135086622,
            133275309,114364328,107043718,101991189,90338345,83257441,80373285,58617616,64444167,46709983,50818468)
ctl <- data.table(CHR = sample(seq_along(CHRLEN), 20000L, TRUE, prob = CHRLEN / sum(CHRLEN)))
ctl[, POS := as.integer(runif(.N, WINDOW + 1, CHRLEN[CHR] - WINDOW))][, id := .I]
CS <- merge(ctl, setnames(win_stats(ctl), "fg_window_best_p", "best_p"), by = "id")[n_fg > 0]
brk <- unique(quantile(CS$n_fg, seq(0, 1, 0.1))); CS[, dens := cut(n_fg, brk, include.lowest = TRUE)]
L[, dens := cut(n_fg, brk, include.lowest = TRUE)]
summ <- list()
for (tr in c("FDR 5%", "suggestive (p < 1e-5)", "all leads")) {
  tier <- tr
  nz <- if (tr == "all leads") L[apoe == FALSE & n_fg > 0] else L[apoe == FALSE & n_fg > 0][tier == tr]
  if (!nrow(nz)) next
  for (thr in c(0.05, 1e-3, 1e-4, 5e-8)) {
    k <- sum(nz$fg_window_best_p < thr)
    cr <- sapply(seq_len(nrow(nz)), function(i) { cc <- CS[dens == nz$dens[i]]; if (!nrow(cc)) NA_real_ else mean(cc$best_p < thr) })
    ek <- sum(cr, na.rm = TRUE)
    summ[[length(summ) + 1]] <- data.table(tier = tier, statistic = sprintf("FinnGen window (250 kb) best p < %g", thr), n_loci = nrow(nz),
      observed = k, expected = round(ek, 2), ratio = round(k / max(ek, 1e-9), 2),
      p_value = if (ek > 0) signif(poisson.test(k, ek, alternative = "greater")$p.value, 3) else NA_real_) }
  cx <- nz[!is.na(concordant)]
  if (nrow(cx)) summ[[length(summ) + 1]] <- data.table(tier = tier, statistic = "direction sign test (exact variant, binomial vs 0.5)",
      n_loci = nrow(cx), observed = sum(cx$concordant), expected = nrow(cx) / 2, ratio = round(mean(cx$concordant), 2),
      p_value = signif(binom.test(sum(cx$concordant), nrow(cx), 0.5, alternative = "greater")$p.value, 3))
}
S8b <- rbindlist(summ); fwrite(S8b, file.path(TAB, "S8b_replication_summary.csv")); print(S8b)

## ---- GWAS Catalog: strict AD/dementia associations near each lead; gene-symbol audit ----
GC <- fread(CATF, sep = "\t", quote = "", fill = TRUE, showProgress = FALSE)
setnames(GC, make.names(names(GC), unique = TRUE))
GC[, chr := suppressWarnings(as.integer(CHR_ID))]; GC[, bp := suppressWarnings(as.integer(CHR_POS))]
GC[, pv := suppressWarnings(as.numeric(P.VALUE))]
tl <- tolower(paste(GC$DISEASE.TRAIT, GC$MAPPED_TRAIT))
isad <- grepl("alzheimer|dementia", tl)
bad  <- grepl("snp x snp|cerebrospinal|csf|p-tau|biomarker|lewy|frontotemporal|polygenic|vascular dementia", tl) | grepl(" x ", GC$SNPS, fixed = TRUE)
GC[, ad_strict := isad & !bad & is.finite(pv) & pv < 5e-8]
GC[, ad_any := isad & !bad]
LOG("GWAS Catalog rows: %d; strict genome-wide AD/dementia associations: %d", nrow(GC), sum(GC$ad_strict))
GCp <- GC[!is.na(chr) & !is.na(bp) & chr <= 22]
near <- function(pts, tab) {
  h <- pts[, .(CHR, start = POS - WINDOW, end = POS + WINDOW, id)]; setkey(h, CHR, start, end)
  t2 <- tab[, .(CHR = chr, start = bp, end = bp, pv, trait = DISEASE.TRAIT, genes = MAPPED_GENE, pmid = PUBMEDID)]; setkey(t2, CHR, start, end)
  ov <- foverlaps(t2, h, type = "within", nomatch = NULL, by.x = c("CHR", "start", "end"), by.y = c("CHR", "start", "end"))
  ov[order(pv), .(n = .N, best_p = min(pv), example = paste(unique(head(paste0(genes, " [", trait, "; PMID ", pmid, "]"), 2)), collapse = " | ")), by = id] }
ca <- near(L[, .(CHR, POS, id)], GCp[ad_strict == TRUE]); cy <- near(L[, .(CHR, POS, id)], GCp)
L <- merge(L, setnames(ca, c("id", "catalog_AD_n", "catalog_AD_best_p", "catalog_AD_example")), by = "id", all.x = TRUE)
L <- merge(L, cy[, .(id, catalog_any_n = n)], by = "id", all.x = TRUE)
L[is.na(catalog_AD_n), catalog_AD_n := 0L][is.na(catalog_any_n), catalog_any_n := 0L]
L[, novelty := fifelse(apoe, "APOE region", fifelse(catalog_AD_n > 0, "AD association within 250 kb in the GWAS Catalog", "no AD association within 250 kb"))]
bg <- { CSi <- copy(CS)[, id := .I]; m <- near(CSi[, .(CHR, POS, id)], GCp[ad_strict == TRUE]); mean(CSi$id %in% m$id) }
LOG("background: %.1f%% of matched control loci carry a strict AD association within %d bp", 100 * bg, WINDOW)
setorder(L, -apoe, pvalue)
S8 <- L[, .(Tier = tier, CHR, POS, REF, ALT, `Locus label` = label, `APOE region` = apoe, MAF, `ADSP P` = signif(pvalue, 3), `ADSP OR (ALT)` = round(exp(Est), 3),
            `FinnGen P (exact variant)` = signif(fg_exact_p, 3), `FinnGen beta (ALT, aligned)` = round(fg_beta_alt, 4), `FinnGen AF (ALT)` = signif(fg_af_alt, 3),
            `Direction concordant` = concordant, `FinnGen variants in window` = n_fg, `FinnGen best P in window` = signif(fg_window_best_p, 3),
            Verdict = verdict, `Catalog AD associations in window` = catalog_AD_n, `Catalog AD best P` = signif(catalog_AD_best_p, 3),
            `Catalog example` = catalog_AD_example, `Catalog associations, any trait` = catalog_any_n, Novelty = novelty)]
fwrite(S8, file.path(TAB, "S8_replication_finngen.csv"))
LOG("S8 written: %d leads; verdicts: %s", nrow(S8), paste(names(table(S8$Verdict)), table(S8$Verdict), sep = " = ", collapse = "; "))

## gene-symbol audit for the gene x mask sets (Supplementary Table 3) and PSEN1
s3 <- fread(file.path(TAB, "S3_fdr_gene_mask_sets.csv"))
genes <- unique(c(s3$`Gene name`, "PSEN1"))
sym <- function(x) toupper(trimws(x))
GC[, mapped := sym(MAPPED_GENE)]; GC[, reported := sym(REPORTED.GENE.S.)]
S9 <- rbindlist(lapply(genes, function(g) {
  pat <- paste0("(^|[^A-Z0-9])", g, "($|[^A-Z0-9])")
  hit <- GC[grepl(pat, mapped) | grepl(pat, reported)]
  data.table(Gene = g, `Masks admitted at FDR 5%` = paste(s3[`Gene name` == g, Category], collapse = ", "),
             `Catalog associations (any trait)` = nrow(hit), `Strict AD/dementia associations` = sum(hit$ad_strict),
             `AD/dementia associations (any p)` = sum(hit$ad_any), `Best AD P` = if (any(hit$ad_strict)) signif(min(hit[ad_strict == TRUE]$pv), 3) else NA_real_,
             `First AD report` = if (any(hit$ad_any)) { h <- hit[ad_any == TRUE][order(pv)][1]; sprintf("%s %s (PMID %s)", h$FIRST.AUTHOR, substr(h$DATE, 1, 4), h$PUBMEDID) } else "",
             `Other traits (examples)` = paste(head(unique(hit[ad_any == FALSE]$DISEASE.TRAIT), 3), collapse = "; "))
}))
fwrite(S9, file.path(TAB, "S9_catalog_genes.csv")); print(S9[, 1:6])
writeLines(log_lines, file.path(TAB, "build_log_replication.txt")); LOG("done.")
