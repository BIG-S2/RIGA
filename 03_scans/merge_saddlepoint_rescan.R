#!/usr/bin/env Rscript
## merge_saddlepoint_rescan.R (2026-10-13, review point 1): merge the saddlepoint gene-centric rescans of the stratified design
## (coding_spa/, noncoding_spa/, ncRNA_spa/), apply the paper's rule (Benjamini-Hochberg FDR 5% within each family) to the
## primary test (saddlepoint burden, STAAR-B), and compare with the deployed STAAR-O discoveries (normal approximation).
## Families: coding (all five masks), noncoding (seven masks), ncRNA. Every set with >= 2 rare variants is tested, as deployed.
## Writes revision/tables/review/spa_rescan_{coding,noncoding,ncRNA}.csv.gz and spa_rescan_fdr_sets.csv.
suppressPackageStartupMessages(library(data.table))
D <- "${ADSP_ROOT}/output/matched_lpa_pc_pair33_projcell"; TAB <- "${ADSP_ROOT}/revision/tables/review"
flat <- function(x) {   # one result object (list of matrices or a matrix of mode list) -> data.table
  if (is.null(x) || !length(x)) return(NULL)
  if (is.matrix(x) || is.data.frame(x)) x <- list(x)
  rbindlist(lapply(x, function(m) { if (is.null(m) || !length(m)) return(NULL); d <- as.data.frame(m, stringsAsFactors = FALSE)
    as.data.table(lapply(d, function(v) unlist(v))) }), fill = TRUE) }
read_dir <- function(sub, pat) {
  fs <- list.files(file.path(D, sub), pattern = pat, full.names = TRUE)
  cat(sub, ":", length(fs), "chunk files\n")
  rbindlist(lapply(fs, function(f) { e <- new.env(); load(f, envir = e); flat(get(ls(e)[1], envir = e)) }), fill = TRUE) }
fam <- list(coding = c("coding_spa", "^ADSP_AD_Coding_\\d+\\.Rdata$", "coding", "^ADSP_AD_Coding_\\d+\\.Rdata$"),
            noncoding = c("noncoding_spa", "^ADSP_AD_Noncoding_\\d+\\.Rdata$", "noncoding", "^ADSP_AD_Noncoding_\\d+\\.Rdata$"),
            ncRNA = c("ncRNA_spa", "^ADSP_AD_ncRNA_\\d+\\.Rdata$", "ncRNA", "^ADSP_AD_ncRNA_\\d+\\.Rdata$"))
sig <- list(); FAMS <- strsplit(Sys.getenv("FAMILIES", "coding,noncoding,ncRNA"), ",")[[1]]
for (f in intersect(names(fam), FAMS)) {
  a <- read_dir(fam[[f]][1], fam[[f]][2]); setnames(a, make.names(names(a), unique = TRUE))
  a[, `:=`(gene = Gene.name, cat = Category, p_spa = as.numeric(STAAR.B), nSNV = as.numeric(X.SNV), cMAC = as.numeric(cMAC))]
  a <- unique(a[is.finite(p_spa), .(gene, Chr, cat, nSNV, cMAC, p_spa, burden_1_25 = as.numeric(Burden.1.25.))], by = c("gene", "cat"))
  o <- read_dir(fam[[f]][3], fam[[f]][4]); setnames(o, make.names(names(o), unique = TRUE))
  o <- unique(o[, .(gene = Gene.name, cat = Category, p_staarO = as.numeric(STAAR.O), p_acatv = as.numeric(ACAT.V.1.25.))], by = c("gene", "cat"))
  m <- merge(a, o, by = c("gene", "cat"), all = TRUE)
  m[, q_spa := p.adjust(p_spa, "BH")]; m[, q_staarO := p.adjust(p_staarO, "BH")]
  cat(sprintf("%s: %d sets with SPA p; %d with STAAR-O p; FDR5%% SPA burden %d; FDR5%% STAAR-O %d; both %d\n", f, sum(is.finite(m$p_spa)), sum(is.finite(m$p_staarO)),
              sum(m$q_spa < 0.05, na.rm = TRUE), sum(m$q_staarO < 0.05, na.rm = TRUE), sum(m$q_spa < 0.05 & m$q_staarO < 0.05, na.rm = TRUE)))
  syn <- m[cat == "synonymous"]; if (nrow(syn)) cat(sprintf("  synonymous: %d sets, SPA FDR5%% %d, SPA p<1e-4 %d (expected %.1f)\n", nrow(syn), sum(syn$q_spa < 0.05, na.rm = TRUE), sum(syn$p_spa < 1e-4, na.rm = TRUE), 1e-4 * sum(is.finite(syn$p_spa))))
  lam <- function(p) { p <- p[is.finite(p) & p > 0]; median(qchisq(p, 1, lower.tail = FALSE)) / qchisq(0.5, 1) }
  cat(sprintf("  lambda SPA burden %.3f | STAAR-O %.3f\n", lam(m$p_spa), lam(m$p_staarO)))
  fwrite(m, file.path(TAB, sprintf("spa_rescan_%s.csv.gz", f)))
  s <- m[(q_spa < 0.05) | (q_staarO < 0.05)][, family := f]; setorder(s, p_spa); sig[[f]] <- s
  print(s[, .(gene, cat, nSNV, cMAC, p_spa = signif(p_spa, 3), q_spa = signif(q_spa, 3), p_staarO = signif(p_staarO, 3), q_staarO = signif(q_staarO, 3), p_acatv = signif(p_acatv, 3))]) }
fwrite(rbindlist(sig, fill = TRUE), file.path(TAB, sprintf("spa_rescan_fdr_sets_%s.csv", paste(FAMS, collapse = "_")))); cat("DONE\n")
