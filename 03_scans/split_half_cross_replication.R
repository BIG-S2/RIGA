#!/usr/bin/env Rscript
## *_qc7.R (2026-10-14): same as the original but reading the QC_label7 rescans (individual_analysis_qc7/) instead of the QC2-scan + QC7-filter tables.
## crosstest_projcell_split.R (2026-10-08): internal replication of the stratified single-variant scan across the two
## held-out halves. Discovery in one half = BH FDR 5% over its QC_label7 variants (MAC >= 50, as in the main scan), pruned
## to the strongest variant per 500-kb window outside the APOE region (chr19:44.4-46.4 Mb, one lead); replication in the
## other half = same direction and P < 0.05 / (number of non-APOE discovery leads). The 102 full-sample FDR variants
## (Supplementary Table 2) are also reported in both halves (descriptive, sign concordance). Writes S13 tables.
suppressPackageStartupMessages({ library(data.table) })
ROOT <- "${ADSP_ROOT}"; OUT <- "${ADSP_ROOT}/output"; TAB <- file.path(ROOT, "revision/tables/supp_v3")
rd <- function(h) { e <- new.env(); load(sprintf("%s/matched_lpa_pc_pair33_projcell_split%s/individual_analysis_qc7/results_individual_analysis_genome.Rdata", OUT, h), envir = e)
  x <- as.data.table(get(ls(e)[1], envir = e)); x[, key := paste(CHR, POS, REF, ALT, sep = ":")]; x }
A <- rd("A"); B <- rd("B"); cat("tested variants A", nrow(A), " B", nrow(B), "\n")
lam <- function(x) { c2 <- qchisq(x$pvalue, 1, lower.tail = FALSE); median(c2, na.rm = TRUE) / qchisq(0.5, 1) }
apoe <- function(x) x$CHR == 19 & x$POS >= 44.4e6 & x$POS <= 46.4e6
prune <- function(d) { d <- d[order(pvalue)]; keep <- logical(nrow(d)); for (i in seq_len(nrow(d))) if (!any(keep & d$CHR == d$CHR[i] & abs(d$POS - d$POS[i]) <= 5e5)) keep[i] <- TRUE; d[keep] }
one_dir <- function(D, R, nd, nr) {
  D[, q := p.adjust(pvalue, "BH")]; disc <- D[q < 0.05]; disc[, apoe := apoe(disc)]
  leads <- rbind(prune(disc[apoe == FALSE]), disc[apoe == TRUE][which.min(pvalue)], fill = TRUE); leads[, apoe := apoe(leads)]
  n_disc <- sum(!leads$apoe); thr <- if (n_disc > 0) 0.05 / n_disc else NA_real_
  cat(sprintf("%s -> %s: FDR variants %d, leads %d (non-APOE %d), replication threshold %s\n", nd, nr, nrow(disc), nrow(leads), n_disc, format(thr, digits = 2)))
  r <- R[leads[, .(key)], on = "key", nomatch = NULL]
  z <- merge(leads[, .(key, CHR, POS, REF, ALT, MAF, N, Est, pvalue, q, apoe)], r[, .(key, MAF_rep = MAF, N_rep = N, Est_rep = Est, pvalue_rep = pvalue)], by = "key", all.x = TRUE)
  z[, same_direction := sign(Est) == sign(Est_rep)]
  z[, verdict := fifelse(apoe, "APOE (reference)", fifelse(is.na(pvalue_rep), "not tested in the other half",
                 fifelse(same_direction & pvalue_rep < thr, sprintf("replicated (same direction, P < 0.05/%d)", n_disc),
                 fifelse(same_direction & pvalue_rep < 0.05, "same direction, nominal only", "not replicated"))))]
  z[, `:=`(discovery_half = nd, replication_half = nr)]; z[order(pvalue)] }
cross <- rbind(one_dir(copy(A), B, "A", "B"), one_dir(copy(B), A, "B", "A"))
print(cross[, .(discovery_half, key, MAF, Est, pvalue, Est_rep, pvalue_rep, verdict)])
fwrite(cross, file.path(TAB, "S14b_split_crosstest.csv"))
## the 102 full-sample FDR variants in both halves
s2 <- fread(file.path(TAB, "S2_fdr_single_variants.csv")); s2[, key := paste(CHR, POS, REF, ALT, sep = ":")]
f <- merge(s2[, .(key, `Locus label (Figure 3a)`, `Beta (log-odds)`, P)], A[, .(key, Est_A = Est, P_A = pvalue, N_A = N)], by = "key", all.x = TRUE)
f <- merge(f, B[, .(key, Est_B = Est, P_B = pvalue, N_B = N)], by = "key", all.x = TRUE)
f[, concordant := sign(Est_A) == sign(Est_B)]; f[, both_with_full := sign(Est_A) == sign(`Beta (log-odds)`) & sign(Est_B) == sign(`Beta (log-odds)`)]
f <- merge(f, s2[, .(key, CHR, POS)], by = "key"); f[, apoe := CHR == 19 & POS >= 44.4e6 & POS <= 46.4e6]   # APOE region by position
cat(sprintf("full-sample FDR variants: %d; tested in both halves %d; A/B sign concordance %d/%d (non-APOE %d/%d); non-APOE with P<0.05 in A: %d, in B: %d\n",
            nrow(f), sum(!is.na(f$Est_A) & !is.na(f$Est_B)), sum(f$concordant, na.rm = TRUE), sum(!is.na(f$concordant)),
            sum(f$concordant & !f$apoe, na.rm = TRUE), sum(!is.na(f$concordant) & !f$apoe), sum(f$P_A < 0.05 & !f$apoe, na.rm = TRUE), sum(f$P_B < 0.05 & !f$apoe, na.rm = TRUE)))
print(f[apoe == FALSE, .(key, `Locus label (Figure 3a)`, `Beta (log-odds)`, P, Est_A, P_A, Est_B, P_B)])
fwrite(f[order(P)], file.path(TAB, "S14_split_fullsample_fdr.csv"))
fwrite(f[apoe == FALSE | key == "19:44909967:TGG:T"][order(P)], file.path(TAB, "S14c_split_discovery_variants.csv"))
sm <- data.table(half = c("A", "B"), tested = c(nrow(A), nrow(B)), N = c(A$N[1], B$N[1]), lambda = c(lam(A), lam(B)),
                 fdr_variants = c(sum(p.adjust(A$pvalue, "BH") < 0.05), sum(p.adjust(B$pvalue, "BH") < 0.05)),
                 fdr_non_apoe = c(sum(p.adjust(A$pvalue, "BH") < 0.05 & !apoe(A)), sum(p.adjust(B$pvalue, "BH") < 0.05 & !apoe(B))))
print(sm); fwrite(sm, file.path(TAB, "S14_split_scan_summary.csv"))
sp <- fread(file.path(TAB, "S14_split_summary.csv")); fwrite(merge(sp, sm[, .(half, tested, lambda, fdr_variants, fdr_non_apoe)], by = "half"), file.path(TAB, "S14a_split_summary.csv")); cat("DONE\n")
