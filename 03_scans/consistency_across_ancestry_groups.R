#!/usr/bin/env Rscript
## consistency_across_ancestry_groups.R (2026-10-14): as race_consistency_projcell_qc7.R, but on the 22 single variants that the full
## QC_label7 rescan reports (APOE lead + the 21 non-APOE FDR 5% variants in 11 loci), taken with their locus labels from
## S14c_split_discovery_variants.csv. (a) each variant inside each self-reported group, with a fixed-effect inverse-variance
## meta-analysis of the group estimates and Cochran's Q; (b) every non-APOE variant below 5e-8 inside a group, with its full-sample
## statistics. Writes S12b_race_discovery_variants.csv, S12b_race_meta.csv and S12c_race_gw_hits.csv (supp_v3).
suppressPackageStartupMessages(library(data.table))
OUT <- "${ADSP_ROOT}/output"; TAB <- "${ADSP_ROOT}/revision/tables/supp_v3"
rd <- function(d) { e <- new.env(); load(file.path(OUT, d, "individual_analysis_qc7/results_individual_analysis_genome.Rdata"), envir = e)
  x <- as.data.table(get(ls(e)[1], envir = e)); x[, key := paste(CHR, POS, REF, ALT, sep = ":")]; x }
apoe <- function(x) x$CHR == 19 & x$POS >= 44.4e6 & x$POS <= 46.4e6
groups <- c(White = 5, `Hispanic/Latino` = 6, `Black or African American` = 4, Asian = 2)
L <- fread(file.path(TAB, "S14c_split_discovery_variants.csv"))[, .(key, Locus = `Locus label (Figure 3a)`, P)][order(P)]
disc <- L$key; lab <- L$Locus; cat("variants:", length(disc), "\n")
full <- rd("matched_lpa_pc_pair33_projcell")
F <- full[match(disc, key), .(key, Locus = lab, N_full = N, MAF_full = MAF, Est_full = Est, SE_full = Est_se, P_full = pvalue)]
a <- list(); h <- list()
for (g in names(groups)) { x <- rd(sprintf("matched_lpa_pc_pair33_projcell_race%d", groups[g]))
  d <- data.table(Group = g, Locus = lab); d[, key := disc]; a[[g]] <- d[, c("N", "MAF", "Est", "SE", "P") := x[match(disc, key), .(N, MAF, Est, Est_se, pvalue)]]
  gw <- x[pvalue < 5e-8 & !apoe(x)]
  if (nrow(gw)) h[[g]] <- merge(gw[, .(Group = g, key, CHR, POS, REF, ALT, N, MAF, Est, SE = Est_se, P = pvalue)],
                                full[, .(key, MAF_full = MAF, Est_full = Est, P_full = pvalue)], by = "key", all.x = TRUE)
  rm(x); gc() }
A <- rbindlist(a); A <- merge(A, F, by = c("key", "Locus"))[order(match(key, disc), match(Group, names(groups)))]
meta <- A[!is.na(Est), { w <- 1 / SE^2; b <- sum(w * Est) / sum(w); se <- sqrt(1 / sum(w)); Q <- sum(w * (Est - b)^2)
  .(groups_tested = .N, Est_meta = b, SE_meta = se, P_meta = 2 * pnorm(-abs(b / se)), Q = Q, P_het = pchisq(Q, .N - 1, lower.tail = FALSE)) }, by = .(key, Locus)]
meta <- merge(meta, F[, .(key, Est_full, P_full)], by = "key")[order(match(key, disc))]
print(A[, .(Group, Locus, N, MAF, Est, SE, P, Est_full, P_full)], digits = 3, nrows = 200); print(meta, digits = 3, nrows = 100)
fwrite(A, file.path(TAB, "S12b_race_discovery_variants.csv")); fwrite(meta, file.path(TAB, "S12b_race_meta.csv"))
H <- rbindlist(h); if (nrow(H)) { print(H[order(Group, P)], digits = 3); fwrite(H[order(Group, P)], file.path(TAB, "S12c_race_gw_hits.csv")) }
cat("DONE\n")
