#!/usr/bin/env Rscript
## supplementary_table8_gene_sets.R (2026-10-13): Supplementary Tables S3 (gene x mask sets at FDR 5%) and the gene-set rows of S7
## (threshold tiers) for the primary statistic SB on sets with cMAC >= 10, from the tables of gene_set_primary_statistic.R.
## The single-variant row of S7 is carried over from the existing S7 until the QC_label7 rescan is summarized.
suppressPackageStartupMessages(library(data.table))
TAB <- "${ADSP_ROOT}/revision/tables"; REV <- file.path(TAB, "review"); OUT <- file.path(TAB, "supp_v3")
fam_label <- c(coding = "Gene-centric coding", noncoding = "Gene-centric noncoding", ncRNA = "Gene-centric ncRNA"); fixed_alpha <- c(coding = 5e-7, noncoding = 3.57e-7, ncRNA = 2.5e-6)
s3 <- list(); tiers <- list()
for (f in names(fam_label)) {
  d <- fread(file.path(REV, sprintf("genesets_primary_%s.csv.gz", f)))[tested == TRUE & is.finite(SB) & SB > 0]
  m <- nrow(d); k <- sum(d$q_SB < 0.05); hit <- d[q_SB < 0.05][order(SB)]
  lam <- qchisq(median(d$SB), df = 1, lower.tail = FALSE) / qchisq(0.5, df = 1)
  cat(sprintf("%s: %d sets with cMAC >= 10, lambda %.3f, FDR 5%% -> %d (critical %.3g, largest p %.3g)\n", f, m, lam, k, 0.05 * k / m, if (k) max(hit$SB) else NA))
  tiers[[f]] <- data.table(Family = fam_label[f], `Tests (m)` = m, Masks = uniqueN(d$Category), lambda_GC = round(lam, 4), `Fixed threshold` = fixed_alpha[f], `Admitted at fixed` = sum(d$SB < fixed_alpha[f]),
    `Bonferroni 0.05/m` = 0.05 / m, `Admitted at Bonferroni` = sum(d$SB < 0.05 / m), `Admitted at FDR 5% (k)` = k, `FDR critical value 0.05 k/m` = 0.05 * k / m, `Largest admitted p` = if (k) max(hit$SB) else NA_real_)
  if (k) s3[[f]] <- data.table(Analysis = fam_label[f], hit[, .(`Gene name`, Chr, Category, `#SNV`, cMAC, `Burden(1,25)`, `SKAT(1,25)`, `ACAT-V(1,25)`, `STAAR-B(1,25)`, `STAAR-S(1,25)`, `STAAR-A(1,25)`, `ACAT-O`, `STAAR-O`, SB, `BH q value` = q_SB)]) }
S3 <- rbindlist(s3); setorder(S3, Analysis, SB); fwrite(S3, file.path(OUT, "S3_fdr_gene_mask_sets.csv")); cat("S3:", nrow(S3), "rows\n")
old <- fread(file.path(OUT, "S7_threshold_tiers.csv")); sv <- old[!Family %in% fam_label]
S7 <- rbind(sv, rbindlist(tiers), fill = TRUE); fwrite(S7, file.path(OUT, "S7_threshold_tiers.csv")); print(S7[, .(Family, `Tests (m)`, `Admitted at FDR 5% (k)`, `FDR critical value 0.05 k/m`)]); cat("DONE\n")
