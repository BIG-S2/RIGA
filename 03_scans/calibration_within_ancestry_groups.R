#!/usr/bin/env Rscript
## *_qc7.R (2026-10-14): same as the original but reading the QC_label7 rescans (individual_analysis_qc7/) instead of the QC2-scan + QC7-filter tables.
## lambda_projcell_race.R <race code> (2026-10-07): genomic-control lambda of the per-ancestry stratified scan
## (QC_label7, MAC >= 50) overall and by minor-allele-frequency bin; appends one block of rows to
## revision/tables/supp_v3/S12_lambda_by_ancestry.csv. Race codes: 5 White, 6 Hispanic/Latino, 4 Black, 2 Asian.
suppressPackageStartupMessages(library(data.table))
r <- commandArgs(TRUE)[1]; lab <- c(`5` = "White", `6` = "Hispanic/Latino", `4` = "Black or African American", `2` = "Asian", all = "All (stratified design)")[r]
f <- if (r == "all") "${ADSP_ROOT}/output/matched_lpa_pc_pair33_projcell/individual_analysis_qc7/results_individual_analysis_genome.Rdata" else
     sprintf("${ADSP_ROOT}/output/matched_lpa_pc_pair33_projcell_race%s/individual_analysis_qc7/results_individual_analysis_genome.Rdata", r)
g <- as.data.table(get(load(f)))[is.finite(pvalue) & pvalue > 0]
n_sub <- if (r == "all") 30894L else nrow(fread(sprintf("${ADSP_ROOT}/output/matched/matched_samples_projcell_race%s.csv", r), select = "IID"))
lam <- function(p) qchisq(median(p), 1, lower.tail = FALSE) / qchisq(0.5, 1)
g[, bin := cut(MAF, c(0, 0.01, 0.05, 0.5), labels = c("MAF < 1%", "1-5%", ">= 5%"), include.lowest = TRUE)]
apoe <- g$CHR == 19 & g$POS >= 44.4e6 & g$POS <= 46.4e6
rows1 <- data.table(Group = lab, N = n_sub, Bin = "all variants (MAC >= 50)", Tests = nrow(g), lambda_GC = round(lam(g$pvalue), 4),
                    `n < 5e-8 outside APOE` = sum(g$pvalue < 5e-8 & !apoe), `n < 1e-5 outside APOE` = sum(g$pvalue < 1e-5 & !apoe))
rows2 <- rbindlist(lapply(levels(g$bin), function(bb) { x <- g[bin == bb]
  data.table(Group = lab, N = n_sub, Bin = bb, Tests = nrow(x), lambda_GC = round(lam(x$pvalue), 4), `n < 5e-8 outside APOE` = NA_integer_, `n < 1e-5 outside APOE` = NA_integer_) }))
rows <- rbind(rows1, rows2)
out <- "${ADSP_ROOT}/revision/tables/supp_v3/S12_lambda_by_ancestry.csv"
old <- if (file.exists(out)) fread(out)[Group != lab] else data.table()
fwrite(rbind(old, rows, fill = TRUE), out); print(rows)
