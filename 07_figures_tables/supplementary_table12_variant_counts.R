#!/usr/bin/env Rscript
## supplementary_table12_variant_counts.R (2026-10-14): replace the single-variant test counts of Supplementary Table 6 (per chromosome and Total)
## with those of the QC_label7 rescan. The builder's merge failed because Table 6 carries a "Total" row (character Chromosome).
suppressPackageStartupMessages(library(data.table))
OUT <- "${ADSP_ROOT}/revision/tables/supp_v3"
g <- as.data.table(get(load("${ADSP_ROOT}/output/matched_lpa_pc_pair33_projcell/individual_analysis_qc7/results_individual_analysis_genome.Rdata")))
cnt <- g[is.finite(pvalue) & pvalue > 0, .N, by = .(Chromosome = as.character(CHR))]; rm(g)
s6 <- fread(file.path(OUT, "S6_variant_counts_by_chromosome.csv"), colClasses = list(character = "Chromosome")); cc <- grep("^Single-variant tests", names(s6), value = TRUE)[1]
s6[, (cc) := cnt$N[match(Chromosome, cnt$Chromosome)]]; s6[Chromosome == "Total", (cc) := sum(cnt$N)]
fwrite(s6, file.path(OUT, "S6_variant_counts_by_chromosome.csv")); print(s6[, c("Chromosome", cc), with = FALSE]); cat("DONE\n")
