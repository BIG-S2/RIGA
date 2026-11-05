#!/usr/bin/env Rscript
## diversity_accounting.R (2026-10-13, review point 7): who the stratified design keeps and how much each group
## contributes. Pool = the participants with an AD case/control label, sex, age and a sequencing project (built with the prep
## of src/prepare_sample_and_one_to_one_matching.R, as ancestry_group_null_models.R does); analytic = the 30,894 of the stratified
## design. Information: a participant in cell c contributes pi_c (1 - pi_c), pi_c the cell case fraction (the per-subject weight
## of the within-cell score variance for a variant of equal frequency); the share is the group total over the overall total.
## Dimensions: genetic ancestry (1000 Genomes super-population projection, adsp_ref_ancestry.tsv), self-reported race/ethnicity,
## age band (age at onset for cases, age at last examination for controls), contributing cohort.
## Writes revision/tables/supp_v3/S1h_diversity_accounting.csv.
suppressPackageStartupMessages({ library(data.table) })
SRC <- "${ADSP_ROOT}/src"; OUT <- "${ADSP_ROOT}/output"
src <- readLines(file.path(SRC, "prepare_sample_and_one_to_one_matching.R")); cut <- grep("^m\\.out <- matchit", src)[1]
eval(parse(text = paste(src[1:(cut - 1)], collapse = "\n"))); pool <- as.data.table(dfx)
qc <- fread("${ADSP_ROOT}/pheno/gcad.r5.wgs.58507.2024.11.03.quality.metrics.ALL.csv", select = c("Sample_name", "Study"), colClasses = "character")
setnames(qc, "Sample_name", "IID"); qc <- unique(qc, by = "IID"); pool <- merge(pool, qc, by = "IID", all.x = TRUE)
pool[is.na(Study) | Study == "", Study := "Unknown"]; pool[, cell_SC := paste(Study, SeqCenter, sep = "|")]
cat("pool:", nrow(pool), "\n")
an <- fread(file.path(OUT, "matched/matched_samples_lpa_pc.pair33_projcell.csv"), select = c("IID", "cell_SC", "AD"), colClasses = list(character = "IID"))
pool[, analytic := IID %in% an$IID]; cat("analytic in pool:", sum(pool$analytic), "of", nrow(an), "\n")
cf <- an[, .(pi = mean(AD)), by = cell_SC]; pool <- merge(pool, cf, by = "cell_SC", all.x = TRUE)
pool[, w := fifelse(analytic, pi * (1 - pi), 0)]
anc <- fread("${ADSP_ROOT}/adsp_ref_ancestry.tsv", select = c("IID", "SuperPop")); pool <- merge(pool, anc, by = "IID", all.x = TRUE)
race_lab <- c(`1` = "American Indian/Alaska Native", `2` = "Asian", `3` = "Native Hawaiian/Pacific Islander", `4` = "Black or African American", `5` = "White", `6` = "Hispanic/Latino")
pool[, Race_lab := fifelse(as.character(Race) %in% names(race_lab), race_lab[as.character(Race)], "Other/unknown")]
pool[, AgeN := as.numeric(gsub("\\+", "", Age))]
pool[, AgeBand := cut(AgeN, c(-Inf, 64, 74, 84, Inf), labels = c("< 65", "65-74", "75-84", ">= 85"))]
tab <- function(dim, col) pool[, .(n_pool = .N, n_analytic = sum(analytic), info = sum(w)), by = .(group = get(col))][, dimension := dim]
D <- rbind(tab("Genetic ancestry (1000 Genomes super-population)", "SuperPop"), tab("Self-reported race/ethnicity", "Race_lab"),
           tab("Age band", "AgeBand"), tab("Contributing cohort", "Cohort"))
D[, `:=`(pct_retained = round(100 * n_analytic / n_pool, 1), pool_share = round(100 * n_pool / nrow(pool), 2),
         analytic_share = round(100 * n_analytic / sum(pool$analytic), 2), info_share = round(100 * info / sum(pool$w), 2))]
setcolorder(D, c("dimension", "group", "n_pool", "n_analytic", "pct_retained", "pool_share", "analytic_share", "info_share")); D[, info := NULL]
D <- D[order(dimension, -n_pool)]
fwrite(D, "${ADSP_ROOT}/revision/tables/supp_v3/S1h_diversity_accounting.csv")
print(D[dimension != "Contributing cohort"]); cat("cohorts with zero retained:", D[dimension == "Contributing cohort" & n_analytic == 0, paste(group, collapse = "; ")], "\n")
