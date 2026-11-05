# Merge the 291 ADSP_AD_Individual_Analysis chunks of ONE arm into results_individual_analysis_genome.Rdata
# (the object name every downstream script expects).  merge_chunks() copied from
# adsp_sv_qq_manhattan_pair33_cohort_3arm.R; refuses to overwrite.  Usage: Rscript merge_single_variant_chunks.R <dir>
suppressPackageStartupMessages({ library(data.table) })
dir <- commandArgs(trailingOnly = TRUE)[1]; stopifnot(!is.na(dir), dir.exists(dir))
GF <- file.path(dir, "results_individual_analysis_genome.Rdata")
if (file.exists(GF)) stop("refusing to overwrite ", GF)
fs <- list.files(dir, "^ADSP_AD_Individual_Analysis_\\d+\\.Rdata$", full.names = TRUE)
stopifnot(length(fs) == 291); cat("merging", length(fs), "chunks in", dir, "\n")
results_individual_analysis_genome <- rbindlist(lapply(fs, function(f) { e <- new.env(); load(f, envir = e)
  as.data.table(get(ls(e)[1], envir = e)) }), fill = TRUE, use.names = TRUE)
cat("rows:", format(nrow(results_individual_analysis_genome), big.mark = ","), " cols:", paste(names(results_individual_analysis_genome), collapse = ","), "\n")
save(results_individual_analysis_genome, file = GF); cat("wrote", GF, "\n")
