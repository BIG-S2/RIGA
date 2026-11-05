# ukb_genome_wide_coding_scan.R (2026-10-19): copy of revision/scripts/ukb_gw_coding_spa_reconstructed.R for the
# root-checked STAAR/STAARpipeline library (ukb_rootcheck_20261019). The ONLY code change is the output directory
# (see the marked line); R_LIBS_USER, log names and resources are set in jobs/03_gw_scan.sbatch.
# Genome-wide gene-centric CODING scan in UK Biobank, run separately for the RIGA and
# non-RIGA arms. Direct port of src/STAARpipeline_Gene_Centric_Coding.lpa_pc.r.
#
# WHY GENOME-WIDE. The targeted 427-test replication cannot answer the RIGA-vs-no-RIGA
# question: the S7 calibration metrics (synonymous-mask O/E at 1e-4, count of study-wide
# significant synonymous genes) require a genome-wide null distribution, and a deliberately
# enriched 38-gene list has none. This scan reproduces the ADSP S7 experiment in UKB.
#
# WHY use_SPA = TRUE here. The companion non-SPA run (ukb_gw_coding.R) produces STAAR-O, which
# is what the ADSP S7 metrics use -- but at UKB's 1:60-1:86 case:control imbalance the non-SPA
# omnibus is severely miscalibrated in BOTH arms (observed synonymous O/E@1e-4 of 127-205 and
# bulk lambda 0.24-0.34: a conservative bulk with an explosive tail). Absolute calibration
# therefore cannot be read off that run.
#
# SPA gives a burden-only statistic (STAAR-B) and no STAAR-O. That is fine for this purpose:
# the negative-control argument does not depend on which statistic is used -- synonymous rare-
# variant burden should be null whatever tests it -- so synonymous STAAR-B behaviour is a valid
# and, at this imbalance, a TRUSTWORTHY calibration metric. This run is the primary; the non-SPA
# run is kept only for comparability with the ADSP STAAR-O numbers.
#
# Usage:  Rscript ukb_gw_coding.R <arrayid> <arm>      arm = noRIGA | RIGA

suppressPackageStartupMessages({
  library(gdsfmt); library(SeqArray); library(SeqVarTools)
  library(STAAR);  library(STAARpipeline)
})

REV <- "${ADSP_ROOT}/revision"
UKB <- file.path(REV, "ukb")

args    <- commandArgs(TRUE)
arrayid <- as.numeric(args[1])
ARM     <- "reconstructed"   ## 2026-10-13: the corrected full cohort (ukb_01c_reconstruct_null.R)
stopifnot(ARM == "reconstructed")

NULLFILE <- file.path(UKB, "ukb_obj_nullmodel_reconstructed.RData"); if (FALSE) switch(ARM,
  noRIGA = file.path(UKB, "ukb_obj_nullmodel.RData"),
  RIGA   = file.path(UKB, "ukb_obj_nullmodel_RIGA.RData"))
## 2026-10-19 rootcheck: output directory only (default gw_coding_spa/rootcheck; the validation job overrides it
## through UKB_ROOTCHECK_OUTDIR). The root-checked library is selected by R_LIBS_USER in the sbatch, not here.
output_path <- Sys.getenv("UKB_ROOTCHECK_OUTDIR", file.path(UKB, "gw_coding_spa", "rootcheck"))
dir.create(output_path, showWarnings = FALSE, recursive = TRUE)

load(file.path(UKB, "agds_dir_ukb.RData"))                                  # agds_dir
Annotation_name_catalog <- get(load("${ADSP_ROOT}/output/Annotation_name_catalog.Rdata"))

nenv <- new.env(); load(NULLFILE, envir = nenv)
obj_nullmodel <- get(if ("obj_nullmodel" %in% ls(nenv)) "obj_nullmodel" else ls(nenv)[1], envir = nenv)
obj_nullmodel$sparse_kins  <- TRUE
obj_nullmodel$use_SPA      <- TRUE       # SPA burden: the correctly calibrated statistic here
obj_nullmodel$relatedness  <- TRUE
cat("arm:", ARM, " null N =", length(obj_nullmodel$id_include), "\n")

QC_label                <- "annotation/info/QC_label2"
variant_type            <- "SNV"
geno_missing_imputation <- "mean"
Annotation_dir          <- "annotation/info/FunctionalAnnotation"
Use_annotation_weights  <- TRUE
Annotation_name <- c("CADD","LINSIGHT","FATHMM.XF","aPC.EpigeneticActive",
                     "aPC.EpigeneticRepressed","aPC.EpigeneticTranscription",
                     "aPC.Conservation","aPC.LocalDiversity","aPC.Mappability",
                     "aPC.TF","aPC.Protein")
output_file_name <- "UKB_AD_Coding"

## ---- chunking ---------------------------------------------------------------------
# 20 genes per task, not ADSP's 50: UKB is 5-11x larger than ADSP, so a 50-gene task would
# risk exceeding the wall clock. More, smaller tasks also make partial results usable sooner.
gene_num_in_array  <- 20
group.num.allchr   <- ceiling(table(genes_info[, 2]) / gene_num_in_array)
chr                <- which.max(arrayid <= cumsum(group.num.allchr))
group.num          <- group.num.allchr[chr]
groupid            <- if (chr == 1) arrayid else arrayid - cumsum(group.num.allchr)[chr - 1]

genes_info_chr <- genes_info[genes_info[, 2] == chr, ]
sub_seq_num    <- dim(genes_info_chr)[1]
sub_seq_id     <- if (groupid < group.num)
                    ((groupid - 1) * gene_num_in_array + 1):(groupid * gene_num_in_array) else
                    ((groupid - 1) * gene_num_in_array + 1):sub_seq_num
cat("arrayid", arrayid, "-> chr", chr, " genes", min(sub_seq_id), "-", max(sub_seq_id), "\n")

outfile <- file.path(output_path, paste0(output_file_name, "_", arrayid, ".Rdata"))
if (file.exists(outfile)) { cat("already done:", outfile, "\n"); quit(save = "no") }

genofile <- seqOpen(agds_dir[chr])
results_coding <- c()
t0 <- Sys.time()
for (kk in sub_seq_id) {
  gene_name <- genes_info_chr[kk, 1]
  r <- tryCatch(
    Gene_Centric_Coding(chr = chr, gene_name = gene_name, genofile = genofile,
                        obj_nullmodel = obj_nullmodel,
                        rare_maf_cutoff = 0.01, rv_num_cutoff = 2,
                        QC_label = QC_label, variant_type = variant_type,
                        geno_missing_imputation = geno_missing_imputation,
                        Annotation_dir = Annotation_dir,
                        Annotation_name_catalog = Annotation_name_catalog,
                        Use_annotation_weights = Use_annotation_weights,
                        Annotation_name = Annotation_name),
    # A gene failing (too few rare variants, degenerate genotype matrix) must not kill the
    # whole chunk -- ADSP's driver hard-codes gene exclusions instead, which does not
    # transfer to a different cohort's variant content.
    error = function(e) { cat("  skip", gene_name, ":", conditionMessage(e), "\n"); NULL })
  if (!is.null(r)) results_coding <- append(results_coding, r)
}
cat("chunk done in", round(difftime(Sys.time(), t0, units = "mins"), 1), "min;",
    length(results_coding), "result blocks\n")

save(results_coding, file = outfile)
seqClose(genofile)
cat("wrote", outfile, "\nDONE\n")
