#!/usr/bin/env Rscript
## ukb_single_variant_direct_test.R (2026-10-19): copy of revision/scripts/ukb_sv_direct_reconstructed.R for the root-checked STAAR/STAARpipeline library
## (ukb_rootcheck_20261019; library selected by R_LIBS_USER in jobs/05_downstream.sbatch or jobs/06_tables.sbatch).
## The only code changes are the output paths/names marked '2026-10-19 rootcheck' below; nothing existing is overwritten.
## ukb_sv_direct_reconstructed.R (2026-10-14): direct test in the corrected UK Biobank cohort (saddlepoint score test) of every
## non-APOE single variant admitted at FDR 5% in the QC_label7 rescan (Supplementary Table 2) plus APOE e4. Writes
## revision/tables/review/ukb_reconstructed_sv.csv (replaces the five-variant version).
suppressPackageStartupMessages({ library(data.table); library(SeqArray); library(STAAR); library(STAARpipeline) })
UKB <- "${ADSP_ROOT}/revision/ukb"; TAB <- "${ADSP_ROOT}/revision/tables/review"
load(file.path(UKB, "agds_dir_ukb.RData")); obj <- get(load(file.path(UKB, "ukb_obj_nullmodel_reconstructed.RData"))); obj$sparse_kins <- TRUE; obj$relatedness <- TRUE; obj$use_SPA <- TRUE
s2 <- fread("${ADSP_ROOT}/revision/tables/supp_v3/S2_fdr_single_variants.csv")
tg <- s2[!(CHR == 19 & POS >= 44.4e6 & POS <= 46.4e6), .(gene = `Nearest protein-coding gene`, chr = as.integer(CHR), pos = as.integer(POS), ref = REF, alt = ALT, adsp_p = P, adsp_maf = MAF)]
tg <- rbind(tg, data.table(gene = "APOE", chr = 19L, pos = 44908684L, ref = "T", alt = "C", adsp_p = NA_real_, adsp_maf = NA_real_)); cat("variants:", nrow(tg), "\n")
sv <- list()
for (i in seq_len(nrow(tg))) { v <- tg[i]; g <- seqOpen(agds_dir[v$chr])
  r <- tryCatch(Individual_Analysis(chr = v$chr, start_loc = v$pos - 1L, end_loc = v$pos + 1L, genofile = g, obj_nullmodel = obj, mac_cutoff = 1, QC_label = "annotation/info/QC_label2",
                                    variant_type = "variant", geno_missing_imputation = "mean", SPA_p_filter = TRUE, p_filter_cutoff = 0.05), error = function(e) NULL); seqClose(g)
  d <- if (is.null(r)) NULL else as.data.table(as.data.frame(r))[POS == v$pos & REF == v$ref & ALT == v$alt]
  sv[[i]] <- if (is.null(d) || !nrow(d)) data.table(v, note = "variant absent or filtered") else data.table(v, d[1, .(ALT_AF, MAF, N, pvalue)], note = "")
  cat(sprintf("%2d/%d %-10s chr%d:%d %s\n", i, nrow(tg), v$gene, v$chr, v$pos, if (is.null(d) || !nrow(d)) "absent" else sprintf("MAF %.5f p %.3g", d$MAF[1], d$pvalue[1]))) }
fwrite(rbindlist(sv, fill = TRUE), file.path(TAB, "ukb_rootcheck_sv.csv")); cat("DONE\n")   # 2026-10-19 rootcheck: new file
