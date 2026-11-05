#!/usr/bin/env Rscript
## single_variant_saddlepoint_retest.R (2026-10-13): the 102 FDR-significant single variants of the stratified design (Supplementary
## Table 2) re-tested on the deployed null with the normal-approximation score test and with the saddlepoint score test, so the
## two tests can be compared on the discoveries themselves (the permuted-null panel has only 93 sites below MAF 2%).
## Read-only on inputs. Writes revision/tables/review/sv_fdr_spa_retest.csv.
suppressPackageStartupMessages({ library(data.table); library(SeqArray); library(STAAR); library(STAARpipeline) })
OUT <- "${ADSP_ROOT}/output"; TAB <- "${ADSP_ROOT}/revision/tables/review"
obj <- get(load(file.path(OUT, "matched/obj_null_stratum_all_matched_lpa_pc.pair33_projcell.RData"))); obj$sparse_kins <- TRUE; obj$relatedness <- TRUE
agds_dir <- get(load(file.path(OUT, "agds_dir.Rdata"))); QC <- "annotation/info/QC_label7"
s2 <- fread("${ADSP_ROOT}/revision/tables/supp_v3/S2_fdr_single_variants.csv")[, .(CHR = as.integer(CHR), POS = as.integer(POS), REF, ALT, P_published = P, MAF_published = MAF)]
run <- function(chr, from, to, spa) { obj$use_SPA <- spa; g <- seqOpen(agds_dir[chr]); on.exit(seqClose(g))
  r <- if (spa) Individual_Analysis(chr = chr, start_loc = from, end_loc = to, genofile = g, obj_nullmodel = obj, mac_cutoff = 50, QC_label = QC,
                                    variant_type = "variant", geno_missing_imputation = "mean", SPA_p_filter = TRUE, p_filter_cutoff = 0.05)
       else Individual_Analysis(chr = chr, start_loc = from, end_loc = to, genofile = g, obj_nullmodel = obj, mac_cutoff = 50, QC_label = QC,
                                variant_type = "variant", geno_missing_imputation = "mean")
  as.data.table(as.data.frame(r)) }
out <- list()
for (k in seq_len(nrow(s2))) { v <- s2[k]; chr <- v$CHR; from <- v$POS - 1L; to <- v$POS + 1L   # v2: one window per variant (the chr19 window timed out)
  t0 <- proc.time()[3]
  a <- run(chr, from, to, FALSE); b <- run(chr, from, to, TRUE)
  setnames(a, "pvalue", "p_normal"); setnames(b, "pvalue", "p_spa")
  m <- merge(a[, .(CHR, POS, REF, ALT, ALT_AF, MAF, N, p_normal, Score, Score_se = if ("Score_se" %in% names(a)) Score_se else NA_real_)],
             b[, .(CHR, POS, REF, ALT, p_spa)], by = c("CHR", "POS", "REF", "ALT"), all.x = TRUE)
  m <- merge(v, m, by = c("CHR", "POS", "REF", "ALT"), all.x = TRUE); out[[length(out) + 1]] <- m
  cat(sprintf("%d/%d chr%d:%d matched %d (%.1f s)\n", k, nrow(s2), chr, v$POS, sum(is.finite(m$p_normal)), proc.time()[3] - t0)) }
res <- rbindlist(out); res[, MAC := round(2 * N * MAF)]; setorder(res, p_normal)
fwrite(res, file.path(TAB, "sv_fdr_spa_retest.csv"))
cat(sprintf("\n%d variants; normal p reproduces published p (max |log10 ratio|) %.3f\n", nrow(res), max(abs(log10(res$p_normal / res$P_published)), na.rm = TRUE)))
print(res[CHR != 19 | MAC < 500, .(CHR, POS, MAC, MAF = signif(MAF, 3), p_normal = signif(p_normal, 3), p_spa = signif(p_spa, 3))])
cat(sprintf("FDR-significant at the published boundary (3.9e-7): normal %d, SPA %d; rare (MAC<300): normal %d, SPA %d\n",
            sum(res$p_normal < 3.9e-7, na.rm = TRUE), sum(res$p_spa < 3.9e-7, na.rm = TRUE), sum(res$MAC < 300 & res$p_normal < 3.9e-7, na.rm = TRUE), sum(res$MAC < 300 & res$p_spa < 3.9e-7, na.rm = TRUE)))
cat("DONE\n")
