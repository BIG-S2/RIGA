# GENERIC version (arms, dirs, tag on the command line). Post-filter single-variant arms to the corrected QC label (QC_label7 =
# per-variant AN>3636 & HWE_SLP_I > -4 & F_MISSING < 0.05; equivalent to the PRIMARY
# QC_label5 per REPORT_pair33_primary_arm.md section 12.9) and re-tabulate calibration.
#
# WHY POST-FILTERING IS EXACT HERE (report section 11.3). The Individual_Analysis drivers
# (src/STAARpipeline_Individual_Analysis.*.r) all ran with QC_label = QC_label2, the
# defective position-inherited label. A single-variant score test is computed per variant
# against ONE fixed null model, so each variant's p-value does not depend on which other
# variants are in the scan. Dropping the variants that fail QC_label7 therefore reproduces
# EXACTLY what a rerun with QC_label = QC_label7 would produce (this is NOT true of the
# gene-centric burden tests, which aggregate variants).
#
# Join logic is reused from revision/scripts/adsp_sv_qclabel5.R / adsp_sv_plots_qc5.R:
# per chromosome, key = (POS, REF, ALT) with REF/ALT parsed from the aGDS "allele" node,
# deduplicated on the key (first occurrence), right-joined onto the results rows.
# The only change is that the aGDS is read through a variant filter restricted to the
# positions present in the results (identical answer, far less memory than materialising
# 47M allele strings for chr1).
#
# Inputs : <arm>/individual_analysis/results_individual_analysis_genome.Rdata (as run, QC_label2)
# Outputs: <arm>/individual_analysis/results_individual_analysis_genome.QC7.Rdata (NEW, same object name)
#          revision/tables/qc7_sv_calibration.csv     (long: arm x label x maf_bin)
#          revision/tables/qc7_sv_summary.csv         (wide: QC2 vs QC7 per arm, "all" rows)
#          revision/tables/qc7_sv_label_overlap.csv   (QC2/QC5/QC7 membership among tested variants)
suppressPackageStartupMessages({ library(data.table); library(SeqArray) })
setDTthreads(4)
OUT <- "${ADSP_ROOT}/output"; REV <- "${ADSP_ROOT}/revision"
TAB <- file.path(REV, "tables")
ALPHA <- 5e-8
args <- commandArgs(trailingOnly = TRUE); stopifnot(length(args) == 3)   # arms,dirs,tag
ARMS <- data.table(arm = strsplit(args[1], ",")[[1]], dir = file.path(OUT, strsplit(args[2], ",")[[1]], "individual_analysis"))
TAG <- args[3]
ARMS[, src := file.path(dir, "results_individual_analysis_genome.Rdata")]
ARMS[, dst := file.path(dir, "results_individual_analysis_genome.QC7.Rdata")]
stopifnot(all(file.exists(ARMS$src)))
if (any(file.exists(ARMS$dst))) stop("refusing to overwrite existing QC7 file(s): ",
                                     paste(ARMS$dst[file.exists(ARMS$dst)], collapse=", "))
ad <- get(load(file.path(OUT, "agds_dir.Rdata"))); stopifnot(length(ad) == 22, all(file.exists(ad)))
t0 <- Sys.time(); el <- function() sprintf("[%.1f min]", as.numeric(difftime(Sys.time(), t0, units="mins")))

## ---- 1. load the as-run tables; keep the original object untouched ------------------
FULL <- list(); OBJ <- list(); K <- list()
for (i in seq_len(nrow(ARMS))) { a <- ARMS$arm[i]
  e <- new.env(); load(ARMS$src[i], envir=e); nm <- ls(e)[1]; x <- get(nm, envir=e); rm(e)
  OBJ[[a]] <- nm; FULL[[a]] <- x
  stopifnot(all(c("CHR","POS","REF","ALT","MAF","pvalue") %in% names(x)))
  K[[a]] <- data.table(row = seq_len(nrow(x)), CHR = as.integer(x$CHR), POS = as.numeric(x$POS),
                       REF = as.character(x$REF), ALT = as.character(x$ALT))
  cat(sprintf("%s loaded %-9s object=%s class=%s rows=%s cols=%s  CHR range %s-%s\n", el(), a, nm,
              paste(class(x), collapse="/"), format(nrow(x), big.mark=","), ncol(x),
              min(K[[a]]$CHR, na.rm=TRUE), max(K[[a]]$CHR, na.rm=TRUE)))
  print(table(K[[a]]$CHR %in% 1:22, useNA="ifany"))
}

## ---- 2. per-chromosome QC label lookup from the aGDS (QC_label2 / 5 / 7) ------------
FLAG <- setNames(vector("list", nrow(ARMS)), ARMS$arm); DUP <- integer(22)
for (cc in 1:22) {
  posset <- unique(unlist(lapply(K, function(k) k[CHR == cc]$POS)))
  g <- seqOpen(ad[cc]); pos <- as.numeric(seqGetData(g, "position"))
  sel <- pos %in% posset
  seqSetFilter(g, variant.sel = sel, verbose = FALSE)
  al <- seqGetData(g, "allele")
  q <- data.table(POS = pos[sel], REF = sub(",.*$", "", al), ALT = sub("^[^,]*,", "", al),
                  q2 = as.character(seqGetData(g, "annotation/info/QC_label2")) == "PASS",
                  q5 = as.character(seqGetData(g, "annotation/info/QC_label5")) == "PASS",
                  q7 = as.character(seqGetData(g, "annotation/info/QC_label7")) == "PASS")
  rm(al); seqClose(g)
  ndup <- nrow(q) - uniqueN(q, by = c("POS","REF","ALT")); DUP[cc] <- ndup
  ndisc <- if (ndup > 0) q[, .(d = uniqueN(q7) > 1L), by = .(POS,REF,ALT)][d == TRUE, .N] else 0L
  q <- unique(q, by = c("POS","REF","ALT")); setkey(q, POS, REF, ALT)
  for (a in ARMS$arm) { k <- K[[a]][CHR == cc]; if (!nrow(k)) next
    m <- q[k, on = .(POS, REF, ALT)]            # right join: every results row kept, NA if absent from aGDS
    stopifnot(nrow(m) == nrow(k), !anyDuplicated(m$row))   # dedup above guarantees no row multiplication
    FLAG[[a]][[length(FLAG[[a]]) + 1L]] <- m[, .(row, q2, q5, q7)] }
  cat(sprintf("%s chr%-2d aGDS variants %s | at results positions %s | dup keys %d (QC7-discordant %d) | results rows %s\n",
              el(), cc, format(length(pos), big.mark=","), format(sum(sel), big.mark=","), ndup, ndisc,
              paste(sapply(K, function(k) format(k[CHR == cc, .N], big.mark=",")), collapse="/")))
  rm(q, pos, sel, posset); gc(verbose = FALSE)
}

## ---- 3. apply the filter, write the QC7 files, tabulate label membership -------------
lam  <- function(p) median(qchisq(p, 1, lower.tail = FALSE)) / qchisq(0.5, 1)
lamq <- function(p, q) unname(quantile(qchisq(p, 1, lower.tail = FALSE), q) / qchisq(1 - q, 1, lower.tail = FALSE))
is_apoe <- function(CHR, POS) CHR == 19 & POS > 44.4e6 & POS < 46.4e6
met <- function(z) {                       # z: CHR, POS, p (finite, 0<p<=1)
  apo <- is_apoe(z$CHR, z$POS); e <- z[!apo]; s <- z[p < ALPHA]
  s[, bin := paste0(CHR, "_", floor(POS / 1e6))]
  data.table(n = nrow(z), lambda = round(lam(z$p), 4), n_exAPOE = nrow(e), lambda_exAPOE = round(lam(e$p), 4),
             n_lt1e4_exAPOE = sum(e$p < 1e-4), oe1e4_exAPOE = round(sum(e$p < 1e-4) / (nrow(e) * 1e-4), 3),
             n_lt1e5_exAPOE = sum(e$p < 1e-5), oe1e5_exAPOE = round(sum(e$p < 1e-5) / (nrow(e) * 1e-5), 3),
             n_lt1e6_exAPOE = sum(e$p < 1e-6), oe1e6_exAPOE = round(sum(e$p < 1e-6) / (nrow(e) * 1e-6), 3),
             q99_exAPOE = round(lamq(e$p, .99), 4), q999_exAPOE = round(lamq(e$p, .999), 4),
             oe1e4_all = round(sum(z$p < 1e-4) / (nrow(z) * 1e-4), 3),
             nsig_5e8 = nrow(s), loci_1Mb = uniqueN(s$bin),
             loci_nonAPOE = uniqueN(s[!is_apoe(CHR, POS)]$bin),
             apoe_peak = if (any(apo)) round(max(-log10(z[apo]$p)), 1) else NA_real_)
}
CAL <- list(); OV <- list()
for (a in ARMS$arm) {
  x <- FULL[[a]]; Fg <- rbindlist(FLAG[[a]]); setorder(Fg, row)
  q2 <- q5 <- q7 <- rep(NA, nrow(x)); q2[Fg$row] <- Fg$q2; q5[Fg$row] <- Fg$q5; q7[Fg$row] <- Fg$q7
  p <- as.numeric(x$pvalue); pv <- is.finite(p) & p > 0 & p <= 1
  keep <- !is.na(q7) & q7
  OV[[a]] <- data.table(arm = a, n_rows = nrow(x), n_p_valid = sum(pv), n_unmatched_in_aGDS = sum(is.na(q7)),
    n_p_valid_q2_PASS = sum(pv & q2 %in% TRUE), n_p_valid_q2_FAIL = sum(pv & q2 %in% FALSE),
    n_p_valid_q5_PASS = sum(pv & q5 %in% TRUE), n_p_valid_q7_PASS = sum(pv & q7 %in% TRUE),
    n_p_valid_q5_not_q7 = sum(pv & q5 %in% TRUE & q7 %in% FALSE), n_p_valid_q7_not_q5 = sum(pv & q7 %in% TRUE & q5 %in% FALSE),
    pct_kept_q5 = round(100 * sum(pv & q5 %in% TRUE) / sum(pv), 2), pct_kept_q7 = round(100 * sum(pv & q7 %in% TRUE) / sum(pv), 2),
    n_rows_written_QC7 = sum(keep))
  print(OV[[a]])
  # write the filtered table under the ORIGINAL object name, next to the original
  y <- if (is.data.table(x)) x[which(keep)] else x[which(keep), , drop = FALSE]
  assign(OBJ[[a]], y); save(list = OBJ[[a]], file = ARMS[arm == a]$dst); rm(list = OBJ[[a]]); rm(y)
  cat(sprintf("%s wrote %s (%s rows)\n", el(), ARMS[arm == a]$dst, format(sum(keep), big.mark=",")))
  # calibration, QC2 (as run) vs QC7, overall and per MAF bin
  z0 <- data.table(CHR = as.integer(x$CHR), POS = as.numeric(x$POS), p = p, MAF = as.numeric(x$MAF), q7 = q7)[pv]
  z0[, maf := pmin(MAF, 1 - MAF)]; z0[, maf_bin := as.character(cut(maf, c(0, 0.001, 0.005, 0.01, 0.05, 0.5), include.lowest = TRUE))]
  for (L in c("QC2", "QC7")) { z <- if (L == "QC2") z0 else z0[q7 %in% TRUE]
    CAL[[length(CAL) + 1L]] <- cbind(data.table(arm = a, label = L, maf_bin = "all"), met(z))
    for (b in c("[0,0.001]", "(0.001,0.005]", "(0.005,0.01]", "(0.01,0.05]", "(0.05,0.5]")) { zb <- z[maf_bin == b]
      if (nrow(zb)) CAL[[length(CAL) + 1L]] <- cbind(data.table(arm = a, label = L, maf_bin = b), met(zb)) } }
  rm(x, z0, Fg, q2, q5, q7, p, pv, keep); FULL[[a]] <- NULL; FLAG[[a]] <- NULL; gc(verbose = FALSE)
}
cal <- rbindlist(CAL); ov <- rbindlist(OV)
cal[, frac_kept_vs_QC2 := round(n / n[label == "QC2"], 4), by = .(arm, maf_bin)]
fwrite(cal, file.path(TAB, paste0("qc7_sv_", TAG, "_calibration.csv")))
fwrite(ov,  file.path(TAB, paste0("qc7_sv_", TAG, "_label_overlap.csv")))
w <- dcast(cal[maf_bin == "all"], arm ~ label,
           value.var = c("n", "lambda", "lambda_exAPOE", "oe1e4_exAPOE", "oe1e5_exAPOE", "oe1e6_exAPOE",
                         "q999_exAPOE", "nsig_5e8", "loci_1Mb", "loci_nonAPOE", "apoe_peak"))
w <- merge(w, ov[, .(arm, pct_kept_q5, pct_kept_q7)], by = "arm")
w[, frac_kept_QC7 := round(n_QC7 / n_QC2, 4)]
fwrite(w, file.path(TAB, paste0("qc7_sv_", TAG, "_summary.csv")))
cat("\n=== label membership among tested variants ===\n"); print(ov)
cat("\n=== calibration, all variants ===\n"); print(cal[maf_bin == "all"], nrows = 100)
cat("\n=== per-MAF-bin O/E@1e-4 (ex-APOE) ===\n")
print(dcast(cal, arm + maf_bin ~ label, value.var = "oe1e4_exAPOE"))
cat("\naGDS duplicate (POS,REF,ALT) keys per chromosome:\n"); print(DUP)
cat(sprintf("%s DONE\n", el()))
