#!/usr/bin/env Rscript
## conditional_analysis_apoe_trem2.R (2026-10-13, review point 6; v3: normal approximation only; saves every STAAR column of the
## unconditional and conditional rows; adds the primary statistic SB = Cauchy combination of all SKAT and burden p values (STAAR-O without ACAT-V)): conditional gene-centric tests on the deployed stratified
## null (obj_null_stratum_all_matched_lpa_pc.pair33_projcell.RData, QC_label7, deployed annotation weights).
##   A  APOE-region sets | APOE e4 (rs429358) + e2 (rs7412)
##   B  APOE-region sets | e4 + e2 + the regional lead indel chr19:44909967 TGG>T
##   C  TREM2 promoter_DHS and missense | TREM2 coding risk variants R47H (and R62H when present and QC_label7 PASS)
## Each set is tested unconditionally and conditionally, with STAAR-O (use_SPA FALSE, deployed) and the saddlepoint burden
## (use_SPA TRUE). Also writes the per-variant score tests of the rare NECTIN2 missense variants (point 5/6: which alleles
## drive the coding signal and in which direction). Writes revision/tables/review/conditional_*.csv.
suppressPackageStartupMessages({ library(data.table); library(SeqArray); library(STAAR); library(STAARpipeline) })
OUT <- "${ADSP_ROOT}/output"; TAB <- "${ADSP_ROOT}/revision/tables/review"; dir.create(TAB, showWarnings = FALSE)
obj <- get(load(file.path(OUT, "matched/obj_null_stratum_all_matched_lpa_pc.pair33_projcell.RData"))); obj$sparse_kins <- TRUE; obj$relatedness <- TRUE
agds_dir <- get(load(file.path(OUT, "agds_dir.Rdata"))); catalog <- get(load(file.path(OUT, "Annotation_name_catalog.Rdata")))
ANN <- c("CADD", "LINSIGHT", "FATHMM.XF", "aPC.EpigeneticActive", "aPC.EpigeneticRepressed", "aPC.EpigeneticTranscription",
         "aPC.Conservation", "aPC.LocalDiversity", "aPC.Mappability", "aPC.TF", "aPC.Protein")
QC <- "annotation/info/QC_label7"
present <- function(k) {  # keep known loci present in the aGDS with QC_label7 PASS and matching alleles
  keep <- logical(nrow(k))
  for (i in seq_len(nrow(k))) { g <- seqOpen(agds_dir[k$CHR[i]]); seqSetFilterPos(g, chr = as.character(k$CHR[i]), pos = k$POS[i], verbose = FALSE)
    if (length(seqGetData(g, "variant.id"))) { al <- seqGetData(g, "allele"); q <- seqGetData(g, QC); keep[i] <- any(al == paste(k$REF[i], k$ALT[i], sep = ",") & q == "PASS") }
    seqClose(g) }
  cat("known loci kept:", paste(k$POS[keep], collapse = ","), " dropped:", paste(k$POS[!keep], collapse = ","), "\n"); k[keep, ] }
apoe <- present(data.frame(CHR = c(19, 19), POS = c(44908684, 44908822), REF = c("T", "C"), ALT = c("C", "T")))
apoe_lead <- rbind(apoe, present(data.frame(CHR = 19, POS = 44909967, REF = "TGG", ALT = "T")))
trem2 <- present(data.frame(CHR = c(6, 6), POS = c(41161514, 41161469), REF = c("C", "C"), ALT = c("T", "T")))
sets <- data.table(gene = c("NECTIN2", "NECTIN2", "TOMM40", "TOMM40", "TOMM40", "APOE", "CLASRP", "PPP1R37", "TREM2", "TREM2"), chr = c(19, 19, 19, 19, 19, 19, 19, 19, 6, 6),
                   cat = c("missense", "enhancer_DHS", "promoter_DHS", "UTR", "downstream", "enhancer_DHS", "downstream", "enhancer_DHS", "promoter_DHS", "missense"))
flat <- function(r) {   # v2: a single-category call returns a matrix (often of mode list); an all-category call returns a list
  if (is.null(r) || !length(r)) return(NULL)
  m <- if (is.matrix(r) || is.data.frame(r)) r else r[[1]]
  if (is.null(m) || !length(m)) return(NULL)
  d <- as.data.frame(m, stringsAsFactors = FALSE, check.names = FALSE)
  for (cc in names(d)) if (is.list(d[[cc]])) d[[cc]] <- unlist(lapply(d[[cc]], function(z) if (length(z)) z[1] else NA))
  d }
getp <- function(d, cols) { for (cl in cols) if (!is.null(d) && cl %in% names(d)) return(as.numeric(d[[cl]][1])); NA_real_ }
cct <- function(p) { p <- p[is.finite(p)]; p <- pmin(pmax(p, 1e-300), 1 - 1e-16); 0.5 - atan(mean(tan((0.5 - p) * pi))) / pi }
sb_of <- function(d) { if (is.null(d)) return(NA_real_); cols <- grep("^(SKAT|Burden)\\(1,(25|1)\\)", names(d), value = TRUE); if (!length(cols)) return(NA_real_); cct(as.numeric(unlist(d[1, cols]))) }
res <- list(); rows_full <- list()
for (i in seq_len(nrow(sets))) {
  s <- sets[i]; coding <- s$cat == "missense"; known_list <- if (s$gene == "TREM2") list(C = trem2) else list(A = apoe, B = apoe_lead)
  for (spa in c(FALSE)) {
    obj$use_SPA <- spa; g <- seqOpen(agds_dir[s$chr])
    common <- list(chr = s$chr, gene_name = s$gene, category = s$cat, genofile = g, obj_nullmodel = obj, rare_maf_cutoff = 0.01, rv_num_cutoff = 2,
                   QC_label = QC, variant_type = "SNV", geno_missing_imputation = "mean", Annotation_dir = "annotation/info/FunctionalAnnotation",
                   Annotation_name_catalog = catalog, Use_annotation_weights = TRUE, Annotation_name = ANN)
    f0 <- if (coding) Gene_Centric_Coding else Gene_Centric_Noncoding
    u <- flat(tryCatch(do.call(f0, common), error = function(e) { cat("ERR uncond", s$gene, s$cat, conditionMessage(e), "\n"); NULL }))
    pu <- getp(u, "STAAR-O"); sbu <- sb_of(u); if (!is.null(u)) rows_full[[length(rows_full) + 1]] <- cbind(data.table(gene = s$gene, category = s$cat, conditioning = "none"), as.data.table(u))
    for (kn in names(known_list)) {
      fc <- if (coding) { if (spa) Gene_Centric_Coding_cond_spa else Gene_Centric_Coding_cond } else { if (spa) Gene_Centric_Noncoding_cond_spa else Gene_Centric_Noncoding_cond }
      args <- c(common, list(known_loci = known_list[[kn]])); if (!spa) args$method_cond <- "optimal"
      cr <- flat(tryCatch(do.call(fc, args), error = function(e) { cat("ERR cond", s$gene, s$cat, kn, conditionMessage(e), "\n"); NULL }))
      pc <- getp(cr, c("STAAR-O_cond", "STAAR-O")); sbc <- sb_of(cr); if (!is.null(cr)) rows_full[[length(rows_full) + 1]] <- cbind(data.table(gene = s$gene, category = s$cat, conditioning = kn), as.data.table(cr))
      res[[length(res) + 1]] <- data.table(gene = s$gene, category = s$cat, conditioning = kn, n_known = nrow(known_list[[kn]]),
                                           SB_unconditional = sbu, SB_conditional = sbc, STAARO_unconditional = pu, STAARO_conditional = pc)
      cat(sprintf("%-8s %-13s %s  SB %.3g -> %.3g | STAAR-O %.3g -> %.3g\n", s$gene, s$cat, kn, sbu, sbc, pu, pc)) }
    seqClose(g) } }
fwrite(rbindlist(res), file.path(TAB, "conditional_apoe_trem2_v3.csv")); fwrite(rbindlist(rows_full, fill = TRUE), file.path(TAB, "conditional_apoe_trem2_v3_fullrows.csv"))
cat("DONE\n")
