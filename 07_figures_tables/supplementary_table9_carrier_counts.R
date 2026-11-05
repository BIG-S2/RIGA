###########################################################
# supplementary_table9_carrier_counts.R (2026-10-04)
# Supplementary Table S4 -- rare-variant carrier counts in cases and
# controls for the gene x mask sets of the STRATIFIED design D_S that
# pass FDR 5% (read from supp_v3/S3_fdr_gene_mask_sets.csv) plus the
# five PSEN1 coding masks discussed in the text. Derived from
# revision/scripts/build_S4_carriers.R (matched design, 2026-08-06):
# same variant selection as the deployed scan (QC_label7, SNV,
# rare_maf_cutoff 0.01) through the STAARpipelineSummary *_Info functions.
###########################################################
suppressPackageStartupMessages({
  library(data.table); library(SeqArray); library(gdsfmt); library(SeqVarTools)
  library(STAAR); library(STAARpipeline); library(STAARpipelineSummary)
})
OUT   <- "${ADSP_ROOT}/revision/tables/supp_v3"
MATCH <- "${ADSP_ROOT}/output/matched/matched_samples_lpa_pc.pair33_projcell.csv"
log_lines <- character(0)
LOG <- function(...) { s <- sprintf(...); message(s); log_lines <<- c(log_lines, s) }

matched <- fread(MATCH)
LOG("analytic set rows: %d ; cases %d ; controls %d", nrow(matched), sum(matched$AD == 1), sum(matched$AD == 0))
obj_nullmodel <- get(load("${ADSP_ROOT}/output/matched/obj_null_stratum_all_matched_lpa_pc.pair33_projcell.RData"))
obj_nullmodel$use_SPA <- FALSE      # only the variant lists are used here
Annotation_name_catalog <- get(load("${ADSP_ROOT}/output/Annotation_name_catalog.Rdata"))
agds_dir <- get(load("${ADSP_ROOT}/output/agds_dir.Rdata"))
Annotation_name <- c("CADD","LINSIGHT","FATHMM.XF","aPC.EpigeneticActive","aPC.EpigeneticRepressed",
                     "aPC.EpigeneticTranscription","aPC.Conservation","aPC.LocalDiversity",
                     "aPC.Mappability","aPC.TF","aPC.Protein")
QC_label <- "annotation/info/QC_label7"

pheno_id <- as.character(obj_nullmodel$id_include)
y <- matched$AD[match(pheno_id, as.character(matched$IID))]
LOG("null-model N = %d ; matched to phenotype: %d ; cases %d controls %d",
    length(pheno_id), sum(!is.na(y)), sum(y == 1, na.rm = TRUE), sum(y == 0, na.rm = TRUE))
stopifnot(!any(is.na(y)))

s3 <- fread(file.path(OUT, "S3_fdr_gene_mask_sets.csv"))
kind_of <- function(a) ifelse(grepl("coding$", a) & !grepl("noncoding", a), "coding", ifelse(grepl("ncRNA", a), "ncRNA", "noncoding"))
targets <- unique(rbind(
  s3[, .(gene = `Gene name`, chr = as.integer(Chr), category = Category, kind = kind_of(Analysis), source = "FDR 5% set")],
  data.table(gene = "PSEN1", chr = 14L, category = c("plof", "plof_ds", "missense", "disruptive_missense", "synonymous"),
             kind = "coding", source = "Discussed in text (not significant)")), by = c("gene", "chr", "category"))
targets <- targets[!is.na(chr)]
LOG("targets: %d sets on chromosomes %s", nrow(targets), paste(sort(unique(targets$chr)), collapse = ","))

s4 <- data.table(); var_detail <- data.table()
for (chrN in sort(unique(targets$chr))) {
  gp <- agds_dir[chrN]
  if (is.na(gp) || !file.exists(gp)) { LOG("chr%d aGDS missing", chrN); next }
  genofile <- seqOpen(gp); LOG("opened %s", basename(gp))
  for (r in which(targets$chr == chrN)) {
    gn <- targets$gene[r]; ct <- targets$category[r]; kd <- targets$kind[r]
    seqResetFilter(genofile, verbose = FALSE)
    info <- NULL
    res <- try({
      if (kd == "coding") {
        info <- Gene_Centric_Coding_Info(category = ct, chr = chrN, genofile = genofile, obj_nullmodel = obj_nullmodel,
                  gene_name = gn, rare_maf_cutoff = 0.01, method_cond = "optimal", QC_label = QC_label,
                  variant_type = "SNV", geno_missing_imputation = "mean",
                  Annotation_dir = "annotation/info/FunctionalAnnotation",
                  Annotation_name_catalog = Annotation_name_catalog, Annotation_name = Annotation_name)
      } else if (kd == "noncoding") {
        info <- Gene_Centric_Noncoding_Info(category = ct, chr = chrN, genofile = genofile, obj_nullmodel = obj_nullmodel,
                  gene_name = gn, rare_maf_cutoff = 0.01, method_cond = "optimal", QC_label = QC_label,
                  variant_type = "SNV", geno_missing_imputation = "mean",
                  Annotation_dir = "annotation/info/FunctionalAnnotation",
                  Annotation_name_catalog = Annotation_name_catalog, Annotation_name = Annotation_name)
      } else {
        info <- info_ncRNA(chr = chrN, genofile = genofile, obj_nullmodel = obj_nullmodel, gene_name = gn,
                  known_loci = NULL, rare_maf_cutoff = 0.01, method_cond = "optimal", QC_label = QC_label,
                  variant_type = "SNV", geno_missing_imputation = "mean",
                  Annotation_dir = "annotation/info/FunctionalAnnotation",
                  Annotation_name_catalog = Annotation_name_catalog, Annotation_name = Annotation_name)
      }
    }, silent = TRUE)
    if (inherits(res, "try-error") || is.null(info)) {
      LOG("%s / %s : Info() FAILED (%s)", gn, ct, paste(attr(res, "condition")$message, collapse = " ")); next }
    info <- as.data.frame(info)
    LOG("%s / %s : %d variants; cols = %s", gn, ct, nrow(info), paste(names(info), collapse = ", "))
    if (nrow(info) == 0) next
    poscol <- intersect(c("POS", "Position", "position", "pos"), names(info))[1]
    refcol <- intersect(c("REF", "Ref", "ref"), names(info))[1]
    altcol <- intersect(c("ALT", "Alt", "alt"), names(info))[1]
    if (is.na(poscol)) { LOG("%s / %s : no POS column", gn, ct); next }
    want_pos <- as.integer(info[[poscol]])
    seqResetFilter(genofile, verbose = FALSE)
    all_pos <- seqGetData(genofile, "position")
    if (!is.na(refcol) && !is.na(altcol)) {
      key_gds <- paste(all_pos, seqGetData(genofile, "$ref"), seqGetData(genofile, "$alt"), sep = ":")
      key_want <- paste(want_pos, info[[refcol]], info[[altcol]], sep = ":")
      vsel <- which(key_gds %in% key_want)
      LOG("%s / %s : matched on POS:REF:ALT (%d of %d)", gn, ct, length(vsel), length(key_want))
    } else vsel <- which(all_pos %in% want_pos)
    if (!length(vsel)) { LOG("%s / %s : no GDS match", gn, ct); next }
    seqSetFilter(genofile, variant.sel = vsel, sample.id = pheno_id, verbose = FALSE)
    sid_f <- seqGetData(genofile, "sample.id"); ord <- match(pheno_id, sid_f)
    G <- seqGetData(genofile, "$dosage_alt")[ord, , drop = FALSE]
    vpos <- seqGetData(genofile, "position"); vref <- seqGetData(genofile, "$ref"); valt <- seqGetData(genofile, "$alt")
    G[is.na(G)] <- 0
    af <- colMeans(G) / 2; flip <- which(af > 0.5); if (length(flip)) G[, flip] <- 2 - G[, flip]
    maf <- colMeans(G) / 2; keep <- which(maf > 0 & maf < 0.01)
    if (!length(keep)) { LOG("%s / %s : no rare variants after MAF filter", gn, ct); next }
    G <- G[, keep, drop = FALSE]; vpos <- vpos[keep]; vref <- vref[keep]; valt <- valt[keep]; maf <- maf[keep]
    carr <- rowSums(G > 0) > 0
    nc <- sum(y == 1); nk <- sum(y == 0); cc <- sum(carr & y == 1); ck <- sum(carr & y == 0)
    s4 <- rbind(s4, data.table(
      `Gene name` = gn, Chr = chrN, Category = ct, Source = targets$source[r],
      `#SNV in set` = ncol(G), cMAC = sum(G),
      `MAC in cases` = sum(G[y == 1, , drop = FALSE]), `MAC in controls` = sum(G[y == 0, , drop = FALSE]),
      `Carriers (cases)` = cc, `N cases` = nc, `Carriers (controls)` = ck, `N controls` = nk,
      `Carrier freq. (cases)` = round(cc / nc, 5), `Carrier freq. (controls)` = round(ck / nk, 5),
      `Carrier OR (unadjusted)` = round((cc / (nc - cc)) / (ck / (nk - ck)), 3),
      `Fisher P (carrier, unadjusted)` = signif(fisher.test(matrix(c(cc, nc - cc, ck, nk - ck), 2))$p.value, 4)), fill = TRUE)
    LOG("%s / %s : %d SNV, carriers %d/%d cases vs %d/%d controls", gn, ct, ncol(G), cc, nc, ck, nk)
    if (kd == "coding" || ncol(G) <= 400)
      var_detail <- rbind(var_detail, data.table(
        `Gene name` = gn, Category = ct, Chr = chrN, POS = vpos, REF = vref, ALT = valt, MAF = round(maf, 8),
        `MAC total` = colSums(G), `MAC cases` = colSums(G[y == 1, , drop = FALSE]),
        `MAC controls` = colSums(G[y == 0, , drop = FALSE]),
        `Carriers cases` = colSums(G[y == 1, , drop = FALSE] > 0),
        `Carriers controls` = colSums(G[y == 0, , drop = FALSE] > 0)), fill = TRUE)
    seqResetFilter(genofile, verbose = FALSE)
  }
  seqClose(genofile)
}
if (nrow(s4)) { fwrite(s4, file.path(OUT, "S4_carrier_counts.csv")); LOG("S4 written: %d rows", nrow(s4)) }
if (nrow(var_detail)) { fwrite(var_detail, file.path(OUT, "S4b_variant_level_carriers.csv")); LOG("S4b written: %d variants", nrow(var_detail)) }
writeLines(log_lines, file.path(OUT, "build_log_S4.txt"))
LOG("done.")
