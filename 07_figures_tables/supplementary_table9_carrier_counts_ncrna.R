###########################################################
# supplementary_table9_carrier_counts_ncrna.R (2026-10-04) -- carrier counts for the ncRNA
# set(s) of Supplementary Table 4 (info_ncRNA() of STAARpipelineSummary
# fails on this aGDS: "argument 1 is not a vector"). The variant selection
# replicates STAARpipeline::ncRNA exactly: QC_label7 PASS & SNV; GENCODE
# category ncRNA_exonic / ncRNA_exonic;splicing / ncRNA_splicing; gene name
# among the first three comma-separated GENCODE.Info genes; MAF < 1% in the
# analytic set. Appends to supp_v3/S4_carrier_counts.csv (+ S4b).
###########################################################
suppressPackageStartupMessages({ library(data.table); library(SeqArray); library(gdsfmt); library(SeqVarTools) })
OUT   <- "${ADSP_ROOT}/revision/tables/supp_v3"
MATCH <- "${ADSP_ROOT}/output/matched/matched_samples_lpa_pc.pair33_projcell.csv"
log_lines <- character(0); LOG <- function(...) { s <- sprintf(...); message(s); log_lines <<- c(log_lines, s) }
matched <- fread(MATCH)
obj_nullmodel <- get(load("${ADSP_ROOT}/output/matched/obj_null_stratum_all_matched_lpa_pc.pair33_projcell.RData"))
cat_tbl <- get(load("${ADSP_ROOT}/output/Annotation_name_catalog.Rdata"))
agds_dir <- get(load("${ADSP_ROOT}/output/agds_dir.Rdata"))
A <- "annotation/info/FunctionalAnnotation"
pheno_id <- as.character(obj_nullmodel$id_include)
y <- matched$AD[match(pheno_id, as.character(matched$IID))]; stopifnot(!any(is.na(y)))
s3 <- fread(file.path(OUT, "S3_fdr_gene_mask_sets.csv"))
targets <- s3[grepl("ncRNA", Analysis), .(gene = `Gene name`, chr = as.integer(Chr))]
LOG("ncRNA targets: %s", paste(sprintf("%s (chr%d)", targets$gene, targets$chr), collapse = ", "))
s4 <- data.table(); var_detail <- data.table()
for (chrN in unique(targets$chr)) {
  genofile <- seqOpen(agds_dir[chrN])
  filt <- seqGetData(genofile, "annotation/info/QC_label7")
  SNVlist <- (filt == "PASS") & isSNV(genofile)
  vid <- seqGetData(genofile, "variant.id")
  gcat <- seqGetData(genofile, paste0(A, cat_tbl$dir[cat_tbl$name == "GENCODE.Category"]))
  is.in <- (gcat %in% c("ncRNA_exonic", "ncRNA_exonic;splicing", "ncRNA_splicing")) & SNVlist
  vid_nc <- vid[is.in]; rm(gcat, filt, SNVlist)
  seqSetFilter(genofile, variant.id = vid_nc, verbose = FALSE)
  ginfo <- seqGetData(genofile, paste0(A, cat_tbl$dir[cat_tbl$name == "GENCODE.Info"]))
  Gene <- as.character(sapply(strsplit(ginfo, split = "[;]"), function(z) gsub("\\(.*\\)", "", z[1])))
  gl <- strsplit(Gene, ","); g1 <- sapply(gl, "[", 1); g2 <- sapply(gl, "[", 2); g3 <- sapply(gl, "[", 3)
  vid_nc <- seqGetData(genofile, "variant.id"); seqResetFilter(genofile, verbose = FALSE)
  for (r in which(targets$chr == chrN)) {
    gn <- targets$gene[r]
    sel <- union(union(which(g1 == gn), which(g2 == gn)), which(g3 == gn))
    LOG("%s: %d ncRNA SNVs selected on chr%d", gn, length(sel), chrN)
    if (!length(sel)) next
    seqSetFilter(genofile, variant.id = vid_nc[sel], sample.id = pheno_id, verbose = FALSE)
    sid_f <- seqGetData(genofile, "sample.id"); ord <- match(pheno_id, sid_f)
    G <- seqGetData(genofile, "$dosage_alt")[ord, , drop = FALSE]
    vpos <- seqGetData(genofile, "position"); vref <- seqGetData(genofile, "$ref"); valt <- seqGetData(genofile, "$alt")
    G[is.na(G)] <- 0
    af <- colMeans(G) / 2; flip <- which(af > 0.5); if (length(flip)) G[, flip] <- 2 - G[, flip]
    maf <- colMeans(G) / 2; keep <- which(maf > 0 & maf < 0.01)
    G <- G[, keep, drop = FALSE]; vpos <- vpos[keep]; vref <- vref[keep]; valt <- valt[keep]; maf <- maf[keep]
    carr <- rowSums(G > 0) > 0
    nc <- sum(y == 1); nk <- sum(y == 0); cc <- sum(carr & y == 1); ck <- sum(carr & y == 0)
    s4 <- rbind(s4, data.table(
      `Gene name` = gn, Chr = chrN, Category = "ncRNA", Source = "FDR 5% set",
      `#SNV in set` = ncol(G), cMAC = sum(G),
      `MAC in cases` = sum(G[y == 1, , drop = FALSE]), `MAC in controls` = sum(G[y == 0, , drop = FALSE]),
      `Carriers (cases)` = cc, `N cases` = nc, `Carriers (controls)` = ck, `N controls` = nk,
      `Carrier freq. (cases)` = round(cc / nc, 5), `Carrier freq. (controls)` = round(ck / nk, 5),
      `Carrier OR (unadjusted)` = round((cc / (nc - cc)) / (ck / (nk - ck)), 3),
      `Fisher P (carrier, unadjusted)` = signif(fisher.test(matrix(c(cc, nc - cc, ck, nk - ck), 2))$p.value, 4)))
    var_detail <- rbind(var_detail, data.table(
      `Gene name` = gn, Category = "ncRNA", Chr = chrN, POS = vpos, REF = vref, ALT = valt, MAF = round(maf, 8),
      `MAC total` = colSums(G), `MAC cases` = colSums(G[y == 1, , drop = FALSE]), `MAC controls` = colSums(G[y == 0, , drop = FALSE]),
      `Carriers cases` = colSums(G[y == 1, , drop = FALSE] > 0), `Carriers controls` = colSums(G[y == 0, , drop = FALSE] > 0)))
    LOG("%s / ncRNA : %d SNV (MAF<1%%), cMAC %d, carriers %d/%d cases vs %d/%d controls", gn, ncol(G), sum(G), cc, nc, ck, nk)
    seqResetFilter(genofile, verbose = FALSE)
  }
  seqClose(genofile)
}
f4 <- file.path(OUT, "S4_carrier_counts.csv"); f4b <- file.path(OUT, "S4b_variant_level_carriers.csv")
if (nrow(s4)) {
  old <- if (file.exists(f4)) fread(f4) else data.table()
  old <- old[!(Category == "ncRNA" & `Gene name` %in% s4$`Gene name`)]
  fwrite(rbind(old, s4, fill = TRUE), f4); LOG("S4 appended: %d ncRNA rows (total %d)", nrow(s4), nrow(old) + nrow(s4))
  oldb <- if (file.exists(f4b)) fread(f4b) else data.table()
  oldb <- oldb[!(Category == "ncRNA" & `Gene name` %in% s4$`Gene name`)]
  fwrite(rbind(oldb, var_detail, fill = TRUE), f4b)
}
writeLines(log_lines, file.path(OUT, "build_log_S4_ncRNA.txt")); LOG("done.")
