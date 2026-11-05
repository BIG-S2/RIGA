#!/usr/bin/env Rscript
## benchmark_burden_test_by_chromosome.R -- plasmode_burden_chr22_maf001.R for ANY chromosome (2026-09-27, chr21 replication of the
## chr22 burden channel): CHR from the environment selects processed/r5_chr<CHR>.split.gds and outputs/plasmode/burden_chr<CHR>_maf001/;
## window units are named chr<CHR>_win<k>.  Everything else is identical to the chr22 builder.
###############################################################################
## plasmode_burden_chr22_maf001.R -- MAF < 0.1 % variant of plasmode_burden_chr22.R (2026-09-22, user rule:
## rare-variant tests include only variants with MAF < 0.1 %).  Identical SeqArray pass, QC_label7 == PASS,
## pool MAC with the same minor-allele flip; only the MASKS table differs:
##   gene_maf001 / coding_maf001 / syn_maf001 / win100k_maf001 : MAC 1-54  (MAF < 0.001 in the 27,306-subject pool)
##   gene_mac6_54 / coding_mac6_54                              : MAC 6-54  (the same without singletons..quintuplets)
## Outputs -> outputs/plasmode/burden_chr22_maf001/ (same file layout).  The per-sample consistency check of the
## original (rare_chr22_per_sample.csv bins split at MAC 50, not 54) is replaced by a nesting check.
## Original header follows.
###############################################################################
###############################################################################
## plasmode_burden_chr22.R -- DELIVERABLE B rare-call channel (run ONCE via
## plasmode_burden_chr22.sbatch; NOT submitted by the build step).
##
## SeqArray pass over the chr22 aGDS (processed/r5_chr22.split.gds == agds_dir[22]),
## QC_label7 == "PASS" (node from revision/scripts/variant_qc_config.R), restricted to the
## plasmode pool (outputs/plasmode/plasmode_pool.csv), pool MAC with the minor-allele
## flip EXACTLY as adsp_center_batch_rare_chr22.R L83-91 (ac vs 2*nn - ac on the
## non-missing calls; missing calls count 0 minor alleles), and writes per-unit
## minor-allele-count burden matrices (units x N, float32 .npy, columns in
## plasmode_pool.csv row order) for these masks:
##
##   unit_type gene    : gene symbol = first token of annotation/info/FunctionalAnnotation/
##                       genecode_comprehensive_info (GENCODE.Info in Annotation_name_catalog),
##                       genic categories only (GENCODE.Category not intergenic/upstream/downstream)
##                       -- all PASS variants (SNV + indel)          masks MAC 1-5, 6-50, 51-500
##   unit_type win100k : 100-kb windows floor(position / 1e5), ALL PASS variants (intergenic
##                       included)                                     masks MAC 6-50, 51-500
##   unit_type coding  : gene, PASS SNV with STAARpipeline's lof.in.coding categories
##                       (GENCODE.EXONIC.Category stopgain/stoploss/nonsynonymous SNV/synonymous SNV
##                       or GENCODE.Category splicing/exonic;splicing/ncRNA_splicing/
##                       ncRNA_exonic;splicing)                        masks MAC 1-5, 6-50, 51-500
##   unit_type syn     : gene, PASS SNV with GENCODE.EXONIC.Category == "synonymous SNV"
##                                                                     masks MAC 6-50, 51-500
##
## Outputs (outputs/plasmode/burden_chr22/):
##   burden_<mask>.npy          units x N float32 (C order), rows = gene_table_<mask>.csv rows
##   gene_table_<mask>.csv      unit_index, unit, unit_type, n_variants, cMAC, cMAF (= cMAC / 2N),
##                              pos_min, pos_max
##   masks_manifest.csv         mask, unit_type, mac_lo, mac_hi, n_units, n_variants, files
##   annotation_categories.csv  GENCODE.Category / EXONIC.Category tables over PASS variants
##   burden_build.json          nodes used, counts, checks, timing
## Consistency check (full run only): column sums of win100k MAC 6-50 / 51-500 must equal
## MAC6_10 + MAC11_50 / MAC51_500 of rare_chr22_per_sample.csv (same pool, same flip rule).
## Env: NCPU (8), BSIZE (2000), MAXVAR (0 = all; >0 = first MAXVAR PASS variants, smoke), TAG.
###############################################################################
suppressPackageStartupMessages({ library(SeqArray); library(data.table); library(parallel); library(jsonlite) })
REPO  <- "${ADSP_ROOT}"
CHR   <- Sys.getenv("CHR", "22")                                   # per-chromosome build (2026-09-27): CHR from the environment
GDS   <- file.path(REPO, sprintf("processed/r5_chr%s.split.gds", CHR))
POOL  <- Sys.getenv("POOL", file.path(REPO, "outputs/plasmode/plasmode_pool.csv"))
PERS  <- file.path(REPO, sprintf("revision/tables/adsp_center_batch/rare_chr%s_per_sample.csv", CHR))   # not used by this script
OUT   <- Sys.getenv("OUTDIR", file.path(REPO, sprintf("outputs/plasmode/burden_chr%s_maf001", CHR)))
NCPU  <- as.integer(Sys.getenv("NCPU", "8"))
BSIZE <- as.integer(Sys.getenv("BSIZE", "2000"))
MAXVAR <- as.integer(Sys.getenv("MAXVAR", "0"))
TAG   <- Sys.getenv("TAG", "")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
source(file.path(REPO, "revision/scripts/variant_qc_config.R"))
QC_NODE <- qc_label()
ANNO <- "annotation/info/FunctionalAnnotation"
NODE_CAT  <- paste0(ANNO, "/genecode_comprehensive_category")          # GENCODE.Category
NODE_EXC  <- paste0(ANNO, "/genecode_comprehensive_exonic_category")   # GENCODE.EXONIC.Category
NODE_INFO <- paste0(ANNO, "/genecode_comprehensive_info")              # GENCODE.Info
t0 <- Sys.time()
msg <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), sprintf(...), "\n", sep = "")
info <- list(gds = GDS, qc_node = QC_NODE, nodes = c(NODE_CAT, NODE_EXC, NODE_INFO), ncpu = NCPU, bsize = BSIZE, maxvar = MAXVAR)

## ---------- .npy writer (float32, C order) ----------
write_npy_f32 <- function(path, M) {
  stopifnot(is.matrix(M))
  hdr <- sprintf("{'descr': '<f4', 'fortran_order': False, 'shape': (%d, %d), }", nrow(M), ncol(M))
  pad <- 64 - ((10 + nchar(hdr) + 1) %% 64); if (pad == 64) pad <- 0
  hdr <- paste0(hdr, strrep(" ", pad), "\n")
  con <- file(path, "wb"); on.exit(close(con))
  writeBin(as.raw(c(0x93, charToRaw("NUMPY"), 0x01, 0x00)), con)
  writeBin(as.integer(nchar(hdr)), con, size = 2, endian = "little")
  writeBin(charToRaw(hdr), con)
  v <- as.numeric(t(M))                       # C order: row after row
  writeBin(v, con, size = 4, endian = "little")
  invisible(path)
}

## ---------- pool ----------
pool <- fread(POOL, colClasses = list(character = "IID"))
N <- nrow(pool)
g <- seqOpen(GDS, readonly = TRUE)
sid <- seqGetData(g, "sample.id")
stopifnot(all(pool$IID %in% sid))
seqSetFilter(g, sample.id = pool$IID, verbose = FALSE)
sid_sel <- seqGetData(g, "sample.id")
stopifnot(length(sid_sel) == N)
col_order <- match(pool$IID, sid_sel)          # burden column j (GDS order) -> pool row
stopifnot(!anyNA(col_order))
msg("pool N = %d, all present in GDS", N)

## ---------- variants: QC7 PASS ----------
q7 <- seqGetData(g, QC_NODE)
pass <- !is.na(q7) & as.character(q7) == "PASS"
msg("%s PASS: %d of %d (%.2f%%)", QC_NODE, sum(pass), length(pass), 100 * mean(pass))
vid_all <- seqGetData(g, "variant.id")
vsel <- vid_all[pass]
if (MAXVAR > 0 && length(vsel) > MAXVAR) { vsel <- vsel[seq_len(MAXVAR)]; msg("SMOKE: first %d PASS variants", MAXVAR) }
rm(q7)
seqSetFilter(g, variant.id = vsel, verbose = FALSE)
nvar <- length(vsel)
## absolute index (position in the whole file) -> position in the selection
vid_sel <- seqGetData(g, "variant.id")
stopifnot(identical(vid_sel, vsel))
sel_pos <- integer(length(vid_all)); sel_pos[vid_all %in% vsel] <- seq_len(nvar)  # variant.id == 1:n in these files
stopifnot(all(vid_all == seq_along(vid_all)))
info$n_pass <- sum(pass); info$n_selected <- nvar

## ---------- annotation -> unit codes ----------
position <- as.integer(seqGetData(g, "position"))
snv <- isSNV(g)
CA <- as.character(seqGetData(g, NODE_CAT))
EX <- as.character(seqGetData(g, NODE_EXC))
IN <- as.character(seqGetData(g, NODE_INFO))
msg("annotation read: %d variants (%.1f min)", nvar, as.numeric(difftime(Sys.time(), t0, units = "mins")))
cat_tab <- as.data.table(table(GENCODE.Category = CA, useNA = "ifany"))[order(-N)]
exc_tab <- as.data.table(table(GENCODE.EXONIC.Category = EX, useNA = "ifany"))[order(-N)]
fwrite(rbind(cat_tab[, .(node = "GENCODE.Category", level = GENCODE.Category, N)],
             exc_tab[, .(node = "GENCODE.EXONIC.Category", level = GENCODE.EXONIC.Category, N)]),
       file.path(OUT, sprintf("annotation_categories%s.csv", TAG)))
print(cat_tab); print(exc_tab)
NONGENIC <- c("intergenic", "upstream", "downstream", "upstream;downstream", "")
genic <- !is.na(CA) & !(CA %in% NONGENIC)
gene <- sub("\\(.*$", "", sub("[,;].*$", "", IN))          # first token, strip "(dist=..)"
gene[!genic | is.na(gene) | gene %in% c("", "NONE", "UNKNOWN", "NA")] <- NA_character_
lof.in.coding <- (EX == "stopgain") | (EX == "stoploss") | (CA == "splicing") |
  (CA == "exonic;splicing") | (CA == "ncRNA_splicing") | (CA == "ncRNA_exonic;splicing") |
  (EX == "nonsynonymous SNV") | (EX == "synonymous SNV")
lof.in.coding[is.na(lof.in.coding)] <- FALSE
is_coding <- snv & lof.in.coding & !is.na(gene)
is_syn <- snv & !is.na(EX) & EX == "synonymous SNV" & !is.na(gene)
rm(CA, EX, IN); invisible(gc())

unit_of <- function(lab, keep) { f <- factor(ifelse(keep, lab, NA_character_)); list(code = as.integer(f), levels = levels(f)) }
U_gene   <- unit_of(gene, !is.na(gene))
U_win    <- unit_of(sprintf("chr%s_win%05d", CHR, position %/% 100000L), rep(TRUE, nvar))
U_coding <- unit_of(gene, is_coding)
U_syn    <- unit_of(gene, is_syn)
UNITS <- list(gene = U_gene, win100k = U_win, coding = U_coding, syn = U_syn)
for (u in names(UNITS)) msg("unit_type %-8s: %d units, %d variants assigned", u, length(UNITS[[u]]$levels), sum(!is.na(UNITS[[u]]$code)))
info$n_variants_assigned <- sapply(UNITS, function(u) sum(!is.na(u$code)))
info$n_units <- sapply(UNITS, function(u) length(u$levels))
info$n_snv <- sum(snv); info$n_genic <- sum(!is.na(gene)); info$n_coding <- sum(is_coding); info$n_syn <- sum(is_syn)

MASKS <- rbindlist(list(
  ## MAF < 0.1 % in the 27,306-subject pool  <=>  MAC / (2 x 27,306) < 0.001  <=>  MAC <= 54 (user rule, 2026-09-22)
  data.table(mask = "gene_maf001",      unit_type = "gene",    lo = 1L,  hi = 54L),
  data.table(mask = "gene_mac6_54",     unit_type = "gene",    lo = 6L,  hi = 54L),
  data.table(mask = "coding_maf001",    unit_type = "coding",  lo = 1L,  hi = 54L),
  data.table(mask = "coding_mac6_54",   unit_type = "coding",  lo = 6L,  hi = 54L),
  data.table(mask = "syn_maf001",       unit_type = "syn",     lo = 1L,  hi = 54L),
  data.table(mask = "win100k_maf001",   unit_type = "win100k", lo = 1L,  hi = 54L)))
CODES <- lapply(UNITS, function(u) u$code)     # integer vectors, length nvar (NA = not assigned)
NU <- sapply(UNITS, function(u) length(u$levels))
rm(gene, is_coding, is_syn, lof.in.coding, snv); invisible(gc())

## ---------- worker ----------
worker <- function(f) {
  nm <- MASKS$mask
  acc  <- lapply(seq_len(nrow(MASKS)), function(i) matrix(0L, nrow = NU[[MASKS$unit_type[i]]], ncol = N))
  cmac <- lapply(seq_len(nrow(MASKS)), function(i) numeric(NU[[MASKS$unit_type[i]]]))
  nvr  <- lapply(seq_len(nrow(MASKS)), function(i) integer(NU[[MASKS$unit_type[i]]]))
  pmn  <- lapply(seq_len(nrow(MASKS)), function(i) rep(NA_integer_, NU[[MASKS$unit_type[i]]]))
  pmx  <- lapply(seq_len(nrow(MASKS)), function(i) rep(NA_integer_, NU[[MASKS$unit_type[i]]]))
  nblk <- 0L
  ## seqBlockApply passes a SCALAR block-start index (verified on the SeqArray example file, 2026-09-20),
  ## not one index per variant; the PASS selection is non-contiguous in absolute index space, so the
  ## block's variants are start + 0:(ncol(x)-1) only in RELATIVE terms. Under seqParallel each worker
  ## holds a sub-selection, so map this worker's selection to positions in the GLOBAL selection first.
  vid_w  <- seqGetData(f, "variant.id")               # this worker's selected variants, in order
  gpos_w <- sel_pos[vid_w]                            # -> position in the global PASS selection
  stopifnot(length(gpos_w) > 0L, all(gpos_w > 0L))
  seqBlockApply(f, "$dosage_alt", margin = "by.variant", as.is = "none", bsize = BSIZE, var.index = "relative",
    FUN = function(idx, x) {
      ## x: N x nvar_block integer ALT dosage (NA = missing); idx: RELATIVE index of the block's first variant
      gi <- gpos_w[idx + seq_len(ncol(x)) - 1L]       # positions in the global selected PASS set
      stopifnot(length(gi) == ncol(x), !anyNA(gi), nrow(x) == N)
      miss <- is.na(x)
      nn   <- colSums(!miss)
      ac   <- colSums(x, na.rm = TRUE)
      flip <- ac > (2 * nn - ac)                      # minor allele is REF  (rare_chr22 L87-91)
      mac  <- ifelse(flip, 2 * nn - ac, ac)
      if (any(flip)) x[, flip] <- 2L - x[, flip, drop = FALSE]
      x[miss] <- 0L                                    # missing call = 0 minor alleles (rowSums na.rm)
      pos <- position[gi]
      for (i in seq_len(nrow(MASKS))) {
        code <- CODES[[MASKS$unit_type[i]]][gi]
        stopifnot(length(code) == length(mac))
        s <- which(!is.na(code) & mac >= MASKS$lo[i] & mac <= MASKS$hi[i])
        if (!length(s)) next
        u <- code[s]
        rs <- rowsum(t(x[, s, drop = FALSE]), u, reorder = TRUE)
        uu <- as.integer(rownames(rs)); storage.mode(rs) <- "integer"
        acc[[i]][uu, ] <<- acc[[i]][uu, , drop = FALSE] + rs
        cmac[[i]][uu] <<- cmac[[i]][uu] + as.numeric(rowsum(mac[s], u, reorder = TRUE)[, 1])
        nvr[[i]][uu]  <<- nvr[[i]][uu] + as.integer(table(factor(u, levels = uu)))
        pmn[[i]][uu]  <<- pmin(pmn[[i]][uu], as.integer(rowsum(pos[s], u, reorder = TRUE)[, 1] * 0 +
                                                     tapply(pos[s], u, min)), na.rm = TRUE)
        pmx[[i]][uu]  <<- pmax(pmx[[i]][uu], as.integer(tapply(pos[s], u, max)), na.rm = TRUE)
      }
      nblk <<- nblk + 1L
      if (nblk %% 100L == 0L) cat(sprintf("  worker pid %d: %d blocks\n", Sys.getpid(), nblk))
      invisible(NULL)
    })
  list(acc = acc, cmac = cmac, nvr = nvr, pmn = pmn, pmx = pmx, nblk = nblk)
}
combine <- function(a, b) {
  list(acc = mapply(function(x, y) x + y, a$acc, b$acc, SIMPLIFY = FALSE),
       cmac = mapply(`+`, a$cmac, b$cmac, SIMPLIFY = FALSE),
       nvr = mapply(`+`, a$nvr, b$nvr, SIMPLIFY = FALSE),
       pmn = mapply(function(x, y) pmin(x, y, na.rm = TRUE), a$pmn, b$pmn, SIMPLIFY = FALSE),
       pmx = mapply(function(x, y) pmax(x, y, na.rm = TRUE), a$pmx, b$pmx, SIMPLIFY = FALSE),
       nblk = a$nblk + b$nblk)
}
msg("starting genotype pass: %d variants x %d samples, %d workers, bsize %d", nvar, N, NCPU, BSIZE)
if (NCPU <= 1L) {
  ## 2026-09-20: serial path (NCPU=1) -- surfaces the real error message when the cluster
  ## path fails inside seqParallel/.DynamicClusterCall with an empty message (jobs 1643921, 1644112).
  msg("serial genotype pass (NCPU=1)")
  res <- worker(g)
} else {
  cl <- makeCluster(NCPU, outfile = "")            # outfile="" -> worker stderr/stdout reach the job log
  clusterEvalQ(cl, suppressPackageStartupMessages({ library(SeqArray); library(data.table) }))
  clusterExport(cl, c("N", "BSIZE", "MASKS", "CODES", "NU", "sel_pos", "position"))
  res <- seqParallel(cl, g, FUN = worker, split = "by.variant", .combine = combine)
  stopCluster(cl)
}
seqClose(g)
msg("genotype pass done: %d blocks (%.1f min)", res$nblk, as.numeric(difftime(Sys.time(), t0, units = "mins")))

## ---------- write ----------
man <- list(); checks <- list()
for (i in seq_len(nrow(MASKS))) {
  m <- MASKS$mask[i]; ut <- MASKS$unit_type[i]
  keep <- which(res$nvr[[i]] > 0)
  M <- res$acc[[i]][keep, col_order, drop = FALSE]        # columns -> plasmode_pool.csv row order
  gt <- data.table(unit_index = seq_along(keep) - 1L, unit = UNITS[[ut]]$levels[keep], unit_type = ut,
                   n_variants = res$nvr[[i]][keep], cMAC = res$cmac[[i]][keep],
                   cMAF = res$cmac[[i]][keep] / (2 * N), pos_min = res$pmn[[i]][keep], pos_max = res$pmx[[i]][keep])
  stopifnot(all(abs(rowSums(M) - gt$cMAC) < 1e-6))        # matrix row sums == cMAC
  fn_npy <- file.path(OUT, sprintf("burden_%s%s.npy", m, TAG)); fn_gt <- file.path(OUT, sprintf("gene_table_%s%s.csv", m, TAG))
  write_npy_f32(fn_npy, M); fwrite(gt, fn_gt)
  man[[i]] <- data.table(mask = m, unit_type = ut, mac_lo = MASKS$lo[i], mac_hi = MASKS$hi[i], n_units = nrow(gt),
                         n_variants = sum(gt$n_variants), total_cMAC = sum(gt$cMAC), npy = fn_npy, gene_table = fn_gt)
  msg("%-18s: %5d units, %8d variants, total cMAC %.0f -> %s", m, nrow(gt), sum(gt$n_variants), sum(gt$cMAC), basename(fn_npy))
  if (m %in% c("win100k_mac6_50", "win100k_mac51_500")) checks[[m]] <- colSums(M)
}
man <- rbindlist(man); fwrite(man, file.path(OUT, sprintf("masks_manifest%s.csv", TAG)))

## ---------- consistency: the MAF < 0.1 % masks nest (MAC 6-54 subset of MAC 1-54 per unit type) ----------
for (ut in c("gene", "coding")) {
  a <- fread(file.path(OUT, sprintf("gene_table_%s_maf001%s.csv", ut, TAG))); b <- fread(file.path(OUT, sprintf("gene_table_%s_mac6_54%s.csv", ut, TAG)))
  j <- merge(a[, .(unit, n_all = n_variants, c_all = cMAC)], b[, .(unit, n_sub = n_variants, c_sub = cMAC)], by = "unit")
  bad <- j[n_sub > n_all | c_sub > c_all + 1e-6]
  info[[sprintf("check_nesting_%s_bad_units", ut)]] <- nrow(bad)
  msg("nesting check %s: %d units in both masks, %d violate MAC 6-54 <= MAC 1-54", ut, nrow(j), nrow(bad))
}
info$minutes <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
writeLines(toJSON(info, auto_unbox = TRUE, pretty = TRUE), file.path(OUT, sprintf("burden_build%s.json", TAG)))
msg("DONE (%.1f min)", info$minutes)
