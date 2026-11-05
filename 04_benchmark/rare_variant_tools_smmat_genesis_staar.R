#!/usr/bin/env Rscript
###############################################################################
## rare_variant_tools_smmat_genesis_staar.R -- RARE-VARIANT tool comparison on the plasmode permuted null (2026-10-01):
## GMMAT::SMMAT (B/S/O/E), GENESIS::assocTestAggregate (Burden/SKAT/SMMAT) and STAAR::STAAR (STAAR-O without annotations)
## on the SAME gene x variant sets as the harness's coding_maf001 mask (every coding SNV with pool MAC 1-54, grouped by gene;
## plasmode_rv_groups.py gene_variants.csv.gz) and the SAME permuted phenotype Y* of one replicate, under one covariate set:
##   pooled_project  N5  Age + Sex + PC1-20 + centre + Study_cov (project collapsed as the harness), all 27,306 subjects
##   strat_project   S1  Age + Sex + PC1-20 + project x centre cell dummies, subjects of MIXED cells only (cases and controls)
##   saturated_cell  N0  Age + Sex + cohort x project x centre cell dummies, mixed cells only (the harness's oracle)
## Null model = logistic GLM (the sparse-GRM variance component is not identifiable on the KING-unrelated pool: both GMMAT
## and GENESIS fell back to the GLM in the common-variant tool runs, gmmat_genesis_single_variant.R), full-rank design columns.
## Args: --ystar <csv IID,Y> --covset <set> --gds <SeqArray GDS of the merged coding bed> --groups <gene_variants.csv.gz>
##       --out <prefix> [--which smmat,genesis,staar] [--threads 8] [--max-genes 0] [--mask coding|syn|nonsyn]
## Outputs: <out>_genes.tsv.gz (gene, chr, n_var, cMAC, p_<tool>_<test> ...), <out>_meta.json
###############################################################################
suppressPackageStartupMessages({ library(data.table); library(Matrix); library(jsonlite); library(SeqArray) })
getarg <- function(flag, default = NULL) { a <- commandArgs(TRUE); i <- match(flag, a); if (is.na(i) || i == length(a)) default else a[i + 1] }
REPO   <- "${ADSP_ROOT}"
POOL   <- file.path(REPO, "outputs/plasmode/plasmode_pool.csv")
YCSV <- getarg("--ystar"); COVSET <- getarg("--covset", "pooled_project"); GDS <- getarg("--gds"); GROUPS <- getarg("--groups")
OUT <- getarg("--out"); WHICH <- strsplit(getarg("--which", "smmat,genesis,staar"), ",")[[1]]
THREADS <- as.integer(getarg("--threads", "8")); MAXG <- as.integer(getarg("--max-genes", "0")); MASK <- getarg("--mask", "coding")
stopifnot(!is.null(YCSV), !is.null(GDS), !is.null(GROUPS), !is.null(OUT), COVSET %in% c("pooled_project", "strat_project", "saturated_cell"))
t_all <- proc.time()[3]
msg <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), sprintf(...), "\n", sep = "")
lam <- function(p) { p <- p[is.finite(p) & p > 0 & p <= 1]; if (!length(p)) NA_real_ else median(qchisq(1 - p, 1)) / qchisq(0.5, 1) }

## ---------------------------------------------------------------- phenotype + design --
pool <- fread(POOL, colClasses = list(character = "IID"))
y <- fread(YCSV, colClasses = list(character = "IID"))
stopifnot(identical(y$IID, pool$IID))
pool[, Y := as.integer(y$Y)]
PCS <- paste0("PC", 1:20)
cell_lab <- switch(COVSET, pooled_project = NULL, strat_project = paste(pool$Study, pool$center, sep = "|"), saturated_cell = as.character(pool$cell_A))
an <- !is.na(pool$Y)
if (!is.null(cell_lab)) {                      # mixed cells only, as the harness's mixed_rows()
  ct <- data.table(l = cell_lab, Y = pool$Y)[!is.na(Y), .(n1 = sum(Y), N = .N), by = l]
  mixed <- ct[n1 > 0 & n1 < N, l]
  an <- an & cell_lab %in% mixed
  pool[, cell_f := factor(ifelse(an, cell_lab, NA_character_))]
}
zs <- function(v) { v <- as.numeric(v); v[is.na(v)] <- mean(v[an], na.rm = TRUE); (v - mean(v[an])) / sd(v[an]) }
pool[, Age_f := zs(Age)]
pool[, Sex_f := as.integer(as.numeric(Sex))]
for (p in PCS) set(pool, j = p, value = zs(pool[[p]]))
if (COVSET == "pooled_project") {              # harness collapse rule on the analyzed rows under Y*: mono or N < 50 -> Other
  ct <- pool[an, .(n0 = sum(Y == 0), n1 = sum(Y == 1), N = .N), by = Study]
  other <- ct[n0 == 0 | n1 == 0 | N < 50, Study]
  pool[, Study_cov := ifelse(Study %in% other, "Other", Study)]
  pool[, center_f := factor(center)]; pool[, Study_f := factor(Study_cov)]
  rhs <- as.formula(paste("~ Age_f + Sex_f +", paste(PCS, collapse = " + "), "+ center_f + Study_f"))
} else if (COVSET == "strat_project") {
  rhs <- as.formula(paste("~ Age_f + Sex_f +", paste(PCS, collapse = " + "), "+ cell_f"))
} else rhs <- ~ Age_f + Sex_f + cell_f
Xall <- model.matrix(rhs, data = as.data.frame(pool[an]))
qrx <- qr(Xall); keep_cols <- sort(qrx$pivot[seq_len(qrx$rank)])
dropped <- setdiff(colnames(Xall), colnames(Xall)[keep_cols])
Xr <- Xall[, keep_cols, drop = FALSE]; Xr <- Xr[, colnames(Xr) != "(Intercept)", drop = FALSE]
cn <- make.names(colnames(Xr), unique = TRUE); colnames(Xr) <- cn
covdf <- as.data.frame(Xr, check.names = FALSE)
ph <- cbind(data.frame(IID = pool$IID[an], Y = pool$Y[an], stringsAsFactors = FALSE), covdf)
form <- as.formula(paste("Y ~", paste(cn, collapse = " + ")))
msg("covset %s: analyzed %d of %d subjects, %d cases; design %d columns (rank %d, dropped %d: %s)", COVSET, nrow(ph), nrow(pool), sum(ph$Y),
    ncol(Xall), qrx$rank, length(dropped), paste(head(dropped, 8), collapse = ","))
meta <- list(covset = COVSET, n_analyzed = nrow(ph), n_cases = sum(ph$Y), n_design_cols = ncol(Xall), design_rank = qrx$rank,
             dropped_columns = dropped, mask = MASK, threads = THREADS)

## ---------------------------------------------------------------- gene x variant sets --
gv <- fread(GROUPS, colClasses = list(character = c("gene", "chr", "snp", "cat", "ref", "alt")))
if (MASK == "syn") gv <- gv[cat == "syn"] else if (MASK == "nonsyn") gv <- gv[cat == "nonsyn"]
gv[, variant.id := bim_index + 1L]            # seqBED2GDS numbers variants 1..M in .bim order (checked below)
gt <- gv[, .(chr = chr[1], n_var = .N, cMAC = sum(mac), pos_min = min(pos), pos_max = max(pos)), by = gene]
if (MAXG > 0) { gt <- gt[seq_len(min(MAXG, nrow(gt)))]; gv <- gv[gene %in% gt$gene] }
msg("mask %s: %d genes, %d variants", MASK, nrow(gt), nrow(gv))
g <- seqOpen(GDS)
sid <- seqGetData(g, "sample.id"); sid <- sub("^0_", "", sid)          # seqBED2GDS may compose FID_IID; the fam FID is 0
stopifnot(all(ph$IID %in% sid))
vid_all <- seqGetData(g, "variant.id"); pos_all <- seqGetData(g, "position"); chr_all <- as.character(seqGetData(g, "chromosome"))
stopifnot(identical(as.integer(vid_all), seq_along(vid_all)))
chk <- gv[sample(.N, min(.N, 2000L))]
stopifnot(all(pos_all[chk$variant.id] == chk$pos), all(chr_all[chk$variant.id] == chk$chr))
seqClose(g)
res <- copy(gt)
SCR <- Sys.getenv("SIM_SCRATCH", tempdir())
dir.create(SCR, showWarnings = FALSE, recursive = TRUE)

## ---------------------------------------------------------------- GMMAT::SMMAT --
if ("smmat" %in% WHICH) {
  suppressPackageStartupMessages(library(GMMAT))
  t0 <- proc.time()[3]
  obj <- glmmkin(fixed = form, data = ph, kins = NULL, id = "IID", family = binomial(link = "logit"))
  msg("GMMAT glmmkin (GLM, no kinship): converged %s (%.0f s)", obj$converged, proc.time()[3] - t0)
  grp <- file.path(SCR, sprintf("smmat_%s.group", MASK))
  fwrite(gv[, .(gene, chr, pos, ref, alt, weight = 1)], grp, sep = "\t", col.names = FALSE)
  t1 <- proc.time()[3]
  sm <- tryCatch(SMMAT(obj, geno.file = GDS, group.file = grp, MAF.range = c(1e-7, 0.5), MAF.weights.beta = c(1, 25), miss.cutoff = 1,
                       method = "davies", tests = c("B", "S", "O", "E"), auto.flip = TRUE, ncores = THREADS, verbose = FALSE),
                 error = function(e) { msg("SMMAT ERROR: %s", conditionMessage(e)); NULL })
  if (!is.null(sm)) {
    sm <- as.data.table(sm)
    pc <- grep("^(B\\.pval|S\\.pval|O\\.pval|E\\.pval)$", names(sm), value = TRUE)
    res <- merge(res, sm[, c("group", "n.variants", pc), with = FALSE], by.x = "gene", by.y = "group", all.x = TRUE, sort = FALSE)
    setnames(res, c("n.variants", pc), c("n_var_smmat", paste0("p_smmat_", sub("\\.pval", "", pc))))
    meta$smmat <- list(seconds = proc.time()[3] - t1, n_genes = sum(is.finite(res$p_smmat_O)),
                       lambda = sapply(paste0("p_smmat_", sub("\\.pval", "", pc)), function(cc) lam(res[[cc]])))
    msg("SMMAT: %d genes; lambda B %.3f S %.3f O %.3f E %.3f (%.0f s)", meta$smmat$n_genes, lam(res$p_smmat_B), lam(res$p_smmat_S), lam(res$p_smmat_O), lam(res$p_smmat_E), proc.time()[3] - t1)
  }
  rm(obj); invisible(gc())
  fwrite(res, paste0(OUT, "_genes.tsv.gz"), sep = "\t")
}

## ---------------------------------------------------------------- GENESIS::assocTestAggregate --
if ("genesis" %in% WHICH) {
  suppressPackageStartupMessages({ library(GENESIS); library(SeqVarTools); library(Biobase); library(GenomicRanges); library(BiocParallel) })
  t0 <- proc.time()[3]
  g <- seqOpen(GDS)
  sid_g <- seqGetData(g, "sample.id")
  sdf <- data.frame(sample.id = sid_g, stringsAsFactors = FALSE)
  m <- match(sub("^0_", "", sid_g), ph$IID)
  sdf$Y <- ph$Y[m]
  for (cc in cn) sdf[[cc]] <- covdf[[cc]][m]
  seqData <- SeqVarData(g, sampleData = AnnotatedDataFrame(sdf))
  nullmod <- tryCatch(fitNullModel(seqData, outcome = "Y", covars = cn, cov.mat = NULL, family = "binomial", sample.id = sid_g[!is.na(m)], verbose = FALSE),
                      error = function(e) { msg("GENESIS fitNullModel ERROR: %s", conditionMessage(e)); NULL })
  if (!is.null(nullmod)) {
    msg("GENESIS fitNullModel (GLM): converged %s (%.0f s)", isTRUE(nullmod$converged), proc.time()[3] - t0)
    t1 <- proc.time()[3]
    gr <- GRanges(seqnames = gv$chr, ranges = IRanges(start = gv$pos, width = 1L))
    grl <- split(gr, factor(gv$gene, levels = gt$gene))
    it <- SeqVarListIterator(seqData, variantRanges = grl, verbose = FALSE)
    nw <- max(1L, min(THREADS, 8L))
    assoc <- tryCatch(assocTestAggregate(it, null.model = nullmod, AF.max = 1, weight.beta = c(1, 25), test = "SMMAT", sparse = TRUE,
                                         BPPARAM = if (nw > 1) MulticoreParam(nw) else SerialParam(), verbose = FALSE),
                      error = function(e) { msg("GENESIS assocTestAggregate ERROR: %s", conditionMessage(e)); NULL })
    if (!is.null(assoc)) {
      ar <- as.data.table(assoc$results); ar[, gene := names(assoc$variantInfo)]
      pc <- grep("^pval", names(ar), value = TRUE); pc <- pc[!grepl("method", pc)]   # pval_theta.method is a character column
      keep <- c("gene", "n.site", pc)
      ar <- ar[, keep, with = FALSE]
      setnames(ar, c("n.site", pc), c("n_var_genesis", paste0("p_genesis_", sub("^pval_?", "", pc))))
      res <- merge(res, ar, by = "gene", all.x = TRUE, sort = FALSE)
      pcols <- grep("^p_genesis_", names(res), value = TRUE)
      meta$genesis <- list(seconds = proc.time()[3] - t1, n_genes = sum(is.finite(res[[pcols[1]]])), lambda = sapply(pcols, function(cc) lam(res[[cc]])))
      msg("GENESIS SMMAT: %d genes; lambda %s (%.0f s)", meta$genesis$n_genes, paste(sprintf("%s %.3f", sub("p_genesis_", "", pcols), meta$genesis$lambda), collapse = ", "), proc.time()[3] - t1)
    }
  }
  seqClose(g); rm(nullmod); invisible(gc())
  fwrite(res, paste0(OUT, "_genes.tsv.gz"), sep = "\t")
}

## ---------------------------------------------------------------- STAAR-O (no annotation) --
if ("staar" %in% WHICH) {
  suppressPackageStartupMessages({ library(STAAR); library(parallel) })
  t0 <- proc.time()[3]
  nullglm <- fit_null_glm(form, data = ph, family = binomial(link = "logit"))
  msg("STAAR fit_null_glm: converged %s (%.0f s)", nullglm$converged, proc.time()[3] - t0)
  t1 <- proc.time()[3]
  gl <- split(gv$variant.id, factor(gv$gene, levels = gt$gene))
  genes <- gt$gene
  chunks <- split(seq_along(genes), (seq_along(genes) - 1L) %% max(1L, THREADS))
  getp <- function(x) { x <- unlist(x); if (length(x)) as.numeric(x[1]) else NA_real_ }
  worker <- function(ix) {
    gw <- seqOpen(GDS)
    on.exit(seqClose(gw))
    seqSetFilter(gw, sample.id = seqGetData(gw, "sample.id")[match(ph$IID, sub("^0_", "", seqGetData(gw, "sample.id")))], verbose = FALSE)
    sid_w <- sub("^0_", "", seqGetData(gw, "sample.id")); ord <- match(ph$IID, sid_w); stopifnot(!anyNA(ord))
    out <- vector("list", length(ix))
    for (k in seq_along(ix)) {
      vids <- gl[[genes[ix[k]]]]
      r <- c(p_staar_O = NA_real_, p_staar_S_1_25 = NA_real_, p_staar_B_1_25 = NA_real_, p_staar_A_1_25 = NA_real_, p_staar_S_1_1 = NA_real_, p_staar_B_1_1 = NA_real_, p_staar_A_1_1 = NA_real_, n_var_staar = length(vids))
      if (length(vids) >= 2) {
        seqSetFilter(gw, variant.id = vids, verbose = FALSE)
        G <- seqGetData(gw, "$dosage_alt")[ord, , drop = FALSE]
        G[is.na(G)] <- 0
        st <- tryCatch(STAAR(G, nullglm, rare_maf_cutoff = 1, rv_num_cutoff = 2), error = function(e) NULL)
        if (!is.null(st)) {   # STAAR returns the per-test p as 1-row data frames (SKAT(1,25), STAAR-S(1,25) identical w/o annotation)
          r[1:7] <- c(getp(st$results_STAAR_O), getp(st$results_STAAR_S_1_25), getp(st$results_STAAR_B_1_25), getp(st$results_STAAR_A_1_25),
                      getp(st$results_STAAR_S_1_1), getp(st$results_STAAR_B_1_1), getp(st$results_STAAR_A_1_1))
          r["n_var_staar"] <- as.numeric(st$num_variant)
        }
      }
      out[[k]] <- r
    }
    data.table(gene = genes[ix], do.call(rbind, out))
  }
  sr <- mclapply(chunks, worker, mc.cores = max(1L, THREADS), mc.preschedule = FALSE)
  bad <- sapply(sr, function(x) inherits(x, "try-error") || is.null(x))
  if (any(bad)) msg("STAAR: %d worker chunk(s) failed", sum(bad))
  sr <- rbindlist(sr[!bad])
  res <- merge(res, sr, by = "gene", all.x = TRUE, sort = FALSE)
  meta$staar <- list(seconds = proc.time()[3] - t1, n_genes = sum(is.finite(res$p_staar_O)),
                     lambda = sapply(grep("^p_staar_", names(res), value = TRUE), function(cc) lam(res[[cc]])))
  msg("STAAR: %d genes; lambda O %.3f S(1,25) %.3f B(1,25) %.3f A(1,25) %.3f (%.0f s)", meta$staar$n_genes, lam(res$p_staar_O), lam(res$p_staar_S_1_25), lam(res$p_staar_B_1_25), lam(res$p_staar_A_1_25), proc.time()[3] - t1)
}
setorder(res, chr, pos_min)
fwrite(res, paste0(OUT, "_genes.tsv.gz"), sep = "\t")
meta$seconds <- proc.time()[3] - t_all
writeLines(toJSON(meta, auto_unbox = TRUE, pretty = TRUE, digits = NA), paste0(OUT, "_meta.json"))
msg("DONE in %.0f s -> %s_genes.tsv.gz (%d genes, %d columns)", meta$seconds, OUT, nrow(res), ncol(res))
