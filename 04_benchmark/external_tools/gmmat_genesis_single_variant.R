#!/usr/bin/env Rscript
###############################################################################
## gmmat_genesis_single_variant.R -- GMMAT and GENESIS on the REAL plasmode genotypes for one replicate's phenotype
## (large common-variant simulation, 2026-09-23).  "As deployed": the null model is a logistic GLMM with the
## fastSparseGRM kinship matrix restricted to the pool (27,306 x 27,306 with 746,494 off-diagonal entries: the KING-unrelated
## pool still holds pairs above the 4th-degree fastSparseGRM threshold), fixed effects Age + Sex + PC1-20 + centre + Study_cov
## (sequencing project collapsed exactly as the harness's pooled_project arm: single-label or N < 50 -> Other, on
## the analyzed rows under this replicate's Y).
##   GMMAT 1.4.2  : glmmkin(binomial) -> glmm.score on the PLINK bed (plain score p; GMMAT has no SPA)
##   GENESIS 2.36 : fitNullModel(cov.mat = sparse GRM) -> assocTestSingle(test = "Score.SPA") on the GDS
## Args: --ystar <csv: IID,Y; Y = NA for subjects dropped by the generator> --out <prefix> [--which gmmat,genesis]
##       [--threads 8] [--maxiter 50]
## Outputs: <out>_pvals.tsv.gz (snp, p_gmmat, p_genesis_spa, p_genesis_score, in .bim order), <out>_meta.json
###############################################################################
suppressPackageStartupMessages({ library(data.table); library(Matrix); library(jsonlite) })
getarg <- function(flag, default = NULL) { a <- commandArgs(TRUE); i <- match(flag, a); if (is.na(i) || i == length(a)) default else a[i + 1] }
REPO   <- "${ADSP_ROOT}"
POOL   <- file.path(REPO, "outputs/plasmode/plasmode_pool.csv")
BFILE  <- file.path(REPO, "outputs/plasmode/pool_281k")
GDS    <- file.path(REPO, "outputs/plasmode/pool_281k.gds")          # SNPRelate GDS (GENESIS)
SEQGDS <- file.path(REPO, "outputs/plasmode/pool_281k.seq.gds")      # SeqArray GDS (GMMAT glmm.score, parallel reader)
GRM    <- file.path(REPO, "outputs/plasmode/pool_281k_sparseGRM.rds")
YCSV   <- getarg("--ystar"); OUT <- getarg("--out"); WHICH <- strsplit(getarg("--which", "gmmat,genesis"), ",")[[1]]
THREADS <- as.integer(getarg("--threads", "8")); MAXITER <- as.integer(getarg("--maxiter", "100"))
MAXITER_GMMAT <- as.integer(getarg("--maxiter-gmmat", "30"))   # the AI-REML never converges on this pool; fail fast, then the GLM fallback
stopifnot(!is.null(YCSV), !is.null(OUT))
t_all <- proc.time()[3]
msg <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), sprintf(...), "\n", sep = "")

pool <- fread(POOL, colClasses = list(character = "IID"))
y <- fread(YCSV, colClasses = list(character = "IID"))
stopifnot(identical(y$IID, pool$IID))
pool[, Y := as.integer(y$Y)]
bim <- fread(paste0(BFILE, ".bim"), header = FALSE, colClasses = list(character = "V2"))
M <- nrow(bim)
PCS <- paste0("PC", 1:20)
an <- !is.na(pool$Y)
## Study_cov: the harness's collapse rule on the analyzed rows under Y (plasmode_null_harness.collapse_levels, min_n 50)
ct <- pool[an, .(n0 = sum(Y == 0), n1 = sum(Y == 1), N = .N), by = Study]
other <- ct[n0 == 0 | n1 == 0 | N < 50, Study]
pool[, Study_cov := ifelse(Study %in% other, "Other", Study)]
## quantitative covariates standardized on the analyzed rows (the raw PCs are ~1e-3 in scale: coefficients ~150 and a
## reciprocal condition number ~1e-18 in the tools' working-weight solves, 2026-09-23 smoke); the fitted model is unchanged
zs <- function(v) { v <- as.numeric(v); v[is.na(v)] <- mean(v[an], na.rm = TRUE); (v - mean(v[an])) / sd(v[an]) }
pool[, Age_f := zs(Age)]
pool[, Sex_f := as.integer(as.numeric(Sex))]
for (p in PCS) set(pool, j = p, value = zs(pool[[p]]))
pool[, center_f := factor(center)]
pool[, Study_f := factor(Study_cov)]
ph <- as.data.frame(pool[an, c("IID", "Y", "Age_f", "Sex_f", PCS, "center_f", "Study_f"), with = FALSE])
## Full-rank design: centre and the collapsed project are nested for some levels, so the two dummy sets are collinear.
## The harness and REGENIE drop the dependent columns; GMMAT / GENESIS refuse a rank-deficient design.  Build the model
## matrix on the WHOLE pool (covariates have no NA), keep the QR-selected columns, and hand both tools numeric covariates.
rhs <- as.formula(paste("~ Age_f + Sex_f +", paste(PCS, collapse = " + "), "+ center_f + Study_f"))
Xall <- model.matrix(rhs, data = as.data.frame(pool))
stopifnot(nrow(Xall) == nrow(pool))
qrx <- qr(Xall[an, , drop = FALSE])
keep_cols <- sort(qrx$pivot[seq_len(qrx$rank)])
dropped <- setdiff(colnames(Xall), colnames(Xall)[keep_cols])
Xr <- Xall[, keep_cols, drop = FALSE]
Xr <- Xr[, colnames(Xr) != "(Intercept)", drop = FALSE]
cn <- make.names(colnames(Xr), unique = TRUE)
colnames(Xr) <- cn
msg("design: %d columns, rank %d, dropped %d collinear column(s): %s", ncol(Xall), qrx$rank, length(dropped),
    paste(head(dropped, 12), collapse = ", "))
covdf <- as.data.frame(Xr, check.names = FALSE)
ph <- cbind(ph[, c("IID", "Y")], covdf[an, , drop = FALSE])
form <- as.formula(paste("Y ~", paste(cn, collapse = " + ")))
K <- readRDS(GRM)
K <- forceSymmetric(K[ph$IID, ph$IID])          # dsCMatrix: the class GMMAT routes to its sparse solver (a dgCMatrix hits "Not a matrix.")
msg("analyzed %d of %d subjects, %d cases; Study_cov %d levels (%d -> Other); GRM %d x %d, %d off-diagonal nonzeros",
    nrow(ph), nrow(pool), sum(ph$Y), nlevels(pool$Study_f), length(other), nrow(K), ncol(K), length(K@x) - nrow(K))
meta <- list(n_analyzed = nrow(ph), n_cases = sum(ph$Y), n_study_levels = nlevels(pool$Study_f), n_study_other = length(other),
             grm_offdiag_nnz = length(K@x) - nrow(K), covariates = "Age+Sex+PC1-20+centre+Study_cov (full-rank columns)",
             n_design_cols = ncol(Xall), design_rank = qrx$rank, dropped_columns = dropped, maxiter = MAXITER)
p_gmmat <- rep(NA_real_, M); p_gspa <- rep(NA_real_, M); p_gsc <- rep(NA_real_, M)
lam <- function(p) { p <- p[is.finite(p) & p > 0 & p <= 1]; if (!length(p)) NA_real_ else median(qchisq(1 - p, 1)) / qchisq(0.5, 1) }

## ---------------------------------------------------------------- GMMAT --
if ("gmmat" %in% WHICH) {
  suppressPackageStartupMessages(library(GMMAT))
  t0 <- proc.time()[3]
  fallback <- FALSE
  obj <- tryCatch(glmmkin(fixed = form, data = ph, kins = K, id = "IID", family = binomial(link = "logit"), maxiter = MAXITER_GMMAT),
                  error = function(e) { msg("GMMAT glmmkin ERROR: %s", conditionMessage(e)); NULL })
  if (is.null(obj)) {
    ## AI-REML did not converge and GMMAT's Brent fallback supports dense kinship only ("Not a matrix." on a sparse GRM):
    ## refit the null WITHOUT the kinship term (what glmmkin returns when theta -> 0) and record the fallback.
    fallback <- TRUE
    obj <- tryCatch(glmmkin(fixed = form, data = ph, kins = NULL, id = "IID", family = binomial(link = "logit"), maxiter = MAXITER),
                    error = function(e) { msg("GMMAT glmmkin (kins = NULL) ERROR: %s", conditionMessage(e)); NULL })
    if (!is.null(obj)) msg("GMMAT: refitted without the kinship term (fallback_glm)")
  }
  if (!is.null(obj)) {
    meta$gmmat <- list(converged = isTRUE(obj$converged), fallback_glm = fallback, theta = unname(as.numeric(obj$theta)), null_fit_s = proc.time()[3] - t0,
                       max_abs_coef = max(abs(obj$coefficients)))
    msg("GMMAT glmmkin: converged %s, theta %s, max|coef| %.2f (%.0f s)", obj$converged, paste(signif(obj$theta, 4), collapse = ","),
        max(abs(obj$coefficients)), proc.time()[3] - t0)
    t1 <- proc.time()[3]
    outf <- paste0(OUT, "_gmmat_score.txt")
    ok <- tryCatch({ glmm.score(obj, infile = SEQGDS, outfile = outf, MAF.range = c(1e-7, 0.5), miss.cutoff = 1,
                                nperbatch = 2000, ncores = THREADS); TRUE },   # SeqArray GDS: parallel reader (PLINK path is single-threaded, > 1 h)
                   error = function(e) { msg("GMMAT glmm.score ERROR: %s", conditionMessage(e)); FALSE })
    if (ok) {
      sc <- fread(outf, colClasses = list(character = "SNP"))
      i <- match(sc$SNP, bim$V2); k <- !is.na(i)
      p_gmmat[i[k]] <- as.numeric(sc$PVAL[k])
      meta$gmmat$score_s <- proc.time()[3] - t1; meta$gmmat$n_tested <- sum(is.finite(p_gmmat)); meta$gmmat$lambda_gc <- lam(p_gmmat)
      msg("GMMAT score: %d SNPs tested, lambda %.3f (%.0f s)", sum(is.finite(p_gmmat)), lam(p_gmmat), proc.time()[3] - t1)
      unlink(outf)
    }
  } else meta$gmmat <- list(converged = FALSE, error = "glmmkin failed")
  ## checkpoint: GMMAT's p-values survive a failure in the GENESIS block
  fwrite(data.table(snp = bim$V2, p_gmmat = p_gmmat, p_genesis_spa = p_gspa, p_genesis_score = p_gsc), paste0(OUT, "_pvals.tsv.gz"), sep = "\t")
  writeLines(toJSON(meta, auto_unbox = TRUE, pretty = TRUE, digits = NA), paste0(OUT, "_meta.json"))
}

## -------------------------------------------------------------- GENESIS --
if ("genesis" %in% WHICH) {
  suppressPackageStartupMessages({ library(GENESIS); library(GWASTools); library(BiocParallel) })
  t0 <- proc.time()[3]
  gr <- GdsGenotypeReader(GDS)
  sid <- as.character(getScanID(gr))
  stopifnot(identical(sid, pool$IID))
  sa <- cbind(data.frame(scanID = sid, Y = pool$Y, stringsAsFactors = FALSE), covdf)
  sad <- ScanAnnotationDataFrame(sa)
  gdata <- GenotypeData(gr, scanAnnot = sad)
  gfallback <- FALSE
  nullmod <- tryCatch(fitNullModel(sad, outcome = "Y", covars = cn,
                                   cov.mat = K, sample.id = ph$IID, family = "binomial", max.iter = MAXITER, verbose = FALSE),
                      error = function(e) { msg("GENESIS fitNullModel ERROR: %s", conditionMessage(e)); NULL })
  if (is.null(nullmod)) {
    ## the variance component of the sparse GRM is not identifiable on this unrelated pool (AI-REML oscillates or hits a
    ## non-PD matrix): refit without the kinship term, as for GMMAT, and record the fallback
    gfallback <- TRUE
    nullmod <- tryCatch(fitNullModel(sad, outcome = "Y", covars = cn, cov.mat = NULL, sample.id = ph$IID, family = "binomial",
                                     max.iter = MAXITER, verbose = FALSE),
                        error = function(e) { msg("GENESIS fitNullModel (cov.mat = NULL) ERROR: %s", conditionMessage(e)); NULL })
    if (!is.null(nullmod)) msg("GENESIS: refitted without the kinship term (fallback_glm)")
  }
  if (!is.null(nullmod)) {
    vc <- if (is.null(nullmod$varComp)) NA_real_ else unname(as.numeric(nullmod$varComp))   # NULL for the GLM refit
    meta$genesis <- list(converged = isTRUE(nullmod$converged), fallback_glm = gfallback, varComp = vc,
                         null_fit_s = proc.time()[3] - t0, max_abs_coef = max(abs(nullmod$fixef$Est)))
    msg("GENESIS fitNullModel: converged %s, fallback_glm %s, varComp %s, max|coef| %.2f (%.0f s)", isTRUE(nullmod$converged), gfallback,
        paste(signif(vc, 4), collapse = ","), max(abs(nullmod$fixef$Est)), proc.time()[3] - t0)
    t1 <- proc.time()[3]
    it <- GenotypeBlockIterator(gdata, snpBlock = 2000)      # 2026-09-23: 10,000-SNP blocks x 8 forked workers exceeded 24 GB
    nw <- min(THREADS, 4L)
    assoc <- tryCatch(assocTestSingle(it, null.model = nullmod, test = "Score.SPA",
                                      BPPARAM = if (nw > 1) MulticoreParam(nw) else SerialParam(), verbose = FALSE),
                      error = function(e) { msg("GENESIS assocTestSingle ERROR: %s", conditionMessage(e)); NULL })
    if (!is.null(assoc)) {
      i <- match(as.character(assoc$variant.id), bim$V2); k <- !is.na(i)      # variant.id = the .bim SNP names (snpgdsBED2GDS keeps them)
      p_gspa[i[k]] <- as.numeric(assoc$SPA.pval[k])
      p_gsc[i[k]] <- if ("Score.pval" %in% names(assoc)) as.numeric(assoc$Score.pval[k]) else pchisq(as.numeric(assoc$Score.Stat[k])^2, 1, lower.tail = FALSE)
      meta$genesis$score_s <- proc.time()[3] - t1; meta$genesis$n_tested <- sum(is.finite(p_gspa))
      meta$genesis$lambda_gc_spa <- lam(p_gspa); meta$genesis$lambda_gc_score <- lam(p_gsc)
      meta$genesis$spa_converged_share <- mean(assoc$SPA.converged[k], na.rm = TRUE)
      msg("GENESIS Score.SPA: %d SNPs tested, lambda SPA %.3f / score %.3f (%.0f s)", sum(is.finite(p_gspa)), lam(p_gspa), lam(p_gsc), proc.time()[3] - t1)
    }
  } else meta$genesis <- list(converged = FALSE, error = "fitNullModel failed")
  close(gr)
}
fwrite(data.table(snp = bim$V2, p_gmmat = p_gmmat, p_genesis_spa = p_gspa, p_genesis_score = p_gsc), paste0(OUT, "_pvals.tsv.gz"), sep = "\t")
meta$seconds <- proc.time()[3] - t_all
writeLines(toJSON(meta, auto_unbox = TRUE, pretty = TRUE, digits = NA), paste0(OUT, "_meta.json"))
msg("DONE in %.0f s -> %s_pvals.tsv.gz", meta$seconds, OUT)
