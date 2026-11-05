#!/usr/bin/env Rscript
###############################################################################
## gmmat_genesis_single_variant.R -- GMMAT and GENESIS on the REAL plasmode genotypes for one replicate's phenotype
## (isolated convergence repair, 2026-10-15).  "As deployed": the null model is a logistic GLMM with the
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
THREADS <- as.integer(getarg("--threads", "8")); MAXITER <- as.integer(getarg("--maxiter", "500"))
MAXITER_GMMAT <- as.integer(getarg("--maxiter-gmmat", "500"))
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
  # Failed mixed models are not replaced by an unverified no-kinship model.
  if (!is.null(obj) && (!isTRUE(obj$converged) ||
      any(!is.finite(obj$coefficients)) || any(!is.finite(obj$theta)))) obj <- NULL
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
  # Fixed recovery schedule. Every attempt fits the original sparse-GRM model.
  # Increasing iterations and EM warm-up improves numerical stability without
  # changing participants, covariates, tolerance, or the variance-component model.
  schedules <- list(list(name = "AI500", start = NULL, EM.iter = 0L),
                    list(name = "EM25_AI500_start0.01", start = 0.01, EM.iter = 25L),
                    list(name = "EM50_AI500_start1", start = 1, EM.iter = 50L))
  attempts <- list()
  nullmod <- NULL
  chosen <- NULL
  for (schedule in schedules) {
    fit_warnings <- character()
    fit_error <- NULL
    fit_t0 <- proc.time()[3]
    candidate <- tryCatch(withCallingHandlers(
      fitNullModel(sad, outcome = "Y", covars = cn, cov.mat = K,
                   sample.id = ph$IID, family = "binomial", max.iter = MAXITER,
                   AIREML.tol = 1e-4, start = schedule$start,
                   EM.iter = schedule$EM.iter, drop.zeros = TRUE, verbose = FALSE),
      warning = function(w) { fit_warnings <<- c(fit_warnings, conditionMessage(w)); invokeRestart("muffleWarning") }),
      error = function(e) { fit_error <<- conditionMessage(e); NULL })
    valid <- !is.null(candidate) && isTRUE(candidate$converged)
    boundary <- !is.null(candidate) && isTRUE(candidate$zeroFLAG) && is.null(candidate$varComp)
    boundary_score <- boundary_tolerance <- NA_real_
    if (valid) {
      beta <- as.numeric(candidate$fixef[[1]])
      # GENESIS mixed-model fitted.values are on the working-response scale.
      # Its conditional linear.predictor must be inverse-linked for this check.
      mu <- if (isTRUE(candidate$model$family$mixedmodel))
        as.numeric(plogis(candidate$fit$linear.predictor)) else
        as.numeric(candidate$fit$fitted.values)
      vc <- as.numeric(candidate$varComp)
      valid <- length(mu) == nrow(ph) && all(is.finite(beta)) &&
        all(is.finite(mu)) && all(mu > 0 & mu < 1) &&
        all(is.finite(vc)) && all(vc >= 0) &&
        (length(vc) > 0 || boundary)
      if (valid && boundary) {
        # Independent local score check at the package-estimated zero boundary.
        # This is diagnostic support for a constrained REML/PQL boundary fit,
        # not a manual GLM substitution or a global likelihood proof.
        X <- as.matrix(candidate$model.matrix)
        w <- mu * (1 - mu)
        wx <- X * w
        resid <- as.numeric(candidate$fit$outcome) - mu
        term1 <- as.numeric(crossprod(resid, K %*% resid))
        term2 <- sum(diag(K) * w) - sum(diag(solve(crossprod(X, wx), crossprod(wx, K %*% wx))))
        boundary_score <- 0.5 * (term1 - term2)
        boundary_tolerance <- 1e-6 * max(1, abs(term1), abs(term2))
        valid <- is.finite(boundary_score) && boundary_score <= boundary_tolerance
      }
    }
    attempt <- list(schedule = schedule$name, maxiter = MAXITER,
                    EM_iter = schedule$EM.iter, start = schedule$start,
                    tolerance = 1e-4, elapsed_s = proc.time()[3] - fit_t0,
                    converged = !is.null(candidate) && isTRUE(candidate$converged),
                    valid = valid, in_model_zero_boundary = boundary,
                    boundary_score = boundary_score, boundary_tolerance = boundary_tolerance,
                    varComp = if (is.null(candidate)) NULL else unname(as.numeric(candidate$varComp)),
                    warnings = fit_warnings, error = fit_error)
    attempts[[length(attempts) + 1L]] <- attempt
    msg("GENESIS attempt %s: valid %s, converged %s, boundary %s (%.1f s)",
        schedule$name, valid, attempt$converged, boundary, attempt$elapsed_s)
    if (valid) { nullmod <- candidate; chosen <- attempt; break }
  }
  meta$genesis <- list(converged = !is.null(nullmod), accepted = !is.null(nullmod),
                       fallback_glm = FALSE, attempts = attempts,
                       null_fit_s = proc.time()[3] - t0,
                       package_version = as.character(packageVersion("GENESIS")))
  if (!is.null(nullmod)) {
    meta$genesis$varComp <- unname(as.numeric(nullmod$varComp))
    meta$genesis$in_model_zero_boundary <- chosen$in_model_zero_boundary
    meta$genesis$selected_schedule <- chosen$schedule
    meta$genesis$max_abs_coef <- max(abs(as.numeric(nullmod$fixef[[1]])))
    # Preserve the accepted model for independent inspection before aggregation.
    if (identical(getarg("--save-null", "FALSE"), "TRUE")) saveRDS(nullmod, paste0(OUT, "_genesis_null.rds"))
    t1 <- proc.time()[3]
    it <- GenotypeBlockIterator(gdata, snpBlock = 2000)      # 2026-09-23: 10,000-SNP blocks x 8 forked workers exceeded 24 GB
    nw <- min(THREADS, 4L)
    assoc <- tryCatch(assocTestSingle(it, null.model = nullmod, test = "Score.SPA", recalc.pval.thresh = 0.05,
                                      BPPARAM = if (nw > 1) MulticoreParam(nw) else SerialParam(), verbose = FALSE),
                      error = function(e) { msg("GENESIS assocTestSingle ERROR: %s", conditionMessage(e)); NULL })
    if (!is.null(assoc)) {
      i <- match(as.character(assoc$variant.id), bim$V2); k <- !is.na(i)      # variant.id = the .bim SNP names (snpgdsBED2GDS keeps them)
      # Package default applies SPA only to normal-score p <= 0.05.
      # Above that threshold, SPA.converged is NA and SPA.pval retains the score p.
      score_p <- pchisq(as.numeric(assoc$Score.Stat)^2, 1, lower.tail = FALSE)
      spa_required <- is.finite(score_p) & score_p <= 0.05
      spa_converged <- !is.na(assoc$SPA.converged) & assoc$SPA.converged
      score_retained <- is.na(assoc$SPA.converged) & is.finite(score_p) & score_p > 0.05
      spa_valid <- k & (spa_converged | score_retained) &
        is.finite(assoc$SPA.pval) & assoc$SPA.pval >= 0 & assoc$SPA.pval <= 1
      p_gspa[i[spa_valid]] <- as.numeric(assoc$SPA.pval[spa_valid])
      meta$genesis$spa_recalculation_threshold <- 0.05
      meta$genesis$n_spa_nonconverged <- sum(k & spa_required & !spa_converged)
      meta$genesis$n_score_retained_above_threshold <- sum(k & score_retained)
      p_gsc[i[k]] <- if ("Score.pval" %in% names(assoc)) as.numeric(assoc$Score.pval[k]) else pchisq(as.numeric(assoc$Score.Stat[k])^2, 1, lower.tail = FALSE)
      meta$genesis$score_s <- proc.time()[3] - t1; meta$genesis$n_tested <- sum(is.finite(p_gspa))
      meta$genesis$lambda_gc_spa <- lam(p_gspa); meta$genesis$lambda_gc_score <- lam(p_gsc)
      meta$genesis$spa_converged_share <- mean(assoc$SPA.converged[k], na.rm = TRUE)
      msg("GENESIS Score.SPA: %d SNPs tested, lambda SPA %.3f / score %.3f (%.0f s)", sum(is.finite(p_gspa)), lam(p_gspa), lam(p_gsc), proc.time()[3] - t1)
    }
  } else meta$genesis$error <- "No converged same-model fit passed the acceptance checks"
  close(gr)
}
fwrite(data.table(snp = bim$V2, p_gmmat = p_gmmat, p_genesis_spa = p_gspa, p_genesis_score = p_gsc), paste0(OUT, "_pvals.tsv.gz"), sep = "\t")
meta$repair_policy <- "same sparse GRM, no error-triggered GLM fallback, explicit convergence gate"
meta$seconds <- proc.time()[3] - t_all
writeLines(toJSON(meta, auto_unbox = TRUE, pretty = TRUE, digits = NA), paste0(OUT, "_meta.json"))
msg("DONE in %.0f s -> %s_pvals.tsv.gz", meta$seconds, OUT)
