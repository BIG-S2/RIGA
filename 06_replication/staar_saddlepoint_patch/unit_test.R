#!/usr/bin/env Rscript
## unit_test.R (2026-10-19, rev 2): item 2(d). Synthetic cases that force the UNPATCHED Newton solver to return an
## unconverged iterate with a finite p value, and the demonstration that the patched solver returns NA and that the
## callers then return the bisection value (which agrees with an independent R root finder); plus (rev 2, after review)
## the same forced failure through ALL FOUR C++ callers (STAAR_B_Binary_SPA, STAAR_B_SPA, STAAR_B_SPA_SMMAT,
## STAAR_B_SPA_SMMAT_sparse), a no-root construction that reaches every caller and must yield NA in the component column,
## and -- when a null-model file is given -- the scan's own path: STAAR_Binary_SPA() with the REAL reconstructed sparse-kins
## null (relatedness = TRUE, sparse_kins = TRUE, SPA_p_filter = TRUE), forced Newton failure and forced NA.
##
## Usage:  Rscript unit_test.R original <outdir> [nullmodel.RData]   # DEFAULT library (records the defect)
##         Rscript unit_test.R patched  <outdir> [nullmodel.RData]   # root-checked library (the test proper)
## The patched run reads the original run's CSV when it exists and prints both side by side. Exit status 0 = PASS.
## Checks in the ORIGINAL mode document the defect; the real-null checks of that mode are informative (recorded, pass = TRUE,
## because the original library's partial safeguard may or may not catch a given construction).
suppressPackageStartupMessages({ library(STAAR); library(STAARpipeline); library(Matrix) })
args <- commandArgs(TRUE); mode <- args[1]; outdir <- args[2]; nullfile <- if (length(args) >= 3 && nzchar(args[3])) args[3] else ""
stopifnot(mode %in% c("original", "patched"), !is.na(outdir)); dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
desc <- packageDescription("STAAR"); is_patched <- !is.null(desc$RootcheckPatch)
cat("STAAR loaded from", system.file(package = "STAAR"), " patched:", is_patched, "\n")
if (mode == "patched" && !is_patched) stop("mode=patched but the loaded STAAR has no RootcheckPatch field (R_LIBS_USER?)")
if (mode == "original" && is_patched) stop("mode=original but the loaded STAAR is the patched one")
tol <- .Machine$double.eps^0.25; mi <- 1000L
logf <- file.path(outdir, "unit_test_rootlog.tsv")
if (mode == "patched") { if (file.exists(logf)) file.remove(logf); Sys.setenv(STAAR_ROOTCHECK_LOG = logf) }
read_log <- function() if (file.exists(logf) && file.size(logf) > 0) read.delim(logf, header = FALSE, quote = "", na.strings = c("NA", "nan", "-nan"), stringsAsFactors = FALSE) else data.frame()
n_log <- function() nrow(read_log())
LOGCOL <- list(call_seq = 2, n_variants = 3, cMAC = 4, n_weights = 5, component = 6, p_normal = 7, tail = 8, event = 9, context = 16, n_spa = 17)

## independent reference: bracketed root + Lugannani-Rice tail, same centred-Bernoulli CGF target as STAAR (cf. guarded_spa.R)
ref_tail <- function(q, g, mu, lower) {
  eta <- qlogis(mu); k1 <- function(t) sum(g * (plogis(eta + t * g) - mu)) - q
  lo <- -1; hi <- 1; for (k in 1:200) { if (k1(lo) <= 0 && k1(hi) >= 0) break; if (k1(lo) > 0) lo <- lo * 2; if (k1(hi) < 0) hi <- hi * 2 }
  z <- uniroot(k1, c(lo, hi), tol = 1e-13, maxiter = 5000)$root
  sp <- function(x) pmax(x, 0) + log1p(exp(-abs(x)))
  kap <- sum(sp(eta + z * g) - sp(eta) - z * g * mu); k2 <- sum(g^2 * plogis(eta + z * g) * (1 - plogis(eta + z * g)))
  w <- sign(z) * sqrt(2 * (z * q - kap)); v <- z * sqrt(k2)
  list(root = z, resid = k1(z), p = pnorm(w + log(v / w) / w, lower.tail = lower)) }
ref_two_sided <- function(q, g, mu) ref_tail(abs(q), g, mu, FALSE)$p + ref_tail(-abs(q), g, mu, TRUE)$p

rows <- list(); checks <- list()
chk <- function(name, ok, detail = "") { checks[[length(checks) + 1]] <<- data.frame(check = name, pass = isTRUE(ok), detail = detail); cat(sprintf("  [%s] %s %s\n", if (isTRUE(ok)) "PASS" else "FAIL", name, detail)) }
info <- function(name, detail) { checks[[length(checks) + 1]] <<- data.frame(check = paste0(name, " (informative)"), pass = TRUE, detail = detail); cat(sprintf("  [INFO] %s %s\n", name, detail)) }
relerr <- function(a, b) abs(a - b) / abs(b)

## ---- case A: solver level. One dominant weight (15.968) and 900 small ones (1.006), heterogeneous mu, q = 11.559 ------------
set.seed(153); muA <- runif(901, 0.003, 0.03); GA <- c(15.968, rep(1.006, 900)); qA <- 11.559
xN <- STAAR:::NR_Binary_SPA(muA, GA, qA, 0, tol, mi); rN <- if (is.na(xN)) NA else STAAR:::K1_Binary_SPA(xN, muA, GA, qA)
pN <- STAAR:::Saddle_Binary_SPA(qA, muA, GA, tol, mi, FALSE)
xB <- STAAR:::Bisection_Binary_SPA(muA, GA, qA, -100, 100, tol); rB <- if (is.na(xB)) NA else STAAR:::K1_Binary_SPA(xB, muA, GA, qA)
pB <- STAAR:::Saddle_Binary_SPA_Bisection(qA, muA, GA, tol, mi, -100, 100, FALSE)
refA <- ref_tail(qA, GA, muA, FALSE)
rows[["A"]] <- data.frame(case = "A_solver", mode = mode, newton_root = xN, newton_resid = rN, newton_p = pN, bisect_root = xB, bisect_resid = rB, bisect_p = pB, ref_root = refA$root, ref_p = refA$p)
cat(sprintf("case A: Newton root %s resid %s p %s | bisection root %s resid %s p %s | reference root %.8g p %.8g\n",
            format(xN), format(rN), format(pN), format(xB), format(rB), format(pB), refA$root, refA$p))
if (mode == "patched") {
  chk("A: patched NR_Binary_SPA returns NA on the non-converging case", is.na(xN))
  chk("A: patched Saddle_Binary_SPA returns NA (callers then bisect)", is.na(pN))
  chk("A: patched bisection root satisfies |K1-q| <= 1e-6*max(1,|q|)", is.finite(xB) && abs(rB) <= 1e-6 * max(1, abs(qA)), sprintf("resid %.3g", rB))
  chk("A: patched bisection p equals the independent reference to 1e-6 relative", relerr(pB, refA$p) <= 1e-6, sprintf("%.10g vs %.10g", pB, refA$p))
} else {
  chk("A: ORIGINAL Newton returns a finite root with a large residual (the defect)", is.finite(xN) && abs(rN) > 1, sprintf("root %.4g resid %.4g", xN, rN))
  chk("A: ORIGINAL Saddle_Binary_SPA passes a finite, tiny p through", is.finite(pN) && pN < 1e-10, sprintf("p %.3g vs reference %.4g", pN, refA$p))
}

## ---- case B: caller level (STAAR_B_Binary_SPA). Fitted logistic null (intercept + covariate), 5,000 subjects, 901
## single-carrier variants; variant 1 (weight 22.9, as a MAC-2 LINSIGHT-weighted variant) has one case and one control
## carrier, the 900 others (weight 1.006) one carrier each, five of them cases. This reproduces the SNX22 pathology:
## the unpatched caller returns p ~ 2e-5 for a set whose correct two-sided saddlepoint p is ~ 0.028. ------------------
n <- 5000L; m <- 901L
Gm <- matrix(0, n, m); Gm[1, 1] <- 1; Gm[2, 1] <- 1; for (j in 2:m) Gm[j + 1, j] <- 1
set.seed(7); z <- rnorm(n, 0, 0.8); y <- rbinom(n, 1, plogis(-4.7 + z)); y[1] <- 1L; y[2] <- 0L; y[3:7] <- 1L; y[8:902] <- 0L
fit <- glm(y ~ z, family = binomial); mu <- as.numeric(fitted(fit)); X <- model.matrix(fit); W <- mu * (1 - mu)
XW <- t(X) * rep(W, each = ncol(X)); XXWX_inv <- X %*% solve(t(X) %*% (X * W)); resid <- y - mu     # as in fit_null_glm_Binary_SPA
wB <- cbind(c(22.9, rep(1.006, m - 1)))
Gc <- as.vector((Gm - XXWX_inv %*% (XW %*% Gm)) %*% wB); qB <- sum(resid * Gc)
stopifnot(abs(qB - sum(resid * (Gm %*% wB))) < 1e-8)   # fitted null: centred and uncentred scores coincide (the caller uses the latter)
xN2 <- STAAR:::NR_Binary_SPA(mu, Gc, qB, 0, tol, mi); rN2 <- if (is.na(xN2)) NA else STAAR:::K1_Binary_SPA(xN2, mu, Gc, qB)
p_caller <- as.numeric(STAAR:::STAAR_B_Binary_SPA(Gm, XW, XXWX_inv, resid, mu, wB, tol, mi))
p_bis_up <- STAAR:::Saddle_Binary_SPA_Bisection(abs(qB), mu, Gc, tol, mi, -100, 100, FALSE)
p_bis_lo <- STAAR:::Saddle_Binary_SPA_Bisection(-abs(qB), mu, Gc, tol, mi, -100, 100, TRUE)
refB <- ref_two_sided(qB, Gc, mu)
rows[["B"]] <- data.frame(case = "B_caller", mode = mode, newton_root = xN2, newton_resid = rN2, newton_p = NA, bisect_root = NA, bisect_resid = NA, bisect_p = p_bis_up + p_bis_lo, ref_root = NA, ref_p = refB, caller_p = p_caller)
cat(sprintf("case B: q %.5f | Newton root %s resid %s | caller p %s | bisection two-sided %s | reference %.8g\n", qB, format(xN2), format(rN2), format(p_caller), format(p_bis_up + p_bis_lo), refB))
## the same null expressed as GMMAT/glmmkin quantities (kinship = 0, so Sigma^-1 = diag(W)): inputs of the three other callers
G_sp <- as(Gm, "dgCMatrix"); Sigma_i <- sparseMatrix(i = 1:n, j = 1:n, x = W); Sigma_iX <- X * W; covm <- solve(t(X) %*% (X * W))
P <- diag(W) - Sigma_iX %*% covm %*% t(Sigma_iX)
call4 <- function(residuals) list(
  Binary_SPA   = as.numeric(STAAR:::STAAR_B_Binary_SPA(Gm, XW, XXWX_inv, residuals, mu, wB, tol, mi)),
  SPA          = as.numeric(STAAR:::STAAR_B_SPA(G = Gm, XW = XW, XXWX_inv = XXWX_inv, residuals = residuals, muhat = mu, weights_B = wB, tol = tol, max_iter = mi, p_filter_cutoff = 0.05, G_sp = G_sp, X = X, working = fit$weights, sigma = 1, fam = 1)),
  SPA_SMMAT    = as.numeric(STAAR:::STAAR_B_SPA_SMMAT(G = Gm, XW = XW, XXWX_inv = XXWX_inv, residuals = residuals, muhat = mu, weights_B = wB, tol = tol, max_iter = mi, p_filter_cutoff = 0.05, G_sp = G_sp, P = P)),
  SPA_SMMAT_sparse = as.numeric(STAAR:::STAAR_B_SPA_SMMAT_sparse(G = Gm, XW = XW, XXWX_inv = XXWX_inv, residuals = residuals, muhat = mu, weights_B = wB, tol = tol, max_iter = mi, p_filter_cutoff = 0.05, G_sp = G_sp, Sigma_i = Sigma_i, Sigma_iX = Sigma_iX, cov = covm)))
p4 <- call4(resid)
cat("case B through the four callers:", paste(names(p4), format(unlist(p4), digits = 10), collapse = ", "), "\n")
rows[["B4"]] <- data.frame(case = "B_four_callers", mode = mode, newton_root = NA, newton_resid = NA, newton_p = NA, bisect_root = NA, bisect_resid = NA, bisect_p = NA, ref_root = NA, ref_p = refB, caller_p = NA,
                           p_SPA = p4$SPA, p_SPA_SMMAT = p4$SPA_SMMAT, p_SPA_SMMAT_sparse = p4$SPA_SMMAT_sparse)
if (mode == "patched") {
  chk("B: patched Newton returns NA for the upper-tail root", is.na(xN2))
  chk("B: patched STAAR_B_Binary_SPA equals the bisection two-sided value (1e-8 relative; the caller forms G_cumu in C++)", relerr(p_caller, p_bis_up + p_bis_lo) <= 1e-8, sprintf("%.15g vs %.15g", p_caller, p_bis_up + p_bis_lo))
  chk("B: patched caller p equals the independent reference to 1e-6 relative", relerr(p_caller, refB) <= 1e-6, sprintf("%.10g vs %.10g", p_caller, refB))
  chk("B4: STAAR_B_SPA (GLM, filtered) returns the same bisection value (1e-8 relative)", relerr(p4$SPA, p_caller) <= 1e-8, sprintf("%.12g", p4$SPA))
  chk("B4: STAAR_B_SPA_SMMAT (dense P, filtered) returns the same bisection value (1e-8 relative)", relerr(p4$SPA_SMMAT, p_caller) <= 1e-8, sprintf("%.12g", p4$SPA_SMMAT))
  chk("B4: STAAR_B_SPA_SMMAT_sparse (the scan's caller, filtered) returns the same bisection value (1e-8 relative)", relerr(p4$SPA_SMMAT_sparse, p_caller) <= 1e-8, sprintf("%.12g", p4$SPA_SMMAT_sparse))
  lg <- read_log()
  chk("B/A: diagnostic log records the Newton failures and the bisection fallbacks", nrow(lg) >= 4 && any(grepl("^newton_", lg[[LOGCOL$event]])) && any(grepl("^bisection_ok", lg[[LOGCOL$event]])), sprintf("%d lines: %s", nrow(lg), paste(names(table(lg[[LOGCOL$event]])), table(lg[[LOGCOL$event]]), collapse = ", ")))
  cd <- lg[lg[[LOGCOL$event]] == "call_done", ]
  ## five burden calls so far: the direct STAAR_B_Binary_SPA call above plus the four of call4()
  chk("B4: one call_done line per burden call (5: direct Binary_SPA + the four callers), each with n_weights 1 and n_spa 1 (the component entered SPA)", nrow(cd) == 5 && all(cd[[LOGCOL$n_weights]] == 1) && all(cd[[LOGCOL$n_spa]] == 1) && all(cd[[LOGCOL$n_variants]] == m), sprintf("%d call_done lines", nrow(cd)))
  chk("B4: filtered callers log p_normal < 0.05 for the SPA component; the unfiltered caller logs NA", { e <- lg[grepl("^newton_", lg[[LOGCOL$event]]) & lg[[LOGCOL$call_seq]] >= 1, ]; pn <- e[[LOGCOL$p_normal]]; sum(is.na(pn)) >= 1 && all(pn[!is.na(pn)] < 0.05) }, "")
} else {
  chk("B: ORIGINAL caller passes a p from an unverified root through (differs from the reference)", is.finite(p_caller) && relerr(p_caller, refB) > 1e-3, sprintf("caller %.5g vs reference %.5g", p_caller, refB))
  chk("B4: ORIGINAL STAAR_B_SPA / _SMMAT / _SMMAT_sparse pass the same wrong value through", all(sapply(p4[-1], function(p) is.finite(p) && relerr(p, refB) > 1e-3)), paste(format(unlist(p4[-1]), digits = 5), collapse = ", "))
}

## ---- case C: a target outside the support (no root exists) must surface as NA at the solver level, and (case G) at every
## caller: an artificial residual vector (3 for every subject) makes the score q exceed sup K1, so neither solver can find a
## root; the component must come back as NA (not 1, not a number) from all four callers. -----------------------------------
if (mode == "patched") {
  pC <- STAAR:::Saddle_Binary_SPA_Bisection(sum((1 - muA) * GA) * 1.5, muA, GA, tol, mi, -100, 100, FALSE)
  chk("C: patched bisection reports NA when no root can be bracketed", is.na(pC), format(pC))
}
res_fake <- rep(3, n); n0 <- if (mode == "patched") n_log() else 0
pG <- call4(res_fake)
cat("case G (no root reachable) through the four callers:", paste(names(pG), format(unlist(pG)), collapse = ", "), "\n")
rows[["G"]] <- data.frame(case = "G_no_root_callers", mode = mode, newton_root = NA, newton_resid = NA, newton_p = NA, bisect_root = NA, bisect_resid = NA, bisect_p = NA, ref_root = NA, ref_p = NA, caller_p = pG$Binary_SPA, p_SPA = pG$SPA, p_SPA_SMMAT = pG$SPA_SMMAT, p_SPA_SMMAT_sparse = pG$SPA_SMMAT_sparse)
if (mode == "patched") {
  chk("G: no-root construction: all four callers return NA in the component column (not 1, not a number)", all(is.na(unlist(pG))), paste(format(unlist(pG)), collapse = ","))
  lg <- read_log(); lgG <- lg[seq_len(nrow(lg)) > n0, ]
  chk("G: the log shows the Newton failure and the bisection failure (no bracket / non-finite) for each caller", sum(grepl("^newton_", lgG[[LOGCOL$event]])) >= 4 && sum(lgG[[LOGCOL$event]] %in% c("bisection_nobracket", "bisection_nonfinite", "bisection_unverified")) >= 4, paste(names(table(lgG[[LOGCOL$event]])), table(lgG[[LOGCOL$event]]), collapse = ", "))
} else {
  info("G: ORIGINAL callers on the no-root construction (NA would be correct)", paste(names(pG), format(unlist(pG)), collapse = ", "))
}

## ---- case D: an ordinary well-conditioned problem is untouched (Newton converges, no event logged) ---------------------
if (mode == "patched") {
  set.seed(5); GD <- runif(400, 0.1, 1); muD <- runif(400, 0.01, 0.05); qD <- 0.8
  n_before <- n_log(); xD <- STAAR:::NR_Binary_SPA(muD, GD, qD, 0, tol, mi); n_after <- n_log()
  chk("D: well-conditioned case: Newton root returned, residual within tolerance, no event logged", is.finite(xD) && abs(STAAR:::K1_Binary_SPA(xD, muD, GD, qD)) <= 1e-6 && n_after == n_before, sprintf("root %.6g", xD))
}

## ---- case E: the pipeline's private copy of the solver (STAARpipeline::Individual_Score_Test_SPA) on the case-B construction
p_pipe <- as.numeric(STAARpipeline:::Individual_Score_Test_SPA(Gm %*% wB, XW, XXWX_inv, resid, mu, tol, mi))
cat(sprintf("STAARpipeline Individual_Score_Test_SPA on the weighted burden as a single variant: p %s (reference %.8g)\n", format(p_pipe), refB))
if (mode == "patched") chk("E: patched STAARpipeline single-variant SPA equals the reference to 1e-6 relative", relerr(p_pipe, refB) <= 1e-6, sprintf("%.10g vs %.10g", p_pipe, refB)) else chk("E: ORIGINAL STAARpipeline single-variant SPA differs from the reference", relerr(p_pipe, refB) > 1e-3, sprintf("%.5g vs %.5g", p_pipe, refB))

## ---- case H: the scan's own path. STAAR_Binary_SPA() with the REAL reconstructed null (sparse kins, SPA_p_filter = TRUE,
## p_filter_cutoff = 0.05), 301 variants: variant 1 (annotation PHRED 40 -> weight ~25) with one case and one control carrier,
## 300 singletons (PHRED 0.2 -> weight ~1.1), three of them in cases. Found by scanning constructions on the null (19 Oct):
## both Newton tails hit max_iter for the annotation-weighted columns; reference two-sided p ~ 0.0046. ---------------------
rm(P, G_sp, Gm); invisible(gc())
if (nzchar(nullfile)) {
  cat("\ncase H: real null", nullfile, "\n")
  hres <- tryCatch({
    obj <- get(load(nullfile)); obj$sparse_kins <- TRUE; obj$relatedness <- TRUE; obj$use_SPA <- TRUE
    muH <- obj$fitted.values; rH <- obj$scaled.residuals; nH <- length(muH); yH <- round(muH + rH)
    stopifnot(nH == length(obj$id_include), sum(yH == 1) > 1000)
    XS <- as.matrix(obj$XSigma_i); XXI <- as.matrix(obj$XXSigma_iX_inv)
    cases <- which(yH == 1); ctrls <- which(yH == 0); mH <- 301L; ncH <- 3L; phH <- 0.2
    set.seed(1); c1 <- sample(cases, 1); k1 <- sample(ctrls, 1); sc <- sample(setdiff(cases, c1), ncH); sk <- sample(setdiff(ctrls, k1), mH - 1 - ncH)
    Gs <- sparseMatrix(i = c(c1, k1, sc, sk), j = c(1, 1, 2:mH), x = 1, dims = c(nH, mH)); GH <- as.matrix(Gs)
    phred <- c(40, rep(phH, mH - 1)); A <- data.frame(FAKE = phred)
    maf <- colSums(GH) / (2 * nH); rank <- 1 - 10^(-phred / 10)
    WH <- cbind("Burden(1,25)" = dbeta(maf, 1, 25), "Burden(1,25)-FAKE" = dbeta(maf, 1, 25) * rank, "Burden(1,1)" = dbeta(maf, 1, 1), "Burden(1,1)-FAKE" = dbeta(maf, 1, 1) * rank)
    refH <- t(sapply(colnames(WH), function(k) { gw <- as.numeric(Gs %*% WH[, k]); gt <- gw - as.numeric(XXI %*% (XS %*% gw)); q <- sum(rH * gw)
      c(q = q, q_centred = sum(rH * gt), ref_p = ref_two_sided(q, gt, muH), newton_root = STAAR:::NR_Binary_SPA(muH, gt, abs(q), 0, tol, mi)) }))
    n0 <- if (mode == "patched") n_log() else 0
    r1 <- STAAR_Binary_SPA(GH, obj, A, rare_maf_cutoff = 0.01, rv_num_cutoff = 2, SPA_p_filter = TRUE, p_filter_cutoff = 0.05)
    pH <- c(unlist(r1$results_STAAR_B_1_25)[c("Burden(1,25)", "Burden(1,25)-FAKE")], unlist(r1$results_STAAR_B_1_1)[c("Burden(1,1)", "Burden(1,1)-FAKE")])
    ## forced NA: a copy of the null whose scaled residuals are 3 for every subject (q >> sup K1 for every column)
    obj2 <- obj; obj2$scaled.residuals <- rep(3, nH)
    n1 <- if (mode == "patched") n_log() else 0
    r2 <- STAAR_Binary_SPA(GH, obj2, A, rare_maf_cutoff = 0.01, rv_num_cutoff = 2, SPA_p_filter = TRUE, p_filter_cutoff = 0.05)
    pH2 <- c(unlist(r2$results_STAAR_B_1_25)[c("Burden(1,25)", "Burden(1,25)-FAKE")], unlist(r2$results_STAAR_B_1_1)[c("Burden(1,1)", "Burden(1,1)-FAKE")])
    list(refH = refH, pH = pH, r1 = r1, pH2 = pH2, r2 = r2, n0 = n0, n1 = n1, cMAC = r1$cMAC, num_variant = r1$num_variant) }, error = function(e) e)
  if (inherits(hres, "error")) { chk("H: real-null case ran", FALSE, conditionMessage(hres)) } else with(hres, {
    cat(sprintf("case H: %d variants, cMAC %.0f\n", num_variant, cMAC)); print(cbind(refH, caller_p = pH, caller_p_forcedNA = pH2))
    rows[["H"]] <<- data.frame(case = "H_real_null_sparse_kins", mode = mode, newton_root = refH["Burden(1,25)-FAKE", "newton_root"], newton_resid = NA, newton_p = NA, bisect_root = NA, bisect_resid = NA, bisect_p = NA, ref_root = NA, ref_p = refH["Burden(1,25)-FAKE", "ref_p"], caller_p = pH[["Burden(1,25)-FAKE"]],
                               p_SPA = NA, p_SPA_SMMAT = NA, p_SPA_SMMAT_sparse = pH[["Burden(1,1)-FAKE"]])
    spa_cols <- names(pH)[!is.na(pH) & pH < 0.05]   # stored p < 0.05 => the component went through SPA (filtered caller)
    if (mode == "patched") {
      chk("H: real null, scan path: the annotation-weighted components entered SPA (p < 0.05)", all(c("Burden(1,25)-FAKE", "Burden(1,1)-FAKE") %in% spa_cols), paste(spa_cols, collapse = ","))
      chk("H: real null: patched Newton returns NA for the (1,25)-FAKE upper-tail root (forced failure)", is.na(refH["Burden(1,25)-FAKE", "newton_root"]))
      chk("H: real null: every SPA component equals the independent reference to 1e-6 relative", all(sapply(spa_cols, function(k) relerr(pH[[k]], refH[k, "ref_p"]) <= 1e-6)), paste(sapply(spa_cols, function(k) sprintf("%s %.8g vs %.8g", k, pH[[k]], refH[k, "ref_p"])), collapse = "; "))
      lg <- read_log(); lgH <- lg[seq_len(nrow(lg)) > n0 & seq_len(nrow(lg)) <= n1, ]
      chk("H: real null: the log records a Newton failure resolved by bisection for the forced component (component 1 = Burden(1,25)-FAKE)", any(grepl("^newton_(maxiter|nonfinite|unverified)", lgH[[LOGCOL$event]]) & lgH[[LOGCOL$component]] == 1) && any(grepl("^bisection_ok", lgH[[LOGCOL$event]]) & lgH[[LOGCOL$component]] == 1), paste(names(table(lgH[[LOGCOL$event]])), table(lgH[[LOGCOL$event]]), collapse = ", "))
      cdH <- lgH[lgH[[LOGCOL$event]] == "call_done", ]
      chk("H: real null: call_done line carries n_variants 301, cMAC 302, n_weights 4 and n_spa = number of SPA components", nrow(cdH) == 1 && cdH[[LOGCOL$n_variants]] == 301 && cdH[[LOGCOL$cMAC]] == 302 && cdH[[LOGCOL$n_weights]] == 4 && cdH[[LOGCOL$n_spa]] == length(spa_cols), sprintf("n_spa %s vs %d", if (nrow(cdH)) cdH[[LOGCOL$n_spa]] else NA, length(spa_cols)))
      chk("H: real null, forced no-root: every component is NA in the component column and STAAR-B falls to 1 (all NA)", all(is.na(pH2)) && isTRUE(r2$results_STAAR_B == 1), paste(format(pH2), collapse = ","))
      lgN <- lg[seq_len(nrow(lg)) > n1, ]
      chk("H: real null, forced no-root: the log shows bisection failures (no bracket / non-finite) for the components", sum(lgN[[LOGCOL$event]] %in% c("bisection_nobracket", "bisection_nonfinite", "bisection_unverified")) >= 4, paste(names(table(lgN[[LOGCOL$event]])), table(lgN[[LOGCOL$event]]), collapse = ", "))
    } else {
      info("H: ORIGINAL library on the real-null construction: caller vs reference for the SPA components", paste(sapply(spa_cols, function(k) sprintf("%s caller %.6g ref %.6g (rel %.2g)", k, pH[[k]], refH[k, "ref_p"], relerr(pH[[k]], refH[k, "ref_p"]))), collapse = "; "))
      info("H: ORIGINAL library on the forced no-root construction (NA would be correct)", paste(names(pH2), format(pH2), collapse = ", "))
    }
  })
} else cat("\ncase H skipped (no null-model file given)\n")

## ---- outputs ------------------------------------------------------------------------------------------------------------
cols <- c("case", "mode", "newton_root", "newton_resid", "newton_p", "bisect_root", "bisect_resid", "bisect_p", "ref_root", "ref_p", "caller_p", "p_SPA", "p_SPA_SMMAT", "p_SPA_SMMAT_sparse")
res <- do.call(rbind, lapply(rows, function(r) { for (cc in setdiff(cols, names(r))) r[[cc]] <- NA; r[cols] }))
write.csv(res, file.path(outdir, paste0("unit_test_values_", mode, ".csv")), row.names = FALSE)
ck <- do.call(rbind, checks); write.csv(ck, file.path(outdir, paste0("unit_test_checks_", mode, ".csv")), row.names = FALSE)
if (mode == "patched" && file.exists(file.path(outdir, "unit_test_values_original.csv"))) {
  o <- read.csv(file.path(outdir, "unit_test_values_original.csv")); cat("\nside by side (original library vs patched):\n")
  print(rbind(o, res)[, c("case", "mode", "newton_root", "newton_resid", "newton_p", "bisect_p", "caller_p", "p_SPA_SMMAT_sparse", "ref_p")])
}
cat(sprintf("\nUNIT TEST (%s): %d/%d checks passed\n", mode, sum(ck$pass), nrow(ck)))
if (!all(ck$pass)) quit(status = 1)
