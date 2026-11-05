#!/usr/bin/env Rscript
## split_half_null_models.R (2026-10-08): the D_S null model refitted inside one held-out half (HALF=A or B;
## subjects from split_sample_halves.R). Same specification as D_S: AD ~ Age + Sex + cell + latent class + posterior +
## entropy + PC1-10, sparse GRM restricted to the half, SPA null. Classes with < 30 subjects are folded into the largest
## class and aliased class dummies are folded too; if the full set still fails, the fit falls back to COVSET=nocls
## (Age + Sex + cell + PC1-10, as in the per-ancestry refits). Derived from ancestry_group_null_models.R.
suppressPackageStartupMessages({ library(data.table); library(STAAR); library(GMMAT); library(STAARpipeline) })
H <- Sys.getenv("HALF"); stopifnot(H %in% c("A", "B")); OUT <- "${ADSP_ROOT}/output/matched"
m <- fread(file.path(OUT, sprintf("matched_samples_projcell_split%s.csv", H)), colClasses = list(character = c("IID", "lpa_cluster")))
cat("half", H, ": N", nrow(m), "cases", sum(m$AD == 1), "cells", uniqueN(m$cell_SC), "\n")
m[, Sex := factor(Sex)]; m[, cell_SC := factor(cell_SC)]
tab <- sort(table(m$lpa_cluster), decreasing = TRUE); big <- names(tab)[1]; m[lpa_cluster %in% names(tab)[tab < 30], lpa_cluster := big]
PCS <- paste0("PC", 1:10)
form_of <- function(cs) if (cs == "nocls") as.formula(paste("AD ~ Age + Sex + cell_SC +", paste(PCS, collapse = " + "))) else
  as.formula(paste("AD ~ Age + Sex + cell_SC + lpa_cluster + cluster_prob + entropy_norm +", paste(PCS, collapse = " + ")))
sgrm <- get(load("${ADSP_ROOT}/processed/output.sparseGRM.sGRM.RData"))
gid <- rownames(sgrm); if (is.null(gid)) gid <- colnames(sgrm); gid_s <- sub("^[^_]*_", "", gid)
idx <- match(m$IID, gid_s); stopifnot(!any(is.na(idx))); K <- sgrm[idx, idx]; dimnames(K) <- list(m$IID, m$IID); rm(sgrm)
OPT <- Sys.getenv("OPTIM", "AI")
fit_one <- function(COVSET) {
  FORM <- form_of(COVSET); cat("covariate set:", COVSET, " optimizer:", OPT, "\n")
  repeat { m[, lpa_cluster := factor(lpa_cluster)]
    X <- model.matrix(FORM, data = as.data.frame(m)); q <- qr(X); if (q$rank == ncol(X)) break
    al <- colnames(X)[q$pivot[(q$rank + 1):ncol(X)]]; cat("aliased:", paste(al, collapse = " | "), "\n")
    lv <- sub("^lpa_cluster", "", al[grepl("^lpa_cluster", al)]); if (!length(lv)) stop("aliased non-class column")
    m[, lpa_cluster := as.character(lpa_cluster)]; m[lpa_cluster %in% lv, lpa_cluster := big] }
  cat("design columns:", ncol(X), " classes kept:", nlevels(m$lpa_cluster), "\n")
  df <- as.data.frame(m); rownames(df) <- df$IID; t0 <- proc.time()[3]
  obj <- fit_nullmodel(FORM, data = df, kins = K, id = "IID", use_sparse = TRUE, family = binomial(link = "logit"), verbose = TRUE, use_SPA = TRUE, method.optim = OPT)
  cf <- obj$coefficients; nonpc <- cf[!grepl("^PC", names(cf))]
  cat(sprintf("fit %.0f s; N %d; converged %s; theta %s; max |non-PC coef| %.2f\n", proc.time()[3] - t0, length(obj$id_include), obj$converged, paste(round(obj$theta, 3), collapse = ","), max(abs(nonpc))))
  if (max(abs(nonpc)) > 15) stop("separation on a factor column"); if (!isTRUE(obj$converged)) stop("not converged"); obj }
COVSET <- Sys.getenv("COVSET", "full")
obj_null <- tryCatch(fit_one(COVSET), error = function(e) { cat("FULL SET FAILED:", conditionMessage(e), "\n"); if (COVSET == "nocls") stop(e); fit_one("nocls") })
save(obj_null, file = file.path(OUT, sprintf("obj_null_projcell_split%s.RData", H))); cat("DONE\n")
