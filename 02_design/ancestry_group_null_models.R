#!/usr/bin/env Rscript
## ancestry_group_null_models.R <race code> (2026-10-07): the stratified design D_S fitted INSIDE one self-reported
## ancestry group (project x center cells formed within the group, mixed cells only; same covariates as D_S:
## Age + Sex + cell + latent class + posterior + entropy + PC1-10; sparse GRM; SPA null). Derived from
## riga_study_units_null_model.R. Race codes: 5 White, 6 Hispanic/Latino, 4 Black, 2 Asian.
suppressPackageStartupMessages({ library(data.table); library(STAAR); library(GMMAT) })
r <- as.integer(commandArgs(TRUE)[1]); stopifnot(r %in% c(2L, 4L, 5L, 6L))
SRC <- "${ADSP_ROOT}/src"; OUT <- "${ADSP_ROOT}/output/matched"
src <- readLines(file.path(SRC, "prepare_sample_and_one_to_one_matching.R")); cut <- grep("^m\\.out <- matchit", src)[1]
eval(parse(text = paste(src[1:(cut - 1)], collapse = "\n"))); dfx <- as.data.table(dfx)
qc <- fread("${ADSP_ROOT}/pheno/gcad.r5.wgs.58507.2024.11.03.quality.metrics.ALL.csv", select = c("Sample_name", "Study"), colClasses = "character")
setnames(qc, "Sample_name", "IID"); qc <- unique(qc, by = "IID"); dfx <- merge(dfx, qc, by = "IID", all.x = TRUE); dfx[is.na(Study) | Study == "", Study := "Unknown"]
dfx <- dfx[Race == r]; cat("race", r, ": pool", nrow(dfx), "cases", sum(dfx$AD == 1), "\n")
dfx[, cell_SC := paste(as.character(Study), as.character(SeqCenter), sep = "|")]
ct <- dfx[, .(N = .N, cases = sum(AD == 1)), by = cell_SC]; ct[, mixed := cases > 0 & cases < N]
m <- dfx[cell_SC %in% ct[mixed == TRUE]$cell_SC]
cat("cells", nrow(ct), "mixed", sum(ct$mixed), "subjects kept", nrow(m), "(cases", sum(m$AD == 1), ")\n")
m[, Sex := factor(Sex)]; m[, cell_SC := factor(cell_SC)]; m[, lpa_cluster := as.character(lpa_cluster)]
tab <- sort(table(m$lpa_cluster), decreasing = TRUE); big <- names(tab)[1]
m[lpa_cluster %in% names(tab)[tab < 30], lpa_cluster := big]
PCS <- paste0("PC", 1:10)
## COVSET=nocls drops the latent-class terms: inside one ancestry group the class dummies are nearly collinear with the
## cells and PCs (dpotrf: leading principal minor not positive in the White, Black and Asian fits), so the per-ancestry
## null uses Age + Sex + cell + PC1-10.
COVSET <- Sys.getenv("COVSET", "full"); cat("covariate set:", COVSET, "\n")
FORM <- if (COVSET == "nocls") as.formula(paste("AD ~ Age + Sex + cell_SC +", paste(PCS, collapse = " + "))) else
        as.formula(paste("AD ~ Age + Sex + cell_SC + lpa_cluster + cluster_prob + entropy_norm +", paste(PCS, collapse = " + ")))
repeat {
  m[, lpa_cluster := factor(lpa_cluster)]
  X <- model.matrix(FORM, data = as.data.frame(m)); q <- qr(X); if (q$rank == ncol(X)) break
  al <- colnames(X)[q$pivot[(q$rank + 1):ncol(X)]]; cat("aliased:", paste(al, collapse = " | "), "\n")
  lv <- sub("^lpa_cluster", "", al[grepl("^lpa_cluster", al)])
  if (!length(lv)) stop("aliased non-class column; not saving")
  m[, lpa_cluster := as.character(lpa_cluster)]; m[lpa_cluster %in% lv, lpa_cluster := big] }
cat("design columns:", ncol(X), " classes kept:", nlevels(m$lpa_cluster), "\n")
sgrm <- get(load("${ADSP_ROOT}/processed/output.sparseGRM.sGRM.RData"))
gid <- rownames(sgrm); if (is.null(gid)) gid <- colnames(sgrm); gid_s <- sub("^[^_]*_", "", gid)
idx <- match(m$IID, gid_s); stopifnot(!any(is.na(idx))); K <- sgrm[idx, idx]; dimnames(K) <- list(m$IID, m$IID)
df <- as.data.frame(m); rownames(df) <- df$IID; t0 <- proc.time()[3]
OPT <- Sys.getenv("OPTIM", "AI"); cat("optimizer:", OPT, "\n")
obj_null <- fit_nullmodel(FORM, data = df, kins = K, id = "IID", use_sparse = TRUE, family = binomial(link = "logit"), verbose = TRUE, use_SPA = TRUE, method.optim = OPT)
cf <- obj_null$coefficients; nonpc <- cf[!grepl("^PC", names(cf))]
cat(sprintf("fit %.0f s; N %d; converged %s; theta %s; max |non-PC coef| %.2f\n", proc.time()[3] - t0, length(obj_null$id_include), obj_null$converged, paste(round(obj_null$theta, 3), collapse = ","), max(abs(nonpc))))
if (max(abs(nonpc)) > 15) stop("separation on a factor column; not saving")
fwrite(m, file.path(OUT, sprintf("matched_samples_projcell_race%d.csv", r)))
save(obj_null, file = file.path(OUT, sprintf("obj_null_projcell_race%d.RData", r))); cat("DONE\n")
