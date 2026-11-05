## D_S (2026-09-21): the ALL-SUBJECT STRATIFIED design = plasmode arm S1 (strat_project) on the real labels: every
## subject of a project x centre cell that holds both cases and controls, NO matching, the cell as a factor in the
## null (nothing collapsed) plus the pair33 covariates (SeqCenter dropped: it is nested in the cell).  In the power
## run this design kept 95% of the oracle's information (1:1 project key 65%) with a calibrated null.
## Built from matched_design_null_model.R; the original header follows.
## D_P (2026-09-21): the GCAD SEQUENCING PROJECT (QC-file column `Study`, 65 levels in the phenotyped
## pool; nests in SeqCenter except ADSP-Extension-CC and Amyloid-Pitt) added to the exact matching KEY,
## i.e. the plasmode arm M3 (riga_projectkey) on the real AD labels.  Everything else is the published
## pair33 recipe (src/prepare_sample_and_one_to_one_matching.R prep, 1:1 nearest Mahalanobis on PC1-10 inside
## exact blocks, no replacement, no confidence filter, the 37-column null, sparse GRM, use_SPA=TRUE),
## exactly as build_null_pair33_cohortkey.R (D3) did for Cohort.  Blocks: Sex + age_band + SeqCenter +
## lpa_cluster + Study.  Null formula UNCHANGED (for a 1:1 exact-matched set every key level holds equal
## cases and controls, so a key covariate adds nothing: harness D3n == D3); the QR guard that folds an
## aliased SeqCenter level into "Other" is kept.
## Motivation: on the real 281k genotypes with Y* permuted inside Cohort x project x centre cells the
## project-keyed 1:1 design is calibrated (FPR@1e-3 1.01x, 1e-5 0.9x) where the centre-keyed pair33 is
## not (6.2x / 195x); the power run shows the cost.  This arm asks what the real labels give.
suppressPackageStartupMessages({ library(data.table); library(MatchIt); library(STAAR); library(GMMAT); library(STAARpipeline) })
SRC <- "${ADSP_ROOT}/src"; OUT <- "${ADSP_ROOT}/output/matched"
TAB <- "${ADSP_ROOT}/revision/tables"
src <- readLines(file.path(SRC, "prepare_sample_and_one_to_one_matching.R")); cut <- grep("^m\\.out <- matchit", src)[1]
eval(parse(text = paste(src[1:(cut-1)], collapse = "\n")))          # builds dfx exactly as pair33
dfx <- as.data.table(dfx)
if (!"Cohort" %in% names(dfx)) { ph <- fread("${ADSP_ROOT}/pheno/adsp_combined.csv", select = c("IID","Cohort"))
  dfx <- merge(dfx, ph, by = "IID", all.x = TRUE) }
stopifnot(!any(is.na(dfx$Cohort)))
## the GCAD project label (QC file `Study`); a subject without a QC row gets its own "Unknown" level
qc <- fread("${ADSP_ROOT}/pheno/gcad.r5.wgs.58507.2024.11.03.quality.metrics.ALL.csv",
            select = c("Sample_name", "Study"), colClasses = "character"); setnames(qc, "Sample_name", "IID")
qc <- unique(qc, by = "IID"); dfx <- merge(dfx, qc, by = "IID", all.x = TRUE)
dfx[is.na(Study) | Study == "", Study := "Unknown"]; dfx[, Study := factor(Study)]
cat("pre-match pool:", nrow(dfx), " projects:", uniqueN(dfx$Study), " (Unknown:", sum(dfx$Study == "Unknown"), ")\n")
pt <- dfx[, .(N = .N, cases = sum(AD == 1)), by = Study][order(-N)]; pt[, mono := cases == 0 | cases == N]
cat("projects with cases only or controls only:", sum(pt$mono), "holding", sum(pt[mono == TRUE]$N), "subjects\n")
PCS <- paste0("PC", 1:10); set.seed(0)
dfx[, cell_SC := paste(as.character(Study), as.character(SeqCenter), sep = "|")]
ct <- dfx[, .(N = .N, cases = sum(AD == 1)), by = cell_SC]; ct[, mixed := cases > 0 & cases < N]
cat("project x centre cells:", nrow(ct), " mixed:", sum(ct$mixed), " subjects in mixed cells:", sum(ct[mixed == TRUE]$N), "\n")
matched <- dfx[cell_SC %in% ct[mixed == TRUE]$cell_SC]
matched[, subclass := cell_SC]
cat("D_S rows:", nrow(matched), "subjects (cases", sum(matched$AD == 1), "/ controls", sum(matched$AD == 0), "),",
    uniqueN(matched$cell_SC), "cells,", uniqueN(matched$Study), "projects,", uniqueN(matched$Cohort), "cohorts\n")
csv <- file.path(OUT, "matched_samples_lpa_pc.pair33_projcell.csv"); stopifnot(!file.exists(csv)); fwrite(matched, csv)
for (kk in c("Study", "Cohort")) {
  ret <- merge(dfx[, .(N_pool = .N, cases_pool = sum(AD == 1)), by = kk], matched[, .(N_matched = .N, cases = sum(AD==1)), by = kk], all.x = TRUE)
  ret[is.na(N_matched), `:=`(N_matched = 0L, cases = 0L)]; ret[, retention := round(N_matched / N_pool, 3)]
  fwrite(ret[order(-N_pool)], file.path(TAB, sprintf("pair33_projcell_retention_by_%s.csv", ifelse(kk == "Study", "project", "cohort")))) }

matched[, Sex := factor(Sex)]; matched[, SeqCenter := factor(SeqCenter)]; matched[, lpa_cluster := factor(lpa_cluster)]
matched[, cell_SC := factor(cell_SC)]
FORM <- as.formula(paste("AD ~ Age + Sex + cell_SC + lpa_cluster + cluster_prob + entropy_norm +", paste(PCS, collapse = " + ")))
repeat {                                                              # aliasing guard
  X <- model.matrix(FORM, data = as.data.frame(matched)); q <- qr(X); rk <- q$rank
  if (rk == ncol(X)) break
  al <- colnames(X)[q$pivot[(rk+1):ncol(X)]]; cat("rank", rk, "of", ncol(X), "aliased:", paste(al, collapse=" | "), "\n")
  stop("aliased column(s) in the cell-factor design: ", paste(al, collapse = ", ")) }
cat("design columns:", ncol(X), "\n")
sgrm <- get(load("${ADSP_ROOT}/processed/output.sparseGRM.sGRM.RData"))
gid <- rownames(sgrm); if (is.null(gid)) gid <- colnames(sgrm); gid_s <- sub("^[^_]*_", "", gid)
idx <- match(matched$IID, gid_s); stopifnot(!any(is.na(idx))); K <- sgrm[idx, idx]; dimnames(K) <- list(matched$IID, matched$IID)
df <- as.data.frame(matched); rownames(df) <- df$IID; t0 <- proc.time()[3]
obj_null <- fit_nullmodel(FORM, data = df, kins = K, id = "IID", use_sparse = TRUE, family = binomial(link = "logit"), verbose = TRUE, use_SPA = TRUE)
cat(sprintf("fit time %.1f s\n", proc.time()[3] - t0))
cf <- obj_null$coefficients; nonpc <- cf[!grepl("^PC", names(cf))]
cat("N:", length(obj_null$id_include), " converged:", obj_null$converged, " theta:", round(obj_null$theta, 4), " max |non-PC coef|:", round(max(abs(nonpc)), 3), "\n")
if (max(abs(nonpc)) > 15) stop("separation on a factor column; not saving")
nf <- file.path(OUT, "obj_null_stratum_all_matched_lpa_pc.pair33_projcell.RData"); stopifnot(!file.exists(nf)); save(obj_null, file = nf)
cat("wrote", nf, "\nDONE\n")
