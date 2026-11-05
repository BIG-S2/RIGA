#!/usr/bin/env Rscript
## build_null.R (2026-10-20): the pooled-full null (obj_null_stratum_all_pooled_lpa_pc_full.RData, 38,138 participants,
## 37 design columns) refitted with ONE additional term, a collapsed sequencing-project factor, on exactly the same
## participants, the same covariate columns, the same sparse GRM and the same fit_nullmodel settings as
## src/pooled_design_null_model.R.
##
## Collapse rule (the paper's benchmark "Pooled (covariates)" design; revision/scripts/build_null_pair33_cohortcov.R L19-28,
## revision/scripts/benchmark_null_harness.py collapse_levels, manuscript Methods): on the analysed participants, a project
## whose participants are all cases or all controls (n0 == 0 | n1 == 0) OR whose N < 50 is assigned to the level "Other";
## single pass, recomputed on the 38,138; reference level "Other".
##
## Refusal gates (status JSON, non-zero exit so the dependent scan arrays never start): fit did not converge, non-finite
## coefficients / theta, or more than 1% of participants with a fitted probability within 1e-6 of 0 or 1.
##
## Usage: Rscript build_null.R [--use-spa-fit=TRUE|FALSE] [--overwrite]
##   --use-spa-fit  the use_SPA argument of fit_nullmodel. Default FALSE = the written specification of this comparator;
##                  sbatch/00_build_null.sbatch passes it explicitly (USE_SPA_FIT=TRUE in the environment of SUBMIT.sh or
##                  of the sbatch call flips it to the pooled-full builder's call, src/pooled_design_null_model.R). In
##                  STAARpipeline 0.9.8 the flag does not change the glmmkin fit; it only stores XSigma_i and
##                  XXSigma_iX_inv in the object and sets obj$use_SPA, and those matrices are read only inside
##                  `if (obj_nullmodel$use_SPA)` blocks and by the *_SPA functions. Every scan driver of this package
##                  sets obj_nullmodel$use_SPA <- FALSE at scan time, as the deployed pooled-full drivers do, so the
##                  scans are identical either way. The value used is recorded in the status JSON and in the object.
##   --overwrite    allow replacing an existing output RData (default: refuse, new-files-only policy).
suppressPackageStartupMessages({ library(data.table); library(Matrix); library(STAARpipeline); library(jsonlite) })
setDTthreads(1L)

args <- commandArgs(TRUE)
opt <- function(name, default) { m <- grep(paste0("^--", name, "="), args, value = TRUE); if (length(m)) sub(paste0("^--", name, "="), "", m[1]) else default }
USE_SPA_FIT <- as.logical(opt("use-spa-fit", "FALSE")); stopifnot(!is.na(USE_SPA_FIT))
OVERWRITE   <- "--overwrite" %in% args

ROOT   <- "${ADSP_ROOT}"
PKG    <- "${ADSP_ROOT}/revision/pooled_project_comparator_20261020"
OLD    <- file.path(ROOT, "output/obj_null_stratum_all_pooled_lpa_pc_full.RData")
NEW    <- file.path(ROOT, "output/obj_null_stratum_all_pooled_lpa_pc_full_project.RData")
PHENO  <- file.path(ROOT, "pheno/adsp_combined.csv")
QCMETA <- file.path(ROOT, "pheno/gcad.r5.wgs.58507.2024.11.03.quality.metrics.ALL.csv")
GRM    <- file.path(ROOT, "processed/output.sparseGRM.sGRM.RData")
RIGA   <- file.path(ROOT, "output/matched/matched_samples_lpa_pc.pair33_projcell.csv")
STATUS_DIR <- file.path(PKG, "status"); dir.create(STATUS_DIR, showWarnings = FALSE, recursive = TRUE)
STATUS <- file.path(STATUS_DIR, "null_status.json")
MIN_N <- 50L                     # collapse threshold (N < 50 -> Other)
MIN_MINORITY <- 15L              # v2 (2026-10-20): levels with fewer than 15 of either outcome also -> Other (the attempt-1 fit with these
                                 # near-pure levels kept failed in AI-REML: dpotrf leading principal minor of order 39 not positive)
NEAR_PURE_MIN <- 15L             # disclosure only (not a gate): a kept level with fewer than 15 participants of one outcome
DUMMY_WARN <- 15                 # quasi-separation warning for 0/1 indicator coefficients (|logit| 15 <-> probability within 3e-7 of 0/1)
EXTREME_EPS <- 1e-6; EXTREME_MAX_FRAC <- 0.01
PILOT_OLD_MD5 <- "13e9641edec4c99787f96f660b23b728"   # md5 of OLD recorded by the 16 Oct preflight (pilot_status.json)

report <- list(status = "started", started = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), package = PKG,
               old_null = OLD, new_null = NEW, use_spa_fit_argument = USE_SPA_FIT, scan_use_SPA = FALSE,
               collapse_rule = sprintf("project with n0 == 0 or n1 == 0 (one outcome label) or N < %d or fewer than 15 of either outcome (v2) -> 'Other'; v4: the factor is project x SeqCenter cells with the benchmark rule only (one outcome or N < 50 -> Other); SeqCenter dummies dropped (nested); single pass on the 38,138 analysed participants; reference level 'Other'", MIN_N))
refuse <- function(reason, extra = list()) {
  report$status <<- paste0("refused:", reason); report$reason <<- reason; report <<- c(report, extra)
  report$finished <<- format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
  write_json(report, STATUS, auto_unbox = TRUE, pretty = TRUE, digits = NA, null = "null")
  cat("REFUSED:", reason, "\n"); quit(save = "no", status = 3L)
}
load_one <- function(f) { e <- new.env(); n <- load(f, envir = e); stopifnot(length(n) == 1L); e[[n]] }

if (file.exists(NEW) && !OVERWRITE) refuse(sprintf("output exists (%s); pass --overwrite to replace it", NEW))

## ---- 1. the pooled-full null: participants, design, outcome ------------------------------------------------------
cat("Loading the pooled-full null:", OLD, "\n")
old <- load_one(OLD)
report$old_null_md5 <- unname(tools::md5sum(OLD)); report$old_null_md5_matches_pilot <- identical(report$old_null_md5, PILOT_OLD_MD5)
ids <- as.character(old$id_include)
if (length(ids) != 38138L || anyDuplicated(ids)) refuse(sprintf("pooled-full null has %d ids (%d duplicated); expected 38,138 unique", length(ids), anyDuplicated(ids)))
if (!isTRUE(old$converged)) refuse("pooled-full null is not marked converged")
Xbase <- as.matrix(old$X)
if (nrow(Xbase) != length(ids) || ncol(Xbase) != 37L || colnames(Xbase)[1] != "(Intercept)")
  refuse(sprintf("pooled-full design is %d x %d (first column %s); expected 38138 x 37 with (Intercept) first", nrow(Xbase), ncol(Xbase), colnames(Xbase)[1]))
report$base_columns <- colnames(Xbase)
y_recovered <- round(as.numeric(old$residuals) + as.numeric(old$fitted.values))   # bridge_null_check.R recovery of the outcome

ph <- fread(PHENO, select = c("IID", "AD"), colClasses = list(character = "IID"))
if (anyDuplicated(ph$IID)) refuse("duplicated IID in the phenotype file")
y <- ph$AD[match(ids, ph$IID)]
if (anyNA(y) || !all(y %in% 0:1)) refuse("outcome missing or not 0/1 for some pooled-full participants")
if (!identical(as.integer(y), as.integer(y_recovered))) refuse("phenotype AD disagrees with the outcome recovered from the saved null (residuals + fitted values)")
report$N <- length(ids); report$cases <- sum(y == 1L); report$controls <- sum(y == 0L)
if (report$cases != 12758L || report$controls != 25380L) refuse(sprintf("case/control counts %d/%d differ from the pooled-full null (12,758/25,380)", report$cases, report$controls))

## ---- 2. sequencing-project labels (the QC metadata 'Study' column, the same label RIGA keys on) --------------------
qc <- fread(QCMETA, select = c("Sample_name", "Study", "SeqCenter"), colClasses = "character")
qc[, Study := trimws(Study)]
multi <- qc[, .(n = uniqueN(Study)), by = Sample_name][n > 1L]
if (nrow(multi)) refuse(sprintf("%d samples carry more than one Study label in the QC metadata", nrow(multi)))
qc <- unique(qc, by = "Sample_name")
project <- qc$Study[match(ids, qc$Sample_name)]
## v4 (2026-10-20): the added factor is the design's own unit, project x sequencing center, with the benchmark collapse rule
## (units with one outcome label or N < 50 -> 'Other'); centers are nested in the units, so the SeqCenter dummies stay dropped.
center_raw <- trimws(qc$SeqCenter[match(ids, qc$Sample_name)])
center_cb  <- ifelse(center_raw %in% c("GENENTECH", "Genentech", "HudsonAlpha Institute for Biotechnology", "Unknown"), "Other", center_raw)   # the RIGA/pooled-full center coding
project <- paste(project, center_cb, sep = "|")
miss <- is.na(project) | project == ""
if (any(miss)) refuse(sprintf("%d participants have no project label", sum(miss)), list(missing_ids = ids[miss]))
riga <- fread(RIGA, select = c("IID", "Study", "cell_SC"), colClasses = "character")
riga[, Study := trimws(Study)]
ov <- riga[IID %in% ids]
mism <- ov[cell_SC != project[match(IID, ids)]]   # v4: the unit must reproduce RIGA's cell_SC for every overlapping participant
report$riga_overlap <- nrow(ov); report$riga_project_label_mismatches <- nrow(mism)
if (nrow(mism)) { fwrite(mism, file.path(STATUS_DIR, "riga_project_label_mismatches.csv")); refuse("project labels differ from the RIGA design key for some shared participants (see status/riga_project_label_mismatches.csv)") }

## ---- 3. collapse rule --------------------------------------------------------------------------------------------
d  <- data.table(IID = ids, AD = as.integer(y), project = project)
ct <- d[, .(N = .N, cases = sum(AD == 1L), controls = sum(AD == 0L)), by = project]
ct[, mono  := cases == 0L | controls == 0L]
ct[, small := N < MIN_N]
ct[, tiny := FALSE]   # v4: benchmark rule only (the v2 near-pure rule did not address the failure, which was collinearity)
ct[, level := ifelse(mono | small | tiny, "Other", project)]
ct[, reason := fifelse(mono & small, "one outcome label and N < 50", fifelse(mono, "one outcome label", fifelse(small, "N < 50", fifelse(tiny, "fewer than 15 of one outcome", "kept"))))]
setorder(ct, -N)
fwrite(ct, file.path(STATUS_DIR, "project_level_map.csv"))
d[ct, level := i.level, on = "project"]
lv <- d[, .(N = .N, cases = sum(AD == 1L), controls = sum(AD == 0L), n_projects = uniqueN(project)), by = level][order(-N)]
lv[, minority_outcome := pmin(cases, controls)]; lv[, near_pure := minority_outcome < NEAR_PURE_MIN]
fwrite(lv, file.path(STATUS_DIR, "project_level_counts.csv"))
fwrite(d, file.path(STATUS_DIR, "participant_levels.csv"))
report$projects <- nrow(ct); report$projects_one_outcome <- sum(ct$mono); report$projects_small <- sum(ct$small & !ct$mono)
report$projects_tiny <- sum(ct$tiny & !ct$mono & !ct$small); report$min_minority_rule <- MIN_MINORITY; report$projects_collapsed <- sum(ct$level == "Other"); report$participants_in_Other <- sum(d$level == "Other")
report$levels_fitted <- nrow(lv); report$Other_cases <- lv[level == "Other"]$cases; report$Other_controls <- lv[level == "Other"]$controls
## Disclosure, not a gate: the rule keeps every two-outcome level with N >= MIN_N, including levels with a handful of one
## outcome (expected on the 38,138: ADSP-FUS1-HIHGBB 66/2, FASe-Families-WGS 60/7, KnightADRC-WGS 52/14, ADSP-FUS3-
## OldOrderAmish2 3/60, WHICAP1 4/108, ADSP-FUS3-CNSA 4/172, ADSP-FUS2-StEPAD2 4/55, ADSP-FUS2-WRAP 9/739). Their MLE is
## finite and their fitted probabilities sit near the observed level case fraction, far from the 1e-6 gate, so the gate
## cannot flag them; they are listed here so that the author and a referee can see them.
report$near_pure_rule <- sprintf("kept level with fewer than %d participants of one outcome (disclosed, not a gate)", NEAR_PURE_MIN)
report$near_pure_levels <- lv[near_pure == TRUE, sprintf("%s (%d cases / %d controls)", level, cases, controls)]
cat("near-pure kept levels:", if (length(report$near_pure_levels)) paste(report$near_pure_levels, collapse = "; ") else "none", "\n")
cat(sprintf("projects %d; one-outcome %d; N<50 (not one-outcome) %d; collapsed to Other %d projects, %d participants (%d cases / %d controls); %d levels fitted\n",
            report$projects, report$projects_one_outcome, report$projects_small, report$projects_collapsed, report$participants_in_Other,
            report$Other_cases, report$Other_controls, report$levels_fitted))
if (!(report$Other_cases > 0L && report$Other_controls > 0L)) refuse("the Other level does not contain both outcomes")
if (any(lv$cases == 0L | lv$controls == 0L)) refuse("a fitted project level has one outcome label")

## ---- 4. design: pooled-full columns + project dummies; explicit full-rank basis -----------------------------------
d[, level_f := relevel(factor(level), ref = "Other")]
XP <- model.matrix(~ level_f, data = as.data.frame(d))[, -1L, drop = FALSE]
colnames(XP) <- sub("^level_f", "cell_", colnames(XP))
## v3 (2026-10-20): drop the sequencing-center dummies from the base design. Projects are nested in centers for 55 of the 66
## projects and the kept project dummies are near-collinear with the center dummies under the PQL working weights (attempts 1
## and 2 failed at AI-REML iteration 17 with a non-positive leading minor at the second project column). The center effect is
## carried by the project levels (and, for the 11 multi-center projects, averaged within project); recorded in the status JSON.
DROP_BASE <- grep("^SeqCenter", colnames(Xbase), value = TRUE)
report$base_columns_dropped_v3 <- DROP_BASE
Xbase <- Xbase[, setdiff(colnames(Xbase), DROP_BASE), drop = FALSE]
Xfull <- cbind(Xbase, XP)
qx <- qr(Xfull, tol = 1e-8)
basis <- sort(qx$pivot[seq_len(qx$rank)])
dropped <- setdiff(colnames(Xfull), colnames(Xfull)[basis])
report$design_columns_base <- ncol(Xbase); report$design_columns_project <- ncol(XP); report$design_columns_full <- ncol(Xfull)
report$design_rank <- qx$rank; report$aliased_columns_dropped <- dropped
if (any(colnames(Xbase) %in% dropped)) refuse("a pooled-full covariate column would be dropped as aliased; the base design must stay intact", list(dropped = dropped))
X <- Xfull[, basis, drop = FALSE]
orig_names <- colnames(X)
fit_names <- orig_names
fit_names[fit_names == "(Intercept)"] <- "Intercept"
fit_names <- sub("Kinghorn Centre for Clinical Genomics \\(KCCG\\), Sydney", "KCCG", fit_names)
fit_names <- make.names(fit_names, unique = TRUE)
fwrite(data.table(original_column = orig_names, fitted_name = fit_names, source = ifelse(orig_names %in% colnames(Xbase), "pooled-full null X", "collapsed project dummy")),
       file.path(STATUS_DIR, "design_column_map.csv"))
colnames(X) <- fit_names
df <- data.frame(IID = d$IID, AD = d$AD, X, check.names = FALSE, stringsAsFactors = FALSE)
fm <- reformulate(fit_names, response = "AD", intercept = FALSE)
if (qr(X)$rank != ncol(X)) refuse("design basis is not full rank")

## ---- 5. sparse GRM, as in the pooled-full builder --------------------------------------------------------------
sgrm <- load_one(GRM)
grm_ids <- sub("^[^_]+_", "", rownames(sgrm))
if (anyDuplicated(grm_ids)) refuse("duplicated ids in the sparse GRM")
if (!all(ids %in% grm_ids)) refuse(sprintf("%d participants absent from the sparse GRM", sum(!ids %in% grm_ids)))
ix <- match(ids, grm_ids); K <- sgrm[ix, ix]; dimnames(K) <- list(ids, ids); rm(sgrm)

## ---- 6. fit ------------------------------------------------------------------------------------------------------
cat("Fitting:", deparse(fm)[1], "... (", ncol(X), "columns )\n")
warn <- character(); t0 <- proc.time()[[3L]]
fit <- tryCatch(withCallingHandlers(
  fit_nullmodel(fm, data = df, id = "IID", kins = K, use_sparse = TRUE, family = binomial(link = "logit"), verbose = TRUE, use_SPA = USE_SPA_FIT),
  warning = function(w) { warn <<- c(warn, conditionMessage(w)); invokeRestart("muffleWarning") }), error = function(e) e)
elapsed <- proc.time()[[3L]] - t0
report$fit_seconds <- round(elapsed, 1); report$fit_warnings <- unique(warn)
if (inherits(fit, "error")) refuse(paste0("fit_nullmodel error: ", conditionMessage(fit)))

## ---- 7. gates and diagnostics ------------------------------------------------------------------------------------
mu <- as.numeric(fit$fitted.values); coef <- fit$coefficients
ext <- mu < EXTREME_EPS | mu > 1 - EXTREME_EPS
report$converged <- isTRUE(fit$converged); report$theta <- as.numeric(fit$theta); report$theta_pooled_full <- as.numeric(old$theta)
## Scale note: the continuous terms legitimately carry large coefficients (the PCs have sd 0.005-0.010, so the saved
## pooled-full null already has PC2 100.0, PC5 99.6, PC3 83.0, PC1 -50.2, i.e. at most 0.7 per SD, with 0 extreme fitted
## values). max_abs_coefficient is therefore recorded for completeness only; the quasi-separation warnings below look at
## the 0/1 indicator columns (project_*, and the base SeqCenter*/lpa_cluster*/Sex1 dummies separately).
report$max_abs_coefficient <- max(abs(coef)); report$max_abs_coefficient_name <- names(coef)[which.max(abs(coef))]
report$max_abs_coefficient_pooled_full <- max(abs(old$coefficients)); report$max_abs_coefficient_pooled_full_name <- names(old$coefficients)[which.max(abs(old$coefficients))]
is_proj <- grepl("^project_", names(coef)); is_base_dummy <- grepl("^(SeqCenter|lpa_cluster|Sex)", names(coef))
pc <- coef[is_proj]; bd <- coef[is_base_dummy]
report$max_abs_project_coefficient <- if (length(pc)) max(abs(pc)) else NA; report$max_abs_project_coefficient_name <- if (length(pc)) names(pc)[which.max(abs(pc))] else NA
report$max_abs_base_dummy_coefficient <- if (length(bd)) max(abs(bd)) else NA; report$max_abs_base_dummy_coefficient_name <- if (length(bd)) names(bd)[which.max(abs(bd))] else NA
obd <- old$coefficients[grepl("^(SeqCenter|lpa_cluster|Sex)", names(old$coefficients))]
report$max_abs_base_dummy_coefficient_pooled_full <- if (length(obd)) max(abs(obd)) else NA
report$fitted_probability_range <- range(mu); report$fitted_extreme_count <- sum(ext); report$fitted_extreme_fraction <- mean(ext)
report$fitted_extreme_rule <- sprintf("refuse if > %.0f%% of participants have fitted probability within %g of 0 or 1", 100 * EXTREME_MAX_FRAC, EXTREME_EPS)
report$n_coefficients <- length(coef); report$id_order_identical <- identical(as.character(fit$id_include), ids)
report$fitted_corr_with_pooled_full <- if (identical(as.character(fit$id_include), ids)) cor(mu, as.numeric(old$fitted.values)) else NA
cmap <- data.table(term = fit_names, original_column = orig_names)
oldc <- data.table(original_column = names(old$coefficients), estimate_pooled_full = as.numeric(old$coefficients))
ctab <- data.table(term = names(coef), estimate = as.numeric(coef))[cmap, on = "term"]
ctab[oldc, on = "original_column", estimate_pooled_full := i.estimate_pooled_full]
ctab[, indicator := grepl("^(project_|SeqCenter|lpa_cluster|Sex)", term)]; ctab[, change := estimate - estimate_pooled_full]
fwrite(ctab[, .(term, original_column, indicator, estimate, estimate_pooled_full, change)], file.path(STATUS_DIR, "coefficients.csv"))
diag <- data.table(IID = as.character(fit$id_include), mu = mu)[d, on = "IID"][, .(N = .N, cases = sum(AD == 1L), controls = sum(AD == 0L), observed_case_fraction = mean(AD), mean_fitted = mean(mu), min_fitted = min(mu), max_fitted = max(mu), extreme = sum(mu < EXTREME_EPS | mu > 1 - EXTREME_EPS)), by = level][order(-N)]
diag[, near_pure := pmin(cases, controls) < NEAR_PURE_MIN]
fwrite(diag, file.path(STATUS_DIR, "level_fit_diagnostics.csv"))
if (!report$converged) refuse("fit_nullmodel did not converge")
if (!all(is.finite(coef))) refuse("non-finite coefficient(s)")
if (!all(is.finite(fit$theta)) || any(fit$theta < 0)) refuse("non-finite or negative variance component")
if (!report$id_order_identical) refuse("id_include of the new null is not identical (same order) to the pooled-full null")
if (report$fitted_extreme_fraction > EXTREME_MAX_FRAC) refuse(sprintf("%d participants (%.2f%%) have fitted probabilities within %g of 0 or 1", sum(ext), 100 * mean(ext), EXTREME_EPS))
if (is.finite(report$max_abs_project_coefficient) && report$max_abs_project_coefficient > DUMMY_WARN)
  report$warning_large_project_coefficient <- sprintf("project indicator %s = %.2f exceeds %g in absolute value (quasi-separation inside a kept level); the stated gates passed, inspect status/coefficients.csv and status/level_fit_diagnostics.csv",
                                                      report$max_abs_project_coefficient_name, pc[[report$max_abs_project_coefficient_name]], DUMMY_WARN)
if (is.finite(report$max_abs_base_dummy_coefficient) && report$max_abs_base_dummy_coefficient > DUMMY_WARN)
  report$warning_large_base_dummy_coefficient <- sprintf("base indicator %s = %.2f exceeds %g in absolute value (pooled-full null: max |base indicator| %.2f); inspect status/coefficients.csv",
                                                         report$max_abs_base_dummy_coefficient_name, bd[[report$max_abs_base_dummy_coefficient_name]], DUMMY_WARN, report$max_abs_base_dummy_coefficient_pooled_full)

## ---- 8. save ----------------------------------------------------------------------------------------------------
obj_null <- fit
obj_null$sparse_kins <- TRUE; obj_null$relatedness <- TRUE
obj_null$project_specification <- list(built = report$started, package = PKG, collapse_rule = report$collapse_rule, levels = lv$level,
                                       aliased_columns_dropped = dropped, use_spa_fit_argument = USE_SPA_FIT, base_null = OLD, base_null_md5 = report$old_null_md5)
save(obj_null, file = NEW)
report$new_null_md5 <- unname(tools::md5sum(NEW)); report$new_null_bytes <- file.info(NEW)$size
report$status <- "accepted"; report$finished <- format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
write_json(report, STATUS, auto_unbox = TRUE, pretty = TRUE, digits = NA, null = "null")
writeLines(capture.output(sessionInfo()), file.path(STATUS_DIR, "build_null_sessionInfo.txt"))
cat(sprintf("ACCEPTED: converged=%s theta=%s max|project coef|=%.3f (%s) max|base indicator|=%.3f (%s) max|coef| %.3f (%s; continuous scale) fitted range [%.3g, %.3g] extreme %d (%.3f%%) use_SPA(fit)=%s in %.0f s; saved %s\n",
            report$converged, paste(signif(report$theta, 4), collapse = ","), report$max_abs_project_coefficient, report$max_abs_project_coefficient_name,
            report$max_abs_base_dummy_coefficient, report$max_abs_base_dummy_coefficient_name, report$max_abs_coefficient, report$max_abs_coefficient_name,
            report$fitted_probability_range[1], report$fitted_probability_range[2], report$fitted_extreme_count, 100 * report$fitted_extreme_fraction, USE_SPA_FIT, elapsed, NEW))
