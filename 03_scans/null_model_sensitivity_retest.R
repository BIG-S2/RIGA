#!/usr/bin/env Rscript
## null_model_sensitivity_retest.R (2026-10-13; v3 = v2 plus a fallback when the mixed-model fit is singular: cells with fewer than 20
## participants after the exclusion are dropped and the fit is repeated once; used for loo:27, whose v2 fit failed in the Cholesky step): refit the stratified-design null under one specification
## and re-test every signal reported in the paper. SPEC (env):
##   full        the deployed specification (sanity check: must reproduce the deployed p values)
##   nocls       no latent-class terms (Age + Sex + cell + PC1-10)            [point 7 ablation]
##   pc20        deployed covariates with PC1-20                              [point 7 PC sensitivity]
##   noAPOEsel   drop the APOE-genotype-selected projects                     [point 4]
##   noSelected  drop APOE-genotype- and onset-selected projects (EOAD1/2)    [point 4]
##   loo:<k>     drop the k-th sequencing project (alphabetical, 1-based)     [point 4 leave-one-project-out]
## Null fitting follows split_half_null_models.R (classes < 30 folded, aliased class dummies folded, sparse GRM,
## use_SPA = TRUE at fit time). Tests follow the deployed drivers (QC_label7 for gene-centric, rare_maf_cutoff 0.01,
## rv_num_cutoff 2, annotation weights). Each unit is tested with use_SPA FALSE (deployed STAAR-O / score test) and TRUE
## (saddlepoint burden / SPA score test). Writes revision/tables/review/targeted_<tag>.csv.
suppressPackageStartupMessages({ library(data.table); library(SeqArray); library(STAAR); library(STAARpipeline); library(GMMAT) })
SPEC <- Sys.getenv("SPEC", "full"); OUT <- "${ADSP_ROOT}/output"; TAB <- "${ADSP_ROOT}/revision/tables/review"
dir.create(TAB, showWarnings = FALSE, recursive = TRUE)
m <- fread(file.path(OUT, "matched/matched_samples_lpa_pc.pair33_projcell.csv"), colClasses = list(character = c("IID", "lpa_cluster")))
m[, project := sub("\\|.*$", "", cell_SC)]; projects <- sort(unique(m$project))
excl <- character(0); covset <- "full"; npc <- 10L; tag <- SPEC
if (SPEC == "nocls") covset <- "nocls"
if (SPEC == "pc20") npc <- 20L
if (SPEC == "noAPOEsel") excl <- c("ADSP-FUS1-APOEextremes", "APOEExtremes1")
if (SPEC == "noSelected") excl <- c("ADSP-FUS1-APOEextremes", "APOEExtremes1", "EOAD1", "EOAD2")
if (startsWith(SPEC, "loo:")) { k <- as.integer(sub("loo:", "", SPEC)); excl <- projects[k]; tag <- sprintf("loo_%02d", k) }
m <- m[!project %in% excl]
ct <- m[, .(N = .N, cases = sum(AD == 1)), by = cell_SC]; m <- m[cell_SC %in% ct[cases > 0 & cases < N]$cell_SC]
cat(sprintf("SPEC %s | excluded %s | N %d (cases %d) in %d cells\n", SPEC, paste(excl, collapse = ","), nrow(m), sum(m$AD == 1), uniqueN(m$cell_SC)))
m[, Sex := factor(Sex)]; m[, cell_SC := factor(cell_SC)]
tab <- sort(table(m$lpa_cluster), decreasing = TRUE); big <- names(tab)[1]; m[lpa_cluster %in% names(tab)[tab < 30], lpa_cluster := big]
PCS <- paste0("PC", seq_len(npc))
FORM <- if (covset == "nocls") as.formula(paste("AD ~ Age + Sex + cell_SC +", paste(PCS, collapse = " + "))) else
  as.formula(paste("AD ~ Age + Sex + cell_SC + lpa_cluster + cluster_prob + entropy_norm +", paste(PCS, collapse = " + ")))
repeat { m[, lpa_cluster := factor(lpa_cluster)]; X <- model.matrix(FORM, data = as.data.frame(m)); q <- qr(X); if (q$rank == ncol(X)) break
  al <- colnames(X)[q$pivot[(q$rank + 1):ncol(X)]]; lv <- sub("^lpa_cluster", "", al[grepl("^lpa_cluster", al)])
  if (!length(lv)) stop("aliased non-class column: ", paste(al, collapse = ","))
  m[, lpa_cluster := as.character(lpa_cluster)]; m[lpa_cluster %in% lv, lpa_cluster := big] }
sgrm <- get(load("${ADSP_ROOT}/processed/output.sparseGRM.sGRM.RData"))
gid <- rownames(sgrm); if (is.null(gid)) gid <- colnames(sgrm); gid_s <- sub("^[^_]*_", "", gid)
idx <- match(m$IID, gid_s); stopifnot(!any(is.na(idx))); K <- sgrm[idx, idx]; dimnames(K) <- list(m$IID, m$IID); rm(sgrm)
df <- as.data.frame(m); rownames(df) <- df$IID
fitk <- function(df, K) fit_nullmodel(FORM, data = df, kins = K, id = "IID", use_sparse = TRUE, family = binomial(link = "logit"), verbose = FALSE, use_SPA = TRUE)
obj <- tryCatch(fitk(df, K), error = function(e) { cat("fit failed:", conditionMessage(e), "\n"); NULL })
if (is.null(obj)) { small <- names(which(table(m$cell_SC) < 20)); cat("dropping", length(small), "cells with < 20 participants:", paste(small, collapse = ","), "\n")
  m <- m[!cell_SC %in% small]; m[, cell_SC := factor(as.character(cell_SC))]; idx <- match(m$IID, gid_s); K <- K[m$IID, m$IID]; df <- as.data.frame(m); rownames(df) <- df$IID
  X <- model.matrix(FORM, data = df); obj <- tryCatch(fitk(df, K), error = function(e) { cat("fit failed again:", conditionMessage(e), "\n"); NULL }) }
if (is.null(obj) && covset != "nocls") {   # second fallback: fold latent classes with < 100 members (a class concentrated in the excluded project is near-aliased)
  m[, lpa_cluster := as.character(lpa_cluster)]; tab <- table(m$lpa_cluster); small <- names(tab)[tab < 100]; cat("folding", length(small), "classes with < 100 members into", big, "\n")
  m[lpa_cluster %in% small, lpa_cluster := big]; m[, lpa_cluster := factor(lpa_cluster)]; df <- as.data.frame(m); rownames(df) <- df$IID; X <- model.matrix(FORM, data = df)
  obj <- tryCatch(fitk(df, K), error = function(e) { cat("fit failed after folding:", conditionMessage(e), "\n"); NULL }) }
if (is.null(obj)) {   # last resort: drop the latent-class terms for this specification and say so
  FORM <- as.formula(paste("AD ~ Age + Sex + cell_SC +", paste(PCS, collapse = " + "))); cat("refitting WITHOUT latent-class terms\n"); X <- model.matrix(FORM, data = df); obj <- fitk(df, K) }
cat(sprintf("null: converged %s, theta %s, design columns %d\n", obj$converged, paste(round(obj$theta, 3), collapse = ","), ncol(X)))
obj$sparse_kins <- TRUE; obj$relatedness <- TRUE
dir.create("${ADSP_ROOT}/output/matched/retest_nulls", showWarnings = FALSE)
save(obj, file = sprintf("${ADSP_ROOT}/output/matched/retest_nulls/obj_%s.RData", gsub("[^A-Za-z0-9_]", "_", tag)))
agds_dir <- get(load(file.path(OUT, "agds_dir.Rdata"))); catalog <- get(load(file.path(OUT, "Annotation_name_catalog.Rdata")))
ANN <- c("CADD", "LINSIGHT", "FATHMM.XF", "aPC.EpigeneticActive", "aPC.EpigeneticRepressed", "aPC.EpigeneticTranscription",
         "aPC.Conservation", "aPC.LocalDiversity", "aPC.Mappability", "aPC.TF", "aPC.Protein")
res <- list()
add <- function(...) res[[length(res) + 1]] <<- data.table(tag = tag, spec = SPEC, N = nrow(m), cases = sum(m$AD == 1), ...)
pick <- function(r, col) if (!is.null(r) && col %in% names(r)) as.numeric(r[[col]][1]) else NA_real_
flat <- function(r) {   # v2: a single-category call returns a matrix (often of mode list); an all-category call returns a list
  if (is.null(r) || !length(r)) return(NULL)
  m <- if (is.matrix(r) || is.data.frame(r)) r else r[[1]]
  if (is.null(m) || !length(m)) return(NULL)
  d <- as.data.frame(m, stringsAsFactors = FALSE, check.names = FALSE)
  for (cc in names(d)) if (is.list(d[[cc]])) d[[cc]] <- unlist(lapply(d[[cc]], function(z) if (length(z)) z[1] else NA))
  d }
## single variants (deployed scan: QC_label2 then QC7 post-filter; exact alleles are matched below)
sv <- data.table(locus = c("APOE", "BIN1", "SHARPIN", "KCNV1", "CELF4", "BAIAP2L1"), chr = c(19, 2, 8, 8, 18, 7),
                 pos = c(44909967, 127133887, 144103704, 110888061, 38835190, 98414763), ref = c("TGG", "T", "G", "CAT", "G", "GGAAT"), alt = c("T", "TGGG", "A", "C", "A", "G"))
for (i in seq_len(nrow(sv))) for (spa in c(FALSE, TRUE)) {
  obj$use_SPA <- spa; g <- seqOpen(agds_dir[sv$chr[i]])
  r <- tryCatch(Individual_Analysis(chr = sv$chr[i], start_loc = sv$pos[i] - 1, end_loc = sv$pos[i] + 1, genofile = g, obj_nullmodel = obj, mac_cutoff = 1,
                QC_label = "annotation/info/QC_label7", variant_type = "variant", geno_missing_imputation = "mean"), error = function(e) NULL); seqClose(g)
  r <- if (is.null(r)) NULL else as.data.table(r)[POS == sv$pos[i] & REF == sv$ref[i] & ALT == sv$alt[i]]
  add(type = "single variant", unit = sv$locus[i], test = if (spa) "SPA score" else "score", p = if (is.null(r) || !nrow(r)) NA_real_ else r$pvalue[1],
      est = if (!is.null(r) && nrow(r) && "Est" %in% names(r)) r$Est[1] else NA_real_, n_var = NA_integer_, cmac = NA_real_) }
## gene-centric sets reported at FDR 5% (coding and noncoding)
gc <- data.table(gene = c("TREM2", "NECTIN2", "NECTIN2", "TOMM40", "TOMM40", "TOMM40", "APOE", "CLASRP", "PPP1R37", "TREM2", "MMP17", "MMP17", "PGAM5", "MIPOL1"),
                 chr = c(6, 19, 19, 19, 19, 19, 19, 19, 19, 6, 12, 12, 12, 14),
                 cat = c("missense", "missense", "enhancer_DHS", "promoter_DHS", "UTR", "downstream", "enhancer_DHS", "downstream", "enhancer_DHS", "promoter_DHS", "upstream", "promoter_DHS", "promoter_DHS", "downstream"))
for (i in seq_len(nrow(gc))) for (spa in c(FALSE, TRUE)) {
  obj$use_SPA <- spa; g <- seqOpen(agds_dir[gc$chr[i]]); coding <- gc$cat[i] == "missense"
  f <- if (coding) Gene_Centric_Coding else Gene_Centric_Noncoding
  r <- tryCatch(f(chr = gc$chr[i], gene_name = gc$gene[i], category = gc$cat[i], genofile = g, obj_nullmodel = obj, rare_maf_cutoff = 0.01, rv_num_cutoff = 2,
                  QC_label = "annotation/info/QC_label7", variant_type = "SNV", geno_missing_imputation = "mean", Annotation_dir = "annotation/info/FunctionalAnnotation",
                  Annotation_name_catalog = catalog, Use_annotation_weights = TRUE, Annotation_name = ANN), error = function(e) { cat("ERR", gc$gene[i], gc$cat[i], conditionMessage(e), "\n"); NULL })
  seqClose(g)
  r1 <- flat(r)
  add(type = if (coding) "coding" else "noncoding", unit = paste(gc$gene[i], gc$cat[i]), test = if (spa) "SPA burden (STAAR-B)" else "STAAR-O",
      p = if (is.null(r1)) NA_real_ else if (spa) pick(r1, "STAAR-B") else pick(r1, "STAAR-O"), est = NA_real_,
      n_var = if (is.null(r1)) NA_integer_ else as.integer(pick(r1, "#SNV")), cmac = if (is.null(r1)) NA_real_ else pick(r1, "cMAC")) }
out <- rbindlist(res, fill = TRUE); print(out[, .(type, unit, test, p = signif(p, 3), n_var, cmac)])
fwrite(out, file.path(TAB, sprintf("targeted_%s.csv", gsub("[^A-Za-z0-9_]", "_", tag)))); cat("DONE\n")
