# UKB 200k WGS replication — step 1: assemble the analysis set.
#
# Phenotype : phecode 290.11 "Alzheimer's disease" (column 373 of the phecode label table).
#             Controls are PheWAS exclusion-range controls: anyone carrying a code in 290-292.99
#             who is not a case is set to NA, so the control set is conservative.
# Covariates: sex, age, age^2, PC1-PC10 (UKB whole-cohort array PCs, field 22009), genotyping batch.
# Relatedness: the pre-computed UKB sparse GRM, subset to the analysis set.
#
# NOTE ON PARSING: read.table/read.csv with proper quoting is used throughout. An earlier
# attempt at this project used a naive comma split and silently shifted every column after a
# quoted field containing a comma; do not reintroduce that.

suppressPackageStartupMessages({ library(data.table); library(Matrix) })

REV  <- "${ADSP_ROOT}/revision"
UKB  <- file.path(REV, "ukb")
dir.create(UKB, showWarnings = FALSE, recursive = TRUE)

PHECODE <- "${UKB_PHECODE_DIR}/phecode/pheno_covar/phecode_labels_exclude_feb8_2024.txt"
COVAR   <- "${UKB_PHECODE_DIR}/phecode/pheno_covar/2age_sex_inter_gbatch_40gpc.txt"
GRM     <- "${HOME_DIR}/overflow/UKB_sparseGRM/output.sparseGRM.sGRM.RData"
WGSIDS  <- file.path(REV, "logs/ukb_wgs_sample_ids.txt")

wgs <- readLines(WGSIDS)
cat("aGDS sample IDs:", length(wgs), "\n")

## ---- phenotype: FID, IID, and phecode 290.11 (column 373) -------------------
hdr <- scan(PHECODE, what = "", nlines = 1, sep = "\t", quiet = TRUE)
col_ad <- which(hdr == "290.11")
stopifnot(length(col_ad) == 1L)
cat("phecode 290.11 is column", col_ad, "of", length(hdr), "\n")

ph <- fread(PHECODE, select = c(1, 2, col_ad), sep = "\t",
            colClasses = list(character = 1:2))
setnames(ph, c("FID", "IID", "AD"))
ph[, IID := as.character(IID)]
cat("phecode rows:", nrow(ph), "\n")
print(table(ph$AD, useNA = "ifany"))

## ---- covariates -------------------------------------------------------------
cv <- fread(COVAR, colClasses = list(character = 1:2))
setnames(cv, 1:2, c("FID", "IID"))
cv[, IID := as.character(IID)]
keep_cov <- c("IID", "sex", "age", "age2", paste0("PC", 1:10), "gbatch")
missing_cov <- setdiff(keep_cov, names(cv))
if (length(missing_cov)) stop("missing covariate columns: ", paste(missing_cov, collapse = ", "))
cv <- cv[, ..keep_cov]
cat("covariate rows:", nrow(cv), "\n")

## ---- merge, restrict to WGS, complete cases ---------------------------------
dat <- merge(ph[, .(IID, AD)], cv, by = "IID")
dat <- dat[IID %in% wgs]
cat("after WGS restriction:", nrow(dat), "\n")
dat <- dat[!is.na(AD)]
cat("after dropping AD == NA (PheWAS exclusion range):", nrow(dat), "\n")
dat <- dat[complete.cases(dat)]
cat("after complete-case on covariates:", nrow(dat), "\n")

# Drop withdrawn-placeholder negative IDs if any leaked through
neg <- grepl("^-", dat$IID)
if (any(neg)) { cat("dropping", sum(neg), "negative (withdrawn-placeholder) IDs\n"); dat <- dat[!neg] }

dat[, AD := as.integer(AD)]
cat("\nFINAL ANALYSIS SET\n")
cat("  n        =", nrow(dat), "\n")
cat("  cases    =", sum(dat$AD == 1), "\n")
cat("  controls =", sum(dat$AD == 0), "\n")
cat("  case:control = 1 :", round(sum(dat$AD == 0) / sum(dat$AD == 1), 1), "\n")
cat("  age  mean", round(mean(dat$age), 1), " sd", round(sd(dat$age), 1), "\n")
cat("  sex  table:\n"); print(table(dat$sex))
cat("  gbatch levels:", length(unique(dat$gbatch)), "\n")

## ---- sparse GRM subset ------------------------------------------------------
e <- new.env(); load(GRM, envir = e)
sg <- get(ls(e)[1], envir = e)
ids <- sub("^[^_]+_", "", rownames(sg))          # "FID_IID" -> "IID"
rownames(sg) <- ids; colnames(sg) <- ids
have <- intersect(dat$IID, ids)
cat("\nGRM covers", length(have), "of", nrow(dat), "analysis subjects\n")
dat <- dat[IID %in% have]
sg  <- sg[dat$IID, dat$IID]
cat("after restricting to GRM coverage: n =", nrow(dat),
    " cases =", sum(dat$AD == 1), " controls =", sum(dat$AD == 0), "\n")

dat <- as.data.frame(dat)
saveRDS(list(pheno = dat, grm = sg), file.path(UKB, "ukb_ad_analysis_set.rds"))
write.csv(dat, file.path(UKB, "ukb_ad_pheno.csv"), row.names = FALSE, quote = TRUE)
cat("\nSaved", file.path(UKB, "ukb_ad_analysis_set.rds"), "\n")

## ---- agds_dir for UKB -------------------------------------------------------
agds_dir <- sprintf("${HOME_DIR}/overflow/wgs_agds/ukb.200k.wgs.chr%d.pass.annotated.gds", 1:22)
stopifnot(all(file.exists(agds_dir)))
save(agds_dir, file = file.path(UKB, "agds_dir_ukb.RData"))
cat("Saved UKB agds_dir (22 autosomes, all present)\n")
cat("\nDONE\n")
