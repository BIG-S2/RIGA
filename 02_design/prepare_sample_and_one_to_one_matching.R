# NO-FILTER + TRUE 1:1 PAIR MATCHING  -- the missing cell of the Supp. Fig. S7 design.
# ============================================================================
# The three published S7 arms are
#     noRIGA  no filter, no matching                 N = 38,138
#     filter  confidence filter, no matching         N = 15,888
#     RIGA    confidence filter + exact matching     N = 14,419
# so "matching, but WITHOUT the confidence filter" has never been run.
#
# WHAT THIS DOES DIFFERENTLY FROM src/adsp_match.lpa_pc.R
#   1. the confidence filter (cluster_prob >= 0.7 & entropy_norm <= 0.5) is NOT applied
#   2. matching is genuine 1:1 PAIR matching, not stratification.
#
# Point 2 needs saying plainly, because every existing RIGA script is misleading about it:
# they call matchit(..., method = "exact", ratio = 1), but method="exact" STRATIFIES and
# silently IGNORES ratio=. The published RIGA arm is therefore 5,149:9,270 = 1:1.80, not 1:1.
# Here each case is paired with exactly one control, nearest in Mahalanobis distance over
# PC1..PC10, within blocks that are exact on Sex + age_band + SeqCenter + lpa_cluster --
# i.e. RIGA's own blocking variables, so the arm stays comparable to it.
#
# Age is kept as an exact age_band constraint rather than a continuous caliper. Measured
# (revision/tables/S15_adsp_pairmatch_counts.csv): the caliper version retains 21,258 vs
# 19,408, a 9.5% gain, which does not justify deviating from RIGA's cut(Age, seq(30,100,3)).
#
# TWO NULL MODELS ARE FITTED ON THE IDENTICAL COHORT, deliberately.
# In S7 the covariate set is a deterministic function of the matching factor -- the unmatched
# arms carry Sex + SeqCenter + lpa_cluster as covariates, the matched arm drops them because
# matching absorbed them. So the matching effect is perfectly aliased with a 19-column change
# to the null and is identified in no existing cell. Fitting both formulas on one cohort
# breaks that alias:
#     pair14  AD ~ Age + cluster_prob + entropy_norm + PC1..PC10          (RIGA's formula)
#     pair33  AD ~ Age + Sex + SeqCenter + lpa_cluster + cluster_prob
#                + entropy_norm + PC1..PC10                              (noRIGA's formula)
#
# use_SPA = TRUE at fit time, exactly as the other arms: it only ADDS the SPA fields, the
# glmmkin fit is identical, and the coding driver sets obj_nullmodel$use_SPA <- FALSE so the
# scan returns STAAR-O. Fitting with use_SPA = FALSE would permanently omit those fields.

suppressPackageStartupMessages({
  library(data.table); library(dplyr); library(MatchIt); library(readr); library(tidyr)
  library(gdsfmt); library(SeqArray); library(SeqVarTools)
  library(STAAR);  library(STAARpipeline)
})
set.seed(0)

cov_wgs <- read_csv("${ADSP_ROOT}/pheno/adsp_combined.csv",
                    show_col_types = FALSE) %>% drop_na(AD, Sex, Age)
cov_wgs$Age <- as.numeric(gsub("\\+", "", cov_wgs$Age))

## SAME DA-EM cluster file as the RIGA winner -- only the filter line is absent.
clu <- fread("${ADSP_ROOT}/adsp_lpa_clusters_20PC.tsv"); setDT(clu)
K <- length(grep("^post_k[0-9]+$", names(clu)))
if (!"entropy_norm" %in% names(clu)) clu[, entropy_norm := entropy / log(K)]
cat(sprintf("NO confidence filter: keeping all %d clustered subjects (the filter would keep %d).\n",
            nrow(clu), sum(clu$cluster_prob >= 0.7 & clu$entropy_norm <= 0.5)))

seqmeta <- read_csv("pheno/gcad.r5.wgs.58507.2024.11.03.quality.metrics.ALL.csv",
                    show_col_types = FALSE) %>%
  transmute(IID = Sample_name, SeqCenter = as.character(SeqCenter))

dat0 <- merge(cov_wgs, clu, by = "IID") %>%
  left_join(seqmeta, by = "IID") %>%
  mutate(Sex = as.factor(Sex), AD = as.factor(AD), SeqCenter = as.character(SeqCenter))
dat0$SeqCenter[is.na(dat0$SeqCenter)] <- "Unknown"
dat0$SeqCenter[dat0$SeqCenter %in% c(
  "MEDGENOME","Unknown","Genentech","GENENTECH",
  "HudsonAlpha Institute for Biotechnology","Macrogen","MACROGEN")] <- "Other"
dat0$SeqCenter <- droplevels(as.factor(dat0$SeqCenter))
dat0 <- dat0 %>% tidyr::drop_na(Age, Sex, AD, SeqCenter) %>%
  mutate(age_band = cut(Age, breaks = seq(30, 100, by = 3), right = TRUE)) %>%
  tidyr::drop_na(age_band)
dat0$lpa_cluster  <- as.factor(dat0$lpa_cluster)
dat0$entropy_norm <- if (K > 1) dat0$entropy / log(K) else 0
cat(sprintf("Pool entering matching: %d  (cases %d / controls %d)\n",
            nrow(dat0), sum(dat0$AD == 1), sum(dat0$AD == 0)))

## ---- 1:1 pair matching ----------------------------------------------------
# MatchIt stores `data` as an unevaluated expression; match.data() cannot resolve an inline
# as.data.frame(...) later ("A valid dataset could not be found"). Bind a NAMED object and
# pass it back explicitly.
dfx <- as.data.frame(dat0)
PCS <- paste0("PC", 1:10)
cat("\nPair matching: nearest neighbour, Mahalanobis over PC1..PC10,\n",
    "  exact blocks on Sex + age_band + SeqCenter + lpa_cluster, ratio = 1, no replacement\n")
m.out <- matchit(as.formula(paste("AD ~", paste(PCS, collapse = " + "))),
                 data     = dfx,
                 method   = "nearest",
                 distance = "mahalanobis",
                 exact    = ~ Sex + age_band + SeqCenter + lpa_cluster,
                 replace  = FALSE, ratio = 1)
matched <- match.data(m.out, data = dfx)
matched$stratum <- matched$subclass
setDT(matched)

# With ratio = 1 every subclass is one case + one control by construction; assert it rather
# than re-running RIGA's good_strata filter, which would silently mask a failed match.
chk <- matched[, .(nc = sum(AD == 1), nk = sum(AD == 0)), by = stratum]
stopifnot(all(chk$nc == 1L), all(chk$nk == 1L))
cat(sprintf("Matched: %d subjects = %d pairs (cases %d / controls %d, ratio %.2f)\n",
            nrow(matched), nrow(chk), sum(matched$AD == 1), sum(matched$AD == 0),
            sum(matched$AD == 0) / sum(matched$AD == 1)))
cat(sprintf("  low-confidence subjects retained (would fail RIGA's filter): %d (%.1f%%)\n",
            sum(matched$cluster_prob < 0.7 | matched$entropy_norm > 0.5),
            100 * mean(matched$cluster_prob < 0.7 | matched$entropy_norm > 0.5)))
cat(sprintf("  LPA clusters represented: %d (RIGA's filtered cohort has 9)\n",
            uniqueN(matched$lpa_cluster)))

write.csv(matched, "${ADSP_ROOT}/output/matched/matched_samples_lpa_pc.nofilter_pair.csv",
          row.names = FALSE)

## ---- sparse GRM aligned to the matched sample -----------------------------
sgrm <- get(load("${ADSP_ROOT}/processed/output.sparseGRM.sGRM.RData"))
grm_iids <- sub("^[^_]+_", "", sgrm@Dimnames[[1]])
sgrm@Dimnames[[1]] <- grm_iids; sgrm@Dimnames[[2]] <- grm_iids
matched_ids <- as.character(matched$IID)
stopifnot(all(matched_ids %in% rownames(sgrm)))
matched_cohort <- sgrm[matched_ids, matched_ids]
rm(sgrm); gc(verbose = FALSE)

## ---- two null models, ONE cohort ------------------------------------------
fit_and_save <- function(fmla, path, tag) {
  cat("\n>>>", tag, ":", deparse(fmla, width.cutoff = 500), "\n")
  obj_null <- fit_nullmodel(fmla, data = as.data.frame(matched), id = "IID",
                            kins = matched_cohort, use_sparse = TRUE,
                            family = binomial(link = "logit"), verbose = TRUE,
                            use_SPA = TRUE)
  save(obj_null, file = path)
  cat("saved", path, " (N =", length(obj_null$id_include),
      ", design cols =", ncol(obj_null$X), ")\n")
}
fit_and_save(as.formula(paste("AD ~ Age + cluster_prob + entropy_norm +",
                              paste(PCS, collapse = " + "))),
             "${ADSP_ROOT}/output/matched/obj_null_stratum_all_matched_lpa_pc.nofilter_pair14.RData",
             "pair14 (RIGA's formula)")
fit_and_save(as.formula(paste("AD ~ Age + Sex + SeqCenter + lpa_cluster + cluster_prob + entropy_norm +",
                              paste(PCS, collapse = " + "))),
             "${ADSP_ROOT}/output/matched/obj_null_stratum_all_matched_lpa_pc.nofilter_pair33.RData",
             "pair33 (noRIGA's formula)")
cat("\nDONE\n")
