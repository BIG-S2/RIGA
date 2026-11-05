# lpa_pc, POOLED (FULL) — same as adsp_pooled.lpa_pc.R but WITHOUT the
# cluster-confidence filter. Uses all 38,138 demographic-QC-passing subjects,
# including those whose LPA cluster assignment is uncertain.
#
# This pairs with adsp_pooled.lpa_pc.R (confident-only, ~28,919) to test
# whether the confidence filter materially changes results in the pooled setting.

suppressPackageStartupMessages({
  library(data.table); library(dplyr); library(readr); library(tidyr)
  library(gdsfmt); library(SeqArray); library(SeqVarTools)
  library(STAAR);  library(STAARpipeline)
})
set.seed(0)

cov_wgs <- read_csv("${ADSP_ROOT}/pheno/adsp_combined.csv",
                    show_col_types = FALSE) %>% drop_na(AD, Sex, Age)
cov_wgs$Age <- as.numeric(gsub("\\+", "", cov_wgs$Age))

# ALL subjects with LPA assignment, no confidence filter
clu <- fread("${ADSP_ROOT}/adsp_lpa_clusters_20PC.tsv")
setDT(clu)
K <- length(grep("^post_k[0-9]+$", names(clu)))
if (!"entropy_norm" %in% names(clu)) clu[, entropy_norm := entropy / log(K)]
cat(sprintf("LPA cluster file rows: %d (no confidence filter applied)\n", nrow(clu)))

seqmeta <- read_csv("pheno/gcad.r5.wgs.58507.2024.11.03.quality.metrics.ALL.csv",
                    show_col_types = FALSE) %>%
  transmute(IID = Sample_name, SeqCenter = as.character(SeqCenter))

dat <- merge(cov_wgs, clu, by = "IID") %>%
  left_join(seqmeta, by = "IID") %>%
  mutate(Sex = as.factor(Sex), AD = as.factor(AD),
         SeqCenter = as.character(SeqCenter))

dat$SeqCenter[is.na(dat$SeqCenter)] <- "Unknown"
dat$SeqCenter[dat$SeqCenter %in% c(
  "MEDGENOME","Unknown","Genentech","GENENTECH",
  "HudsonAlpha Institute for Biotechnology","Macrogen","MACROGEN")] <- "Other"
dat$SeqCenter <- droplevels(as.factor(dat$SeqCenter))
dat$lpa_cluster <- as.factor(dat$lpa_cluster)

dat <- dat %>% tidyr::drop_na(Age, Sex, AD, SeqCenter)
dat$entropy_norm <- if (K > 1) dat$entropy / log(K) else 0

cat("\nPooled (FULL) cohort size:", nrow(dat), "\n")
cat("Cases / controls:", sum(dat$AD == 1), "/", sum(dat$AD == 0), "\n")
cat("SeqCenter levels:\n"); print(table(dat$SeqCenter))
cat("lpa_cluster levels:\n"); print(table(dat$lpa_cluster))

sgrm <- get(load("${ADSP_ROOT}/processed/output.sparseGRM.sGRM.RData"))
grm_iids <- sub("^[^_]+_", "", sgrm@Dimnames[[1]])
sgrm@Dimnames[[1]] <- grm_iids; sgrm@Dimnames[[2]] <- grm_iids
ids <- as.character(dat$IID)
ids <- ids[ids %in% grm_iids]
dat <- dat[match(ids, dat$IID), ]
sgrm_sub <- sgrm[ids, ids]

pc_cov_cols <- paste0("PC", 1:10)
fmla <- as.formula(paste(
  "AD ~ Age + Sex + SeqCenter + lpa_cluster + cluster_prob + entropy_norm +",
  paste(pc_cov_cols, collapse = " + ")
))
cat("\nFitting null model:", deparse(fmla), "\n")
obj_null <- fit_nullmodel(
  fmla,
  data       = dat,
  id         = "IID",
  kins       = sgrm_sub,
  use_sparse = TRUE,
  family     = binomial(link = "logit"),
  verbose    = TRUE,
  use_SPA    = TRUE
)

save(obj_null,
     file = "${ADSP_ROOT}/output/obj_null_stratum_all_pooled_lpa_pc_full.RData")
cat("\nSaved: ${ADSP_ROOT}/output/obj_null_stratum_all_pooled_lpa_pc_full.RData\n")
