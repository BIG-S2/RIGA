# UKB 200k WGS replication — step 2: fit the STAAR null model.
#
# use_SPA = TRUE is MANDATORY here. The case:control ratio is roughly 1:86, and without
# saddlepoint correction the burden/score statistics are badly anticonservative in the tail —
# which is precisely the failure mode this manuscript is about.

suppressPackageStartupMessages({
  library(Matrix); library(STAAR); library(STAARpipeline)
})

REV <- "${ADSP_ROOT}/revision"
UKB <- file.path(REV, "ukb")

obj <- readRDS(file.path(UKB, "ukb_ad_analysis_set.rds"))
pheno <- obj$pheno; grm <- obj$grm
cat("n =", nrow(pheno), " cases =", sum(pheno$AD == 1),
    " controls =", sum(pheno$AD == 0), "\n")
cat("GRM dim:", paste(dim(grm), collapse = " x "),
    " nnz:", length(grm@x), "\n")

pheno$sex <- as.factor(pheno$sex)

fmla <- as.formula(paste("AD ~ age + age2 + sex +",
                         paste(paste0("PC", 1:10), collapse = " + ")))
cat("\nnull formula:", deparse(fmla), "\n")
cat("fitting with use_sparse=TRUE, use_SPA=TRUE ...\n")
t0 <- Sys.time()

obj_nullmodel <- fit_nullmodel(
  fmla,
  data       = pheno,
  id         = "IID",
  kins       = grm,
  use_sparse = TRUE,
  family     = binomial(link = "logit"),
  verbose    = TRUE,
  use_SPA    = TRUE
)

cat("\nelapsed:", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min\n")
cat("converged:", isTRUE(obj_nullmodel$converged), "\n")
cat("use_SPA in object:", isTRUE(obj_nullmodel$use_SPA), "\n")

save(obj_nullmodel, file = file.path(UKB, "ukb_obj_nullmodel.RData"))
cat("Saved", file.path(UKB, "ukb_obj_nullmodel.RData"), "\n")
cat("DONE\n")
