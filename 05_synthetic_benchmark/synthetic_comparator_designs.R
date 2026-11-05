#!/usr/bin/env Rscript
# synthetic_comparator_designs.R -- GMMAT and GENESIS arms of the RIGA simulation benchmark.
#
# Invoked once per (scenario, config) by synthetic_comparator_designs.py::run_r_comparators.
# Both tools are fitted against the SAME empirical GCTA-style GRM written by
# Python (grm.f64, N x N float64, column-major == row-major since symmetric),
# and the SAME hard-call genotypes (geno.i8, M x N int8, C order).
#
# Outputs <which>_pvals.txt: one p-value per line, M lines, in variant order,
# NA where the tool dropped the variant.

suppressWarnings(suppressMessages({
  library(GMMAT)
  library(GENESIS)
  library(GWASTools)
}))

args <- commandArgs(trailingOnly = TRUE)
getarg <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (is.na(i)) return(default)
  args[i + 1L]
}

workdir    <- getarg("--workdir")
bfile      <- getarg("--bfile")
geno_path  <- getarg("--geno")
grm_path   <- getarg("--grm")
pheno_path <- getarg("--pheno")
M          <- as.integer(getarg("--m"))
N          <- as.integer(getarg("--n"))
block_size <- as.integer(getarg("--block-size", "2000"))
which      <- strsplit(getarg("--which", "gmmat,genesis"), ",")[[1]]

t_start <- proc.time()[3]

## ---------------------------------------------------------------- phenotype
ph <- read.csv(pheno_path, stringsAsFactors = FALSE)
stopifnot(nrow(ph) == N)
ph$id     <- as.character(ph$id)
ph$center <- factor(ph$center)
ph$Y      <- as.integer(ph$Y)

pc_cols <- grep("^PC[0-9]+$", names(ph), value = TRUE)
pc_cols <- pc_cols[order(as.integer(sub("PC", "", pc_cols)))]
covars  <- c("center", "age", "gender", pc_cols)
form    <- as.formula(paste("Y ~", paste(covars, collapse = " + ")))
cat("Covariates:", paste(covars, collapse = " + "), "\n")
cat("N =", N, " M =", M, " cases =", sum(ph$Y), " controls =", sum(ph$Y == 0), "\n")

## ---------------------------------------------------------------------- GRM
con <- file(grm_path, "rb")
G <- matrix(readBin(con, "double", n = as.numeric(N) * N, size = 8),
            nrow = N, ncol = N)
close(con)
rownames(G) <- colnames(G) <- ph$id
cat(sprintf("GRM: diag mean %.4f sd %.4f | offdiag mean %.5f sd %.5f\n",
            mean(diag(G)), sd(diag(G)),
            mean(G[lower.tri(G)]), sd(G[lower.tri(G)])))

write_pvals <- function(name, p) {
  stopifnot(length(p) == M)
  # Write non-finite entries as "nan", not R's "NA".  numpy.loadtxt on the
  # Python side cannot parse "NA" and raises, which run_comparators catches by
  # marking the tool failed — so a single variant dropped by a MAF or
  # missingness filter would silently void BOTH R comparators for that unit.
  s <- format(p, digits = 10, scientific = TRUE)
  s[!is.finite(p)] <- "nan"
  write.table(s, file.path(workdir, paste0(name, "_pvals.txt")),
              row.names = FALSE, col.names = FALSE, quote = FALSE)
}

## -------------------------------------------------------------------- GMMAT
if ("gmmat" %in% which) {
  t0 <- proc.time()[3]
  obj <- glmmkin(fixed = form, data = ph, kins = G, id = "id",
                 family = binomial(link = "logit"))
  cat(sprintf("GMMAT glmmkin: theta = %s  (%.1fs)\n",
              paste(sprintf("%.5g", obj$theta), collapse = ", "),
              proc.time()[3] - t0))
  outf <- file.path(workdir, "gmmat_score.txt")
  glmm.score(obj, infile = bfile, outfile = outf,
             MAF.range = c(1e-7, 0.5), miss.cutoff = 1, nperbatch = 2000)
  sc <- read.table(outf, header = TRUE, stringsAsFactors = FALSE)
  p <- rep(NA_real_, M)
  idx <- as.integer(sub("^snp", "", sc$SNP)) + 1L
  keep <- !is.na(idx) & idx >= 1L & idx <= M
  p[idx[keep]] <- sc$PVAL[keep]
  write_pvals("gmmat", p)
  cat(sprintf("GMMAT done: %d/%d variants, lambda_GC = %.4f (%.1fs)\n",
              sum(!is.na(p)), M,
              median(qchisq(p[!is.na(p)], df = 1, lower.tail = FALSE),
                     na.rm = TRUE) / qchisq(0.5, 1),
              proc.time()[3] - t0))
}

## ------------------------------------------------------------------ GENESIS
if ("genesis" %in% which) {
  t0 <- proc.time()[3]
  con <- file(geno_path, "rb")
  gv <- readBin(con, "integer", n = as.numeric(M) * N, size = 1, signed = TRUE)
  close(con)
  # Python wrote C order (M x N, variant-major); R fills column-major, so read
  # into an (N x M) matrix and transpose to get variants in rows.
  geno <- matrix(as.integer(gv), nrow = N, ncol = M)
  geno <- t(geno)
  rm(gv); gc(verbose = FALSE)
  stopifnot(dim(geno) == c(M, N))

  snpID <- seq_len(M)
  chrom <- as.integer((snpID - 1L) %/% block_size + 1L)
  pos   <- as.integer((snpID - 1L) %% block_size + 1L)

  mgr <- MatrixGenotypeReader(genotype = geno, snpID = snpID,
                              chromosome = chrom, position = pos,
                              scanID = seq_len(N))
  scanAnnot <- ScanAnnotationDataFrame(
    data.frame(scanID = seq_len(N), ph[, c("Y", "center", "age", "gender", pc_cols)],
               stringsAsFactors = FALSE))
  gdata <- GenotypeData(mgr, scanAnnot = scanAnnot)

  Gg <- G
  rownames(Gg) <- colnames(Gg) <- as.character(seq_len(N))
  nullmod <- fitNullModel(scanAnnot, outcome = "Y", covars = covars,
                          cov.mat = Gg, family = "binomial", verbose = FALSE)
  cat(sprintf("GENESIS fitNullModel: varComp = %s converged = %s (%.1fs)\n",
              paste(sprintf("%.5g", nullmod$varComp), collapse = ", "),
              nullmod$converged, proc.time()[3] - t0))

  it <- GenotypeBlockIterator(gdata, snpBlock = 5000)
  assoc <- assocTestSingle(it, null.model = nullmod, test = "Score",
                           verbose = FALSE)
  p <- rep(NA_real_, M)
  idx <- as.integer(assoc$variant.id)
  keep <- !is.na(idx) & idx >= 1L & idx <= M
  p[idx[keep]] <- assoc$Score.pval[keep]
  write_pvals("genesis", p)
  cat(sprintf("GENESIS done: %d/%d variants, lambda_GC = %.4f (%.1fs)\n",
              sum(!is.na(p)), M,
              median(qchisq(p[!is.na(p)], df = 1, lower.tail = FALSE),
                     na.rm = TRUE) / qchisq(0.5, 1),
              proc.time()[3] - t0))
}

cat(sprintf("synthetic_comparator_designs.R total %.1fs\n", proc.time()[3] - t_start))
