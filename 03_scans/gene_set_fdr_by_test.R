#!/usr/bin/env Rscript
## gene_set_fdr_by_test.R (2026-10-13; FAMILY=coding/noncoding/ncRNA; v2 adds SB = STAAR-O without ACAT-V, i.e. the Cauchy
## combination of every SKAT and burden p value, as in Wang and Heard-Costa 2026): on the deployed ADSP gene-set results (normal approximation, coding/) and the
## saddlepoint rescan (coding_spa/), apply the paper's rule (BH FDR 5% within the coding family) to each candidate primary
## test, with and without a cumulative-MAC floor of 10, and list the discoveries. Read-only. Writes
## revision/tables/review/coding_fdr_by_test.csv.
suppressPackageStartupMessages(library(data.table))
D <- "${ADSP_ROOT}/output/matched_lpa_pc_pair33_projcell"; TAB <- "${ADSP_ROOT}/revision/tables/review"
flat <- function(x) { if (is.null(x) || !length(x)) return(NULL); if (is.matrix(x) || is.data.frame(x)) x <- list(x)
  rbindlist(lapply(x, function(m) { if (is.null(m) || !length(m)) return(NULL); d <- as.data.frame(m, stringsAsFactors = FALSE)
    as.data.table(lapply(d, function(v) unlist(v))) }), fill = TRUE) }
FAM <- Sys.getenv("FAMILY", "coding"); PAT <- c(coding = "^ADSP_AD_Coding_\\d+\\.Rdata$", noncoding = "^ADSP_AD_Noncoding_\\d+\\.Rdata$", ncRNA = "^ADSP_AD_ncRNA_\\d+\\.Rdata$")[[FAM]]
read_dir <- function(sub) { fs <- list.files(file.path(D, sub), pattern = PAT, full.names = TRUE)
  r <- rbindlist(lapply(fs, function(f) { e <- new.env(); load(f, envir = e); flat(get(ls(e)[1], envir = e)) }), fill = TRUE)
  setnames(r, make.names(names(r), unique = TRUE)); r }
o <- read_dir(FAM); s <- read_dir(paste0(FAM, "_spa"))
for (cc in grep("^Burden\\.|^STAAR\\.", names(o), value = TRUE)) set(o, j = cc, value = as.numeric(o[[cc]]))
cct <- function(P) { P <- pmin(pmax(P, 1e-300), 1 - 1e-16); T <- rowMeans(tan((0.5 - P) * pi), na.rm = TRUE); 0.5 - atan(T) / pi }   # Cauchy combination, as STAAR-B does
bcols <- grep("^Burden\\.1\\.(25|1)\\.", names(o), value = TRUE); cat("burden columns combined:", length(bcols), "\n")
sbcols <- grep("^(SKAT|Burden)\\.1\\.(25|1)\\.", names(o), value = TRUE); cat("SKAT+burden columns combined:", length(sbcols), "\n")
o[, SB := cct(as.matrix(o[, ..sbcols]))]
o[, B := cct(as.matrix(o[, ..bcols])[, , drop = FALSE] * 1)]
o <- unique(o[, .(gene = Gene.name, cat = Category, cMAC = as.numeric(cMAC), nSNV = as.numeric(X.SNV), O = as.numeric(STAAR.O), SB, B, B125 = as.numeric(STAAR.B.1.25.))], by = c("gene", "cat"))
s <- unique(s[, .(gene = Gene.name, cat = Category, SPAB = as.numeric(STAAR.B))], by = c("gene", "cat"))
m <- merge(o, s, by = c("gene", "cat"), all.x = TRUE)
cat(sprintf("%s sets %d; cMAC<10: %d (%.1f%%)\n", FAM, nrow(m), sum(m$cMAC < 10, na.rm = TRUE), 100 * mean(m$cMAC < 10, na.rm = TRUE)))
out <- list()
for (floor in c(0, 10)) for (t in c("O", "SB", "B", "B125", "SPAB")) {
  x <- m[cMAC >= floor & is.finite(get(t))]; q <- p.adjust(x[[t]], "BH"); hit <- x[q < 0.05]
  syn <- x[cat == "synonymous"]; nsyn <- sum(p.adjust(syn[[t]], "BH") < 0.05)
  out[[length(out) + 1]] <- data.table(test = t, cmac_floor = floor, n_tests = nrow(x), n_fdr = nrow(hit), n_syn_fdr = nsyn,
    hits = paste(sprintf("%s/%s(cMAC %.0f, p %.2g)", hit$gene, hit$cat, hit$cMAC, hit[[t]]), collapse = "; "))
  cat(sprintf("%-5s floor %2d: %6d tests, FDR5%% %d (synonymous %d): %s\n", t, floor, nrow(x), nrow(hit), nsyn, out[[length(out)]]$hits)) }
fwrite(rbindlist(out), file.path(TAB, sprintf("%s_fdr_by_test_v2.csv", FAM))); cat("DONE\n")
