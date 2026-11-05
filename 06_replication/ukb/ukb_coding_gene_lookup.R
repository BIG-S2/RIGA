#!/usr/bin/env Rscript
## ukb_coding_gene_lookup.R (2026-10-19): copy of revision/scripts/ukb_coding_lookup_reconstructed.R for the root-checked STAAR/STAARpipeline library
## (ukb_rootcheck_20261019; library selected by R_LIBS_USER in jobs/05_downstream.sbatch or jobs/06_tables.sbatch).
## The only code changes are the output paths/names marked '2026-10-19 rootcheck' below; nothing existing is overwritten.
## ukb_coding_lookup_reconstructed.R (2026-10-13; corrected full cohort, saddlepoint burden, from gw_coding_spa/reconstructed) -- adapted from ukb_coding_lookup_projcell.R (2026-10-06) -- UK Biobank (WES, 152,351 participants, 1,753 AD cases) gene-based
## coding results for the genes named in the stratified-design manuscript: every mask, both arms
## (noRIGA = all UKB subjects, RIGA = the RIGA-matched set), saddlepoint burden p (STAAR-B) from
## revision/ukb/gw_coding_spa/<arm>/UKB_AD_Coding_*.Rdata. Writes revision/tables/supp_v3/S10_ukb_coding_lookup.csv.
suppressPackageStartupMessages(library(data.table))
ROOT <- "${ADSP_ROOT}"; UKB <- file.path(ROOT, "revision/ukb/gw_coding_spa"); TAB <- file.path(ROOT, "revision/tables/supp_v3")
GENES <- c("NECTIN2", "TREM2", "PSEN1", "TOMM40", "APOE", "PGAM5", "MMP17", "SHARPIN", "BIN1", "ABCA7", "SORL1")
one_dt <- function(m) {
  if (is.null(m) || !length(m)) return(NULL)
  x <- as.data.table(as.data.frame(m, stringsAsFactors = FALSE, check.names = FALSE))
  for (cc in names(x)) if (is.list(x[[cc]])) set(x, j = cc, value = unlist(x[[cc]], use.names = FALSE))
  for (cc in names(x)) if (!is.character(x[[cc]])) set(x, j = cc, value = as.character(x[[cc]]))
  x }
to_dt <- function(obj) {
  if (is.list(obj) && !is.data.frame(obj) && !is.matrix(obj))
    return(rbindlist(lapply(Filter(function(z) is.matrix(z) || is.data.frame(z), obj), one_dt), fill = TRUE, use.names = TRUE))
  one_dt(obj) }
out <- list(); n_genes <- list()
for (arm in c("rootcheck")) {   # 2026-10-19 rootcheck: gw_coding_spa/rootcheck
  files <- list.files(file.path(UKB, arm), pattern = "^UKB_AD_Coding_\\d+\\.Rdata$", full.names = TRUE)
  cat(arm, ":", length(files), "chunks\n"); tot <- 0L
  for (f in files) {
    e <- new.env(); nm <- load(f, envir = e)
    for (n in nm) { d <- to_dt(get(n, envir = e)); if (is.null(d) || !"Gene name" %in% names(d)) next
      tot <- tot + nrow(d)
      h <- d[`Gene name` %in% GENES]; if (nrow(h)) { h[, arm := arm]; out[[length(out) + 1]] <- h } } }
  n_genes[[arm]] <- tot
}
R <- rbindlist(out, fill = TRUE, use.names = TRUE)
pcol <- if ("STAAR-B" %in% names(R)) "STAAR-B" else grep("^Burden", names(R), value = TRUE)[1]
keep <- intersect(c("arm", "Gene name", "Chr", "Category", "#SNV", "cMAC", "STAAR-B", "Burden(1,25)", "Burden(1,1)", "STAAR-O"), names(R))
R <- R[, ..keep]; setorder(R, arm, `Gene name`, Category)
fwrite(R, file.path(TAB, "S10_ukb_coding_lookup_rootcheck.csv"))   # 2026-10-19 rootcheck: new file
cat("tests per arm:", paste(names(n_genes), unlist(n_genes), collapse = "; "), "\n"); print(R[, c("arm", "Gene name", "Category", "#SNV", "cMAC", pcol), with = FALSE])
