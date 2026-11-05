#!/usr/bin/env Rscript
## gene_set_primary_statistic.R (2026-10-13): the primary gene-set results of the stratified design from the deployed STAAR
## outputs (normal approximation). Primary statistic SB = Cauchy combination of every SKAT and burden p value (both weight
## families, unweighted and annotation-weighted), i.e. STAAR-O without its ACAT-V component. Sets with cMAC < 10 are
## removed before testing; Benjamini-Hochberg FDR 5% within each family (coding, noncoding, ncRNA). STAAR-O and the
## component summaries are kept as secondary columns. Writes revision/tables/review/genesets_primary_<family>.csv.gz and
## genesets_primary_fdr_sets.csv (all families) plus genesets_primary_acatv_only.csv (deployed STAAR-O discoveries lost).
suppressPackageStartupMessages(library(data.table))
D <- "${ADSP_ROOT}/output/matched_lpa_pc_pair33_projcell"; TAB <- "${ADSP_ROOT}/revision/tables/review"
PAT <- c(coding = "^ADSP_AD_Coding_\\d+\\.Rdata$", noncoding = "^ADSP_AD_Noncoding_\\d+\\.Rdata$", ncRNA = "^ADSP_AD_ncRNA_\\d+\\.Rdata$")
flat <- function(x) { if (is.null(x) || !length(x)) return(NULL); if (is.matrix(x) || is.data.frame(x)) x <- list(x)
  rbindlist(lapply(x, function(m) { if (is.null(m) || !length(m)) return(NULL); d <- as.data.frame(m, stringsAsFactors = FALSE)
    as.data.table(lapply(d, function(v) unlist(v))) }), fill = TRUE) }
cct <- function(P) { P <- pmin(pmax(P, 1e-300), 1 - 1e-16); 0.5 - atan(rowMeans(tan((0.5 - P) * pi), na.rm = TRUE)) / pi }
sig <- list(); lost <- list()
for (fam in names(PAT)) {
  fs <- list.files(file.path(D, fam), pattern = PAT[[fam]], full.names = TRUE)
  o <- rbindlist(lapply(fs, function(f) { e <- new.env(); load(f, envir = e); flat(get(ls(e)[1], envir = e)) }), fill = TRUE)
  num <- setdiff(names(o), c("Gene name", "Chr", "Category")); for (cc in num) set(o, j = cc, value = as.numeric(o[[cc]]))
  sb <- grep("^(SKAT|Burden)\\(1,(25|1)\\)", names(o), value = TRUE); o[, SB := cct(as.matrix(o[, ..sb]))]
  o <- unique(o[is.finite(SB)], by = c("Gene name", "Category"))
  o[, tested := cMAC >= 10]; o[, q_SB := NA_real_]; o[tested == TRUE, q_SB := p.adjust(SB, "BH")]
  o[, q_STAARO_deployed := p.adjust(`STAAR-O`, "BH")]                    # the deployed rule: STAAR-O, no floor
  setcolorder(o, c("Gene name", "Chr", "Category", "#SNV", "cMAC", "tested", "SB", "q_SB", "STAAR-O", "q_STAARO_deployed",
                   intersect(c("STAAR-S(1,25)", "STAAR-B(1,25)", "STAAR-A(1,25)", "STAAR-S(1,1)", "STAAR-B(1,1)", "STAAR-A(1,1)", "ACAT-O"), names(o))))
  cat(sprintf("%s: %d sets, %d with cMAC >= 10 (%d combined p values); FDR 5%% primary %d; deployed STAAR-O %d\n", fam, nrow(o), sum(o$tested), length(sb),
              sum(o$q_SB < 0.05, na.rm = TRUE), sum(o$q_STAARO_deployed < 0.05, na.rm = TRUE)))
  fwrite(o, file.path(TAB, sprintf("genesets_primary_%s.csv.gz", fam)))
  s <- o[q_SB < 0.05 & !is.na(q_SB)][, family := fam]; setorder(s, SB); sig[[fam]] <- s
  l <- o[q_STAARO_deployed < 0.05 & (is.na(q_SB) | q_SB >= 0.05)][, family := fam]; lost[[fam]] <- l
  print(s[, .(`Gene name`, Category, `#SNV`, cMAC = round(cMAC), SB = signif(SB, 3), q_SB = signif(q_SB, 3), `STAAR-O` = signif(`STAAR-O`, 3))])
  if (nrow(l)) { cat("  deployed STAAR-O discoveries not significant under the primary test:\n"); print(l[, .(`Gene name`, Category, cMAC = round(cMAC), SB = signif(SB, 3), `STAAR-O` = signif(`STAAR-O`, 3))]) } }
fwrite(rbindlist(sig, fill = TRUE), file.path(TAB, "genesets_primary_fdr_sets.csv")); fwrite(rbindlist(lost, fill = TRUE), file.path(TAB, "genesets_primary_acatv_only.csv")); cat("DONE\n")
