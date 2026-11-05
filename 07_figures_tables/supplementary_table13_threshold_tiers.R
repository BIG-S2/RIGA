## supplementary_table13_threshold_tiers.R (2026-09-29) -- what each rigorously defined threshold admits in the D_S (stratified) results:
## fixed STAARpipelineSummary alpha, Bonferroni over the tests actually performed (FWER 0.05), Benjamini-Hochberg FDR 0.05.
suppressPackageStartupMessages(library(data.table))
OUT <- "${ADSP_ROOT}/output/matched_lpa_pc_pair33_projcell"; TAB <- "${ADSP_ROOT}/revision/tables"
load_cat <- function(dir, cats) rbindlist(lapply(cats, function(cc) { x <- as.data.table(as.data.frame(get(load(file.path(dir, paste0(cc, ".Rdata")))), stringsAsFactors = FALSE))
  x[, .(gene = as.character(`Gene name`), chr = as.integer(Chr), category = cc, p = as.numeric(`STAAR-O`))] }))
tiers <- function(d, fam, fixed) { d <- d[is.finite(p) & p > 0]; m <- nrow(d); bonf <- 0.05 / m; d[, q := p.adjust(p, "BH")]
  bh_line <- if (any(d$q < 0.05)) max(d[q < 0.05]$p) else NA_real_
  cat(sprintf("\n== %s: %d tests; fixed %.2e -> %d; Bonferroni 0.05/n = %.2e -> %d; BH q<0.05 -> %d (largest admitted p %.2e)\n",
              fam, m, fixed, sum(d$p < fixed), bonf, sum(d$p < bonf), sum(d$q < 0.05), bh_line))
  print(d[q < 0.05][order(p), .(gene, chr, category, p = signif(p, 3), q = signif(q, 3))])
  data.table(family = fam, n_tests = m, fixed = fixed, n_fixed = sum(d$p < fixed), bonferroni = bonf, n_bonferroni = sum(d$p < bonf), n_fdr05 = sum(d$q < 0.05), bh_boundary = bh_line) }
r1 <- tiers(load_cat(file.path(OUT, "coding"), c("plof", "plof_ds", "missense", "disruptive_missense", "synonymous")), "coding", 5e-7)
r2 <- tiers(load_cat(file.path(OUT, "noncoding"), c("UTR", "upstream", "downstream", "promoter_CAGE", "promoter_DHS", "enhancer_CAGE", "enhancer_DHS")), "noncoding", 3.57e-7)
g <- as.data.table(get(load(file.path(OUT, "individual_analysis", "results_individual_analysis_genome.QC7.Rdata")))); g <- g[is.finite(pvalue) & pvalue > 0]
m <- nrow(g); g[, q := p.adjust(pvalue, "BH")]; bh <- if (any(g$q < 0.05)) max(g[q < 0.05]$pvalue) else NA_real_
cat(sprintf("\n== single variants: %d tests; 5e-8 -> %d; Bonferroni 0.05/n = %.2e -> %d; BH q<0.05 -> %d (largest p %.2e); 1e-6 -> %d; 1e-5 -> %d\n",
            m, sum(g$pvalue < 5e-8), 0.05 / m, sum(g$pvalue < 0.05 / m), sum(g$q < 0.05), bh, sum(g$pvalue < 1e-6), sum(g$pvalue < 1e-5)))
loci <- function(x) { x <- copy(x); x[, bin := paste0(CHR, ":", POS %/% 1e6, "Mb")]; x[order(pvalue)][!duplicated(bin), .(CHR, POS, pvalue = signif(pvalue, 3), bin)] }
cat("-- 1-Mb loci at BH q<0.05:\n"); print(loci(g[q < 0.05]), nrows = 60)
cat("-- 1-Mb loci at p < 1e-6:\n"); print(loci(g[pvalue < 1e-6]), nrows = 60)
r3 <- data.table(family = "single_variant", n_tests = m, fixed = 5e-8, n_fixed = sum(g$pvalue < 5e-8), bonferroni = 0.05 / m, n_bonferroni = sum(g$pvalue < 0.05 / m), n_fdr05 = sum(g$q < 0.05), bh_boundary = bh)
fwrite(rbindlist(list(r1, r2, r3)), file.path(TAB, "stratified_threshold_tiers.csv")); cat("DONE\n")
