#!/usr/bin/env Rscript
## export_gene_centric_results.R -- flatten the D_S gene-centric results to CSV for plot_gene_centric_summary_projcell.py:
##   noncoding (379 chunks, 7 categories) and ncRNA (222 chunks) -> revision/tables/gene_centric_noncoding_projcell_all.csv
##   SKAT-MC (merged windows) with gene positions            -> revision/tables/skatmc_projcell_all.csv
## Columns: gene, category, chr, gene_start (hg38, STAARpipeline::genes_info or the SKAT-MC gene index), n_variants, cmac, p.
suppressPackageStartupMessages({ library(data.table); library(STAARpipeline) })
OUT <- "${ADSP_ROOT}/output"; TAB <- "${ADSP_ROOT}/revision/tables"
ARM <- file.path(OUT, "matched_lpa_pc_pair33_projcell")
flat <- function(d, pat) {
  fs <- list.files(d, pat, full.names = TRUE); cat(d, ":", length(fs), "chunks\n")
  keep <- c("Gene name", "Chr", "Category", "#SNV", "cMAC", "STAAR-O")
  rows <- lapply(fs, function(f) { e <- new.env(); load(f, envir = e); obj <- get(ls(e)[1], envir = e)
    mats <- if (is.matrix(obj)) list(obj) else Filter(function(z) is.matrix(z) && nrow(z) >= 1, obj)   # noncoding: list of 1-row matrices; ncRNA: one genes x columns matrix
    rbindlist(lapply(mats, function(z) { h <- intersect(keep, colnames(z))
      as.data.table(matrix(as.character(z[, h, drop = FALSE]), nrow = nrow(z), dimnames = list(NULL, h))) }), fill = TRUE) })   # the matrices are list-mode
  r <- rbindlist(rows, fill = TRUE)
  setnames(r, c("Gene name", "Chr", "Category", "#SNV", "cMAC", "STAAR-O"), c("gene", "chr", "category", "n_variants", "cmac", "p"), skip_absent = TRUE)
  if (!"category" %in% names(r)) r[, category := "ncRNA"]
  for (cc in c("chr", "n_variants", "cmac")) if (!cc %in% names(r)) r[, (cc) := NA]      # the ncRNA matrices carry no Chr column
  r[, `:=`(chr = suppressWarnings(as.integer(chr)), n_variants = suppressWarnings(as.integer(n_variants)), cmac = suppressWarnings(as.numeric(cmac)), p = suppressWarnings(as.numeric(p)))]
  unique(r[is.finite(p) & p > 0 & p <= 1], by = c("gene", "category"))
}
gi <- as.data.table(genes_info)
setnames(gi, c("hgnc_symbol", "chromosome_name", "start_position"), c("gene", "chr", "gene_start"))
gi <- unique(gi[chr %in% as.character(1:22), .(gene, chr = as.integer(chr), gene_start = as.numeric(gene_start))], by = "gene")

nc <- flat(file.path(ARM, "noncoding"), "^ADSP_AD_Noncoding_\\d+\\.Rdata$")
nr <- flat(file.path(ARM, "ncRNA"), "^ADSP_AD_ncRNA_\\d+\\.Rdata$")
g <- rbindlist(list(nc, nr), fill = TRUE)
g <- merge(g, gi[, .(gene, chr_gi = chr, gene_start)], by = "gene", all.x = TRUE)
g[, chr := fifelse(is.na(chr), chr_gi, chr)]; g[, chr_gi := NULL]
## genes absent from genes_info (the ncRNA symbols): Ensembl 86 / GRCh38 gene starts by symbol
miss <- unique(g[is.na(gene_start), gene])
if (length(miss)) {
  suppressPackageStartupMessages(library(EnsDb.Hsapiens.v86))
  ens <- as.data.table(ensembldb::genes(EnsDb.Hsapiens.v86, filter = AnnotationFilter::GeneNameFilter(miss),
                                        columns = c("gene_name", "seq_name", "gene_seq_start"), return.type = "data.frame"))
  ens <- ens[seq_name %in% as.character(1:22), .(gene = gene_name, chr_e = as.integer(seq_name), start_e = as.numeric(gene_seq_start))]
  ens <- unique(ens[order(gene, start_e)], by = "gene")
  g <- merge(g, ens, by = "gene", all.x = TRUE)
  g[is.na(gene_start), `:=`(gene_start = start_e, chr = fifelse(is.na(chr), chr_e, chr))]; g[, c("chr_e", "start_e") := NULL]
  cat(sprintf("EnsDb lookup: %d of %d missing symbols placed\n", sum(miss %in% ens$gene), length(miss)))
}
cat(sprintf("gene-centric rows %d (%d without a position)\n", nrow(g), sum(is.na(g$gene_start))))
fwrite(g[order(category, chr, gene_start)], file.path(TAB, "gene_centric_noncoding_projcell_all.csv"))

sk <- fread(file.path(ARM, "skatmc_by_window/skatmc_all_merged.tsv.gz"))
sk <- sk[is.finite(p_value) & p_value > 0 & p_value <= 1, .(gene = set_id, category = "SKAT-MC", chr = as.integer(chr), n_variants = n_var_test, p = p_value)]
idx <- rbindlist(lapply(list.files(file.path(OUT, "skatmc_gene_meta_rare"), "^gene_index_chr[0-9]+_maf01\\.tsv\\.gz$", full.names = TRUE),
                        function(f) { d <- fread(f); d[, chr := as.integer(sub(".*gene_index_chr([0-9]+)_.*", "\\1", basename(f)))]; d[, .(gene = set_id, chr, gene_start = as.numeric(gene_min_pos))] }))
sk <- merge(sk, idx, by = c("gene", "chr"), all.x = TRUE)
cat(sprintf("SKAT-MC rows %d (%d without a position)\n", nrow(sk), sum(is.na(sk$gene_start))))
fwrite(sk[order(chr, gene_start)], file.path(TAB, "skatmc_projcell_all.csv"))
cat("DONE\n")
