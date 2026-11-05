#!/usr/bin/env Rscript
# Extract existing Figure 4 results once; no association analysis is performed.
# Rscript figure4_extract_gene_set_results.R [repository-root] [cache-directory]
args <- commandArgs(trailingOnly=TRUE)
root <- if (length(args)) args[1] else "${ADSP_ROOT}"
cache <- if (length(args)>1) args[2] else "/tmp/fig4_ng_20261014"
dir.create(cache, recursive=TRUE, showWarnings=FALSE)
suppressPackageStartupMessages(library(data.table))
src <- file.path(root,"output/matched_lpa_pc_pair33_projcell/individual_analysis_qc7/results_individual_analysis_genome.Rdata")
cached <- file.path(cache,c("sv_chr.i32","sv_pos.i32","sv_logp.f32","sv_qq_logp.f32","sv_stats.csv","chromosomes.csv"))
if (!all(file.exists(cached)) || any(file.info(cached)$mtime < file.info(src)$mtime)) {
cat("Reading existing single-variant results:",src,"\n"); flush.console()
e <- new.env(); nm <- load(src,envir=e); g <- e[[nm[1]]]
stopifnot(all(c("CHR","POS","pvalue") %in% names(g)))
p <- as.numeric(g$pvalue)
stopifnot(all(is.finite(p)),all(p>0),all(p<=1))
n <- length(p)
writeBin(as.integer(g$CHR),file.path(cache,"sv_chr.i32"),size=4,endian="little")
writeBin(as.integer(g$POS),file.path(cache,"sv_pos.i32"),size=4,endian="little")
writeBin(-log10(p),file.path(cache,"sv_logp.f32"),size=4,endian="little")
chr <- data.table(chr=as.integer(g$CHR),pos=as.numeric(g$POS))[,.(min=min(pos),max=max(pos)),by=chr][order(chr)]
fwrite(chr,file.path(cache,"chromosomes.csv"))
rm(g,e); gc()
ps <- sort(p); rm(p); gc()
k <- max(c(0L,which(ps <= .05*seq_along(ps)/n)))
lambda <- qchisq(median(ps),1,lower.tail=FALSE)/qchisq(.5,1)
writeBin(-log10(ps),file.path(cache,"sv_qq_logp.f32"),size=4,endian="little")
fwrite(data.table(n=n,k=k,critical=.05*k/n,lambda=lambda,min_p=min(ps),max_logp=-log10(min(ps)),source=src),file.path(cache,"sv_stats.csv"))
cat("Single variants:",n,"tests;",k,"FDR discoveries\n"); flush.console()
rm(ps); gc()
} else cat("Using extraction cache newer than its source\n")
# The exact coordinate table used by the original figure, without loading its dependencies.
pk <- "${R_LIBS_ROOT}/x86_64-pc-linux-gnu-library/4.4/STAARpipelineSummary/R/sysdata"
ge <- new.env(); lazyLoad(pk,ge)
fwrite(ge$genes,file.path(cache,"gene_coordinates.csv"))
for(fam in c("coding","noncoding")) {
  f <- file.path(root,"revision/tables/review",paste0("genesets_primary_",fam,".csv.gz"))
  d <- fread(cmd=paste("gzip -cd",shQuote(f)),select=c("Gene name","Chr","Category","#SNV","cMAC","tested","SB","q_SB"))
  d <- d[tested==TRUE]
  stopifnot(all(is.finite(d$SB)),all(d$SB>0),all(d$SB<=1))
  fwrite(d,file.path(cache,paste0(fam,".csv")))
  cat(fam,":",nrow(d),"tested sets;",sum(d$q_SB<.05),"FDR discoveries\n")
}
cat("Extraction complete:",cache,"\n")
