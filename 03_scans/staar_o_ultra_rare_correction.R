#############################################################################
# Post-hoc Wang/Heard-Costa (HGG Advances 2026, PMID 41560352) fix +
# STAAR-O reporting for the gene-centric CODING results.
#
# Operates on ALREADY-COMPUTED per-gene result matrices (no GDS / null model
# needed) -- purely recombines the stored component p-values.
#
# Two corrections:
#   (1) Drop genes with cMAC < 10  (ultra-rare regime where the STAAR-O ACAT-V
#       component degenerates to the burden test and double-counts it).
#   (2) Recompute STAAR-O WITHOUT the ACAT-V component, i.e. the Cauchy/ACAT
#       combination of the SKAT + Burden omnibus pieces only
#       (STAAR-S(1,25),STAAR-S(1,1),STAAR-B(1,25),STAAR-B(1,1)).
#
# For SPA dirs (use_SPA=TRUE) the stored output is burden-only (STAAR-B); there
# is no STAAR-O/ACAT-V to remove, so only correction (1) is applied to STAAR-B.
#
# Non-destructive: writes new *.WHCfix.csv / *.WHCfix.png alongside the originals.
#############################################################################
suppressPackageStartupMessages({ library(STAAR); library(data.table) })

alpha     <- 5e-7
cmac_min  <- 10
cats      <- c("plof","plof_ds","missense","disruptive_missense","synonymous")
results_name <- "ADSP_AD_Coding"
njobs     <- 379

## STAAR-O omnibus components (per the STAAR package construction)
comp_all <- c("STAAR-S(1,25)","STAAR-S(1,1)",
              "STAAR-B(1,25)","STAAR-B(1,1)",
              "STAAR-A(1,25)","STAAR-A(1,1)")            # full STAAR-O
comp_noA <- c("STAAR-S(1,25)","STAAR-S(1,1)",
              "STAAR-B(1,25)","STAAR-B(1,1)")            # STAAR-O minus ACAT-V

num <- function(x) suppressWarnings(as.numeric(x))

## read every per-gene 1-row matrix from a results dir into one data.frame
load_dir <- function(dir){
  files <- list.files(dir, pattern=paste0("^",results_name,"_[0-9]+\\.Rdata$"), full.names=TRUE)
  rows <- list(); k <- 0
  for(f in files){
    e <- new.env(); ok <- tryCatch({ load(f, envir=e); TRUE }, error=function(...) FALSE)
    if(!ok) next
    obj <- get(ls(e)[1], envir=e)
    if(!is.list(obj)) next
    for(el in obj){
      if(is.null(el) || !is.matrix(el) || nrow(el) < 1) next
      cn <- colnames(el); if(is.null(cn)) next
      v  <- as.list(el[1, , drop=TRUE]); names(v) <- cn
      k <- k+1; rows[[k]] <- v
    }
  }
  if(k==0) return(NULL)
  as.data.frame(rbindlist(rows, fill=TRUE), stringsAsFactors=FALSE, check.names=FALSE)
}

## row-wise Cauchy combination of a set of p-value columns
cct_rows <- function(df, cols){
  cols <- cols[cols %in% names(df)]
  P <- sapply(cols, function(c) num(df[[c]]))
  if(is.null(dim(P))) P <- matrix(P, nrow=nrow(df))
  apply(P, 1, function(p){ p <- p[is.finite(p)]; if(!length(p)) return(NA_real_); CCT(p) })
}

process_dir <- function(dir, tag, log){
  cat("\n=================  ", tag, "  =================\n", file=log, append=TRUE)
  cat("dir:", dir, "\n", file=log, append=TRUE)
  df <- load_dir(dir)
  if(is.null(df)){ cat("  (no results)\n", file=log, append=TRUE); return(invisible()) }

  df$Category <- df[["Category"]]
  df$cMAC     <- num(df[["cMAC"]])
  df$Chr      <- df[["Chr"]]
  df$Gene     <- df[["Gene name"]]
  df$nSNV     <- num(df[["#SNV"]])
  is_spa      <- !("STAAR-O" %in% names(df))

  if(!is_spa){
    ## sanity: does CCT over the 6 omnibus pieces reproduce stored STAAR-O?
    df$STAARO_stored <- num(df[["STAAR-O"]])
    df$STAARO_recomp <- cct_rows(df, comp_all)
    df$STAARO_noACATV<- cct_rows(df, comp_noA)
    ok <- is.finite(df$STAARO_stored) & is.finite(df$STAARO_recomp) & df$STAARO_stored>0 & df$STAARO_recomp>0
    d  <- abs(log10(df$STAARO_recomp[ok]) - log10(df$STAARO_stored[ok]))
    cat(sprintf("  STAAR-O recompute check (n=%d): median|dlog10|=%.4g  max=%.4g\n",
                sum(ok), median(d), max(d)), file=log, append=TRUE)
    report_p <- "STAARO_noACATV"; base_p <- "STAARO_stored"
  } else {
    df$STAARB_stored <- num(df[["STAAR-B"]])
    report_p <- "STAARB_stored"; base_p <- "STAARB_stored"
    cat("  SPA dir: burden-only output; no STAAR-O/ACAT-V. Applying cMAC>=10 filter to STAAR-B.\n",
        file=log, append=TRUE)
  }

  ## ---- before/after counts per category --------------------------------
  cat(sprintf("\n  %-22s %12s %14s %14s\n","category","n_genes","sig_BEFORE","sig_AFTER_fix"),
      file=log, append=TRUE)
  keep_all <- df$Category %in% cats
  out <- df[keep_all,]
  out$kept <- out$cMAC >= cmac_min
  for(ct in cats){
    sub <- out[out$Category==ct,]
    nb  <- sum(num(sub[[base_p]])   < alpha, na.rm=TRUE)                       # before: stored, all cMAC
    na_ <- sum(sub$kept & num(sub[[report_p]]) < alpha, na.rm=TRUE)           # after: corrected metric, cMAC>=10
    cat(sprintf("  %-22s %12d %14d %14d\n", ct, nrow(sub), nb, na_),
        file=log, append=TRUE)
  }

  ## ---- write corrected significant tables ------------------------------
  keepcols <- c("Gene","Chr","Category","nSNV","cMAC")
  if(!is_spa) keepcols <- c(keepcols, "STAARO_stored","STAARO_noACATV")
  else        keepcols <- c(keepcols, "STAARB_stored")
  sig <- out[out$kept & num(out[[report_p]]) < alpha, ]
  sig <- sig[order(num(sig[[report_p]])), keepcols]
  fn  <- file.path(dir, "coding_sig.WHCfix.csv")
  write.csv(sig, fn, row.names=FALSE)
  cat(sprintf("\n  wrote %s  (%d sig genes, cMAC>=%d, alpha=%g)\n",
              fn, nrow(sig), cmac_min, alpha), file=log, append=TRUE)

  ## ---- QQ plot of the reported metric, by category ---------------------
  png(file.path(dir, "gene_centric_coding_qqplot.WHCfix.png"), width=800, height=800, res=110)
  par(mar=c(4.5,4.5,2,1))
  syms <- c(plof=0, plof_ds=1, missense=2, disruptive_missense=3, synonymous=4)
  plotted <- FALSE
  maxlp <- 1
  for(ct in cats){
    p <- num(out[out$Category==ct & out$kept, report_p]); p <- p[is.finite(p) & p>0]
    if(length(p)) maxlp <- max(maxlp, -log10(min(p)))
  }
  for(ct in cats){
    p <- num(out[out$Category==ct & out$kept, report_p]); p <- p[is.finite(p) & p>0]
    if(!length(p)) next
    o <- -log10(sort(p)); ex <- -log10(ppoints(length(p)))
    if(!plotted){
      plot(ex, o, xlim=c(0,max(ex)+0.3), ylim=c(0,maxlp+0.5), pch=syms[[ct]],
           xlab=expression(Expected~~-log[10](p)), ylab=expression(Observed~~-log[10](p)),
           main=paste0("STAAR-O (no ACAT-V), cMAC>=",cmac_min," -- ", tag), cex=0.7)
      abline(0,1,col="red"); plotted <- TRUE
    } else points(ex, o, pch=syms[[ct]], cex=0.7)
  }
  legend("topleft", legend=c("pLoF","pLoF+D","Missense","Disruptive Missense","Synonymous"),
         pch=c(0,1,2,3,4), bty="n")
  dev.off()
  cat("  wrote QQ plot.\n", file=log, append=TRUE)
}

## ---------------------------------------------------------------------------
## args[1] = base output dir (default origLPA). Processes <base>/coding and
##          <base>/coding_SPA if present.
args   <- commandArgs(TRUE)
base   <- ifelse(length(args)>=1, args[1], "${ADSP_ROOT}/output/matched_lpa_pc_origLPA")
logf   <- file.path(base, "WHCfix_report.txt")
tagbase<- basename(base)
cat("Wang/Heard-Costa fix + STAAR-O reporting\n", file=logf)        # truncate/create

process_dir(file.path(base,"coding"),     paste0(tagbase," non-SPA (STAAR-O)"), logf)
process_dir(file.path(base,"coding_SPA"), paste0(tagbase," SPA (STAAR-B only)"), logf)

cat("\nDONE. Report at:", logf, "\n", file=logf, append=TRUE)
cat("\nDONE\n")
