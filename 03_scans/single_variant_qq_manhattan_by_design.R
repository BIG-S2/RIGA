# 5-ARM VERSION (2026-09-21): + D_S (all subjects of mixed project x centre cells, cell factor, no matching).
# 4-ARM VERSION (2026-09-21): adds the project-keyed arm D_P (GCAD sequencing project in the exact matching key,
# matched_design_null_model.R); every output carries the suffix _5arm and nothing of the 3-arm run is overwritten.
# Gene-annotated single-variant QQ + Manhattan for the pair33 arms, plus a table of the
# latest signals (revision/tables/latest_sv_signals.csv).  House style copied from
# adsp_sv_qq_manhattan_pair33_cohort_3arm.R: colours #191970 / #C87814, 95% null band, red y = x,
# FULL y-axis with NO cap, png res 190.  Lambda / n genome-wide significant / 1-Mb loci /
# non-APOE loci / APOE peak computed exactly as in that script.
#
# Usage:  Rscript annot_sv_qq_manhattan.R [auto|qc7|asrun]
#   auto  (default) uses <dir>/results_individual_analysis_genome.QC7.Rdata for ALL arms if all three
#         exist, otherwise the as-run QC_label2 files (titles then say "QC_label2 (as run)").
#   Rerun for other arms by editing the ARMS table below.
#
# Annotation: lead variant (min p) of each 1-Mb bin with p < 5e-8 is annotated with the NEAREST
# protein-coding gene from STAARpipeline::genes_info (GRCh38, the same table the gene-centric coding
# drivers iterate over); fallback = GENCODE.Info channel of the per-chromosome aGDS.  Established AD
# loci (30 GRCh38 lead positions, 150-kb window) are flagged and drawn in dark green.
#
# Outputs (all NEW files):
#   figures/supp_fig_S62_sv_qq_pair33_cohort_5arm_annotated.png
#   figures/supp_fig_S63_sv_manhattan_pair33_cohort_5arm_annotated.png
#   tables/latest_sv_signals.csv            one row per (arm, genome-wide-significant 1-Mb locus)
#   tables/annot_sv_5arm_stats.csv          lambda / nsig / loci / non-APOE loci / APOE peak / file used
#   tables/annot_sv_3arm_known_ad_loci.csv  window-min p of the 30 established loci in each arm
#   tables/annot_sv_3arm_calibration.csv    lambda + O/E@1e-4 by MAF bin (house-style table)
#   tables/annot_sv_run_info.txt            gene source, file choice, timestamps
suppressPackageStartupMessages({ library(data.table) })
setDTthreads(4)
OUT <- "${ADSP_ROOT}/output"; REV <- "${ADSP_ROOT}/revision"
FIG <- file.path(REV,"figures"); TAB <- file.path(REV,"tables")
ALPHA <- 5e-8; SUGG <- 1e-5; KNOWN_WIN <- 150e3
COL_A <- "#191970"; COL_B <- "#C87814"; COL_K <- "#006400"; COL_L <- "black"   # odd/even chr, known AD locus, other lead
t0 <- Sys.time(); el <- function() sprintf("[%.1f min]", as.numeric(difftime(Sys.time(), t0, units="mins")))

## ---------------------------------------------------------------- arms (edit to rerun) ----
NP <- nrow(fread(file.path(OUT, "matched/matched_samples_lpa_pc.pair33_projectkey.csv"), select = "IID"))
NS <- nrow(fread(file.path(OUT, "matched/matched_samples_lpa_pc.pair33_projcell.csv"), select = "IID"))
ARMS <- data.table(arm=c("pair33","pair33coh","pair33key","pair33proj","pair33pcell"),
  lab=c("pair33: 1:1 pairs, 37-column null, N=19,408",
        "pair33 + Cohort covariate (same pairs), N=19,408",
        "Cohort in the matching key (D3), N=13,148",
        sprintf("Sequencing project in the matching key (D_P), N=%s", format(NP, big.mark=",")),
        sprintf("Project x centre cells as strata, all subjects, no matching (D_S), N=%s", format(NS, big.mark=","))),
  dir=file.path(OUT, c("matched_lpa_pc_nofilter_pair33","matched_lpa_pc_pair33_cohortcov",
                       "matched_lpa_pc_pair33_cohortkey","matched_lpa_pc_pair33_projectkey","matched_lpa_pc_pair33_projcell"), "individual_analysis"))
NARM <- nrow(ARMS)
ARMS[, asrun := file.path(dir, "results_individual_analysis_genome.Rdata")]
ARMS[, qc7   := file.path(dir, "results_individual_analysis_genome.QC7.Rdata")]
mode <- commandArgs(TRUE)[1]; if (is.na(mode)) mode <- "auto"
use_qc7 <- switch(mode, auto = all(file.exists(ARMS$qc7)), qc7 = TRUE, asrun = FALSE, stop("mode must be auto|qc7|asrun"))
if (use_qc7 && !all(file.exists(ARMS$qc7))) stop("QC7 file(s) missing: ", paste(ARMS$qc7[!file.exists(ARMS$qc7)], collapse=", "))
ARMS[, path := if (use_qc7) qc7 else asrun]; stopifnot(all(file.exists(ARMS$path)))
QCTAG <- if (use_qc7) "QC_label7" else "QC_label2 (as run)"; QCF <- if (use_qc7) "QC7" else "QC2asrun"
cat(sprintf("%s mode=%s  use_qc7=%s  QCTAG='%s'\n", el(), mode, use_qc7, QCTAG)); print(ARMS[, .(arm, path)])

## ---------------------------------------------------------------- established AD loci (GRCh38) ----
## 2026-09-08: SPI1 (rs1057233), WWOX (rs62039712), IQCK (rs7185636) were GRCh37 positions; corrected to GRCh38 from the aGDS rsID
## lookup (tables/annot_verify_sv_known_rsid_positions.csv) and Ensembl REST for rs7185636 (not in aGDS). ABCA7 -> rs4147929, CD2AP -> rs9349407.
KNOWN <- fread(text="locus,chr,pos
APOE,19,44908684
BIN1,2,127135234
CLU,8,27607002
PICALM,11,86156833
CR1,1,207518704
MS4A6A,11,60172122
ABCA7,19,1063444
TREM2,6,41161514
SORL1,11,121564878
CD2AP,6,47485642
EPHA1,7,143413669
PTK2B,8,27362470
HLA-DRB1,6,32615322
SLC24A4,14,92460608
INPP5D,2,233117202
FERMT2,14,52924962
CASS4,20,56443204
SPI1,11,47354897
ADAM10,15,58753575
NYAP1,7,100406823
ECHDC3,10,11678309
ABI3,17,49219935
PLCG2,16,81908423
APH1B,15,63277703
WWOX,16,79321960
IQCK,16,19796841
ACE,17,63471557
SHARPIN,8,144103704
TSPOAP1,17,58331728
ADAMTS1,21,26775872")
KNOWN[, `:=`(chr=as.integer(chr), pos=as.numeric(pos))]
# nearest known locus (within KNOWN_WIN) for vectors of CHR/POS
known_hit <- function(chr_v, pos_v) { if (!length(chr_v)) return(character(0))
  vapply(seq_along(chr_v), function(i) { ci <- chr_v[i]; pi <- pos_v[i]
    k <- KNOWN[chr == ci & abs(pos - pi) <= KNOWN_WIN]; if (!nrow(k)) NA_character_ else k$locus[which.min(abs(k$pos - pi))] }, character(1)) }

## ---------------------------------------------------------------- gene table ----
GEN <- NULL; GSRC <- NA_character_
g <- tryCatch(STAARpipeline::genes_info, error=function(e) NULL)
if (!is.null(g)) GSRC <- "STAARpipeline::genes_info (exported object; columns hgnc_symbol/chromosome_name/start_position/end_position; GRCh38; 18,445 protein-coding genes, same table the gene-centric coding drivers iterate over)"
if (is.null(g)) { e <- new.env(); r <- tryCatch({ data(genes_info, package="STAARpipeline", envir=e); get("genes_info", envir=e) }, error=function(err) NULL)
  if (!is.null(r)) { g <- r; GSRC <- "data(genes_info, package='STAARpipeline')" } }
if (!is.null(g)) { GEN <- as.data.table(g); setnames(GEN, 1:4, c("gene","CHR","start","end"))
  GEN <- GEN[, .(gene=as.character(gene), CHR=suppressWarnings(as.integer(as.character(CHR))), start=as.numeric(start), end=as.numeric(end))][!is.na(CHR) & !is.na(start)]
  cat(sprintf("%s gene table: %s genes on %d chromosomes\n", el(), format(nrow(GEN), big.mark=","), uniqueN(GEN$CHR))) }
# nearest protein-coding gene; distance 0 if inside the gene body
nearest_gene <- function(chr_v, pos_v) {
  if (!length(chr_v)) return(data.table(gene=character(0), dist=numeric(0)))
  if (!is.null(GEN)) return(rbindlist(lapply(seq_along(chr_v), function(i) { ci <- chr_v[i]; pi <- pos_v[i]; gg <- GEN[CHR == ci]
    if (!nrow(gg)) return(data.table(gene=NA_character_, dist=NA_real_))
    d <- ifelse(pi < gg$start, gg$start - pi, ifelse(pi > gg$end, pi - gg$end, 0))
    j <- which.min(d); data.table(gene=gg$gene[j], dist=d[j]) })))
  # ---- fallback: GENCODE.Info of the exact variant in the per-chromosome aGDS (distance NA: gene at the site) ----
  suppressPackageStartupMessages({ library(gdsfmt); library(SeqArray) })
  ad <- get(load(file.path(OUT, "agds_dir.Rdata"))); cat_ <- get(load(file.path(OUT, "Annotation_name_catalog.Rdata")))
  ch <- paste0("annotation/info/FunctionalAnnotation", cat_$dir[cat_$name == "GENCODE.Info"])
  GSRC <<- paste0("aGDS GENCODE.Info channel (", ch, ") at the lead variant")
  out <- data.table(gene=rep(NA_character_, length(chr_v)), dist=NA_real_)
  for (cc in unique(chr_v)) { ii <- which(chr_v == cc); gf <- seqOpen(ad[cc])
    pos <- as.numeric(seqGetData(gf, "position")); sel <- pos %in% pos_v[ii]
    seqSetFilter(gf, variant.sel=sel, verbose=FALSE); v <- as.character(seqGetData(gf, ch)); pp <- pos[sel]; seqClose(gf)
    for (k in ii) { m <- v[pp == pos_v[k]]; m <- m[nzchar(m)]; if (length(m)) out$gene[k] <- sub("[(,;].*$", "", m[1]) } }
  out }

## ---------------------------------------------------------------- load arms ----
rd <- function(path){ e <- new.env(); load(path, envir=e); nm <- ls(e)[1]; x <- as.data.table(get(nm, envir=e)); rm(e)
  cat(sprintf("%s loaded %s\n   object=%s rows=%s cols: %s\n", el(), path, nm, format(nrow(x), big.mark=","), paste(names(x), collapse=", ")))
  pc <- grep("^pvalue$", names(x), ignore.case=TRUE, value=TRUE)[1]; if (is.na(pc)) pc <- grep("^p(val|value)?(_|\\.)?", names(x), ignore.case=TRUE, value=TRUE)[1]
  mc <- grep("^MAF$", names(x), ignore.case=TRUE, value=TRUE)[1]
  setnames(x, c(pc, mc), c("p", "MAF"), skip_absent=TRUE)
  if (!"REF" %in% names(x)) x[, REF := NA_character_]; if (!"ALT" %in% names(x)) x[, ALT := NA_character_]
  x[is.finite(p) & p > 0 & p <= 1, .(CHR=as.integer(CHR), POS=as.numeric(POS), REF=as.character(REF), ALT=as.character(ALT), p=as.numeric(p), MAF=as.numeric(MAF))] }
lam <- function(p) round(median(qchisq(1-p,1))/qchisq(0.5,1),3)
is_apoe <- function(CHR, POS) CHR == 19 & POS > 44.4e6 & POS < 46.4e6

D <- list(); ST <- list(); CAL <- list(); BM <- list(); LEAD <- list(); KW <- list()
for (i in seq_len(nrow(ARMS))) { a <- ARMS$arm[i]
  x <- rd(ARMS$path[i]); setorder(x, p); n <- nrow(x); x[, exp_p := -log10((seq_len(n)-0.5)/n)]
  cat(sprintf("   MAF range %.3g-%.3g (n MAF>0.5: %s)\n", min(x$MAF, na.rm=TRUE), max(x$MAF, na.rm=TRUE), sum(x$MAF > 0.5, na.rm=TRUE)))
  x[, bin := floor(POS/1e6)]
  # (a) house-style statistics
  s <- x[p < ALPHA]; s[, binid := paste0(CHR, "_", bin)]
  ST[[a]] <- list(n=n, lambda=lam(x$p), nsig=nrow(s), loci=uniqueN(s$binid),
    loci_nonAPOE=uniqueN(s[!is_apoe(CHR, POS)]$binid), apoe=max(-log10(x[is_apoe(CHR, POS)]$p)))
  # min p per 1-Mb bin over ALL variants (cross-arm table + QQ labels)
  BM[[a]] <- x[, .(p_bin_min=min(p), POS_bin_min=POS[which.min(p)], MAF_bin_min=MAF[which.min(p)], exp_bin_min=exp_p[which.min(p)]), by=.(CHR, bin)]
  # (b) lead variant per genome-wide-significant 1-Mb locus
  if (nrow(s)) { L <- s[, .SD[which.min(p)], by=.(CHR, bin)]; cnt <- s[, .(n_sig_in_bin=.N), by=.(CHR, bin)]
    L <- cnt[L, on=.(CHR, bin)]; L[, arm := a]; setorder(L, p); LEAD[[a]] <- L }
  # (c) established loci: 150-kb window minimum in this arm
  KW[[a]] <- rbindlist(lapply(seq_len(nrow(KNOWN)), function(k) { w <- x[CHR == KNOWN$chr[k] & abs(POS - KNOWN$pos[k]) <= KNOWN_WIN]
    if (!nrow(w)) return(data.table(arm=a, locus=KNOWN$locus[k], chr=KNOWN$chr[k], pos=KNOWN$pos[k], n_var_window=0L, p_min=NA_real_, POS_min=NA_real_, REF_min=NA_character_, ALT_min=NA_character_, MAF_min=NA_real_, exp_min=NA_real_))
    j <- which.min(w$p); data.table(arm=a, locus=KNOWN$locus[k], chr=KNOWN$chr[k], pos=KNOWN$pos[k], n_var_window=nrow(w), p_min=w$p[j], POS_min=w$POS[j], REF_min=w$REF[j], ALT_min=w$ALT[j], MAF_min=w$MAF[j], exp_min=w$exp_p[j]) }))
  # calibration by MAF bin (house style)
  x[, maf := pmin(MAF, 1-MAF)]; x[, mbin := cut(maf, c(0,0.001,0.005,0.01,0.05,0.5), include.lowest=TRUE)]
  CAL[[a]] <- rbind(x[, .(arm=a, maf_bin="all", n=.N, lambda=lam(p), n_lt1e4=sum(p<1e-4), exp_lt1e4=.N*1e-4)],
                    x[, .(arm=a, n=.N, lambda=lam(p), n_lt1e4=sum(p<1e-4), exp_lt1e4=.N*1e-4), by=.(maf_bin=as.character(mbin))])
  keep <- c(which(x$p < 1e-3), round(seq(sum(x$p < 1e-3)+1, n, length.out=min(150000, n))))
  D[[a]] <- x[unique(keep), .(CHR, POS, p, MAF, exp_p)]; rm(x, s); gc(verbose=FALSE)
  cat(sprintf("%s %-9s n=%s lambda=%.3f nsig=%s loci=%d nonAPOE=%d APOEpeak=%.1f\n", el(), a, format(n, big.mark=","), ST[[a]]$lambda, format(ST[[a]]$nsig, big.mark=","), ST[[a]]$loci, ST[[a]]$loci_nonAPOE, ST[[a]]$apoe))
}
cal <- rbindlist(CAL); cal[, OE_1e4 := round(n_lt1e4/exp_lt1e4, 3)]; cal[, qc := QCTAG]
fwrite(cal, file.path(TAB, "annot_sv_5arm_calibration.csv")); print(cal)
stats <- rbindlist(lapply(ARMS$arm, function(a) data.table(arm=a, label=ARMS[arm==a]$lab, qc=QCTAG, file=ARMS[arm==a]$path,
  n_variants=ST[[a]]$n, lambda=ST[[a]]$lambda, n_gwsig_variants=ST[[a]]$nsig, n_loci_1Mb=ST[[a]]$loci, n_loci_nonAPOE=ST[[a]]$loci_nonAPOE, apoe_peak_log10p=round(ST[[a]]$apoe, 2))))
fwrite(stats, file.path(TAB, "annot_sv_5arm_stats.csv")); print(stats[, -"file"])

## ---------------------------------------------------------------- annotate leads + known loci ----
lead <- rbindlist(LEAD, use.names=TRUE)
if (!nrow(lead)) lead <- data.table(arm=character(0), CHR=integer(0), bin=numeric(0), POS=numeric(0), p=numeric(0), nearest_gene=character(0), known_locus=character(0), known_ad_locus=logical(0))
if (nrow(lead)) { ng <- nearest_gene(lead$CHR, lead$POS); lead[, `:=`(nearest_gene=ng$gene, dist_bp=ng$dist)]
  lead[, known_locus := known_hit(CHR, POS)]; lead[, known_ad_locus := !is.na(known_locus)] }
kw <- rbindlist(KW, use.names=TRUE)
kk <- kw[!is.na(POS_min)]; if (nrow(kk)) { ng <- nearest_gene(kk$chr, kk$POS_min); kk[, `:=`(nearest_gene_at_min=ng$gene, dist_bp_at_min=ng$dist)]
  kw <- merge(kw, kk[, .(arm, locus, nearest_gene_at_min, dist_bp_at_min)], by=c("arm","locus"), all.x=TRUE) }
kw[, `:=`(gwsig=p_min < ALPHA, suggestive=p_min < SUGG, qc=QCTAG)]
fwrite(kw[order(arm, p_min)], file.path(TAB, "annot_sv_5arm_known_ad_loci.csv"))
cat("\n=== established AD loci with window-min p < 1e-5 ===\n"); print(kw[suggestive %in% TRUE][order(arm, p_min), .(arm, locus, chr, POS_min, MAF_min, p_min=signif(p_min,3), nearest_gene_at_min)], nrows=200)

# (f) latest-signals table: one row per (arm, genome-wide-significant 1-Mb locus), plus bin-min p in the other arms
if (nrow(lead)) {
  for (b in ARMS$arm) { v <- BM[[b]][lead, on=.(CHR, bin)]$p_bin_min; lead[, (paste0("p_min_bin_", b)) := v] }
  sig <- lead[, c("arm","CHR","bin","POS","REF","ALT","MAF","p","n_sig_in_bin","nearest_gene","dist_bp","known_ad_locus","known_locus", paste0("p_min_bin_", ARMS$arm)), with=FALSE]
  setnames(sig, c("CHR","POS","bin","MAF","p"), c("chr","lead_pos","bin_Mb","lead_MAF","lead_p"))
  sig[, lead_maf_minor := pmin(lead_MAF, 1-lead_MAF)]; sig[, apoe_region := is_apoe(chr, lead_pos)]; sig[, qc := QCTAG]
  setorder(sig, arm, lead_p); fwrite(sig, file.path(TAB, "latest_sv_signals_5arm.csv"))
  cat(sprintf("\n%s wrote latest_sv_signals_5arm.csv (%d rows)\n", el(), nrow(sig)))
  for (a in ARMS$arm) { cat(sprintf("\n=== %s: %d genome-wide-significant loci ===\n", a, sum(sig$arm==a)))
    print(sig[arm==a, .(chr, lead_pos, REF, ALT, MAF=signif(lead_MAF,3), p=signif(lead_p,3), n_sig_in_bin, nearest_gene, dist_bp, known_locus)], nrows=60) }
} else { fwrite(data.table(arm=character(), note="no genome-wide-significant variants in any arm"), file.path(TAB, "latest_sv_signals_5arm.csv")) }

## ---------------------------------------------------------------- label helpers ----
# label text: nearest gene, with the established locus in parentheses when the names differ
lab_text <- function(gene, known) ifelse(is.na(known), gene, ifelse(is.na(gene) | gene == known, known, paste0(gene, " (", known, ")")))
# Manhattan labels for one arm: leads (all if <= 25 loci else top 20) + every established locus with window min p < 1e-5 (+ APOE always)
manh_labels <- function(a) {
  L <- lead[arm == a]; if (nrow(L) > 25) L <- L[order(p)][1:20]
  A <- if (nrow(L)) L[, .(CHR, POS, p, text=lab_text(nearest_gene, known_locus), col=ifelse(known_ad_locus, COL_K, COL_L), src="lead")] else NULL
  K <- kw[arm == a & !is.na(p_min) & (p_min < SUGG | locus == "APOE")]
  B <- if (nrow(K)) K[, .(CHR=chr, POS=POS_min, p=p_min, text=lab_text(nearest_gene_at_min, locus), col=COL_K, src="known")] else NULL
  z <- rbindlist(list(B, A), use.names=TRUE); if (!nrow(z)) return(NULL)
  z <- z[!duplicated(z[, .(CHR, POS)])]                                    # same variant: known-locus label wins
  setorder(z, p); z[, key := sub(" .*$", "", text)]
  z <- z[!duplicated(paste(CHR, key))]                                     # same gene twice on a chromosome: keep the strongest
  z[, -"key"] }
# label rows in ABSOLUTE y: each label starts one `step` above its point and is pushed up until it is at least
# `step` away from every already-placed label within `w` (x units); processed left to right
place_labels <- function(x, y, w, step) { o <- order(x); yl <- y + step
  for (k in seq_along(o)) { i <- o[k]; prev <- o[seq_len(k-1)]
    repeat { moved <- FALSE
      for (j in prev) if (abs(x[i] - x[j]) < w && abs(yl[i] - yl[j]) < step*0.98) { yl[i] <- yl[j] + step; moved <- TRUE }
      if (!moved) break } }
  yl }

## ---------------------------------------------------------------- (e) QQ, 3 panels side by side ----
XMAX <- ceiling(max(sapply(D, function(z) max(z$exp_p))))
QCAP <- ceiling(max(sapply(D, function(z) max(-log10(z$p)))) * 1.04); MCAP <- QCAP   # full range, identical across arms
png(file.path(FIG, paste0("supp_fig_S62_sv_qq_pair33_cohort_5arm_", QCF, "_annotated.png")), width=1500*NARM, height=1560, res=190)
par(mfrow=c(1,NARM), mar=c(5.2,5.6,5.6,1.6), cex.lab=1.35, cex.axis=1.15, font.lab=2)
for (a in ARMS$arm) { z <- D[[a]]; s <- ST[[a]]
  plot(NA, xlim=c(0,XMAX), ylim=c(0,QCAP), xlab="Expected -log10(p-value)", ylab="Observed -log10(p-value)", main="", yaxs="i")
  nn <- s$n; k <- round(10^seq(0, log10(nn), length.out=500))
  polygon(c(-log10((k-0.5)/nn), rev(-log10((k-0.5)/nn))), c(-log10(qbeta(0.975,k,nn-k+1)), rev(-log10(qbeta(0.025,k,nn-k+1)))), col="grey85", border=NA)
  o <- -log10(z$p); on <- o > QCAP
  points(z$exp_p[!on], o[!on], pch=1, cex=0.45, col=COL_A)
  if (any(on)) points(z$exp_p[on], rep(QCAP*0.985, sum(on)), pch=17, cex=0.8, col="red3")
  abline(0,1,col="red",lwd=2); title(main=paste0(ARMS$lab[ARMS$arm==a], "  [", QCTAG, "]"), cex.main=1.02, line=3.5)
  mtext(sprintf("lambda = %.3f    %s variants    %s genome-wide sig    %s non-APOE loci", s$lambda,
        format(s$n,big.mark=","), format(s$nsig,big.mark=","), format(s$loci_nonAPOE,big.mark=",")), side=3, line=1.9, cex=0.84)
  mtext(sprintf("full y-axis (no cap); APOE peak = %.0f; labels = lead variant of the strongest 1-Mb loci at p<5e-8, max 10 (nearest gene)", s$apoe), side=3, line=0.6, cex=0.8, col="grey30")
  # labels: lead variant of the strongest genome-wide-significant 1-Mb loci (max 10)
  top <- head(BM[[a]][p_bin_min < ALPHA][order(p_bin_min)], 10)   # 2026-09-08: only loci at p<5e-8 (was: 10 strongest regardless)
  ng <- nearest_gene(top$CHR, top$POS_bin_min); top[, `:=`(gene=ng$gene, known=known_hit(CHR, POS_bin_min))]
  top[, `:=`(x=exp_bin_min, y=-log10(p_bin_min), text=lab_text(gene, known), col=ifelse(is.na(known), COL_L, COL_K))]
  setorder(top, y); gap <- QCAP*0.045; yl <- top$y
  for (j in seq_len(nrow(top))) if (j > 1) yl[j] <- max(yl[j], yl[j-1] + gap)
  xl <- top$x - XMAX*0.10
  segments(top$x, top$y, xl, yl, col="grey55", lwd=0.7)
  text(xl, yl, top$text, pos=2, cex=0.78, col=top$col, font=ifelse(top$col == COL_K, 4, 3), offset=0.25)
  points(top$x, top$y, pch=1, cex=0.9, col=top$col, lwd=1.4)
  legend("topleft", bty="o", bg="white", cex=1.02, legend=c("observed","y = x","95% null band","established AD locus (150 kb)"),
         pch=c(1,NA,15,1), lty=c(NA,1,NA,NA), col=c(COL_A,"red","grey85",COL_K), pt.lwd=c(1,1,1,1.6)) }
dev.off(); cat(el(), "wrote QQ\n")

## ---------------------------------------------------------------- (d) Manhattan, 3 panels stacked ----
CHRLEN <- c(248956422,242193529,198295559,190214555,181538259,170805979,159345973,145138636,138394717,133797422,135086622,
            133275309,114364328,107043718,101991189,90338345,83257441,80373285,58617616,64444167,46709983,50818468)
off <- c(0, cumsum(as.numeric(CHRLEN))[-22]); names(off) <- 1:22; mid <- off + CHRLEN/2; GW <- sum(as.numeric(CHRLEN))
png(file.path(FIG, paste0("supp_fig_S63_sv_manhattan_pair33_cohort_5arm_", QCF, "_annotated.png")), width=2300, height=1180*NARM, res=190)
par(mfrow=c(NARM,1), mar=c(4.8,5.6,5.2,1.6), cex.lab=1.3, cex.axis=1.1, font.lab=2)
for (a in ARMS$arm) { z <- D[[a]][CHR %in% 1:22]; s <- ST[[a]]
  z[, gpos := POS + off[as.character(CHR)]]; z[, col := ifelse(CHR %% 2 == 1, COL_A, COL_B)]
  o <- -log10(z$p); on <- o > MCAP
  plot(z$gpos[!on], o[!on], pch=20, cex=0.36, col=z$col[!on], xaxt="n", yaxs="i", xlab="Chromosome", ylab="-log10(p-value)",
       ylim=c(0,MCAP), xlim=c(0,GW), main="")
  if (any(on)) points(z$gpos[on], rep(MCAP*0.98, sum(on)), pch=17, cex=0.75, col="red3")
  axis(1, at=mid, labels=1:22, tick=FALSE, cex.axis=0.92); abline(h=-log10(ALPHA), col="red", lty=3, lwd=2)
  title(main=paste0(ARMS$lab[ARMS$arm==a], "  [", QCTAG, "]"), cex.main=1.2, line=3.2)
  mtext(sprintf("lambda = %.3f    %s genome-wide significant variants    %s independent loci    %s outside APOE", s$lambda,
        format(s$nsig,big.mark=","), format(s$loci,big.mark=","), format(s$loci_nonAPOE,big.mark=",")), side=3, line=1.7, cex=0.84)
  nl <- sum(lead$arm == a)
  mtext(sprintf("full y-axis (no cap); APOE peak = %.0f; black = nearest gene of %s lead variants; green = established AD loci with window min p < 1e-5", s$apoe,
        if (nl <= 25) sprintf("all %d", nl) else "the 20 strongest"), side=3, line=0.5, cex=0.8, col="grey30")
  lb <- manh_labels(a)
  if (!is.null(lb) && nrow(lb)) { lb <- lb[CHR %in% 1:22]; lb[, `:=`(gpos=POS + off[as.character(CHR)], y=-log10(p))]
    step <- MCAP*0.036; lb[, yl := place_labels(gpos, y, w=0.075*GW, step=step)]
    lb[, right := yl > MCAP*0.93]; lb[right == TRUE, yl := pmin(y, MCAP*0.93)]
    lb[, xl := gpos]
    segments(lb[right==FALSE]$gpos, lb[right==FALSE]$y + MCAP*0.006, lb[right==FALSE]$gpos, lb[right==FALSE]$yl - MCAP*0.004, col="grey55", lwd=0.7)
    text(lb[right==FALSE]$xl, lb[right==FALSE]$yl, lb[right==FALSE]$text, pos=3, offset=0.15, cex=0.74, col=lb[right==FALSE]$col, font=ifelse(lb[right==FALSE]$col == COL_K, 4, 3))
    if (any(lb$right)) text(lb[right==TRUE]$xl, lb[right==TRUE]$yl, lb[right==TRUE]$text, pos=4, offset=0.3, cex=0.95, col=lb[right==TRUE]$col, font=ifelse(lb[right==TRUE]$col == COL_K, 4, 3))
    points(lb$gpos, lb$y, pch=1, cex=0.85, col=lb$col, lwd=1.2) }
  legend("topleft", bty="o", bg="white", cex=1.0, horiz=TRUE, legend=c(sprintf("p = %g", ALPHA), "lead variant (nearest gene)", "established AD locus"),
         lty=c(3,NA,NA), lwd=c(2,NA,NA), pch=c(NA,1,1), col=c("red",COL_L,COL_K)) }
dev.off(); cat(el(), "wrote Manhattan\n")

writeLines(c(sprintf("run: %s", format(Sys.time())), sprintf("mode=%s use_qc7=%s QCTAG=%s", mode, use_qc7, QCTAG),
             paste0("gene_source: ", GSRC), paste0("files: ", paste(ARMS$arm, ARMS$path, sep="=", collapse=" ; ")),
             sprintf("loci per arm (1-Mb bins, p<5e-8): %s", paste(sapply(ARMS$arm, function(a) sprintf("%s=%d", a, ST[[a]]$loci)), collapse=", "))),
           file.path(TAB, "annot_sv_run_info_5arm.txt"))
cat("gene source:", GSRC, "\n"); cat(el(), "DONE\n")
