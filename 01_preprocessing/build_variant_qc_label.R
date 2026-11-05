# Corrected QC label: PER-VARIANT, not per-position.
#
# THE BUG BEING FIXED. Add_RUTH_QC_label.R assigns QC_label2 with foverlaps(var, bed) --
# a variant PASSES if its POSITION falls inside any BED interval. The BED intervals come
# from variants that passed the HWE filter, which are overwhelmingly ultra-rare (the HWE
# test has no power there, so p ~ 1). A common variant therefore inherits PASS from an
# ultra-rare neighbour. Measured on chr22: of 229,088 variants the single-variant scan
# tested, only 1,361 (0.6%) passed the HWE filter themselves; 227,727 (99.4%) inherited it.
#
# THE FIX. Match on the exact variant (POS, REF, ALT) against the RUTH site-only VCFs,
# which carry AN / HWE_SLP_I / FIBC_I / F_MISSING per allele.
#
# THRESHOLD is a parameter, not a decision baked into the code -- see THRESH below. The
# node is written as QC_label4, ADDITIVE: QC_label/2/3 are left untouched so nothing that
# currently runs changes behaviour until a driver is pointed at the new node.
#
# usage: Rscript build_variant_qc_label.R <chr> [thresh_expr] [--write]
#   thresh_expr is R code over columns AN, SLP, FIBC, FMISS. Default is the documented
#   filter. Without --write the script only reports what the label WOULD be.

suppressPackageStartupMessages({ library(data.table); library(SeqArray); library(gdsfmt) })
args   <- commandArgs(TRUE)
chr    <- as.integer(args[1])
THRESH <- if (length(args) >= 2 && !grepl("^--", args[2])) args[2] else "AN > 3636 & abs(SLP) < 4"
WRITE  <- any(args == "--write")
# node name is a parameter: --node=QC_labelN. Nodes are ADDITIVE, so several filter
# definitions can coexist in the aGDS and be compared without rebuilding anything.
NODE   <- { a <- grep("^--node=", args, value=TRUE)
            if (length(a)) sub("^--node=", "", a[1]) else "QC_label4" }

RUTHDIR <- "${ADSP_ROOT}/data/ruth"
GDS     <- sprintf("${ADSP_ROOT}/processed/r5_chr%d.split.gds", chr)
cat(sprintf("chr%d | threshold: %s | write: %s\n", chr, THRESH, WRITE))

## ---- 1. per-variant RUTH table ----------------------------------------------------
fs <- list.files(RUTHDIR, pattern=sprintf("^ruth\\.site_only\\.chr%d\\.[0-9-]+\\.vcf(\\.gz)?$", chr),
                 full.names=TRUE)
if (!length(fs)) stop("no RUTH VCFs for chr", chr)
cat("  reading", length(fs), "RUTH VCFs ...\n")
# Concatenate with a single shell pass. An earlier R-level loop using cat(append=TRUE)
# introduced a blank line between files and fread() stopped early at 572,544 of 8,496,648
# records -- hence the shell redirect here, which cannot produce that boundary.
tmp <- tempfile(fileext=".tsv")
FMT <- "%POS\\t%REF\\t%ALT\\t%INFO/AN\\t%INFO/HWE_SLP_I\\t%INFO/FIBC_I\\t%INFO/F_MISSING\\n"
cmd <- paste0("for f in ", paste(shQuote(fs), collapse=" "),
              "; do bcftools query -f ", shQuote(FMT), " \"$f\"; done > ", shQuote(tmp))
st <- system2("bash", c("-c", shQuote(cmd)))
if (st != 0) stop("bcftools extraction failed (exit ", st, ")")
r <- fread(tmp, header=FALSE, col.names=c("POS","REF","ALT","AN","SLP","FIBC","FMISS"))
unlink(tmp)
r <- r[!is.na(POS)]
cat(sprintf("  RUTH records: %s\n", format(nrow(r), big.mark=",")))
r[, keep := eval(parse(text=THRESH), envir=.SD)]
cat(sprintf("  pass threshold (per variant): %s (%.1f%%)\n",
            format(sum(r$keep), big.mark=","), 100*mean(r$keep)))
# RUTH chunk boundaries can emit a variant twice (8,496,648 records vs 8,496,642 GDS
# variants on chr22). Deduplicate on the variant key before joining, or the join
# multiplies rows and breaks the 1:1 alignment to GDS variant order.
ndup <- nrow(r) - uniqueN(r, by=c("POS","REF","ALT"))
if (ndup > 0) { cat(sprintf("  duplicate RUTH keys removed: %s\n", format(ndup, big.mark=","))) 
                r <- unique(r, by=c("POS","REF","ALT")) }
setkey(r, POS, REF, ALT)

## ---- 2. align to the GDS variant order --------------------------------------------
g <- seqOpen(GDS, readonly = !WRITE)
al  <- seqGetData(g, "allele")
v <- data.table(idx = seq_along(al),
                POS = as.numeric(seqGetData(g, "position")),
                REF = sub(",.*$", "", al),
                ALT = sub("^[^,]*,", "", al))
setkey(v, POS, REF, ALT)
m <- r[v, on=.(POS,REF,ALT)]
setorder(m, idx)
stopifnot(nrow(m) == length(al), !is.na(m$idx), all(m$idx == seq_along(al)))
lab <- ifelse(!is.na(m$keep) & m$keep, "PASS", "FAIL")
cat(sprintf("  GDS variants: %s | matched to RUTH: %s (%.2f%%) | unmatched -> FAIL: %s\n",
            format(length(al), big.mark=","), format(sum(!is.na(m$keep)), big.mark=","),
            100*mean(!is.na(m$keep)), format(sum(is.na(m$keep)), big.mark=",")))

## ---- 3. compare to the existing labels --------------------------------------------
old2 <- seqGetData(g, "annotation/info/QC_label2") == "PASS"
new4 <- lab == "PASS"
cat(sprintf("\n  QC_label2 PASS : %s (%.1f%%)\n", format(sum(old2), big.mark=","), 100*mean(old2)))
cat(sprintf("  %s PASS : %s (%.1f%%)\n", NODE, format(sum(new4), big.mark=","), 100*mean(new4)))
cat(sprintf("  PASS under QC_label2 but FAIL under %s (inherited labels removed): %s\n",
            NODE, format(sum(old2 & !new4), big.mark=",")))
cat(sprintf("  FAIL under QC_label2 but PASS under %s (recovered): %s\n",
            NODE, format(sum(!old2 & new4), big.mark=",")))

if (WRITE) {
  ann <- index.gdsn(g, "annotation/info")
  if (NODE %in% ls.gdsn(ann)) delete.gdsn(index.gdsn(g, paste0("annotation/info/", NODE)))
  add.gdsn(ann, NODE, val=lab, compress="LZMA_RA", closezip=TRUE)
  cat(sprintf("\n  WROTE annotation/info/%s to %s\n", NODE, basename(GDS)))
} else cat(sprintf("\n  (dry run -- pass --write to add annotation/info/%s)\n", NODE))
seqClose(g)
cat("DONE\n")
