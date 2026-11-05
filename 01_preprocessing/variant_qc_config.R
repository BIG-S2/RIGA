# variant_qc_config.R -- SINGLE SOURCE OF TRUTH for the variant-QC node used by every
# STAARpipeline driver in src/.  Written 2026-09-08 as the implementation sketch for
# removing QC-label confusion WITHOUT touching the 22 aGDS files.
#
# How to wire it in (one two-line edit per driver; nothing else changes):
#   src/STAARpipeline_Individual_Analysis.*.r      line 33-34   (## QC_label / QC_label <- "annotation/info/QC_label2")
#   src/STAARpipeline_Gene_Centric_Coding.*.r      line 31-32   (most files; SPA/origLPA variants at 39-40, admixq/pgmm/vbgmm at 11-12)
#   src/STAARpipeline_Gene_Centric_Noncoding.*.r   line 31-32
#   src/STAARpipeline_Gene_Centric_ncRNA.*.r       line 29-30
#   src/STAARpipeline_Sliding_Window.*.r           line 29-30
#   src/Association_Analysis_PreStep.r             line 24      (QC_label3 -> chunking on a label no driver tests)
#   src/STAARpipelineSummary_Individual_Analysis*.r line 45-48  (QC_label3 in 7 files)
# replace the literal with:
#   source("${ADSP_ROOT}/revision/scripts/variant_qc_config.R")
#   QC_label <- qc_label()            # primary, or an explicitly declared legacy label
#
# The exact line of every driver is listed in revision/tables/qclabel_provenance.csv
# (column qc_label_line).

## The primary label.  QC_label7 = per-variant RUTH  AN>3636 & HWE_SLP_I>-4 & F_MISSING<0.05
## (revision/scripts/build_variant_qc_label.R via qclabel7.slurm, 2026-08-28).  QC_label5 differs by
## the FMISS clause only and was shown to be a no-op (REPORT_pair33_primary_arm.md section 12.9).
QC_LABEL_PRIMARY <- "annotation/info/QC_label7"

## Labels that may still be requested EXPLICITLY (for reproducing an archived result dir).
## Anything else -- including the bare placeholder "annotation/info/QC_label" (all PASS) and
## the multiallelic-excluding QC_label3 -- is refused.
QC_LABEL_LEGACY <- c(
  QC_label5 = "annotation/info/QC_label5",   # equivalent to primary minus FMISS
  QC_label2 = "annotation/info/QC_label2"    # DEFECTIVE position-inherited label (section 11.1c)
)

## qc_label(): returns the node path to use.
##   * default                     -> QC_LABEL_PRIMARY
##   * env  ADSP_QC_LABEL=QC_label2 -> the named legacy label, but ONLY if the caller also sets
##                                     ADSP_QC_LABEL_REASON (free text, logged), so a legacy
##                                     run can never happen silently.
qc_label <- function() {
  want <- Sys.getenv("ADSP_QC_LABEL", unset = "")
  if (want == "" || want == "QC_label7" || want == QC_LABEL_PRIMARY) {
    message(sprintf("[qclabel_config] using PRIMARY %s", QC_LABEL_PRIMARY))
    return(QC_LABEL_PRIMARY)
  }
  key <- sub("^annotation/info/", "", want)
  if (!key %in% names(QC_LABEL_LEGACY))
    stop(sprintf("[qclabel_config] refusing QC label '%s': not primary and not in the declared legacy set (%s)",
                 want, paste(names(QC_LABEL_LEGACY), collapse = ", ")))
  reason <- Sys.getenv("ADSP_QC_LABEL_REASON", unset = "")
  if (reason == "")
    stop(sprintf("[qclabel_config] legacy label %s requested without ADSP_QC_LABEL_REASON; refusing", key))
  message(sprintf("[qclabel_config] using LEGACY %s  (reason: %s)", QC_LABEL_LEGACY[[key]], reason))
  QC_LABEL_LEGACY[[key]]
}

## Optional guard for drivers that have NOT yet been converted: call
##   assert_primary_qc(QC_label)
## right after the literal assignment; it aborts the job if the literal is not the primary.
assert_primary_qc <- function(x) {
  if (!identical(x, QC_LABEL_PRIMARY) && Sys.getenv("ADSP_QC_LABEL", unset = "") == "")
    stop(sprintf("[qclabel_config] driver hard-codes %s but primary is %s; set ADSP_QC_LABEL/ADSP_QC_LABEL_REASON to override",
                 x, QC_LABEL_PRIMARY))
  invisible(TRUE)
}
