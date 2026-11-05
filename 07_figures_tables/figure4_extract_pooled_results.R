#!/usr/bin/env Rscript
# Extract existing gene-set statistics only: no association scans or model fits.
# Matches gene_set_primary_statistic.R's component selection, equal weights,
# endpoint clipping and cMAC>=10 rule. Primary RIGA tables are never written.
# Usage: Rscript --vanilla adsp/src/figure4_extract_gene_comparison.R [outdir]
suppressPackageStartupMessages(library(data.table))
setDTthreads(1L)
args <- commandArgs(TRUE)
OUT <- if (length(args)) args[1] else "/tmp/riga_fig_comparison_20261015"
ROOT <- "${ADSP_ROOT}"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
ARMS <- c(pooled = "pooled_lpa_pc_full")
PAT <- c(coding = "^ADSP_AD_Coding_[0-9]+\\.Rdata$",
         noncoding = "^ADSP_AD_Noncoding_[0-9]+\\.Rdata$", ncRNA = "^ADSP_AD_ncRNA_[0-9]+\\.Rdata$")
SB_PATTERN <- "^(SKAT|Burden)\\(1,(25|1)\\)"
P_FLOOR <- 1e-300
P_CEILING <- 1 - 1e-16
clip_p <- function(P) pmin(pmax(P, P_FLOOR), P_CEILING)

# Same algebra as the legacy CCT, evaluated without subtracting p from 0.5
# in the extreme tails and without subtracting atan(T)/pi from 0.5 for T>0.
cct_stable <- function(P) {
  P <- clip_p(P)
  V <- matrix(NA_real_, nrow(P), ncol(P))
  lo <- !is.na(P) & P < 0.5
  hi <- !is.na(P) & P > 0.5
  mid <- !is.na(P) & P == 0.5
  V[lo] <- 1 / tan(pi * P[lo])
  V[hi] <- -1 / tan(pi * (1 - P[hi]))
  V[mid] <- 0
  z <- rowMeans(V, na.rm = TRUE)
  ans <- rep(NA_real_, length(z))
  pos <- !is.na(z) & z > 0
  neg <- !is.na(z) & z < 0
  zero <- !is.na(z) & z == 0
  ans[pos] <- atan(1 / z[pos]) / pi
  ans[neg] <- 1 - atan(-1 / z[neg]) / pi
  ans[zero] <- 0.5
  ans
}
cct_legacy <- function(P) {
  P <- clip_p(P)
  0.5 - atan(rowMeans(tan((0.5 - P) * pi), na.rm = TRUE)) / pi
}
lambda_gc <- function(p) {
  if (!length(p)) return(NA_real_)
  qchisq(median(p), 1, lower.tail = FALSE) / qchisq(0.5, 1)
}
clip_note <- paste0("Same explicit clipping as legacy: p<", P_FLOOR,
                    " -> ", P_FLOOR, "; p>", format(P_CEILING, digits = 17),
                    " -> ", format(P_CEILING, digits = 17),
                    ". Raw zero P values are counted, not discarded; exact sub-floor tails cannot be recovered.")

# A constant vector has CCT equal to its component p. This directly checks
# the numerical problem relevant to inflated pooled tails.
fixtures <- c(0, 1e-300, 1e-100, 1e-30, 1e-16, 1e-10, 0.01, 0.5, 0.9, 1)
M <- matrix(rep(fixtures, 4), nrow = length(fixtures), ncol = 4)
fixture <- data.table(input_p = fixtures, expected_after_clipping = clip_p(fixtures),
                      SB = cct_stable(M), SB_legacy = cct_legacy(M))
stopifnot(all(abs(fixture$SB - fixture$expected_after_clipping) <=
                2e-14 * pmax(fixture$expected_after_clipping, P_FLOOR)))
fwrite(fixture, file.path(OUT, "cct_numerical_fixtures.csv"))

# Both source drivers use ceil(genes per chromosome/50). Retrieve the same
# annotation object from the installed package, without opening genotype data.
expected <- 379L
coverage_basis <- "Fallback 379, independently asserted in revision/scripts/annot_verify_coding_3arm.R"
gene_info <- tryCatch({
  ns <- loadNamespace("STAARpipeline")
  get("genes_info", envir = ns)
}, error = function(e) NULL)
if (!is.null(gene_info)) {
  expected <- as.integer(sum(ceiling(table(gene_info[, 2]) / 50)))
  coverage_basis <- paste0("STAARpipeline::genes_info; ", nrow(gene_info),
                           " genes; sum(ceiling(table(chr)/50))")
}
stopifnot(expected > 0L)
cat("Expected chunks:", expected, "(", coverage_basis, ")\n")

snapshots <- list()
coverage <- list()
manifest <- list()
for (arm in names(ARMS)) for (fam in names(PAT)) {
  expected <- if (fam == "ncRNA") 222L else 379L
  d <- file.path(ROOT, "output", ARMS[[arm]], fam)
  fs <- list.files(d, pattern = PAT[[fam]], full.names = TRUE)
  ids <- as.integer(sub(".*_([0-9]+)\\.Rdata$", "\\1", fs))
  ord <- order(ids); fs <- fs[ord]; ids <- ids[ord]
  fi <- file.info(fs)
  s <- data.table(arm, family = fam, chunk_id = ids, input_file = fs,
                  bytes = fi$size, mtime = as.character(fi$mtime),
                  mtime_numeric = as.numeric(fi$mtime))
  snapshots[[paste(arm, fam, sep = "_")]] <- s
  manifest[[length(manifest) + 1L]] <- s
  missing <- setdiff(seq_len(expected), ids)
  extra <- setdiff(ids, seq_len(expected))
  coverage[[length(coverage) + 1L]] <- data.table(
    arm, family = fam, expected_chunks = expected, observed_chunks = length(ids),
    complete = !length(missing) && !length(extra) && !anyDuplicated(ids),
    missing_ids = paste(missing, collapse = ";"), extra_ids = paste(extra, collapse = ";"),
    coverage_basis = if (fam == "ncRNA") "STAARpipeline::ncRNA_gene; 21104 catalog rows; sum(ceiling(table(chr)/100))" else coverage_basis, extraction = if (fam == "coding") "previously_extracted_complete_family" else "requested_complete_family")
}
cov <- rbindlist(coverage)
fwrite(cov, file.path(OUT, "gene_chunk_coverage.csv"))
fwrite(rbindlist(manifest), file.path(OUT, "gene_input_manifest.csv"))

read_one <- function(file, expected_object) {
  e <- new.env(parent = emptyenv())
  objects <- load(file, envir = e)
  if (!expected_object %in% objects)
    stop("Expected object ", expected_object, " absent from ", file)
  x <- get(expected_object, envir = e)
  if (is.matrix(x) || is.data.frame(x)) x <- list(x)
  rows <- lapply(x, function(m) {
    if (is.null(m) || !length(m)) return(NULL)
    z <- as.data.table(lapply(as.data.frame(m, stringsAsFactors = FALSE,
                                          check.names = FALSE), unlist))
    if (anyDuplicated(names(z))) stop("Duplicated columns in ", file)
    comps <- grep(SB_PATTERN, names(z), value = TRUE)
    if (!length(comps)) stop("No SKAT/burden components in ", file)
    z[, component_schema := paste(comps, collapse = "|")]
    z[, expected_components := length(comps)]
    z
  })
  z <- rbindlist(rows, fill = TRUE)
  if (nrow(z)) z[, input_file := file]
  list(data = z, object_names = paste(objects, collapse = ";"),
       list_entries = length(x), nonempty_entries = sum(vapply(rows, Negate(is.null), logical(1))))
}

all_summaries <- list()
for (arm in names(ARMS)) for (fam in c("noncoding", "ncRNA")) {
  expected <- if (fam == "ncRNA") 222L else 379L
  key <- paste(arm, fam, sep = "_")
  # Avoid depending on data.table's scoping of the external arm variable.
  coverage_row <- which(cov[["arm"]] == arm & cov[["family"]] == fam)
  cv <- cov[coverage_row]
  if (!isTRUE(cv$complete)) stop("Incomplete coding family: ", key)
  snapshot <- snapshots[[key]]
  data_parts <- vector("list", nrow(snapshot))
  file_audit <- vector("list", nrow(snapshot))
  for (i in seq_len(nrow(snapshot))) {
    f <- snapshot$input_file[i]
    before <- file.info(f)
    x <- read_one(f, paste0("results_", fam))
    hash <- unname(tools::md5sum(f))
    after <- file.info(f)
    if (before$size != after$size || as.numeric(before$mtime) != as.numeric(after$mtime) ||
        before$size != snapshot$bytes[i] || as.numeric(before$mtime) != snapshot$mtime_numeric[i])
      stop("Input changed after snapshot/during extraction: ", f)
    data_parts[[i]] <- x$data
    file_audit[[i]] <- data.table(input_file = f, chunk_id = snapshot$chunk_id[i],
                                 bytes = after$size, mtime = as.character(after$mtime), md5 = hash,
                                 object_names = x$object_names, list_entries = x$list_entries,
                                 nonempty_entries = x$nonempty_entries, result_rows = nrow(x$data))
    if (i %% 50L == 0L) cat(key, "read", i, "of", nrow(snapshot), "chunks\n")
  }
  o <- rbindlist(data_parts, fill = TRUE)
  fwrite(rbindlist(file_audit), file.path(OUT, paste0(key, "_input_audit.csv")))
  required <- c("Gene name", "Chr", "Category", "#SNV", "cMAC")
  if (!all(required %in% names(o))) stop("Missing required columns in ", key)
  comps <- grep(SB_PATTERN, names(o), value = TRUE)
  for (cc in c("#SNV", "cMAC", comps)) set(o, j = cc, value = as.numeric(o[[cc]]))
  if (any(!is.finite(o$cMAC)) || any(!is.finite(o[["#SNV"]])))
    stop("Nonfinite count metadata in ", key)
  schemas <- o[, .(rows = .N, expected_components = unique(expected_components)),
               by = .(Category, component_schema)]
  fwrite(schemas, file.path(OUT, paste0(key, "_component_schemas.csv")))

  P <- as.matrix(o[, ..comps])
  invalid <- !is.na(P) & (!is.finite(P) | P < 0 | P > 1)
  o[, invalid_components := rowSums(invalid)]
  o[, valid_components := rowSums(!is.na(P))]
  o[, missing_components := expected_components - valid_components]
  o[, raw_zero_components := rowSums(P == 0, na.rm = TRUE)]
  o[, raw_one_components := rowSums(P == 1, na.rm = TRUE)]
  o[, clipped_low_components := rowSums(P < P_FLOOR, na.rm = TRUE)]
  o[, clipped_high_components := rowSums(P > P_CEILING, na.rm = TRUE)]
  problems <- o[invalid_components > 0 | missing_components > 0 | valid_components == 0,
                c(required, "input_file", "invalid_components", "missing_components", "valid_components"), with = FALSE]
  fwrite(problems, file.path(OUT, paste0(key, "_component_problems.csv")))
  if (any(o$invalid_components > 0)) stop("Invalid component probabilities in ", key)
  if (any(o$missing_components > 0))
    stop("Missing expected components require review before extraction: ", key)
  o[, SB := cct_stable(P)]
  o[, SB_legacy := cct_legacy(P)]
  rm(P, data_parts); gc(FALSE)

  # Match the original gene/category identity, but never silently discard
  # conflicting duplicates. Identical duplicates are logged and collapsed.
  dkeys <- o[, .N, by = .(`Gene name`, Category)][N > 1]
  dup <- o[dkeys, on = .(`Gene name`, Category), nomatch = 0]
  fwrite(dup[, c(required, "SB", "SB_legacy", "input_file"), with = FALSE],
         file.path(OUT, paste0(key, "_duplicate_rows.csv")))
  if (nrow(dup)) {
    check_cols <- c("Chr", "#SNV", "cMAC", comps, "component_schema")
    conflicts <- dup[, .(distinct_results = uniqueN(.SD)),
                     by = .(`Gene name`, Category), .SDcols = check_cols][distinct_results > 1]
    if (nrow(conflicts)) {
      fwrite(conflicts, file.path(OUT, paste0(key, "_duplicate_conflicts.csv")))
      stop("Conflicting duplicated gene/category results in ", key)
    }
  }
  before_dedup <- nrow(o)
  o <- unique(o, by = c("Gene name", "Category"))
  o[, tested := cMAC >= 10 & is.finite(SB)]
  o[, `:=`(q_SB = NA_real_, q_SB_legacy = NA_real_)]
  o[tested == TRUE, q_SB := p.adjust(SB, "BH")]
  o[tested == TRUE, q_SB_legacy := p.adjust(SB_legacy, "BH")]
  o[, legacy_log10_error := abs(log10(SB_legacy) - log10(SB))]
  o[, legacy_discovery_disagrees := tested & ((q_SB < 0.05) != (q_SB_legacy < 0.05))]
  setnames(o, "Gene name", "gene")
  keep <- c("gene", "Chr", "Category", "#SNV", "cMAC", "SB", "q_SB", "tested",
            "SB_legacy", "q_SB_legacy", "legacy_log10_error", "legacy_discovery_disagrees",
            "expected_components", "valid_components", "missing_components", "raw_zero_components",
            "raw_one_components", "clipped_low_components", "clipped_high_components", "input_file")
  result <- o[, ..keep]
  fwrite(result, file.path(OUT, paste0(key, ".csv")))
  fwrite(result[SB_legacy == 0 | legacy_log10_error > 1e-4 | legacy_discovery_disagrees == TRUE],
         file.path(OUT, paste0(key, "_legacy_numerical_audit.csv")))
  p <- o[tested == TRUE]$SB
  pl <- o[tested == TRUE]$SB_legacy
  k <- sum(o$q_SB < 0.05, na.rm = TRUE)
  summary <- data.table(arm, family = fam, expected_chunks = expected, observed_chunks = nrow(snapshot),
                        sets_before_dedup = before_dedup, duplicate_rows_removed = before_dedup - nrow(o),
                        sets = nrow(o), tested = length(p), minP = min(p), lambda_GC = lambda_gc(p),
                        BH_05 = k, BH_critical = if (k) 0.05 * k / length(p) else NA_real_,
                        legacy_minP = min(pl), legacy_lambda_GC = lambda_gc(pl),
                        legacy_BH_05 = sum(o$q_SB_legacy < 0.05, na.rm = TRUE),
                        legacy_zero_p = sum(pl == 0),
                        tested_sets_with_raw_zero_components = sum(o[tested == TRUE]$raw_zero_components > 0),
                        tested_sets_with_missing_components = sum(o[tested == TRUE]$missing_components > 0),
                        numerical_audit_rows = sum(result$SB_legacy == 0 | result$legacy_log10_error > 1e-4 |
                                                   result$legacy_discovery_disagrees, na.rm = TRUE),
                        clipping = clip_note)
  all_summaries[[key]] <- summary
  print(summary[, .(arm, sets, tested, minP, lambda_GC, BH_05, legacy_zero_p, legacy_BH_05)])
  fwrite(rbindlist(all_summaries), file.path(OUT, "noncoding_comparison_summary.csv"))
}
writeLines(c(
  "Extraction of existing outputs only; no scans or model fits performed.",
  "Expected complete pooled chunk coverage: 379 regulatory and 222 ncRNA, verified against installed STAARpipeline annotation catalogs.",
  "Only complete pooled regulatory and ncRNA families were extracted. No partial-family BH results were produced.",
  "SB selects every column matching ^(SKAT|Burden)\\(1,(25|1)\\), exactly as gene_set_primary_statistic.R; equal component weights.",
  "tested requires cMAC>=10 and finite SB. Structural differences in component schemas are retained and logged by category.",
  clip_note,
  "SB uses mathematically equivalent cotangent/upper-tail identities. SB_legacy reproduces the older floating-point expression for audit.",
  "Input file counts, indices, hashes, sizes, timestamps, duplicates and component checks are supplied separately.",
  "Primary RIGA tables and result files were neither overwritten nor recalculated."
), file.path(OUT, "gene_extraction_provenance.txt"))
capture.output(sessionInfo(), file = file.path(OUT, "gene_extraction_sessionInfo.txt"))
cat("DONE: complete pooled regulatory and ncRNA families extracted.\n")
