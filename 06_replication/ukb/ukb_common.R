## ukb_common.R (2026-10-19, rev 2): helpers shared by validate.R, compare_chunks.R, scan_audit_rootcheck.R
## and selftest_attribution.R. Read-only with respect to every existing file; sourced, never run.
##
## rev 2 (after review): solver events are attributed to result rows through the gene|category CONTEXT that the patched
## STAARpipeline writes into every log line (STAAR_ROOTCHECK_CONTEXT, set by coding() before each STAAR_Binary_SPA call),
## not through the (#SNV, cMAC, n_weights) key, which is not unique within a chunk (plof == plof_ds whenever a gene has no
## disruptive missense variant) and which rev 1 had wrong for disruptive_missense (26 weight columns, not 2). The
## (#SNV, cMAC) of every attributed event is still verified against its row. The disruptive_missense call of a gene also
## feeds the missense row's "Burden(1,25)-Disruptive"/"Burden(1,1)-Disruptive" columns (coding.R copies weight columns 0 and
## n_weights/2 of that single call), so its events on those two components are attributed to the missense row as well.
suppressPackageStartupMessages(library(data.table))

PKG_DIR  <- "${ADSP_ROOT}/revision/ukb_rootcheck_20261019"
REV      <- "${ADSP_ROOT}/revision"
UKB      <- file.path(REV, "ukb")
ORIG_DIR <- file.path(UKB, "gw_coding_spa", "reconstructed")      # original scan (never written)
NEW_DIR  <- file.path(UKB, "gw_coding_spa", "rootcheck")          # root-checked scan
NTASK    <- 934L
P_COLS_RE <- "^(Burden\\(|STAAR-B)"                                # every p-value column of a coding result row
## order of the STAAR_Binary_SPA calls inside STAARpipeline::coding() (code order). NB: the OUTPUT list order is
## plof, plof_ds, missense, disruptive_missense, synonymous; the calls happen in this order:
CATEGORY_CALL_ORDER <- c("plof_ds", "plof", "synonymous", "missense", "disruptive_missense")
P_FILTER_CUTOFF <- 0.05                                            # SPA_p_filter cutoff of the scan (Gene_Centric_Coding default)
RTOL_IDENTICAL  <- 1e-8                                            # "identical" (bit-identical in practice)
RTOL_POLISH     <- 1e-3                                            # bound for a column whose only events are Newton polishes

chunk_id <- function(f) as.integer(sub(".*_(\\d+)\\.Rdata$", "\\1", basename(f)))

## one chunk file -> long table: chunk, block (position in the list; 5 blocks per gene in the order plof, plof_ds, missense,
## disruptive_missense, synonymous, NULL when the category was not tested), gene, chr, category, nSNV, cMAC, column, col_pos,
## col_kind (weight / aggregate / disruptive_copy), weight_idx (0-based index of the weight column in the C++ call: the
## stored order Burden(1,25), Burden(1,25)-<anno>..., Burden(1,1), Burden(1,1)-<anno>... is the order of weights_B), p
## (character kept verbatim in `p_chr` so that bit-identity can be tested, numeric in `p`)
read_chunk_long <- function(f) {
  e <- new.env()
  ok <- isTRUE(tryCatch({ load(f, envir = e); TRUE }, error = function(z) FALSE))
  if (!ok) return(NULL)
  obj <- get(ls(e)[1], envir = e)
  if (!is.list(obj) || !length(obj)) return(NULL)
  out <- list()
  for (b in seq_along(obj)) {
    z <- obj[[b]]
    if (is.null(z) || !is.matrix(z) || nrow(z) < 1) next
    cn <- colnames(z)
    pc <- cn[grepl(P_COLS_RE, cn)]
    kind <- ifelse(grepl("^STAAR-B", pc), "aggregate", ifelse(grepl("-Disruptive$", pc), "disruptive_copy", "weight"))
    widx <- rep(NA_integer_, length(pc)); widx[kind == "weight"] <- seq_len(sum(kind == "weight")) - 1L
    for (r in seq_len(nrow(z))) {
      vals <- unlist(lapply(pc, function(cc) { v <- z[r, cc][[1]]; if (length(v)) v[1] else NA }))
      out[[length(out) + 1]] <- data.table(block = b, gene = as.character(z[r, "Gene name"][[1]]), chr = as.integer(z[r, "Chr"][[1]]),
                                         category = as.character(z[r, "Category"][[1]]), nSNV = as.integer(z[r, "#SNV"][[1]]),
                                         cMAC = as.numeric(z[r, "cMAC"][[1]]), column = pc, col_pos = seq_along(pc), col_kind = kind,
                                         weight_idx = widx, p_chr = as.character(vals), p = suppressWarnings(as.numeric(vals)))
    }
  }
  res <- rbindlist(Filter(Negate(is.null), out))
  if (!nrow(res)) return(NULL)
  res[, chunk := chunk_id(f)]
  res[]
}

read_scan_long <- function(dir, chunks = NULL) {
  fs <- list.files(dir, pattern = "^UKB_AD_Coding_\\d+\\.Rdata$", full.names = TRUE)
  if (!is.null(chunks)) fs <- fs[chunk_id(fs) %in% chunks]
  rbindlist(Filter(Negate(is.null), lapply(fs, read_chunk_long)))
}

## relative difference used throughout: |a-b| / max(|a|, |b|, 1e-300); NA handling explicit
rel_diff <- function(a, b) {
  d <- abs(a - b) / pmax(abs(a), abs(b), 1e-300)
  d[is.na(a) & is.na(b)] <- 0
  d[xor(is.na(a), is.na(b))] <- Inf
  d
}

## join an original and a new long table on (chunk, gene, category, column); duplicates (a gene tested twice in a chunk is
## impossible, but block order could differ) are resolved by matching the k-th occurrence
compare_long <- function(old, new) {
  old <- copy(old); new <- copy(new)
  old[, occ := seq_len(.N), by = .(chunk, gene, category, column)]
  new[, occ := seq_len(.N), by = .(chunk, gene, category, column)]
  m <- merge(old[, .(chunk, gene, chr, category, column, occ, col_kind_old = col_kind, weight_idx_old = weight_idx, nSNV_old = nSNV, cMAC_old = cMAC, p_old_chr = p_chr, p_old = p)],
             new[, .(chunk, gene, category, column, occ, col_kind_new = col_kind, weight_idx_new = weight_idx, nSNV_new = nSNV, cMAC_new = cMAC, p_new_chr = p_chr, p_new = p)],
             by = c("chunk", "gene", "category", "column", "occ"), all = TRUE)
  m[, col_kind := fifelse(is.na(col_kind_new), col_kind_old, col_kind_new)]
  m[, weight_idx := fifelse(is.na(weight_idx_new), weight_idx_old, weight_idx_new)]
  m[, c("col_kind_old", "col_kind_new", "weight_idx_old", "weight_idx_new") := NULL]
  m[, identical_chr := !is.na(p_old_chr) & !is.na(p_new_chr) & p_old_chr == p_new_chr]
  m[, rel := rel_diff(p_old, p_new)]
  m[, status := fifelse(is.na(p_old_chr) & !is.na(p_new_chr), "new_only",
               fifelse(!is.na(p_old_chr) & is.na(p_new_chr), "old_only",
               fifelse(identical_chr | rel <= RTOL_IDENTICAL, "identical_1e-8",
               fifelse(rel <= 1e-6, "changed_le_1e-6",
               fifelse(rel <= 1e-3, "changed_le_1e-3", "changed_material")))))]
  m[]
}

## ---- diagnostic event logs written by the patched library (STAAR_ROOTCHECK_LOG) ------------------------------------
## 17 tab-separated columns (rev 2; rev-1 logs had 15 and are read with context = "" and n_spa = NA):
##   pkg, call_seq, n_variants, cMAC, n_weights, component, p_normal, tail, event, q, x_newton, resid_newton, x_final,
##   resid_final, iter, context, n_spa
## event lines: one per solver EVENT (newton_nonfinite/newton_maxiter/newton_unverified/newton_polished,
##   bisection_ok/bisection_ok_expanded, bisection_nonfinite/bisection_nobracket/bisection_unverified);
## call_done lines: one per burden call (every STAAR_B_*SPA* caller writes it before returning) with n_spa = number of
##   components of the call that entered the saddlepoint branch; context = "<gene>|<category>" set by STAARpipeline.
LOG_COLS <- c("pkg", "call_seq", "n_variants", "cMAC", "n_weights", "component", "p_normal", "tail", "event",
              "q", "x_newton", "resid_newton", "x_final", "resid_final", "iter", "context", "n_spa")
LOG_NUM  <- c("call_seq", "n_variants", "cMAC", "n_weights", "component", "p_normal", "q", "x_newton", "resid_newton", "x_final", "resid_final", "iter", "n_spa")
FAILURE_EVENTS <- c("newton_nonfinite", "newton_maxiter", "newton_unverified")
BISECT_FAIL    <- c("bisection_nonfinite", "bisection_nobracket", "bisection_unverified")
BISECT_OK      <- c("bisection_ok", "bisection_ok_expanded")

read_one_log <- function(f, chunk = NA_integer_) {
  if (!file.exists(f) || file.size(f) == 0) return(NULL)
  d <- fread(f, header = FALSE, sep = "\t", fill = TRUE, colClasses = "character", na.strings = character(0), quote = "")
  if (ncol(d) == 15) { d[, V16 := ""]; d[, V17 := "NA"] }           # rev-1 log format (no context, no call_done lines)
  if (ncol(d) != 17) stop("unexpected log format (", ncol(d), " columns): ", f)
  setnames(d, LOG_COLS)
  for (cc in LOG_NUM) set(d, j = cc, value = suppressWarnings(as.numeric(d[[cc]])))
  d[, component := as.integer(component)]
  d[, chunk := as.integer(chunk)]
  d[, log_file := basename(f)]
  d[, log_line := seq_len(.N)]
  d
}

read_root_logs <- function(dir, files = NULL) {
  fs <- if (is.null(files)) list.files(dir, pattern = "^rootcheck_\\d+\\.tsv$", full.names = TRUE) else files
  if (!length(fs)) return(data.table())
  d <- rbindlist(Filter(Negate(is.null), lapply(fs, function(f) read_one_log(f, chunk = suppressWarnings(as.integer(sub("^rootcheck_(\\d+)\\.tsv$", "\\1", basename(f))))))))
  if (!nrow(d)) return(d)
  d[is.na(context), context := ""]
  has_ctx <- nzchar(d$context) & grepl("|", d$context, fixed = TRUE)
  d[, gene := fifelse(has_ctx, sub("\\|.*$", "", context), NA_character_)]
  d[, category := fifelse(has_ctx, sub("^[^|]*\\|", "", context), NA_character_)]
  d[, kind := fifelse(event == "call_done", "call_done",
              fifelse(event %in% FAILURE_EVENTS, "newton_failure",
              fifelse(event == "newton_polished", "newton_polished",
              fifelse(event %in% BISECT_FAIL, "bisection_failure",
              fifelse(event %in% BISECT_OK, "bisection", "other")))))]
  d[]
}

## ---- attribution of events to result rows --------------------------------------------------------------------------
## Exact: by the gene|category context of the log line (one STAAR_Binary_SPA call per (gene, category) in coding()).
## attribution = "set"                       event belongs to the row (gene, category) of the run that wrote the log
##               "disruptive_copy"           copy of a disruptive_missense event (weight column 0 or n_weights/2) onto the
##                                           missense row's Burden(1,25)-Disruptive / Burden(1,1)-Disruptive column
##               "no_context"                log line without a gene|category context (never in a coding() run)
##               "unmatched_gene_absent"     the gene has no row at all in the run (gene skipped by the driver after the call)
##               "unmatched_category_absent" the gene has rows but not this category (should not happen)
##               "column_unresolved"         the component index has no weight column in the row (should not happen)
## key_ok: the event's (n_variants, cMAC) equal the row's (#SNV, cMAC) -- verified for every attributed event.
events_by_set <- function(log, new_long) {
  empty <- data.table(chunk = integer(), gene = character(), category = character(), column = character(), kind = character(), attribution = character(), key_ok = logical(), call_seq = numeric(), component = integer(), tail = character(), event = character(), p_normal = numeric())
  if (!nrow(log)) return(empty)
  ev <- copy(log[pkg == "STAAR" & event != "call_done"])
  if (!nrow(ev)) return(empty)
  sets <- unique(new_long[, .(chunk, gene, category, nSNV_set = nSNV, cMAC_set = cMAC)])
  genes_present <- unique(sets[, .(chunk, gene)])[, gene_present := TRUE]
  wmap <- unique(new_long[col_kind == "weight", .(chunk, gene, category, component = weight_idx, column)])
  ev <- merge(ev, sets, by = c("chunk", "gene", "category"), all.x = TRUE, sort = FALSE)
  ev <- merge(ev, genes_present, by = c("chunk", "gene"), all.x = TRUE, sort = FALSE)
  ev[, attribution := "set"]
  ev[is.na(gene) | is.na(category), attribution := "no_context"]
  ev[attribution == "set" & is.na(nSNV_set) & is.na(gene_present), attribution := "unmatched_gene_absent"]
  ev[attribution == "set" & is.na(nSNV_set) & !is.na(gene_present), attribution := "unmatched_category_absent"]
  ev[, key_ok := NA]
  ev[attribution == "set", key_ok := n_variants == nSNV_set & abs(cMAC - cMAC_set) <= 1e-6 * pmax(1, abs(cMAC_set))]
  ev[, via_disruptive_copy := FALSE]
  ## the disruptive_missense call also feeds the missense row (coding.R: Burden(1,25)-Disruptive <- weight column 0,
  ## Burden(1,1)-Disruptive <- weight column n_weights/2 of that call; NA -> 1 there)
  cp <- ev[attribution == "set" & category == "disruptive_missense" & (component == 0L | component == as.integer(n_weights %/% 2))]
  if (nrow(cp)) {
    cp <- copy(cp)
    cp[, column_copy := fifelse(component == 0L, "Burden(1,25)-Disruptive", "Burden(1,1)-Disruptive")]
    cp[, c("category", "nSNV_set", "cMAC_set") := NULL]
    cp[, category := "missense"]
    cp <- merge(cp, sets[category == "missense"], by = c("chunk", "gene", "category"), sort = FALSE)     # only genes WITH a missense row
    cp[, `:=`(attribution = "disruptive_copy", via_disruptive_copy = TRUE)]
  }
  ev <- merge(ev, wmap, by = c("chunk", "gene", "category", "component"), all.x = TRUE, sort = FALSE)
  ev[attribution == "set" & is.na(column), attribution := "column_unresolved"]
  if (nrow(cp)) { cp[, column := column_copy]; cp[, column_copy := NULL]; ev <- rbind(ev, cp, fill = TRUE) }
  setorder(ev, chunk, call_seq, component, tail, via_disruptive_copy)
  ev[]
}

## ambiguous (#SNV, cMAC) keys within a chunk: the reason rev 1's key-based attribution was not exact (reported only)
ambiguous_keys <- function(new_long) {
  s <- unique(new_long[, .(chunk, gene, category, nSNV, cMAC_key = round(cMAC, 6))])
  k <- s[, .(n_sets = .N, sets = paste(paste(gene, category, sep = ":"), collapse = ";")), by = .(chunk, nSNV, cMAC_key)][n_sets > 1]
  k[order(chunk, nSNV)]
}

## the call_done lines: one per burden call; row presence, key check, SPA component counts and the call-order diagnostic
calls_summary <- function(log, new_long) {
  cd <- if (nrow(log)) copy(log[pkg == "STAAR" & event == "call_done"]) else data.table()
  res <- list(calls = cd, n_calls = nrow(cd), n_spa_components = if (nrow(cd)) sum(cd$n_spa, na.rm = TRUE) else 0L,
              n_calls_no_context = 0L, n_calls_without_row = 0L, n_calls_key_mismatch = 0L, n_call_order_anomalies = 0L, n_rows_without_call = NA_integer_)
  if (!nrow(cd)) return(res)
  sets <- unique(new_long[, .(chunk, gene, category, nSNV_set = nSNV, cMAC_set = cMAC, block)])
  cd <- merge(cd, sets, by = c("chunk", "gene", "category"), all.x = TRUE, sort = FALSE)
  cd[, has_row := !is.na(nSNV_set)]
  cd[, key_ok := has_row & n_variants == nSNV_set & abs(cMAC - cMAC_set) <= 1e-6 * pmax(1, abs(cMAC_set))]
  res$n_calls_no_context <- sum(is.na(cd$gene))
  res$n_calls_without_row <- sum(!cd$has_row & !is.na(cd$gene))
  res$n_calls_key_mismatch <- sum(cd$has_row & !cd$key_ok)
  ## rows of the run that have no call_done line (possible only with a rev-1 log, which has no call_done lines)
  rows <- unique(new_long[, .(chunk, gene, category)])
  res$n_rows_without_call <- nrow(rows[!cd[has_row == TRUE, .(chunk, gene, category)], on = c("chunk", "gene", "category")])
  ## call order: genes in block order (= the driver's gene order), categories in CATEGORY_CALL_ORDER
  gene_rank <- sets[, .(grank = min(block)), by = .(chunk, gene)]
  cd <- merge(cd, gene_rank, by = c("chunk", "gene"), all.x = TRUE, sort = FALSE)
  cd[, crank := match(category, CATEGORY_CALL_ORDER)]
  setorder(cd, chunk, call_seq)
  cd[, expected_pos := grank * 10 + crank]
  cd[, order_anomaly := FALSE]
  cd[has_row == TRUE, order_anomaly := c(FALSE, diff(expected_pos) <= 0), by = chunk]
  res$n_call_order_anomalies <- sum(cd$order_anomaly, na.rm = TRUE)
  res$calls <- cd[]
  res
}

## touched columns: every (chunk, gene, category, column) with at least one attributed event, and the aggregates they feed:
## STAAR-B(1,25) combines the Burden(1,25)* columns (on the missense row incl. Burden(1,25)-Disruptive), STAAR-B(1,1) the
## Burden(1,1)* columns, STAAR-B all of them. touched_by_failure = at least one Newton failure / bisection / bisection
## failure on the column (its value may change freely); otherwise the column was only polished (change <= RTOL_POLISH).
touched_columns <- function(ev) {
  empty <- data.table(chunk = integer(), gene = character(), category = character(), column = character(), n_events = integer(), touched_by_failure = logical(), touch_kinds = character())
  if (!nrow(ev)) return(empty)
  t <- ev[attribution %in% c("set", "disruptive_copy") & !is.na(column),
          .(n_events = .N, touched_by_failure = any(kind %in% c("newton_failure", "bisection", "bisection_failure")), touch_kinds = paste(sort(unique(kind)), collapse = ";")),
          by = .(chunk, gene, category, column)]
  if (!nrow(t)) return(empty)
  aggr <- function(sub, name) if (nrow(sub)) sub[, .(column = name, n_events = sum(n_events), touched_by_failure = any(touched_by_failure), touch_kinds = paste(sort(unique(unlist(strsplit(touch_kinds, ";", fixed = TRUE)))), collapse = ";")), by = .(chunk, gene, category)] else NULL
  agg <- rbindlist(list(aggr(t[startsWith(column, "Burden(1,25)")], "STAAR-B(1,25)"), aggr(t[startsWith(column, "Burden(1,1)")], "STAAR-B(1,1)"), aggr(t, "STAAR-B")))
  rbind(t, agg[, names(t), with = FALSE])[]
}

fmt_p <- function(x) ifelse(is.na(x), "NA", formatC(x, digits = 6, format = "g"))
