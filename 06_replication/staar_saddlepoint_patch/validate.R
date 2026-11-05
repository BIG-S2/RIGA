#!/usr/bin/env Rscript
## validate.R (2026-10-19, rev 2): item 2 of ukb_rootcheck_20261019. Runs with the ROOT-CHECKED library (R_LIBS_USER
## set by jobs/02_validate.sbatch) after the unit test (2d) and the three chunk reruns + compare_chunks.R (2b, 2c) have
## produced their files in <outdir>. This script
##   (a) recomputes the three reference sets from the EXACT masks cached by the 18 October targeted run (same null, same
##       G, same annotations, same disruptive subset, same Cauchy combination) and, independently, through the pipeline
##       path Gene_Centric_Coding(category = ...) on the aGDS; both must agree with the reference values to 1e-6 relative:
##       SNX22 synonymous 0.0518593609546751, TREM2 missense 2.10583109606644e-05, NECTIN2 missense 0.000714815559568427;
##   (e) checks that the chunk reruns reproduce those values for the same sets (chunks 676, 326, 849);
##   then assembles VALIDATION.txt and the machine-readable gate file gate/ROOTCHECK_GATE.json ("pass": true/false).
## Usage: Rscript validate.R <outdir>      (outdir = results/validation)
suppressPackageStartupMessages({ library(data.table); library(jsonlite); library(Matrix); library(SeqArray); library(STAAR); library(STAARpipeline) })
source("${ADSP_ROOT}/revision/ukb_rootcheck_20261019/scripts/ukb_common.R")
args <- commandArgs(TRUE); outdir <- args[1]; stopifnot(!is.na(outdir)); dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
GATE_DIR <- file.path(PKG_DIR, "gate"); dir.create(GATE_DIR, showWarnings = FALSE)
t0 <- Sys.time(); checks <- list()
chk <- function(name, ok, detail = "") { checks[[length(checks) + 1]] <<- list(check = name, pass = isTRUE(ok), detail = detail); cat(sprintf("  [%s] %s %s\n", if (isTRUE(ok)) "PASS" else "FAIL", name, detail)) }
REF <- data.table(gene = c("SNX22", "TREM2", "NECTIN2"), chr = c(15L, 6L, 19L), category = c("synonymous", "missense", "missense"), chunk = c(676L, 326L, 849L),
                  ref_staarb = c(0.0518593609546751, 2.10583109606644e-05, 0.000714815559568427), expected_n = c(68L, 156L, 324L), expected_cMAC = c(944, 9383, 11656))
SAFE <- file.path(REV, "five_analysis_execution_20261016/snx22/targeted_safeguard_20261018")
Sys.setenv(STAAR_ROOTCHECK_LOG = file.path(outdir, "validate_rootlog.tsv")); if (file.exists(Sys.getenv("STAAR_ROOTCHECK_LOG"))) file.remove(Sys.getenv("STAAR_ROOTCHECK_LOG"))

## ---- 0. the library in use ----------------------------------------------------------------------------------------
d1 <- packageDescription("STAAR"); d2 <- packageDescription("STAARpipeline")
chk("library: STAAR loaded from the root-checked library", !is.null(d1$RootcheckPatch) && grepl("staar_rootcheck_4.4", system.file(package = "STAAR")), system.file(package = "STAAR"))
chk("library: STAARpipeline loaded from the root-checked library", !is.null(d2$RootcheckPatch) && grepl("staar_rootcheck_4.4", system.file(package = "STAARpipeline")), system.file(package = "STAARpipeline"))
build_ok <- file.exists(file.path(PKG_DIR, "build/BUILD_OK")); chk("library: BUILD_OK record present", build_ok)
for (g in REF$gene) { s <- fread(file.path(SAFE, g, "summary.csv")); chk(sprintf("reference file for %s carries the expected value", g), isTRUE(all.equal(s$guarded_STAAR_B, REF[gene == g, ref_staarb], tolerance = 1e-12)), sprintf("%.15g", s$guarded_STAAR_B)) }

## ---- (d) unit test results ---------------------------------------------------------------------------------------
ut <- file.path(outdir, "unit_test_checks_patched.csv")
if (file.exists(ut)) { u <- fread(ut); chk("(d) synthetic unit test: every check passed with the patched library", all(u$pass), sprintf("%d/%d", sum(u$pass), nrow(u)))
  chk("(d) unit test covered the four C++ callers and the real sparse-kins null (cases B4, G, H present)", any(grepl("^B4:", u$check)) && any(grepl("^G:", u$check)) && any(grepl("^H: real null", u$check)), sprintf("B4 %d, G %d, H %d checks", sum(grepl("^B4:", u$check)), sum(grepl("^G:", u$check)), sum(grepl("^H:", u$check))))
} else chk("(d) synthetic unit test ran", FALSE, "unit_test_checks_patched.csv missing")
uo <- file.path(outdir, "unit_test_checks_original.csv")
if (file.exists(uo)) { u <- fread(uo); chk("(d) synthetic unit test: the defect is reproduced with the original library", all(u$pass), sprintf("%d/%d", sum(u$pass), nrow(u))) }

## ---- (a) exact masks through STAAR_Binary_SPA, and the pipeline path ----------------------------------------------
obj <- get(load(file.path(UKB, "ukb_obj_nullmodel_reconstructed.RData"))); obj$sparse_kins <- TRUE; obj$relatedness <- TRUE; obj$use_SPA <- TRUE
stopifnot(isTRUE(obj$converged), length(obj$id_include) == 182238L)
cat("null model: N =", length(obj$id_include), " cases =", sum(round(obj$fitted.values + obj$scaled.residuals) == 1), "\n")
combine_exact <- function(p) { stopifnot(!is.null(names(p)), all(is.finite(p))); cb <- function(x) { x <- x[x < 1]; if (!length(x) || sum(x) <= 0) 1 else CCT(x) }
  c("STAAR-B(1,25)" = cb(p[startsWith(names(p), "Burden(1,25)")]), "STAAR-B(1,1)" = cb(p[startsWith(names(p), "Burden(1,1)")]), "STAAR-B" = cb(p)) }
load(file.path(UKB, "agds_dir_ukb.RData")); catalog <- get(load("${ADSP_ROOT}/output/Annotation_name_catalog.Rdata"))
AN <- c("CADD", "LINSIGHT", "FATHMM.XF", "aPC.EpigeneticActive", "aPC.EpigeneticRepressed", "aPC.EpigeneticTranscription", "aPC.Conservation", "aPC.LocalDiversity", "aPC.Mappability", "aPC.TF", "aPC.Protein")
flat <- function(r) { if (is.null(r) || !length(r)) return(NULL); m <- if (is.matrix(r) || is.data.frame(r)) r else r[[1]]; if (is.null(m) || !length(m)) return(NULL)
  d <- as.data.frame(m, stringsAsFactors = FALSE, check.names = FALSE); for (cc in names(d)) if (is.list(d[[cc]])) d[[cc]] <- unlist(lapply(d[[cc]], function(z) if (length(z)) z[1] else NA)); d }
aval <- list()
for (i in seq_len(nrow(REF))) {
  g <- REF$gene[i]; cat("\n== reference set", g, REF$category[i], "==\n")
  if (g == "SNX22") { cache <- readRDS(file.path(REV, "five_analysis_execution_20261016/snx22/exact_mask_private.rds")); G <- cache$G; A <- cache$annotations; disruptive <- NULL
  } else { cache <- readRDS(file.path(SAFE, g, "exact_mask_private.rds")); G <- cache$G; A <- cache$A; disruptive <- cache$disruptive }
  chk(sprintf("(a) %s exact mask: %d variants, cMAC %.0f", g, REF$expected_n[i], REF$expected_cMAC[i]), ncol(G) == REF$expected_n[i] && abs(sum(G) - REF$expected_cMAC[i]) < 1e-6, sprintf("%d / %.3f", ncol(G), sum(G)))
  r <- STAAR_Binary_SPA(G, obj, A, rare_maf_cutoff = 0.01, rv_num_cutoff = 2, SPA_p_filter = TRUE, p_filter_cutoff = 0.05)
  p25 <- unlist(r$results_STAAR_B_1_25); p11 <- unlist(r$results_STAAR_B_1_1)
  p <- c(p25[startsWith(names(p25), "Burden(")], p11[startsWith(names(p11), "Burden(")])
  if (!is.null(disruptive)) { dm <- STAAR_Binary_SPA(G[, disruptive, drop = FALSE], obj, annotation_phred = NULL, rare_maf_cutoff = 0.01, rv_num_cutoff = 2, SPA_p_filter = TRUE, p_filter_cutoff = 0.05)
    p <- c(p, "Burden(1,25)-Disruptive" = as.numeric(dm$results_STAAR_B_1_25[1, "Burden(1,25)"]), "Burden(1,1)-Disruptive" = as.numeric(dm$results_STAAR_B_1_1[1, "Burden(1,1)"])) }
  chk(sprintf("(a) %s: no NA component with the patched solver", g), all(is.finite(p)), paste(names(p)[!is.finite(p)], collapse = ","))
  agg <- combine_exact(p[is.finite(p)]); staarb <- unname(agg["STAAR-B"]); rel <- abs(staarb - REF$ref_staarb[i]) / REF$ref_staarb[i]
  chk(sprintf("(a) %s exact-mask STAAR-B agrees with the reference to 1e-6 relative", g), rel <= 1e-6, sprintf("%.15g vs %.15g (rel %.2e)", staarb, REF$ref_staarb[i], rel))
  ## component-level agreement with the guarded recomputation (informative: the guarded run used a different root finder)
  cg <- fread(file.path(SAFE, g, "components_guarded.csv")); cmpc <- merge(data.table(component = names(p), patched_p = as.numeric(p)), cg[, .(component, original_p, guarded_p, branch)], by = "component")
  cmpc[, rel_guarded := abs(patched_p - guarded_p) / pmax(guarded_p, 1e-300)]; cmpc[, rel_original := abs(patched_p - original_p) / pmax(original_p, 1e-300)]
  fwrite(cmpc, file.path(outdir, sprintf("reference_%s_components.csv", g)))
  chk(sprintf("(a) %s: every component agrees with the guarded recomputation to 1e-5 relative", g), all(cmpc$rel_guarded <= 1e-5), sprintf("max rel %.2e (%s)", max(cmpc$rel_guarded), cmpc$component[which.max(cmpc$rel_guarded)]))
  cat(sprintf("   components changed vs the ORIGINAL run (>1e-6 rel): %s\n", paste(cmpc[rel_original > 1e-6, sprintf("%s %.3g->%.3g", component, original_p, patched_p)], collapse = "; ")))
  ## pipeline path on the aGDS
  gf <- seqOpen(agds_dir[REF$chr[i]])
  pr <- tryCatch(Gene_Centric_Coding(chr = REF$chr[i], gene_name = g, category = REF$category[i], genofile = gf, obj_nullmodel = obj, rare_maf_cutoff = 0.01, rv_num_cutoff = 2,
                                     QC_label = "annotation/info/QC_label2", variant_type = "SNV", geno_missing_imputation = "mean", Annotation_dir = "annotation/info/FunctionalAnnotation",
                                     Annotation_name_catalog = catalog, Use_annotation_weights = TRUE, Annotation_name = AN), error = function(e) { cat("ERR pipeline", conditionMessage(e), "\n"); NULL })
  seqClose(gf); d <- flat(pr); pipe_b <- if (is.null(d)) NA_real_ else as.numeric(d[["STAAR-B"]][1]); rel2 <- abs(pipe_b - REF$ref_staarb[i]) / REF$ref_staarb[i]
  chk(sprintf("(a) %s pipeline path (Gene_Centric_Coding) STAAR-B agrees with the reference to 1e-6 relative", g), isTRUE(rel2 <= 1e-6), sprintf("%.15g (rel %.2e; #SNV %s cMAC %s)", pipe_b, rel2, if (is.null(d)) NA else d[["#SNV"]][1], if (is.null(d)) NA else d$cMAC[1]))
  aval[[g]] <- data.table(gene = g, category = REF$category[i], reference = REF$ref_staarb[i], exact_mask_patched = staarb, pipeline_patched = pipe_b, rel_exact = rel, rel_pipeline = rel2)
}
fwrite(rbindlist(aval), file.path(outdir, "reference_sets_patched.csv"))

## ---- (b), (c): the chunk reruns, compared by compare_chunks.R ------------------------------------------------------
sj <- file.path(outdir, "chunks_status.json")
if (file.exists(sj)) {
  s <- fromJSON(sj)
  chk("(b) three chunks rerun and compared", length(s$chunks) == 3 && all(c(676, 326, 849) %in% s$chunks), paste(s$chunks, collapse = ","))
  ## attribution (rev 2): by the gene|category context of every log line, verified against the row's (#SNV, cMAC)
  chk("(b) attribution: every solver event carries a gene|category context (patched STAARpipeline)", s$n_events_no_context == 0, sprintf("%d of %d events without context", s$n_events_no_context, s$n_events))
  chk("(b) attribution: every event's (gene, category) has a result row in the rerun", s$n_events_unmatched_category_absent == 0, sprintf("%d category rows absent; %d events in genes skipped by the driver (informative)", s$n_events_unmatched_category_absent, s$n_events_unmatched_gene_absent))
  chk("(b) attribution: #SNV and cMAC of every attributed event equal those of its row", s$n_events_key_mismatch == 0 && s$n_calls_key_mismatch == 0, sprintf("%d event / %d call mismatches over %d events, %d calls", s$n_events_key_mismatch, s$n_calls_key_mismatch, s$n_events_attributed, s$n_calls))
  chk("(b) attribution: every event's component index maps to a result column", s$n_events_column_unresolved == 0, sprintf("%d unresolved", s$n_events_column_unresolved))
  chk("(b) attribution: every result row has a call_done line and no call lacks a row", isTRUE(s$n_rows_without_call == 0) && s$n_calls_without_row == 0 && s$n_calls_no_context == 0, sprintf("%s rows without call, %d calls without row, %d calls without context; %d burden calls, %d SPA components", format(s$n_rows_without_call), s$n_calls_without_row, s$n_calls_no_context, s$n_calls, s$n_spa_components))
  cat(sprintf("   informative: %d ambiguous (#SNV, cMAC) keys covering %d sets (not used for attribution); %d call-order anomalies; %d disruptive->missense event copies\n", s$n_ambiguous_keys, s$n_sets_in_ambiguous_keys, s$n_call_order_anomalies, s$n_events_disruptive_copies))
  chk("(b1) every set without a solver event is identical to the original to 1e-8 relative", isTRUE(s$b1_no_event_sets_identical), sprintf("%d violations over %d event-free sets", s$b1_violations, s$n_sets_without_event))
  chk("(b2) sets/columns with Newton polishing only changed by at most 1e-3 relative", isTRUE(s$b2_polish_only_sets_le_1e3) && s$c_polish_columns_exceeding_1e3 == 0, sprintf("%d set violations, %d polish-only columns above 1e-3", s$b2_violations, s$c_polish_columns_exceeding_1e3))
  chk("(b3) every material change is explained by a Newton failure / bisection / NA event on the set", isTRUE(s$b3_material_changes_explained), sprintf("%d material-change sets, %d violations", s$n_sets_material_change, s$b3_violations))
  chk("(c) per component: every column without a solver event of its own is identical to 1e-8 (inside event sets too)", isTRUE(s$c_untouched_columns_identical), sprintf("%d violations over %d untouched columns (%d inside event sets); %d touched columns, %d of them changed", s$c_violations, s$c_untouched_columns_checked, s$c_untouched_columns_in_event_sets, s$n_touched_columns, s$n_touched_columns_changed))
  chk("(c2) every logged event entered the saddlepoint branch through the filter (p_normal < 0.05)", isTRUE(s$c2_events_entered_spa_through_filter), sprintf("%d violations", s$c2_violations))
  chk("(b) no component left NA by the patched solver in the three chunks", s$n_components_new_NA == 0, sprintf("%d NA components", s$n_components_new_NA))
  sets <- fread(file.path(outdir, "chunks_sets.csv"))
  for (i in seq_len(nrow(REF))) { r <- sets[chunk == REF$chunk[i] & gene == REF$gene[i] & category == REF$category[i]]
    rel <- if (nrow(r) == 1) abs(r$staarb_new - REF$ref_staarb[i]) / REF$ref_staarb[i] else NA
    chk(sprintf("(e) chunk %d rerun: %s %s STAAR-B equals the reference to 1e-6 relative", REF$chunk[i], REF$gene[i], REF$category[i]), isTRUE(rel <= 1e-6), sprintf("%s -> %s (ref %.15g, rel %s)", fmt_p(r$staarb_old[1]), fmt_p(r$staarb_new[1]), REF$ref_staarb[i], format(rel))) }
} else chk("(b) chunk comparison present", FALSE, "chunks_status.json missing")

## ---- gate ----------------------------------------------------------------------------------------------------------
ck <- rbindlist(lapply(checks, as.data.table)); pass <- all(ck$pass)
gate <- list(pass = pass, timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), host = Sys.info()[["nodename"]], library = system.file(package = "STAAR"),
             staar_rootcheck = d1$RootcheckPatch, staarpipeline_rootcheck = d2$RootcheckPatch, n_checks = nrow(ck), n_pass = sum(ck$pass), failed = ck[pass == FALSE, check],
             reference_sets = rbindlist(aval), validation_dir = outdir)
write_json(gate, file.path(GATE_DIR, "ROOTCHECK_GATE.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA)
fwrite(ck, file.path(outdir, "validation_checks.csv"))
vt <- c(sprintf("VALIDATION of the root-checked STAAR/STAARpipeline library  (%s)", format(Sys.time())), sprintf("library: %s", system.file(package = "STAAR")),
        sprintf("STAAR: %s", d1$RootcheckPatch), sprintf("STAARpipeline: %s", d2$RootcheckPatch), "",
        sprintf("RESULT: %s  (%d of %d checks passed)", if (pass) "PASS" else "FAIL", sum(ck$pass), nrow(ck)), "",
        sprintf("%-6s %s  %s", ifelse(ck$pass, "PASS", "FAIL"), ck$check, ck$detail), "",
        "Reference sets (STAAR-B): reference / exact-mask recomputation / pipeline recomputation",
        sprintf("  %-8s %-10s %.15g / %.15g / %.15g", rbindlist(aval)$gene, rbindlist(aval)$category, rbindlist(aval)$reference, rbindlist(aval)$exact_mask_patched, rbindlist(aval)$pipeline_patched), "",
        if (file.exists(file.path(outdir, "chunks_summary.txt"))) c("Chunk comparison (compare_chunks.R):", readLines(file.path(outdir, "chunks_summary.txt"))) else "chunk comparison missing",
        "", sprintf("elapsed %.1f min", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
writeLines(vt, file.path(PKG_DIR, "VALIDATION.txt")); writeLines(vt, file.path(outdir, "VALIDATION.txt"))
cat("\n", paste(vt[1:8], collapse = "\n"), "\n\ngate:", file.path(GATE_DIR, "ROOTCHECK_GATE.json"), " pass =", pass, "\n")
if (!pass) quit(status = 1)
