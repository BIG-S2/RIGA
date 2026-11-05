#!/usr/bin/env Rscript
## 2026-09-18 fix: the match-set column is named mkey (data.table(key=) is the constructor's key argument and made setkeyv fail on "KEY0").
###############################################################################
## benchmark_matching.R -- DELIVERABLE B MatchIt pre-step (the RIGA / pair33 recipe on a pool
## whose AD label may be the real AD or a plasmode Y*).
##
## Data prep is copied VERBATIM from riga_design_harness.R L64-87 (itself verbatim from
## src/prepare_sample_and_one_to_one_matching.R), then restricted to the rows of the plasmode
## pool (outputs/plasmode/plasmode_pool.csv); do_match is riga_design_harness.R L96-111
## verbatim (matchit nearest, distance = "mahalanobis" on PC1-10, exact = key,
## replace = FALSE, ratio = 1, set.seed(0)).
##
## Keys: KEY0 = Sex + age_band + SeqCenter + lpa_cluster   (deployed pair33 key)
##       KEY3 = KEY0 + Cohort                               (riga D3)
##       KEY_P = KEY0 + Study                               (GCAD project label)
## NOTE the SeqCenter used in the key is the riga-merged 11-level factor (label of the
## KCCG level is the long Kinghorn string there; identical partition to pool$center).
##
## Usage:
##   Rscript benchmark_matching.R --pool <plasmode_pool.csv> --out <dir> --tag <tag>
##          [--ystar <csv with columns IID,Y>]     # Y* replaces AD; default = real AD
##          [--keys KEY0,KEY3,KEY_P]                # subset of keys
##          [--check-deployed]                      # real-Y reproduction check (V7)
## Output: <out>/match_idx_<tag>.csv   long format, one row per matched subject:
##             key, pair_id, IID, row_idx (0-based row in plasmode_pool.csv), fam_pos, Y
##         <out>/match_summary_<tag>.csv  key, n_rows_in, n_case, n_control, n_pairs, minutes
##         (+ <out>/match_check_<tag>.json when --check-deployed)
###############################################################################
suppressPackageStartupMessages({
  library(data.table); library(dplyr); library(MatchIt); library(readr); library(tidyr)
})
options(stringsAsFactors = FALSE, width = 230, scipen = 6, warn = 1)
REPO <- "${ADSP_ROOT}"
SEED_MATCH <- 0L
tick <- function(x) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), x))

## ---- args --------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
getopt <- function(flag, default = NULL) { i <- match(flag, args); if (is.na(i)) default else args[i + 1] }
POOL  <- getopt("--pool", file.path(REPO, "outputs/plasmode/plasmode_pool.csv"))
OUTD  <- getopt("--out", file.path(REPO, "outputs/plasmode/match"))
TAG   <- getopt("--tag", "realY")
YSTAR <- getopt("--ystar", NULL)
KEYS  <- strsplit(getopt("--keys", "KEY0,KEY3,KEY_P"), ",")[[1]]
CHECK <- "--check-deployed" %in% args
dir.create(OUTD, showWarnings = FALSE, recursive = TRUE)
cat(sprintf("MatchIt %s | R %s | pool=%s | tag=%s | ystar=%s | keys=%s | check=%s\n",
            packageVersion("MatchIt"), getRversion(), POOL, TAG, ifelse(is.null(YSTAR), "<real AD>", YSTAR),
            paste(KEYS, collapse = ","), CHECK))

###############################################################################
## 1. data prep -- VERBATIM from riga_design_harness.R L64-87 (absolute paths only)
###############################################################################
set.seed(0)
cov_wgs <- read_csv("${ADSP_ROOT}/pheno/adsp_combined.csv",
                    show_col_types = FALSE) %>% drop_na(AD, Sex, Age)
cov_wgs$Age <- as.numeric(gsub("\\+", "", cov_wgs$Age))
clu <- fread(file.path(REPO, "adsp_lpa_clusters_20PC.tsv")); setDT(clu)
K <- length(grep("^post_k[0-9]+$", names(clu)))
if (!"entropy_norm" %in% names(clu)) clu[, entropy_norm := entropy / log(K)]
cat(sprintf("NO confidence filter: keeping all %d clustered subjects (the filter would keep %d).\n",
            nrow(clu), sum(clu$cluster_prob >= 0.7 & clu$entropy_norm <= 0.5)))
seqmeta <- read_csv("${ADSP_ROOT}/pheno/gcad.r5.wgs.58507.2024.11.03.quality.metrics.ALL.csv",
                    show_col_types = FALSE) %>%
  transmute(IID = Sample_name, SeqCenter = as.character(SeqCenter))
dat0 <- merge(cov_wgs, clu, by = "IID") %>%
  left_join(seqmeta, by = "IID") %>%
  mutate(Sex = as.factor(Sex), AD = as.factor(AD), SeqCenter = as.character(SeqCenter))
dat0$SeqCenter[is.na(dat0$SeqCenter)] <- "Unknown"
dat0$SeqCenter[dat0$SeqCenter %in% c(
  "MEDGENOME","Unknown","Genentech","GENENTECH",
  "HudsonAlpha Institute for Biotechnology","Macrogen","MACROGEN")] <- "Other"
dat0$SeqCenter <- droplevels(as.factor(dat0$SeqCenter))
dat0 <- dat0 %>% tidyr::drop_na(Age, Sex, AD, SeqCenter) %>%
  mutate(age_band = cut(Age, breaks = seq(30, 100, by = 3), right = TRUE)) %>%
  tidyr::drop_na(age_band)
dat0$lpa_cluster  <- as.factor(dat0$lpa_cluster)
dat0$entropy_norm <- if (K > 1) dat0$entropy / log(K) else 0
cat(sprintf("Pool entering matching: %d  (cases %d / controls %d)\n",
            nrow(dat0), sum(dat0$AD == 1), sum(dat0$AD == 0)))
dfx <- as.data.frame(dat0)
PCS <- paste0("PC", 1:10)
## ---- end of verbatim block ----------------------------------------------
stopifnot(nrow(dfx) == 38129, "Cohort" %in% names(dfx), !anyDuplicated(dfx$IID))
dfx$IID <- as.character(dfx$IID)

###############################################################################
## 2. restrict to the plasmode pool; attach Study / row_idx / fam_pos; optional Y*
###############################################################################
pool <- fread(POOL, colClasses = list(character = c("IID", "Study", "Cohort")))
stopifnot(all(c("IID", "row_idx", "fam_pos", "Study", "AD", "center") %in% names(pool)))
n_pool <- nrow(pool)
dfx <- dfx[dfx$IID %in% pool$IID, ]
mi <- match(dfx$IID, pool$IID)
dfx$Study   <- pool$Study[mi]
dfx$row_idx <- pool$row_idx[mi]
dfx$fam_pos <- pool$fam_pos[mi]
## consistency: the riga-merged SeqCenter partition == pool$center partition
sc_map <- table(as.character(dfx$SeqCenter), pool$center[mi])
stopifnot(all(rowSums(sc_map > 0) == 1), all(colSums(sc_map > 0) == 1))
cat(sprintf("pool rows %d; entering matching after the verbatim prep: %d (dropped %d: age_band NA or unclustered)\n",
            n_pool, nrow(dfx), n_pool - nrow(dfx)))
stopifnot(all(as.integer(as.character(dfx$AD)) == pool$AD[mi]))
if (!is.null(YSTAR)) {
  ys <- fread(YSTAR, colClasses = list(character = "IID"))
  stopifnot(all(c("IID", "Y") %in% names(ys)), all(dfx$IID %in% ys$IID))
  yv <- ys$Y[match(dfx$IID, ys$IID)]
  stopifnot(all(yv %in% c(0L, 1L)))
  dfx$AD <- factor(yv, levels = c(0, 1))
  cat(sprintf("Y* applied: cases %d / controls %d\n", sum(yv == 1), sum(yv == 0)))
}
dfx$Cohort <- factor(dfx$Cohort)
dfx$Study  <- factor(dfx$Study)

###############################################################################
## 3. matching -- riga_design_harness.R L96-111 verbatim (do_match)
###############################################################################
mform <- as.formula(paste("AD ~", paste(PCS, collapse = " + ")))
do_match <- function(key, ratio, tag) {
  set.seed(SEED_MATCH); t0 <- Sys.time()
  m.out <- matchit(mform, data = dfx, method = "nearest", distance = "mahalanobis",
                   exact = reformulate(key), replace = FALSE, ratio = ratio)
  md <- match.data(m.out, data = dfx)
  idx <- match(as.character(md$IID), dfx$IID); stopifnot(!anyNA(idx))
  st <- rep(NA_character_, nrow(dfx)); st[idx] <- as.character(md$subclass)
  sz <- table(table(st))
  tick(sprintf("%s: key=%s ratio=%d -> %d rows, %d strata, stratum sizes {%s} (%.1f min)", tag,
               paste(key, collapse = "+"), ratio, sum(!is.na(st)), uniqueN(na.omit(st)),
               paste(sprintf("%s:%d", names(sz), as.integer(sz)), collapse = ", "),
               as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  st
}
KEY0  <- c("Sex", "age_band", "SeqCenter", "lpa_cluster")
KEY3  <- c(KEY0, "Cohort")
KEY_P <- c(KEY0, "Study")
KEYLIST <- list(KEY0 = KEY0, KEY3 = KEY3, KEY_P = KEY_P)
stopifnot(all(KEYS %in% names(KEYLIST)))

out <- list(); summ <- list()
for (kn in KEYS) {
  t0 <- Sys.time()
  st <- do_match(KEYLIST[[kn]], 1L, kn)
  sel <- !is.na(st)
  out[[kn]] <- data.table(mkey = kn, pair_id = st[sel], IID = dfx$IID[sel], row_idx = dfx$row_idx[sel],
                          fam_pos = dfx$fam_pos[sel], Y = as.integer(as.character(dfx$AD[sel])))
  summ[[kn]] <- data.table(mkey = kn, n_rows_in = nrow(dfx), n_case = sum(dfx$AD == 1), n_control = sum(dfx$AD == 0),
                           n_matched = sum(sel), n_pairs = uniqueN(st[sel]),
                           minutes = round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 2))
}
out <- rbindlist(out); summ <- rbindlist(summ)
setorder(out, mkey, pair_id)
fwrite(out, file.path(OUTD, sprintf("match_idx_%s.csv", TAG)))
fwrite(summ, file.path(OUTD, sprintf("match_summary_%s.csv", TAG)))
print(summ)

###############################################################################
## 4. --check-deployed: reproduction of the deployed pair33 pairs (riga L138-145 style),
##    restricted to the unrelated pool
###############################################################################
if (CHECK) {
  stopifnot(is.null(YSTAR), "KEY0" %in% KEYS)
  prev <- file.path(REPO, "output/matched/matched_samples_lpa_pc.nofilter_pair.csv")
  pv <- fread(prev, select = c("IID", "subclass"), colClasses = list(character = c("IID", "subclass")))
  pv_pool <- pv[IID %in% pool$IID]
  cur <- out[mkey == "KEY0"]
  ov_subj <- length(intersect(pv_pool$IID, cur$IID))
  ## deployed pairs with BOTH members in the pool, and how many of those are reproduced as pairs
  pp <- pv_pool[, .N, by = subclass][N == 2, subclass]
  dep_pairs <- pv_pool[subclass %in% pp][order(subclass, IID)][, .(a = IID[1], b = IID[2]), by = subclass]
  cur_pairs <- cur[order(pair_id, IID)][, .(a = IID[1], b = IID[2]), by = pair_id]
  key_dep <- paste(dep_pairs$a, dep_pairs$b); key_cur <- paste(cur_pairs$a, cur_pairs$b)
  same_pair <- sum(key_dep %in% key_cur)
  chk <- list(deployed_rows = nrow(pv), deployed_rows_in_pool = nrow(pv_pool),
              deployed_pairs_both_in_pool = nrow(dep_pairs),
              rebuilt_rows = nrow(cur), rebuilt_pairs = nrow(cur_pairs),
              overlap_subjects = ov_subj,
              overlap_pct_of_deployed_in_pool = 100 * ov_subj / nrow(pv_pool),
              overlap_pct_of_rebuilt = 100 * ov_subj / nrow(cur),
              deployed_pairs_reproduced_exactly = same_pair,
              deployed_pairs_reproduced_pct = 100 * same_pair / nrow(dep_pairs),
              pass_99pct = (100 * ov_subj / nrow(pv_pool)) >= 99)
  cat(sprintf("D0 check vs deployed pair33 file restricted to the pool: deployed %d (in pool %d, %d complete pairs), rebuilt %d (%d pairs), overlap %d subjects (%.2f%% of deployed-in-pool, %.2f%% of rebuilt), %d deployed pairs reproduced exactly (%.2f%%)\n",
              chk$deployed_rows, chk$deployed_rows_in_pool, chk$deployed_pairs_both_in_pool, chk$rebuilt_rows,
              chk$rebuilt_pairs, ov_subj, chk$overlap_pct_of_deployed_in_pool, chk$overlap_pct_of_rebuilt,
              same_pair, chk$deployed_pairs_reproduced_pct))
  writeLines(jsonlite::toJSON(chk, auto_unbox = TRUE, pretty = TRUE), file.path(OUTD, sprintf("match_check_%s.json", TAG)))
  cat(sprintf("CHECK-DEPLOYED %s\n", ifelse(chk$pass_99pct, "PASS", "BELOW-99PCT (reported, see json)")))
}
tick("done")
