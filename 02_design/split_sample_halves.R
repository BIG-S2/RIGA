#!/usr/bin/env Rscript
## split_sample_halves.R (2026-10-08): held-out split of the 30,894 subjects of the stratified design D_S into two
## halves A and B. Relatives never straddle the split: relatedness clusters are the connected components of (i) the
## nonzero off-diagonal pattern of the sparse GRM that the null model uses and (ii) declared family membership in the
## ADSP family-based phenotype file (FamID within cohort, Father/Mother links). Whole clusters are allocated greedily,
## in random order, to the half that best balances case and control counts within every project x center cell.
## Cells left with a single label in a half are dropped from that half (no within-cell contrast).
suppressPackageStartupMessages({ library(data.table); library(Matrix) })
ROOT <- "${ADSP_ROOT}"; OUT <- "${ADSP_ROOT}/output/matched"
TAB <- file.path(ROOT, "revision/tables/supp_v3"); SEED <- as.integer(Sys.getenv("SEED", "20261008")); set.seed(SEED)
m <- fread(file.path(OUT, "matched_samples_lpa_pc.pair33_projcell.csv"), colClasses = list(character = c("IID", "FID", "SUBJID.x", "SUBJID.y")))
stopifnot(nrow(m) == 30894); n <- nrow(m); par <- seq_len(n)
find <- function(i) { while (par[i] != i) { par[i] <<- par[par[i]]; i <- par[i] }; i }
unite <- function(i, j) { ri <- find(i); rj <- find(j); if (ri != rj) par[rj] <<- ri; invisible(NULL) }
## (i) sparse GRM: pairs related at third degree or closer. The matrix is on the kinship scale (diagonal 0.50, declared
## parent-offspring pairs at 0.25) and keeps near-zero entries inside its blocks, so "any nonzero entry" would chain 3,095
## subjects into one cluster; the edge rule is kinship > KIN_THRESH (default 0.0442, the KING third-degree boundary).
KIN_THRESH <- as.numeric(Sys.getenv("KIN_THRESH", "0.0442"))
sgrm <- get(load("${ADSP_ROOT}/processed/output.sparseGRM.sGRM.RData"))
gid <- rownames(sgrm); if (is.null(gid)) gid <- colnames(sgrm); gid_s <- sub("^[^_]*_", "", gid)
idx <- match(m$IID, gid_s); stopifnot(!any(is.na(idx))); K <- as(as(sgrm[idx, idx], "generalMatrix"), "TsparseMatrix"); rm(sgrm)
ed <- data.table(i = K@i + 1L, j = K@j + 1L, x = K@x)[i < j]
cat(sprintf("GRM off-diagonal pairs: %d; related at > %.4f: %d; at > 0.0221 (fourth degree): %d; > 0.354 (MZ/duplicate): %d, same SUBJID: %d\n", nrow(ed), KIN_THRESH, sum(ed$x > KIN_THRESH), sum(ed$x > 0.0221), sum(ed$x > 0.354), sum(ed$x > 0.354 & m$SUBJID.x[ed$i] == m$SUBJID.x[ed$j])))
ed4 <- ed[x > 0.0221]; ed <- ed[x > KIN_THRESH]; ii <- ed$i; jj <- ed$j
for (k in seq_along(ii)) unite(ii[k], jj[k])
## (ii) declared families
fam <- fread(file.path(ROOT, "pheno/ADSPFamilyBasedPhenotypes_DS_2024.11.22_ALL.csv"), select = c("SUBJID", "Cohort", "FamID", "Father", "Mother"), colClasses = "character")
rowof <- function(ids) { r <- match(ids, m$SUBJID.x); for (col in c("SUBJID.y", "IID")) { r2 <- match(ids, m[[col]]); r[is.na(r)] <- r2[is.na(r)] }; r }
fam[, row := rowof(SUBJID)]; f <- fam[!is.na(row)]; f[, fkey := paste(Cohort, FamID, sep = "|")]
cat("family-file subjects in D_S:", nrow(f), " families:", uniqueN(f$fkey), "\n")
for (rs in split(f$row, f$fkey)) if (length(rs) > 1) for (k in 2:length(rs)) unite(rs[1], rs[k])
for (p in c("Father", "Mother")) { pr <- rowof(f[[p]]); for (k in which(!is.na(pr))) unite(f$row[k], pr[k]) }
## (iii) FID column of the design file, used only if it behaves like a family identifier
fid <- m$FID; fid[is.na(fid) | fid == "" | fid == m$IID | fid == "0"] <- NA; tb <- table(fid)
cat("FID groups of size >1:", sum(tb > 1), " largest FID group:", if (length(tb)) max(tb) else 0L, "\n")
if (length(tb) && max(tb) <= 50) for (rs in split(which(!is.na(fid)), fid[!is.na(fid)])) if (length(rs) > 1) for (k in 2:length(rs)) unite(rs[1], rs[k])
comp <- vapply(seq_len(n), function(i) as.integer(find(i)), integer(1)); m[, cluster := match(comp, unique(comp))]
cs <- m[, .N, by = cluster]; cat("relatedness clusters:", nrow(cs), " singletons:", sum(cs$N == 1), " largest:", max(cs$N), "\n")
## greedy balanced allocation of whole clusters, random order, balancing (cell, AD) counts
m[, stratum := paste(cell_SC, AD)]; S <- sort(unique(m$stratum)); nS <- length(S); sidx <- match(m$stratum, S)
cl <- split(seq_len(n), m$cluster); cntA <- cntB <- integer(nS); half <- character(n)
for (c in sample(length(cl))) { rs <- cl[[c]]; v <- tabulate(sidx[rs], nS)
  dA <- sum(abs(cntA + v - cntB)); dB <- sum(abs(cntA - cntB - v))
  h <- if (dA < dB) "A" else if (dB < dA) "B" else sample(c("A", "B"), 1)
  half[rs] <- h; if (h == "A") cntA <- cntA + v else cntB <- cntB + v }
stopifnot(all(half[ed$i] == half[ed$j])); m[, half := half]
cat("fourth-degree (kinship > 0.0221) pairs straddling the halves:", sum(half[ed4$i] != half[ed4$j]), "of", nrow(ed4), "\n")
for (rs in split(f$row, f$fkey)) stopifnot(uniqueN(half[rs]) == 1)
ct <- m[, .(N = .N, cases = sum(AD == 1)), by = .(half, cell_SC)]; ct[, mixed := cases > 0 & cases < N]
m[ct, mixed := i.mixed, on = c("half", "cell_SC")]
cat("single-label (half, cell) pairs dropped:", sum(!ct$mixed), " subjects dropped:", sum(!m$mixed), "\n")
sm <- m[mixed == TRUE, .(N = .N, cases = sum(AD == 1), controls = sum(AD == 0), cells = uniqueN(cell_SC), clusters = uniqueN(cluster),
                         subjects_in_clusters_gt1 = sum(cluster %in% cs[N > 1]$cluster)), by = half][order(half)]
print(sm); fwrite(sm, file.path(TAB, "S14_split_summary.csv"))
cc <- dcast(m[mixed == TRUE, .(N = .N, cases = sum(AD == 1)), by = .(cell_SC, half)], cell_SC ~ half, value.var = c("N", "cases"))
fwrite(cc, file.path(TAB, "S14_split_cells.csv"))
for (h in c("A", "B")) fwrite(m[half == h & mixed == TRUE], file.path(OUT, sprintf("matched_samples_projcell_split%s.csv", h)))
cat("DONE seed", SEED, "\n")
