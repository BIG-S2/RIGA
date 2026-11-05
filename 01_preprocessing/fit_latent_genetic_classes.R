# LPA / Gaussian‑mixture clustering of ADSP subjects by genetic PCs
# ------------------------------------------------------------------
# Author: Xiaoqi Li
# Date: 2025‑07‑15
#
# Workflow ("train‑on‑unrelated, project‑relatives")
#   1. Identify the mutually unrelated subjects (≈80 %).
#   2. Fit a Gaussian‑mixture/LPA model *only* on that subset.
#   3. Freeze the class‑specific parameters and obtain posterior class
#      probabilities for every subject — related or not — via predict().
#   4. Write out assignments, QC plots, and a race × cluster table.
#   5. ✱ NEW ✱ Harmonise sample IDs in the unrelated list to the WGS IDs
#      (substring matching) and save the corrected list for future runs.
# ------------------------------------------------------------------

# 1 . Setup ────────────────────────────────────────────────────────────

suppressPackageStartupMessages({
  library(data.table)   # fast I/O
  library(mclust)       # Gaussian‑mixture / LPA
  library(dplyr)        # table manipulation
  library(ggplot2)      # quick visual QC
  library(scales)       # for percent_format()
  library(tidyr)        # pivoting tables
  library(readr)        # read_csv()
  library(mclust)
  library(mvtnorm)
})

# ---------- User‑configurable paths ---------------------------------
cov_wgs <- read_csv("${ADSP_ROOT}/pheno/adsp_combined.csv")
unrel_file   <- "${ADSP_ROOT}/processed/output.pruned.unrelated.unrels"  # list of unrelated IIDs (raw)

num_pcs      <- 20                # number of PCs to use
candidate_K  <- 1:20              # mixture sizes to try
lambda_diag  <- 0.1              # diagonal shrinkage λ

# Variable names ----------------------------------------------------
race_var <- "Race"                # update if different
site_var <- "Cohort"                # genotyping site column
cov_wgs <- cov_wgs %>% mutate(Race = as.factor(Race), Cohort = as.factor(Cohort))


# Output files ------------------------------------------------------
cluster_out       <- sprintf("adsp_lpa_clusters_%sPC.tsv", num_pcs)
bic_png           <- sprintf("bic_curve_%sPC_unrelated.png", num_pcs)
icl_png           <- sprintf("icl_curve_%sPC_unrelated.png", num_pcs)
pc_plot_png      <- sprintf("pc1_pc2_clusters_%sPC.png", num_pcs)
race_table_tsv    <- sprintf("cluster_vs_race_%sPC.tsv", num_pcs)
site_plot_png     <- sprintf("cluster_by_site_%sPC.png", num_pcs)
site_table_tsv    <- sprintf("cluster_vs_site_%sPC.tsv", num_pcs)
fixed_unrel_file  <- sub("\\.unrels$", ".matched.unrels", unrel_file)
cluster_prob_png  <- sprintf("cluster_prob_dist_%sPC.png", num_pcs)
entropy_png       <- sprintf("entropy_dist_%sPC.png", num_pcs)
uncertainty_png   <- sprintf("cluster_uncertainty_%sPC.png", num_pcs)
uncertainty_tsv   <- sprintf("cluster_uncertainty_summary_%sPC.tsv", num_pcs)


rowLogSumExp <- function(M) {
  m <- apply(M, 1, max)
  m + log(rowSums(exp(M - m)))
}

free_params_fullcov <- function(K, d) {
  (K - 1) + K * d + K * d * (d + 1) / 2
}
# ------------------------------------------------------------------
# 2 . Load raw unrelated list ──────────────
# Raw unrelated IDs (may use chip IDs, shortened forms, etc.)
unrel_raw <- fread(unrel_file, header = FALSE)$V2  # one IID per line
cat("Loaded", length(unrel_raw), "raw unrelated IIDs from", unrel_file, "\n")

# -------------------------- Harmonisation ----------------------------------------
# 2b. Harmonise: keep WGS IDs whose string appears inside raw list ----
unrel_ids <- cov_wgs$IID[cov_wgs$IID %in% unrel_raw]  # keep only WGS IDs

cat(sprintf("After harmonisation: %d matched unrelated IDs (of %d WGS IDs).\n",
            length(unrel_ids), length(cov_wgs$IID)))

# Save the harmonised list for future runs
writeLines(unrel_ids, fixed_unrel_file)
cat("Corrected unrelated ID list written to:", fixed_unrel_file, "\n")

# ------------------------------------------------------------------
# 3 . Merge phenotype and mark unrelated -----------------------------
cov_wgs <- cov_wgs %>% mutate(unrel = IID %in% unrel_ids)
cat(sprintf("Total subjects: %d  |  Unrelated (matched): %d  |  Related: %d\n",
            nrow(cov_wgs), sum(cov_wgs$unrel), sum(!cov_wgs$unrel)))

# ------------------------------------------------------------------
# 4 . Prepare PC matrix ───────────────────────────────────────────────

pc_cols <- paste0("PC", 1:num_pcs)
if (!all(pc_cols %in% names(cov_wgs))) {
  stop("Requested PC columns not found: ", paste(setdiff(pc_cols, names(cov_wgs)), collapse = ", "))
}

# convert to data.table
setDT(cov_wgs)
# Scale PCs (mean 0, SD 1) – important for mixture fit
cov_wgs[, (pc_cols) := lapply(.SD, scale), .SDcols = pc_cols]

pc_mat       <- as.matrix(cov_wgs[, ..pc_cols])                   # full matrix
pc_mat_unrel <- as.matrix(cov_wgs[unrel == TRUE, ..pc_cols])      # training subset
# ------------------------------------------------------------------
# 5 . Fit Gaussian-mixture (LPA) with shrinkage on unrelated set ----
# ------------------------------------------------------------------

# # Use mclust's priorControl to apply shrinkage Σ_k ← (1-λ)Σ_k + λI
# prior_ctrl <- priorControl(shrinkage = lambda_diag)

# bic_obj <- mclustBIC(
#   pc_mat_unrel, G = candidate_K, modelNames = "VVV",
#   prior = prior_ctrl
# )

# # Fit the chosen model: note **do not** pass `prior` again when `x` (BIC)
# # is supplied; otherwise mclust throws an error.
# best_model <- Mclust(pc_mat_unrel, x = bic_obj)
# K_opt      <- best_model$G
# --- Max-Entropy (Deterministic-Annealing) GMM with model selection by BIC ---
library(mvtnorm)
library(mclust)
library(mvtnorm)

rowLogSumExp <- function(M) {
  m <- apply(M, 1, max)
  m + log(rowSums(exp(M - m)))
}

free_params_fullcov <- function(K, d) {
  (K - 1) + K * d + K * d * (d + 1) / 2
}

# ------------------------------------------------------------
# MBHAC initializer (mclust-style): hc() + hclass()
# ------------------------------------------------------------
init_mbhac_params <- function(X, K,
                              modelName = "VVV",
                              lambda_diag = 0,
                              jitter = 1e-6,
                              hc_obj = NULL) {
  N <- nrow(X); d <- ncol(X)

  # pooled covariance as a safe fallback
  Xc0 <- sweep(X, 2, colMeans(X), "-")
  S0  <- crossprod(Xc0) / N

  if (is.null(hc_obj)) {
    hc_obj <- mclust::hc(X, modelName = modelName)
  }
  cl <- mclust::hclass(hc_obj, G = K)

  # responsibilities (hard)
  r <- matrix(0, nrow = N, ncol = K)
  r[cbind(seq_len(N), cl)] <- 1

  Nk <- pmax(colSums(r), 1e-12)
  pi_k <- as.numeric(Nk / N)

  mu <- matrix(0, nrow = K, ncol = d)
  Sigmas <- vector("list", K)

  for (k in 1:K) {
    w <- r[, k]  # 0/1
    mu[k, ] <- colSums(X * w) / Nk[k]

    # If a cluster is tiny, use pooled covariance for stability
    if (Nk[k] < 2) {
      Sk <- S0
    } else {
      Xc_k <- sweep(X, 2, mu[k, ], "-")
      Sk <- (t(Xc_k) %*% (Xc_k * w)) / Nk[k]
    }

    # optional shrinkage towards identity
    if (lambda_diag > 0) {
      Sk <- (1 - lambda_diag) * Sk + lambda_diag * diag(d)
    }

    Sigmas[[k]] <- Sk + jitter * diag(d)
  }

  list(pi = pi_k, mu = mu, Sigma = Sigmas,
       responsibilities = r,
       classification = cl,
       hc_obj = hc_obj)
}

# ------------------------------------------------------------
# DA-EM GMM with MBHAC initialization
# ------------------------------------------------------------
gmm_da_em <- function(X, K,
                      tau_schedule = seq(5, 1, length.out = 12),
                      max_iter = 200, tol = 1e-6,
                      jitter = 1e-6,
                      lambda_diag = 0,
                      mbhac_modelName = "VVV",
                      hc_obj = NULL) {
  N <- nrow(X); d <- ncol(X)

  init <- init_mbhac_params(
    X = X, K = K,
    modelName = mbhac_modelName,
    lambda_diag = lambda_diag,
    jitter = jitter,
    hc_obj = hc_obj
  )
  mu     <- init$mu
  pi_k   <- init$pi
  Sigmas <- init$Sigma
  r      <- init$responsibilities

  loglik_tau1 <- NA

  for (tau in tau_schedule) {
    last_obj <- -Inf
    for (iter in 1:max_iter) {
      # --- E-step (tempered responsibilities) ---
      L <- sapply(1:K, function(k) {
        mvtnorm::dmvnorm(X, mean = mu[k, ], sigma = Sigmas[[k]], log = TRUE) + log(pi_k[k])
      })
      L_tau <- L / tau
      lse <- rowLogSumExp(L_tau)
      log_r <- sweep(L_tau, 1, lse, "-")
      r <- pmax(exp(log_r), 1e-16)
      r <- r / rowSums(r)

      # --- M-step ---
      Nk <- pmax(colSums(r), 1e-12)
      pi_k <- as.numeric(Nk / N)

      for (k in 1:K) {
        w <- r[, k]
        mu[k, ] <- colSums(X * w) / Nk[k]

        Xc_k <- sweep(X, 2, mu[k, ], "-")
        Sk <- (t(Xc_k) %*% (Xc_k * w)) / Nk[k]

        if (lambda_diag > 0) {
          Sk <- (1 - lambda_diag) * Sk + lambda_diag * diag(d)
        }
        Sigmas[[k]] <- Sk + jitter * diag(d)
      }

      obj <- sum(lse)
      if (abs(obj - last_obj) <= tol * (1 + abs(last_obj))) break
      last_obj <- obj
    }

    # true loglik at tau=1 for BIC
    loglik_tau1 <- sum(rowLogSumExp(L))
  }

  list(G = K, pi = pi_k, mu = mu, Sigma = Sigmas,
       responsibilities = r, loglik = loglik_tau1)
}

# ------------------------------------------------------------
# Model selection wrapper: compute MBHAC tree once and reuse
# ------------------------------------------------------------
select_gmm_da <- function(X, K_grid,
                          tau_schedule = seq(5, 1, length.out = 12),
                          max_iter = 200, tol = 1e-6,
                          jitter = 1e-6,
                          lambda_diag = 0,
                          mbhac_modelName = "VVV") {
  N <- nrow(X); d <- ncol(X)

  hc_obj <- mclust::hc(X, modelName = mbhac_modelName)

  fits <- vector("list", length(K_grid))
  bics <- numeric(length(K_grid))
  icls <- numeric(length(K_grid))

  for (i in seq_along(K_grid)) {
    K <- K_grid[i]
    fit <- gmm_da_em(
      X = X, K = K,
      tau_schedule = tau_schedule,
      max_iter = max_iter, tol = tol,
      jitter = jitter,
      lambda_diag = lambda_diag,
      mbhac_modelName = mbhac_modelName,
      hc_obj = hc_obj
    )

    z <- pmax(fit$responsibilities, 1e-12)
    p <- free_params_fullcov(K, d)
    bic <- 2 * fit$loglik - p * log(N)         # larger is better
    icl <- bic + sum(z * log(z))               # mclust-style ICL

    fits[[i]] <- fit
    bics[i] <- bic
    icls[i] <- icl
  }

  best_idx <- which.max(bics)
  list(best_model = fits[[best_idx]],
       K_opt = fits[[best_idx]]$G,
       bic_table = data.frame(K = K_grid, BIC = bics),
       icl_table = data.frame(K = K_grid, ICL = icls))
}

# -----------------------------
# Usage (same as your pipeline)
# -----------------------------
lambda_diag <- get0("lambda_diag", ifnotfound = 0)

sel <- select_gmm_da(
  X = pc_mat_unrel,
  K_grid = candidate_K,
  tau_schedule = seq(5, 1, length.out = 12),
  max_iter = 200, tol = 1e-6,
  jitter = 1e-6,
  lambda_diag = lambda_diag,
  mbhac_modelName = "VVV"
)

best_model <- sel$best_model
K_opt      <- sel$K_opt
bic_obj    <- sel$bic_table
icl_obj    <- sel$icl_table
cat(sprintf("\nOptimal number of latent classes (by BIC with shrinkage): %d\n", K_opt))

# ------------------------------------------------------------------
# 6 . Project all subjects through the fixed model -------------------
# ------------------------------------------------------------------
# --- Predict function for DA-EM GMM (keeps params frozen) ---
# Uses the same fields we returned earlier: best_model$pi, $mu, $Sigma, $G
# Returns mclust-like list: $classification and $z
if (!exists("rowLogSumExp")) {
  rowLogSumExp <- function(M) {
    m <- apply(M, 1, max)
    m + log(rowSums(exp(M - m)))
  }
}

pred_gmm_da <- function(model, newdata) {
  K <- model$G
  mu <- model$mu
  Sigmas <- model$Sigma
  pi_k <- model$pi

  # log π_k N(x | μ_k, Σ_k)
  L <- sapply(1:K, function(k) {
    mvtnorm::dmvnorm(newdata, mean = mu[k, ], sigma = Sigmas[[k]], log = TRUE) + log(pi_k[k])
  })

  lse   <- rowLogSumExp(L)
  log_r <- sweep(L, 1, lse, "-")
  z     <- pmax(exp(log_r), 1e-16)     # posterior responsibilities
  z     <- z / rowSums(z)              # safety renormalization

  classification <- max.col(z, ties.method = "first")
  list(classification = classification, z = z, loglik = sum(lse))
}

# Predict posteriors for all rows in pc_mat (parameters frozen)
pred_all <- pred_gmm_da(best_model, newdata = pc_mat)

cov_wgs[, lpa_cluster := pred_all$classification]

z_mat <- pred_all$z
K_opt <- ncol(z_mat)  # just to be safe
post_names <- paste0("post_k", seq_len(K_opt))

z_dt <- data.table::as.data.table(z_mat)
data.table::setnames(z_dt, post_names)

# Stable entropy (avoid log(0))
z_safe  <- pmax(z_mat, 1e-12)
entropy <- -rowSums(z_safe * log(z_safe))

cov_wgs[, cluster_prob := apply(z_mat, 1, max)]
cov_wgs <- cbind(cov_wgs, z_dt)
cov_wgs[, entropy := entropy]
# Normalized entropy [0, 1] for interpretability (1 = max uncertainty)
cov_wgs[, entropy_norm := entropy / log(K_opt)]
cov_wgs[, uncertain_cluster := (cluster_prob < 0.7) | (entropy_norm > 0.5)]

cat("\nCluster assignments and posterior probabilities computed for all subjects.\n")

# ------------------------------------------------------------------
# 6b . Summarise and plot cluster assignment uncertainty -------------
# ------------------------------------------------------------------
cat("\n--- Cluster assignment uncertainty ---\n")
cat("cluster_prob: ", paste(names(summary(cov_wgs$cluster_prob)), collapse = "\t"), "\n")
cat("cluster_prob: ", paste(summary(cov_wgs$cluster_prob), collapse = "\t"), "\n")
cat("entropy:      ", paste(names(summary(cov_wgs$entropy)), collapse = "\t"), "\n")
cat("entropy:      ", paste(summary(cov_wgs$entropy), collapse = "\t"), "\n")
cat("entropy_norm: ", paste(names(summary(cov_wgs$entropy_norm)), collapse = "\t"), "\n")
cat("entropy_norm: ", paste(summary(cov_wgs$entropy_norm), collapse = "\t"), "\n")
n_low_prob   <- sum(cov_wgs$cluster_prob < 0.7, na.rm = TRUE)
n_high_ent   <- sum(cov_wgs$entropy_norm > 0.5, na.rm = TRUE)
n_uncertain  <- sum(cov_wgs$uncertain_cluster, na.rm = TRUE)
cat(sprintf("Subjects with cluster_prob < 0.7: %d (%.1f%%)\n", n_low_prob, 100 * n_low_prob / nrow(cov_wgs)))
cat(sprintf("Subjects with entropy_norm > 0.5: %d (%.1f%%)\n", n_high_ent, 100 * n_high_ent / nrow(cov_wgs)))
cat(sprintf("Flagged as uncertain (prob<0.7 or entropy_norm>0.5): %d (%.1f%%)\n", n_uncertain, 100 * n_uncertain / nrow(cov_wgs)))

uncertainty_summary <- data.table(
  metric = c("cluster_prob_min", "cluster_prob_q25", "cluster_prob_median", "cluster_prob_q75", "cluster_prob_max",
             "entropy_min", "entropy_q25", "entropy_median", "entropy_q75", "entropy_max",
             "entropy_norm_median", "n_low_prob", "n_high_entropy", "n_uncertain", "n_total"),
  value = c(
    min(cov_wgs$cluster_prob), quantile(cov_wgs$cluster_prob, 0.25), median(cov_wgs$cluster_prob),
    quantile(cov_wgs$cluster_prob, 0.75), max(cov_wgs$cluster_prob),
    min(cov_wgs$entropy), quantile(cov_wgs$entropy, 0.25), median(cov_wgs$entropy),
    quantile(cov_wgs$entropy, 0.75), max(cov_wgs$entropy),
    median(cov_wgs$entropy_norm), n_low_prob, n_high_ent, n_uncertain, nrow(cov_wgs)
  )
)
fwrite(uncertainty_summary, uncertainty_tsv, sep = "\t")
cat("Uncertainty summary written to:", uncertainty_tsv, "\n")

# Plots: cluster_prob and entropy distributions, and scatter
gg_prob <- ggplot(cov_wgs, aes(x = cluster_prob)) +
  geom_histogram(bins = 40, fill = "steelblue", alpha = 0.8) +
  geom_vline(xintercept = 0.7, linetype = "dashed", colour = "red") +
  labs(title = "Distribution of max cluster probability", subtitle = "Red line: threshold 0.7 for confident assignment",
       x = "cluster_prob", y = "Count") +
  theme_minimal(base_size = 12)
ggsave(filename = cluster_prob_png, plot = gg_prob, width = 5, height = 3.5, dpi = 300)
cat("Cluster prob plot written to:", cluster_prob_png, "\n")

gg_ent <- ggplot(cov_wgs, aes(x = entropy)) +
  geom_histogram(bins = 40, fill = "darkgreen", alpha = 0.6) +
  labs(title = "Distribution of assignment entropy", subtitle = "Higher = more uncertain",
       x = "entropy", y = "Count") +
  theme_minimal(base_size = 12)
ggsave(filename = entropy_png, plot = gg_ent, width = 5, height = 3.5, dpi = 300)
cat("Entropy plot written to:", entropy_png, "\n")

entropy_thresh <- 0.5 * log(K_opt)
gg_unc <- ggplot(cov_wgs, aes(x = cluster_prob, y = entropy)) +
  geom_point(alpha = 0.3, size = 0.8) +
  geom_hline(yintercept = entropy_thresh, linetype = "dashed", colour = "orange") +
  geom_vline(xintercept = 0.7, linetype = "dashed", colour = "red") +
  annotate("rect", xmin = -Inf, xmax = 0.7, ymin = -Inf, ymax = Inf, alpha = 0.08, fill = "red") +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = entropy_thresh, ymax = Inf, alpha = 0.08, fill = "orange") +
  labs(title = "Cluster uncertainty: entropy vs cluster_prob", subtitle = "Flagged region: prob < 0.7 or entropy_norm > 0.5",
       x = "cluster_prob", y = "entropy") +
  theme_minimal(base_size = 12)
ggsave(filename = uncertainty_png, plot = gg_unc, width = 5.5, height = 4, dpi = 300)
cat("Uncertainty scatter plot written to:", uncertainty_png, "\n")

# ------------------------------------------------------------------
# 7 . Save cluster assignments + posteriors
data.table::fwrite(
  cov_wgs[, c("SUBJID","IID","lpa_cluster","cluster_prob","entropy","entropy_norm","uncertain_cluster","unrel", post_names), with = FALSE],
  cluster_out, sep = "\t"
)
cat("\nCluster assignments + posterior probs saved to:", cluster_out, "\n")
# ------------------------------------------------------------------
# 8 . Evaluation vs race / ethnicity --------------------------------
# ------------------------------------------------------------------

race_var <- "Race"        # change to actual column name
if (!race_var %in% names(cov_wgs)) stop("Race column '", race_var, "' not found.")

# Tidy contingency table
race_tab <- cov_wgs[, .N, by = .(lpa_cluster, race = get(race_var))][
  order(lpa_cluster, race)
][
  , `:=`(
    prop_in_cluster = N / sum(N),          # % within cluster
    prop_overall    = N / sum(N) * .N,     # % of total (computed below)
    .by = lpa_cluster
  )
]

# Fix prop_overall = N / grand total
race_tab[, prop_overall := N / sum(N)]

# Write wide-format table for downstream use (unchanged from your version)
race_tab_wide <- dcast(race_tab, lpa_cluster ~ race, value.var = "N", fill = 0)

fwrite(race_tab_wide, race_table_tsv, sep = "\t")
cat("\nCluster × race contingency table written to:", race_table_tsv, "\n")

# ------------------------------------------------------------------
# 9 . Plots ----------------------------------------------------------
# ------------------------------------------------------------------
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
})

# ----- (A) BIC and ICL curve (training set), highlight optimal K ------------

# Our DA-EM selection returned a data.frame with columns K, ICL
bic_dt <- as.data.table(bic_obj)

# Be forgiving with names (e.g., 'G' or 'bic')
if (!"K" %in% names(bic_dt) && "G" %in% names(bic_dt)) setnames(bic_dt, "G", "K")
if (!"BIC" %in% names(bic_dt) && "bic" %in% names(bic_dt)) setnames(bic_dt, "bic", "BIC")

stopifnot(all(c("K","BIC") %in% names(bic_dt)))

bic_dt <- bic_dt[is.finite(BIC)][order(K)]

bic_plot <- ggplot(bic_dt, aes(x = K, y = BIC)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  geom_vline(xintercept = K_opt, linetype = "dashed") +
  annotate(
    "label",
    x = K_opt,
    y = max(bic_dt$BIC, na.rm = TRUE),
    label = paste0("K* = ", K_opt),
    vjust = -0.5, label.size = 0.25
  ) +
  scale_x_continuous(breaks = sort(unique(bic_dt$K))) +
  labs(
    title = "Model Selection by BIC (DA-EM with shrinkage)",
    subtitle = "Dashed line marks the BIC-optimal number of classes",
    x = "Number of classes (K)",
    y = "BIC Score"
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title.position = "plot", panel.grid.minor = element_blank())

ggsave(filename = bic_png, plot = bic_plot, width = 5, height = 3, dpi = 300)
cat("\nICL plot written to:", bic_png, "\n")

# Our DA-EM selection returned a data.frame with columns K, ICL
icl_dt <- as.data.table(icl_obj)

# Be forgiving with names (e.g., 'G' or 'bic')
if (!"K" %in% names(icl_dt) && "G" %in% names(icl_dt)) setnames(icl_dt, "G", "K")
if (!"ICL" %in% names(icl_dt) && "icl" %in% names(icl_dt)) setnames(icl_dt, "icl", "ICL")

stopifnot(all(c("K","ICL") %in% names(icl_dt)))

icl_dt <- icl_dt[is.finite(ICL)][order(K)]

icl_plot <- ggplot(icl_dt, aes(x = K, y = ICL)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  geom_vline(xintercept = K_opt, linetype = "dashed") +
  annotate(
    "label",
    x = K_opt,
    y = max(icl_dt$ICL, na.rm = TRUE),
    label = paste0("K* = ", K_opt),
    vjust = -0.5, label.size = 0.25
  ) +
  scale_x_continuous(breaks = sort(unique(icl_dt$K))) +
  labs(
    title = "Model Selection by ICL (DA-EM with shrinkage)",
    subtitle = "Dashed line marks the BIC-optimal number of classes",
    x = "Number of classes (K)",
    y = "ICL Score"
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title.position = "plot", panel.grid.minor = element_blank())

ggsave(filename = icl_png, plot = icl_plot, width = 5, height = 3, dpi = 300)
cat("\nICL plot written to:", icl_png, "\n")

# ----- (B) Race/Ethnicity – 100% stacked bars by cluster ------------

# If race_tab doesn't exist or lacks needed cols, build it from cov_wgs
need_build_race_tab <- !(exists("race_tab") &&
                         all(c("lpa_cluster","race","N","prop_in_cluster") %in% names(race_tab)))

if (need_build_race_tab) {
  # Try to detect a race/ethnicity column if 'race_var' isn't defined
  race_candidates <- c("race","Race","RACE","ethnicity","Ethnicity","ETHNICITY",
                       "race_ethnicity","RaceEthnicity","RACE_ETHNICITY")
  if (!exists("race_var") || !(race_var %in% names(cov_wgs))) {
    race_var <- intersect(race_candidates, names(cov_wgs))[1]
    if (is.na(race_var)) stop("Couldn't infer a race/ethnicity column in cov_wgs.")
  }

  # Build count table and within-cluster proportions
  race_tab <- cov_wgs[!is.na(get(race_var)),
                      .(N = .N),
                      by = .(lpa_cluster, race = as.character(get(race_var)))]
  race_tab[, prop_in_cluster := N / sum(N), by = lpa_cluster]
}

# Colorblind-friendly palette (Okabe-Ito)
okabe_ito <- c(
  "#000000","#E69F00","#56B4E9","#009E73",
  "#F0E442","#0072B2","#D55E00","#CC79A7","#999999"
)

# Ensure cluster is treated as an ordered factor by size
cluster_sizes <- cov_wgs[, .N, by = lpa_cluster][order(-N)]
race_tab[, lpa_cluster := factor(lpa_cluster, levels = cluster_sizes$lpa_cluster)]

# Ensure 'race' is a factor (stable legend/order)
race_tab[, race := factor(race)]

stacked_plot <- ggplot(race_tab,
                       aes(x = lpa_cluster, y = prop_in_cluster, fill = race)) +
  geom_col(width = 0.8) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  scale_fill_manual(
    values = okabe_ito[seq_len(min(length(okabe_ito), race_tab[, uniqueN(race)]))]
  ) +
  labs(
    title = "Race/Ethnicity Composition Within Each Cluster",
    x = "Cluster",
    y = "Within-cluster proportion",
    fill = "Race/Ethnicity"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title.position = "plot",
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )

race_bar_png <- if (exists("race_bar_png")) race_bar_png else "race_by_cluster_stacked.png"
ggsave(filename = race_bar_png, plot = stacked_plot, width = 5, height = 3, dpi = 300)
cat("\nRace stacked bar plot written to:", race_bar_png, "\n")

# ----- (C) Race/Ethnicity – annotated heatmap of counts & % ---------

heat_dt <- copy(race_tab)[
  , label := paste0(formatC(N, big.mark = ","), "\n", percent(prop_in_cluster, accuracy = 0.1))
]

# Order axes by cluster size and overall race frequency for readability
race_order <- race_tab[, .(N = sum(N)), by = race][order(-N), race]
heat_dt[, race := factor(race, levels = race_order)]

heatmap_plot <- ggplot(heat_dt, aes(x = race, y = lpa_cluster, fill = prop_in_cluster)) +
  geom_tile() +
  geom_text(aes(label = label), size = 3) +
  scale_fill_gradient(low = "white", high = "steelblue",
                      labels = percent_format(accuracy = 1)) +
  labs(
    title = "Cluster × Race/Ethnicity",
    subtitle = "Cell labels show count and within-cluster %",
    x = "Race/Ethnicity",
    y = "Cluster",
    fill = "Within-cluster %"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title.position = "plot",
    axis.text.x = element_text(angle = 30, hjust = 1),
    panel.grid = element_blank()
  )

race_heat_png <- if (exists("race_heat_png")) race_heat_png else "race_by_cluster_heatmap.png"
ggsave(filename = race_heat_png, plot = heatmap_plot, width = 5, height = 3, dpi = 300)
cat("\nRace heatmap written to:", race_heat_png, "\n")

# PC1 vs PC2 coloured by cluster (quick visual)
if (exists("num_pcs") && num_pcs >= 2 && all(c("PC1","PC2") %in% names(cov_wgs))) {
  gg <- ggplot(cov_wgs, aes_string(x = "PC1", y = "PC2", colour = "factor(lpa_cluster)")) +
          geom_point(alpha = 0.6, size = 0.7) +
          labs(colour = "Cluster") +
          theme_bw()
  ggsave(pc_plot_png, gg, width = 7, height = 6, dpi = 300)
}

# ------------------------------------------------------------------
# 10 . Contingency & plot: cluster × site ----------------------------

if (exists("site_var") && site_var %in% names(cov_wgs)) {
  site_tab <- cov_wgs[, .N, by = .(lpa_cluster, site = get(site_var))] |>
              dcast(lpa_cluster ~ site, value.var = "N", fill = 0)
  fwrite(site_tab, site_table_tsv, sep = "\t")

  # Stacked bar: proportion of sites within each cluster
  gg_site <- ggplot(cov_wgs, aes(factor(lpa_cluster), fill = get(site_var))) +
             geom_bar(position = "fill") +
             scale_y_continuous(labels = percent_format()) +
             labs(x = "Cluster", y = "Proportion", fill = "Site") +
             theme_bw()
  ggsave(site_plot_png, gg_site, width = 10, height = 3, dpi = 300)
} else {
  warning(sprintf("Site column '%s' not found – skipping site plot.", if (exists("site_var")) site_var else "<unset site_var>"))
}

# ------------------------------------------------------------------
# 11 . Next steps (manual) ------------------------------------------
# • Inspect the BIC plot to confirm cluster stability (training subset).
# • Review `cluster_prob` distribution; low values (<0.7) indicate uncertain
#   assignment – consider down‑weighting or excluding when matching.
# • Perform case–control matching *within* lpa_cluster strata.  For example,
#   use MatchIt or optmatch with exact = ~ lpa_cluster.
# • Check post‑matching balance on PCs and demographic covariates.
