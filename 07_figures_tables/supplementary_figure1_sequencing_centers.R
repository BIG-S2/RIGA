#!/usr/bin/env Rscript
## Rebuild the "Data Overview (SeqCenter)" supplementary figure as a SINGLE page.
## Panels a, b, c reproduced from ADSP_data_r5_seqcenter.Rmd; panel d is reduced
## from all 20 PCs to PC1-PC6 for readability. Tests cover all 20 PCs.
## Panel b is re-estimated on the current pool. Historical prediction values
## used to calibrate the simulation remain documented separately in Methods.

suppressMessages({
  library(scales)
  library(tidyverse)
})
have_patchwork <- requireNamespace("patchwork", quietly = TRUE)
if (!have_patchwork) library(gridExtra)

## ---- Current analysis pool, retaining all 15 release center labels ----------
## The original one-page figure used 38,757 historical participants and fixed
## prediction accuracies. This version uses the current 38,129 pre-design pool.
## It does not load genotypes or run association scans.
library(data.table)
library(nnet)
BASE <- "${ADSP_ROOT}"
OUT <- file.path(BASE, "revision/figures/seqcenter_v3")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
pcs <- paste0("PC", 1:20)
normalize_center <- function(x) {
  x[x %in% c("Macrogen", "MACROGEN")] <- "MACROGEN"
  x[x %in% c("Genentech", "GENENTECH")] <- "GENENTECH"
  x[x == "HudsonAlpha Institute for Biotechnology"] <- "HudsonAlpha"
  x[x == "Kinghorn Centre for Clinical Genomics (KCCG), Sydney"] <- "KCCG"
  x
}
qc <- fread(file.path(BASE, "pheno/gcad.r5.wgs.58507.2024.11.03.quality.metrics.ALL.csv"),
            select = c("Sample_name", "SeqCenter"))
setnames(qc, "Sample_name", "IID")
stopifnot(!anyDuplicated(qc$IID))
qc[, SeqCenter := normalize_center(SeqCenter)]
center_levels <- sort(unique(qc$SeqCenter))
stopifnot(length(center_levels) == 15L)
ph <- fread(file.path(BASE, "pheno/adsp_combined.csv"))
cl <- fread(file.path(BASE, "adsp_lpa_clusters_20PC.tsv"), select = "IID")
stopifnot(!anyDuplicated(ph$IID), !anyDuplicated(cl$IID))
ph <- ph[!is.na(AD) & !is.na(Sex) & !is.na(Age)]
ph[, Age := as.numeric(gsub("\\+", "", as.character(Age)))]
## Same age-band eligibility as prepare_sample_and_one_to_one_matching.R.
ph <- ph[!is.na(cut(Age, breaks = seq(30, 100, by = 3), right = TRUE))]
cov_qc <- merge(merge(ph, cl, by = "IID"), qc, by = "IID")
stopifnot(nrow(cov_qc) == 38129L, sum(cov_qc$AD == 1) == 12756L,
          uniqueN(cov_qc$SeqCenter) == 14L, uniqueN(cov_qc$Cohort) == 50L,
          !anyNA(cov_qc[, c("SeqCenter", "Cohort", pcs), with = FALSE]))
cov_qc[, SeqCenter := factor(SeqCenter, levels = center_levels)]
cov_qc[, Cohort := factor(Cohort)]
cat("Pool N", nrow(cov_qc), "cases", sum(cov_qc$AD == 1),
    "release centers", length(center_levels), "occupied centers", uniqueN(cov_qc$SeqCenter),
    "cohorts", nlevels(cov_qc$Cohort), "\n")
center_counts <- cov_qc[, .(N = .N, Cases = sum(AD == 1), Controls = sum(AD == 0)), by = SeqCenter]
center_counts <- merge(data.table(SeqCenter = center_levels), center_counts, by = "SeqCenter", all.x = TRUE)
for (nm in c("N", "Cases", "Controls")) set(center_counts, which(is.na(center_counts[[nm]])), nm, 0L)
fwrite(center_counts, file.path(OUT, "center_counts.csv"))
print(center_counts)
cohort_center <- cov_qc[, .(N = .N), by = .(Cohort, SeqCenter)]
fwrite(cohort_center, file.path(OUT, "cohort_center_counts.csv"))
ct <- table(cov_qc$Cohort, droplevels(cov_qc$SeqCenter))
cat("Single-center cohorts", sum(rowSums(ct > 0) == 1), "of", nrow(ct), "\n")
cat("Cohorts with >=90% at one center", sum(apply(ct, 1, max) / rowSums(ct) >= .9), "\n")
set.seed(42)
ct_test <- chisq.test(ct, simulate.p.value = TRUE, B = 2000)
ad_test <- chisq.test(table(droplevels(cov_qc$SeqCenter), cov_qc$AD))
cat("Cohort-center chi-square (Monte Carlo B=2000)\n"); print(ct_test)
cat("AD-center chi-square\n"); print(ad_test)
stat_tests <- data.table(Test = c("Cohort x center (Monte Carlo, B=2000)", "AD x center (Pearson)"),
                        Statistic = c(unname(ct_test$statistic), unname(ad_test$statistic)),
                        df = c(NA_real_, unname(ad_test$parameter)),
                        P = c(ct_test$p.value, ad_test$p.value))
fwrite(stat_tests, file.path(OUT, "categorical_tests.csv"))
kw <- rbindlist(lapply(pcs, function(pc) {
  z <- kruskal.test(cov_qc[[pc]], cov_qc$SeqCenter)
  data.table(PC = pc, statistic = unname(z$statistic), df = unname(z$parameter), P = z$p.value)
}))
kw[, P_BH := p.adjust(P, method = "BH")]
fwrite(kw, file.path(OUT, "pc_center_tests.csv"))
cat("Kruskal-Wallis significant after BH", sum(kw$P_BH < .05), "of", nrow(kw), "\n")

## Re-estimate all three prediction models on the same seeded 70/30 split.
## Preserve the historical multinom formulas. Standardize PCs using the training
## mean/SD to improve optimizer conditioning without changing the model space.
if ("--reuse-predictions" %in% commandArgs(trailingOnly = TRUE)) {
  acc_df <- fread(file.path(OUT, "prediction_accuracy.csv"))
  stopifnot(nrow(acc_df) == 3L, all(acc_df$Convergence == 0L),
            all(acc_df$Train_N == 26690L), all(acc_df$Test_N == 11439L))
} else {
set.seed(42)
train_idx <- sample(seq_len(nrow(cov_qc)), size = .7 * nrow(cov_qc))
train <- as.data.frame(cov_qc[train_idx]); test <- as.data.frame(cov_qc[-train_idx])
for (pc in pcs) {
  mu <- mean(train[[pc]]); sigma <- sd(train[[pc]])
  train[[pc]] <- (train[[pc]] - mu) / sigma
  test[[pc]] <- (test[[pc]] - mu) / sigma
}
formulas <- list(
  as.formula(paste("SeqCenter ~", paste(pcs, collapse = " + "))),
  SeqCenter ~ Cohort,
  as.formula(paste("SeqCenter ~ Cohort +", paste(pcs, collapse = " + ")))
)
model_names <- c("PC1-20 Only", "Cohort Only", "PC1-20 + Cohort")
acc_df <- rbindlist(lapply(seq_along(formulas), function(i) {
  fit <- multinom(formulas[[i]], data = train, MaxNWts = 50000,
                  maxit = 2000, trace = FALSE)
  stopifnot(fit$convergence == 0L)
  pred <- predict(fit, newdata = test)
  acc <- mean(as.character(pred) == as.character(test$SeqCenter))
  cat(model_names[i], "accuracy", acc, "convergence", fit$convergence, "\n")
  fwrite(data.table(Actual = test$SeqCenter, Predicted = pred),
         file.path(OUT, paste0("prediction_", i, ".csv")))
  data.table(Model = model_names[i], Accuracy = acc, Convergence = fit$convergence,
             Train_N = nrow(train), Test_N = nrow(test), Max_iterations = 2000L)
}))
fwrite(acc_df, file.path(OUT, "prediction_accuracy.csv"))

}

## Compact, unique cohort labels keep every x-axis label inside the figure.
cohort_names <- levels(cov_qc$Cohort)
cohort_labels <- sub("^.*\\(([^()]+)\\)\\s*$", "\\1", cohort_names)
manual <- c("Arizona APOE Cohort Study" = "AZ_APOE", "Cache County Study" = "CacheCounty",
            "Case Western Reserve University (CWRU) Autopsy Cohort" = "CWRU_Autopsy",
            "Case Western Reserve University (CWRU) Rapid Decline Cohort" = "CWRU_RDC",
            "IndyAA - Indianapolis African Americans" = "IndyAA", "Mount Sinai Brain Bank" = "MSBB",
            "Mayo Clinic" = "Mayo", "University of Miami" = "UM",
            "University of Pittsburgh ADRC- Kamboh" = "UPitt_Kamboh",
            "University of Washington Families" = "UW_Families", "Vanderbilt University" = "Vanderbilt")
use_manual <- cohort_names %in% names(manual)
cohort_labels[use_manual] <- manual[cohort_names[use_manual]]
long <- nchar(cohort_labels) > 16
cohort_labels[long] <- abbreviate(cohort_labels[long], minlength = 16)
stopifnot(!anyDuplicated(cohort_labels))
fwrite(data.table(Cohort = cohort_names, Figure_label = cohort_labels), file.path(OUT, "cohort_labels.csv"))
cohort_labels <- setNames(cohort_labels, cohort_names)

## ---- Panel (a): AD % by SeqCenter ------------------------------------------
p_a <- cov_qc %>%
  mutate(AD = as.factor(AD)) %>%
  group_by(SeqCenter, AD) %>%
  summarize(n = n(), .groups = "drop") %>%
  ggplot(aes(x = SeqCenter, y = n, fill = AD)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_y_continuous(labels = percent_format()) +
  scale_x_discrete(drop = FALSE) +
  annotate("text", x = which(center_levels == "MACROGEN"), y = 0.45, label = "n = 0", angle = 90, size = 3) +
  scale_fill_manual(values = c("0" = "#4DAF4A", "1" = "#E41A1C"),
                    labels = c("0" = "Control", "1" = "AD Case")) +
  labs(title = "AD status by sequencing center",
       subtitle = "15 centers in the release. MACROGEN has no participants in this pool.",
       x = "Sequencing Center", y = "Percentage", fill = "AD Status", tag = "(a)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.tag = element_text(face = "bold", size = 18))

## ---- Panel (b): SeqCenter prediction accuracy by model ---------------------
p_b <- ggplot(acc_df, aes(x = reorder(Model, Accuracy), y = Accuracy)) +
  geom_segment(aes(xend = Model, y = 0, yend = Accuracy),
               color = "steelblue", linewidth = 1.2) +
  geom_point(size = 5, color = "steelblue") +
  geom_text(aes(label = paste0(round(Accuracy * 100, 1), "%")),
            hjust = -0.25, size = 5) +
  coord_flip() +
  scale_y_continuous(labels = percent_format(), limits = c(0, 1.3),
                     breaks = c(0, 0.3, 0.6, 0.9)) +
  labs(title = "Sequencing-center prediction accuracy",
       subtitle = sprintf("Held-out test set (n = %s)", format(acc_df$Test_N[1], big.mark = ",")),
       x = "", y = "Accuracy", tag = "(b)") +
  theme_minimal(base_size = 14) +
  theme(plot.tag = element_text(face = "bold", size = 18))

## ---- Panel (c): SeqCenter distribution (%) by Cohort -----------------------
n_centers <- length(center_levels)
center_colors <- colorRampPalette(
  c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",
    "#A65628", "#F781BF", "#999999", "#66C2A5", "#FC8D62",
    "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494")
)(n_centers)

p_c <- cov_qc %>%
  group_by(Cohort, SeqCenter) %>%
  summarize(n = n(), .groups = "drop") %>%
  ggplot(aes(x = Cohort, y = n, fill = SeqCenter)) +
  geom_bar(stat = "identity", position = "fill", show.legend = TRUE) +
  scale_y_continuous(labels = percent_format()) +
  scale_x_discrete(labels = cohort_labels) +
  scale_fill_manual(values = setNames(center_colors, center_levels), drop = FALSE) +
  labs(title = "Sequencing-center distribution by cohort",
       subtitle = "Distributions differ across cohorts (Monte Carlo chi-square P < 0.001).",
       x = "Cohort", y = "Percentage", fill = "SeqCenter", tag = "(c)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 8),
        plot.tag = element_text(face = "bold", size = 18))

## ---- Panel (d): PC1-PC6 boxplots by SeqCenter (REDUCED from PC1-PC20) -------
rep_pcs <- paste0("PC", 1:6)
p_d <- cov_qc %>%
  dplyr::select(SeqCenter, all_of(rep_pcs)) %>%
  pivot_longer(cols = all_of(rep_pcs), names_to = "PC", values_to = "value") %>%
  mutate(PC = factor(PC, levels = rep_pcs)) %>%
  ggplot(aes(x = SeqCenter, y = value, fill = SeqCenter)) +
  geom_boxplot(outlier.size = 0.2) +
  scale_x_discrete(drop = FALSE) +
  scale_fill_manual(values = setNames(center_colors, center_levels), drop = FALSE) +
  facet_wrap(~PC, scales = "free_y", ncol = 2) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
        strip.text = element_text(size = 11, face = "bold"),
        legend.position = "none",
        plot.tag = element_text(face = "bold", size = 18)) +
  labs(title = "Genetic PC distributions by sequencing center",
       subtitle = "PC1-PC6 shown. All 20 PCs differ across centers (BH-adjusted P < 0.05).",
       x = "Sequencing Center", y = "PC Value", tag = "(d)")

## ---- Assemble onto one page ------------------------------------------------
out <- file.path(BASE, "revision/manuscript/figures/Data Overview_seqCenter.png")
out <- path.expand(out)

if (have_patchwork) {
  library(patchwork)
  combined <- (p_a + p_b) / p_c / p_d +
    plot_layout(heights = c(3.0, 3.2, 7.8)) +
    plot_annotation(title = "Sequencing-center heterogeneity in the ADSP analysis pool",
                    subtitle = "38,129 participants across 50 study cohorts",
                    theme = theme(plot.title = element_text(face = "bold", size = 17)))
  ggsave(out, combined, width = 13, height = 17, dpi = 300, limitsize = FALSE, bg = "white")
  ggsave(file.path(OUT, "Data Overview_seqCenter.pdf"), combined, width = 13, height = 17, limitsize = FALSE)
} else {
  g <- arrangeGrob(
    arrangeGrob(p_a, p_b, ncol = 2, widths = c(7, 5)),
    p_c, p_d, ncol = 1, heights = c(3.0, 3.2, 7.8)
  )
  ggsave(out, g, width = 13, height = 17, dpi = 300, limitsize = FALSE)
}
cat("Wrote:", out, "\n")
