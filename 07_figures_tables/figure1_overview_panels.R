###########################################################
# figure1_overview_panels.R (2026-10-05) -- regenerate the panels of main
# Figure 1 (cohort heterogeneity, population structure, covariate
# associations) from the CURRENT phenotyped pool (the pre-design pool of
# the stratified design, n = 38,129) and the current 13-class latent-class
# fit. Each panel is written as its own PNG into revision/figures/fig1_v3/;
# fig1_assemble_v3.py composes them in the layout of the original data.png.
###########################################################
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(scales); library(patchwork); library(car)
})
OUT <- "${ADSP_ROOT}/revision/figures/fig1_v3"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
SRC <- "${ADSP_ROOT}/src"
LOGF <- file.path(OUT, "fig1_build_log.txt"); log_lines <- character(0)
`%||%` <- function(a, b) if (is.null(a) || is.na(a)) b else a
LOG <- function(...) { s <- sprintf(...); message(s); log_lines <<- c(log_lines, s) }

## ---- the pre-design pool, exactly as riga_study_units_null_model.R built it ----
src <- readLines(file.path(SRC, "prepare_sample_and_one_to_one_matching.R"))
cut <- grep("^m\\.out <- matchit", src)[1]
eval(parse(text = paste(src[1:(cut - 1)], collapse = "\n")))
dfx <- as.data.table(dfx)
ph_all <- fread("${ADSP_ROOT}/pheno/adsp_combined.csv")
if (!"Cohort" %in% names(dfx)) dfx <- merge(dfx, ph_all[, .(IID, Cohort)], by = "IID", all.x = TRUE)
qc <- fread("${ADSP_ROOT}/pheno/gcad.r5.wgs.58507.2024.11.03.quality.metrics.ALL.csv",
            select = c("Sample_name", "Study"), colClasses = "character")
setnames(qc, "Sample_name", "IID"); qc <- unique(qc, by = "IID")
dfx <- merge(dfx, qc, by = "IID", all.x = TRUE); dfx[is.na(Study) | Study == "", Study := "Unknown"]
LOG("pool: %d subjects, %d cases, %d cohorts, %d projects, %d centers, %d latent classes",
    nrow(dfx), sum(dfx$AD == 1), uniqueN(dfx$Cohort), uniqueN(dfx$Study), uniqueN(dfx$SeqCenter), uniqueN(dfx$lpa_cluster))
short_lab <- function(x) sub("^.*\\(([^()]+)\\)\\s*$", "\\1", x)   # the acronym in parentheses when present
manual <- c("Arizona APOE Cohort Study" = "AZ_APOE", "Cache County Study" = "CacheCounty",
            "Case Western Reserve University (CWRU) Autopsy Cohort" = "CWRU_Autopsy",
            "Case Western Reserve University (CWRU) Rapid Decline Cohort" = "CWRU_RDC",
            "IndyAA - Indianapolis African Americans" = "IndyAA", "Mount Sinai Brain Bank" = "MSBB",
            "Mayo Clinic" = "Mayo", "University of Miami" = "UM", "University of Pittsburgh ADRC- Kamboh" = "UPitt_Kamboh",
            "University of Washington Families" = "UW_Families", "Vanderbilt University" = "Vanderbilt")
dfx[, CohortS := short_lab(as.character(Cohort))]
dfx[as.character(Cohort) %in% names(manual), CohortS := manual[as.character(Cohort)]]
dfx[nchar(CohortS) > 14, CohortS := abbreviate(CohortS, 14)]
LOG("cohort labels: %s", paste(sort(unique(dfx$CohortS)), collapse = ", "))
race_lab <- c(`1` = "American Indian/Alaska Native", `2` = "Asian", `3` = "Native Hawaiian/Pacific Islander",
              `4` = "Black or African American", `5` = "White", `6` = "Hispanic/Latino")
dfx[, RaceF := factor(race_lab[as.character(Race)], levels = race_lab)]
dfx[, ADF := factor(ifelse(AD == 1, "AD", "Control"), levels = c("Control", "AD"))]
dfx[, SexF := factor(Sex)]
dfx[, SeqCenter := as.character(SeqCenter)]   # merge the two label pairs that differ only in capitalization (15 centers)
dfx[SeqCenter %in% c("Macrogen", "MACROGEN"), SeqCenter := "MACROGEN"]; dfx[SeqCenter %in% c("Genentech", "GENENTECH"), SeqCenter := "GENENTECH"]
dfx[, SeqCenter := factor(SeqCenter)]
th <- theme_minimal(base_size = 7) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 5.5),
                                            axis.text.y = element_text(size = 6), axis.title = element_text(size = 7),
                                            plot.title = element_text(size = 7.5, face = "bold"),
                                            legend.title = element_text(size = 6.5), legend.text = element_text(size = 6),
                                            legend.key.size = unit(0.28, "cm"), legend.position = "right",
                                            panel.grid.minor = element_blank(), plot.margin = margin(3, 3, 1, 3))
nox <- theme(axis.text.x = element_blank(), axis.title.x = element_blank(), axis.ticks.x = element_blank())

## ---- panel a: cohort-level distributions (four strips) ----
pa1 <- ggplot(dfx[, .N, by = .(CohortS, ADF)], aes(CohortS, N, fill = ADF)) +
  geom_col(position = "fill") + scale_y_continuous(labels = percent_format(), breaks = c(0, 0.5, 1)) +
  scale_fill_manual(values = c(Control = "#1b2a41", AD = "#1e90ff")) +
  labs(x = NULL, y = "Percentage", fill = "AD status", title = "AD status by cohort") + th + nox
pa2 <- ggplot(dfx[, .N, by = .(CohortS, SexF)], aes(CohortS, N, fill = SexF)) +
  geom_col(position = "fill") + scale_y_continuous(labels = percent_format(), breaks = c(0, 0.5, 1)) +
  scale_fill_manual(values = c("#1b2a41", "#1e90ff", "#9ecae1")) +
  labs(x = NULL, y = "Percentage", fill = "Sex (code)", title = "Sex by cohort") + th + nox
ad <- dfx[, .(lo = min(Age), med = median(Age), hi = max(Age)), by = CohortS]
pa3 <- ggplot(ad, aes(CohortS, med)) + geom_pointrange(aes(ymin = lo, ymax = hi), color = "#1b4f9c", size = 0.12, linewidth = 0.35) +
  scale_y_continuous(breaks = c(40, 60, 80)) + labs(x = NULL, y = "Age", title = "Age by cohort (median, minimum to maximum)") + th + nox
pa4 <- ggplot(dfx[, .N, by = .(CohortS, RaceF)], aes(CohortS, N, fill = RaceF)) +
  geom_col(position = "fill") + scale_y_continuous(labels = percent_format(), breaks = c(0, 0.5, 1)) +
  scale_fill_manual(values = c("#e41a1c", "#b5a642", "#4daf4a", "#2a9d8f", "#4f81f7", "#f032e6"), drop = FALSE) +
  labs(x = NULL, y = "Percentage", fill = "Self-reported race/ethnicity", title = "Race/ethnicity by cohort") + th + nox
LOG("panel a strips built (%d cohorts)", uniqueN(dfx$CohortS))

## ---- panel b: mean PCs by cohort ----
pcs <- paste0("PC", 1:20)
mh <- melt(dfx[, c("CohortS", pcs), with = FALSE], id.vars = "CohortS", variable.name = "PC", value.name = "v")
mh <- mh[, .(m = mean(v)), by = .(CohortS, PC)]
mh[, PC := factor(PC, levels = pcs)]
lim <- max(abs(mh$m))
mh[, PCr := factor(PC, levels = rev(pcs))]
pb <- ggplot(mh, aes(CohortS, PCr, fill = m)) + geom_tile(color = "white", linewidth = 0.2) +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", midpoint = 0, limits = c(-lim, lim)) +
  labs(x = "Cohort", y = NULL, fill = "Mean PC value", title = "Mean genetic PCs by cohort") +
  th + theme(axis.text.y = element_text(size = 5), panel.grid = element_blank())
pab <- (pa1 / pa2 / pa3 / pa4 / pb) + plot_layout(heights = c(1, 1, 1, 1.4, 3.6)) +
  plot_annotation(tag_levels = list(c("a", "", "", "", "b"))) & theme(plot.tag = element_text(face = "bold", size = 9.5))
ggsave(file.path(OUT, "a.png"), pab, width = 7.2, height = 5.1, dpi = 300, bg = "white")
LOG("panels a+b written as one aligned block")
kw <- mh[, .(p = kruskal.test(v ~ CohortS, data = dfx[, .(v = get(as.character(.BY$PC)), CohortS)])$p.value), by = PC]
LOG("panel b written; Kruskal-Wallis PC ~ cohort: %d of 20 PCs with BH-adjusted p < 0.05", sum(p.adjust(kw$p, "BH") < 0.05))

## ---- panel c: BIC curve of the current 13-class fit (the saved curve of the deployed fit) ----
bic_src <- "${ADSP_ROOT}/bic_curve_20PC_unrelated.png"
file.copy(bic_src, file.path(OUT, "c.png"), overwrite = TRUE)
LOG("panel c: copied %s (DA-EM with shrinkage; BIC-optimal K = 13)", basename(bic_src))

## ---- panel d: PC1 vs PC2 colored by latent class ----
dfx[, Class := factor(lpa_cluster, levels = sort(unique(lpa_cluster)))]
pal13 <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f",
           "#bcbd22", "#17becf", "#393b79", "#ad494a", "#e7ba52")
set.seed(1); ord <- sample(nrow(dfx))
pd <- ggplot(dfx[ord], aes(PC1, PC2, color = Class)) + geom_point(size = 0.22, alpha = 0.6, stroke = 0) +
  scale_color_manual(values = pal13) + guides(color = guide_legend(override.aes = list(size = 1.8, alpha = 1), ncol = 2)) +
  labs(color = "Latent class", title = sprintf("Latent genetic classes (K = %d)", nlevels(dfx$Class))) +
  theme_minimal(base_size = 7) + theme(plot.title = element_text(size = 7.5, face = "bold"), axis.text = element_text(size = 6),
                                       axis.title = element_text(size = 7), legend.title = element_text(size = 6.5),
                                       legend.text = element_text(size = 6), legend.key.size = unit(0.25, "cm"),
                                       panel.grid.minor = element_blank(), plot.margin = margin(3, 3, 1, 3))
ggsave(file.path(OUT, "d.png"), pd, width = 3.6, height = 1.55, dpi = 300, bg = "white")
LOG("panel d written: %d classes on %d subjects", nlevels(dfx$Class), nrow(dfx))

## ---- panel e: partial residuals of the covariate-adjusted pooled model ----
## the pooled comparator of the plasmode benchmark: age, sex, center, collapsed project, PC1-20.
## Center is nested in project for all but two projects, so the fit is rank-deficient and
## car::crPlots refuses it; the component + residual plots are computed directly
## (working residuals + the term's fitted contribution, aliased coefficients taken as 0).
dfx[, Study_c := as.character(Study)]
st <- dfx[, .(N = .N, cases = sum(AD == 1)), by = Study_c]
small <- st[N < 50 | cases == 0 | cases == N]$Study_c
dfx[Study_c %in% small, Study_c := "Other"]; dfx[, Study_c := factor(Study_c)]
LOG("pooled model: %d project levels after collapsing %d single-label or small projects to Other", nlevels(dfx$Study_c), length(small))
form <- as.formula(paste("AD ~ Age + SexF + SeqCenter + Study_c +", paste(pcs, collapse = " + ")))
fit <- suppressWarnings(glm(form, data = dfx, family = binomial()))
LOG("pooled glm: converged %s; deviance %.1f on %d df; aliased coefficients %d; max |coef| %.2f",
    fit$converged, fit$deviance, fit$df.residual, sum(is.na(coef(fit))), max(abs(coef(fit)), na.rm = TRUE))
thcr <- theme_bw(base_size = 6.5) + theme(axis.text = element_text(size = 5.2), axis.title = element_text(size = 6.5),
                                           panel.grid.minor = element_blank(), plot.margin = margin(2, 3, 1, 2))
cr_data <- function(fit, term) {
  tt <- suppressWarnings(predict(fit, type = "terms"))
  data.table(x = fit$model[[term]], cr = residuals(fit, type = "working") + tt[, term])
}
cr_panel <- function(fit, term, xlab = term, smooth = "gam", ylim = NULL) {
  d <- cr_data(fit, term); ylab <- "Partial residual"
  if (!is.numeric(d$x)) d[, x := factor(sub("^Kinghorn.*", "KCCG", as.character(x)))]
  if (is.numeric(d$x)) {
    g <- ggplot(d, aes(x, cr)) + geom_point(shape = 16, size = 0.2, alpha = 0.35, stroke = 0) +
      geom_smooth(method = "lm", formula = y ~ x, se = FALSE, color = "#1b4f9c", linetype = "dashed", linewidth = 0.45)
    g <- if (smooth == "gam") g + geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), se = FALSE, color = "#e7298a", linewidth = 0.55)
         else g + geom_smooth(method = "loess", formula = y ~ x, span = 0.75, se = FALSE, color = "#e7298a", linewidth = 0.55)
    g <- g + scale_x_continuous(breaks = scales::breaks_pretty(n = 3)) + labs(x = xlab, y = ylab) + thcr
    if (!is.null(ylim)) g <- g + coord_cartesian(ylim = ylim)
    g
  } else {
    g <- ggplot(d, aes(x, cr)) + geom_boxplot(outlier.size = 0.2, outlier.alpha = 0.35, outlier.stroke = 0, fill = "grey92", linewidth = 0.25) +
      labs(x = xlab, y = ylab) + thcr + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 4.8))
    if (!is.null(ylim)) g <- g + coord_cartesian(ylim = ylim)
    g
  }
}
## the working residuals of near-separated subjects reach +-150; show the central 99.5% (the smooths use all points)
yl_e <- as.numeric(quantile(cr_data(fit, "Age")$cr, c(0.0025, 0.9975)))
LOG("panel e y-range (0.25-99.75%% of the partial residuals): %.1f to %.1f", yl_e[1], yl_e[2])
ab <- cr_data(fit, "Age"); ab[, bin := cut(x, c(0, 50, 60, 70, 80, 90, 130), right = FALSE)]
LOG("age coefficient (log-odds per year) %.4f; mean partial residual by age bin: %s", coef(fit)["Age"],
    paste(sprintf("%s %.2f (n %d)", levels(ab$bin), ab[, mean(cr), by = bin][order(bin)]$V1, ab[, .N, by = bin][order(bin)]$N), collapse = "; "))
pe <- wrap_plots(lapply(c("Age", "SeqCenter", "SexF", paste0("PC", 1:6)), function(tm)
  cr_panel(fit, tm, c(Age = "Age", SeqCenter = "Sequencing center", SexF = "Sex (code)")[tm] %||% tm, ylim = yl_e)), ncol = 3)
ggsave(file.path(OUT, "e.png"), pe, width = 3.6, height = 2.8, dpi = 300, bg = "white")
LOG("panel e written")

## ---- panel f: the same diagnostics inside one prospective cohort (ACT) ----
act <- dfx[grepl("Adult Changes in Thought", Cohort)]
act[, SeqCenter := factor(as.character(SeqCenter))]; act[, RaceF := factor(as.character(Race))]; act[, SexF := factor(Sex)]
LOG("ACT: %d subjects, %d cases, centers %s, race codes %s", nrow(act), sum(act$AD == 1),
    paste(levels(act$SeqCenter), collapse = "/"), paste(levels(act$RaceF), collapse = "/"))
terms_f <- c("Age", if (nlevels(act$SeqCenter) > 1) "SeqCenter", if (nlevels(act$SexF) > 1) "SexF",
             if (nlevels(act$RaceF) > 1) "RaceF", paste0("PC", 1:5))
fit_act <- suppressWarnings(glm(as.formula(paste("AD ~", paste(terms_f, collapse = " + "))), data = act, family = binomial()))
LOG("ACT glm: converged %s; terms %s; aliased %d", fit_act$converged, paste(terms_f, collapse = ", "), sum(is.na(coef(fit_act))))
labs_f <- c(Age = "Age", SeqCenter = "Sequencing center", SexF = "Sex (code)", RaceF = "Race (code)")
pf <- wrap_plots(lapply(terms_f, function(tm) cr_panel(fit_act, tm, if (tm %in% names(labs_f)) labs_f[[tm]] else tm, smooth = "loess")), ncol = 3)
ggsave(file.path(OUT, "f.png"), pf, width = 3.6, height = 2.8, dpi = 300, bg = "white")
LOG("panel f written")
writeLines(log_lines, LOGF); LOG("done.")
