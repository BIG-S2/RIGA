# load libraries 
library(tidyverse)
# library(Matrix)
library(broom)
library(purrr)
library(ResourceSelection)
library(survey)
library(bigmemory)
library("MatchIt")

geno_file = "${ADSP_ROOT}/data/simulated_chr1_v22k_s10k.txt"
cov_file = "${ADSP_ROOT}/data/simulated_y_5cov_5site_h0.5.txt"

variants_df <- read.big.matrix(file = geno_file, sep = " ")
cov_df <- read.table(file = cov_file, header = TRUE, sep = ",", quote = "", comment.char = "")

ipw_analysis <- function(cov_df, variants_df, idx) {
  cov_df$G <- variants_df[idx,]
  #----------------------------------------------------------------------
  # 1) Estimate the propensity scores for AD based on covariates + Site
  #----------------------------------------------------------------------
  # Example formula for logistic model:
  #   AD ~ Cov1 + Cov2 + Cov3 + Cov4 + Cov5 + Site
  # You can modify this to match your simulation setup
  formula_ps <- as.formula("Y ~ Cov1 + Cov2 + Cov3 + Cov4 + Cov5 + Site")
  
  # Fit a logistic regression to get propensity scores
  ps_fit <- glm(
    formula_ps,
    data    = cov_df,
    family  = binomial()
  )
  
  # Extract fitted probability of "treatment" (AD = 1)
  cov_df$ps <- ps_fit$fitted.values
  
  #---------------------------------------------------------
  # 2) Compute the IPW weights for the ATT (or another target)
  #---------------------------------------------------------
  # Average Treatment effect in the Treated (ATT):
  #   w = 1 for treated,  = ps / (1 - ps) for untreated
  cov_df$w_att <- ifelse(cov_df$Y == 1, 1, cov_df$ps / (1 - cov_df$ps))
  
  #------------------------------------------------------------------------
  # 3) Fit a weighted regression of the outcome (Y) on treatment (AD) (plus 
  #    optionally more covariates) using the "survey" package
  #------------------------------------------------------------------------

  # Create survey design with the IPW weights
  design_att <- svydesign(
    ids     = ~1,           # no clustering
    weights = ~w_att,       # IPW weights
    data    = cov_df
  )
  
  # Example: Weighted linear regression of Y on AD
  # You can add other covariates to this model if desired
  fit_ipw <- svyglm(
    G ~ Y,
    design = design_att,
    family = gaussian()
  )
  
  # Return model summary (you could return the fitted model object, too)
  return(fit_ipw)
}
match_analysis <- function(cov_df, variants_df, idx) {
  cov_df$G <- variants_df[idx,]
  #----------------------------------------------------------------------
  # 1) Estimate the propensity scores for AD based on covariates + Site
  #----------------------------------------------------------------------
  # Example formula for logistic model:
  #   AD ~ Cov1 + Cov2 + Cov3 + Cov4 + Cov5 + Site
  # You can modify this to match your simulation setup
  formula_ps <- as.formula("Y ~ Cov1 + Cov2 + Cov3 + Cov4 + Cov5 + Site")
  
  # Fit a logistic regression to get propensity scores
  m.out <- matchit(formula_ps,
                    data = cov_df,
                    method = "nearest",
                    distance = "glm")
  
  
  # Extract matched samples
  m.data <- match.data(m.out) 
  
  lm.test <- glm(Y ~ Cov1 + Cov2 + Cov3 + Cov4 + Cov5 + Site + G,
                 data = m.data, family = binomial(link = "logit"))
  summary(lm.test)
  
  # Return model summary (you could return the fitted model object, too)
  return(lm.test)
}
