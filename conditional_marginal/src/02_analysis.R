# This script estimates marginal and conditional risk differences and
# odds ratios for the effect of treatment on the data in
# conditional_marginal/data/toy_data.rds.

###############################################
##### Load libraries, functions, and data #####
###############################################

# load libraries
library(tidyverse)
# remotes::install_github('openpharma/RobinCar2', ref = '70_output') # install RobinCar2 if not installed
library(RobinCar2)
library(sandwich)
library(lmtest)

# load conditional_summary function for conditional treatment effect estimates
source(file.path("conditional_marginal/funs/conditional_summary.R"))

# read in data
dat <- readRDS(file = file.path("conditional_marginal/data/toy_data.rds"))

###########################
##### Risk difference #####
###########################

### robin_glm marginal approach for risk difference ###

# Marginal risk difference estimate, no covariate adjustment
res0_marg_RD <- robin_glm(
  Y ~ trt,
  data = dat,
  treatment = trt ~ sp(1),
  contrast = "difference",
  family = binomial(link = "logit")
)

# Marginal risk difference estimate, adjusting for X1
res1_marg_RD <- robin_glm(
  Y ~ trt + X1,
  data = dat,
  treatment = trt ~ sp(1),
  contrast = "difference",
  family = binomial(link = "logit")
)

# Marginal risk difference estimate, adjusting for X1 and X2
res2_marg_RD <- robin_glm(
  Y ~ trt + X1 + X2,
  data = dat,
  treatment = trt ~ sp(1),
  contrast = "difference",
  family = binomial(link = "logit")
)

### Robust SE approach for conditional risk difference ###

# Conditional risk difference estimate, no covariate adjustment
res0_cond_RD <- conditional_summary(data = dat, formula = "Y ~ trt", 
                                    type = "RD")

# Conditional risk difference estimate, adjusting for X1
res1_cond_RD <- conditional_summary(data = dat, formula = "Y ~ trt + X1", 
                                    type = "RD")

# Conditional risk difference estimate, adjusting for X1 and X2
res2_cond_RD <- conditional_summary(data = dat, formula = "Y ~ trt + X1 + X2", 
                                    type = "RD")

# Combine results. Conditional and marginal RD estimates are similar 
# (makes sense since RD is collapsible).
# Conditional and marginal p-values are also similar.
RD_results_df <- data.frame(
  formula = c("Y ~ trt", "Y ~ trt + X1", "Y ~ trt + X1 + X2"),
  marg_RD = c(res0_marg_RD$estimate, res1_marg_RD$estimate, 
              res2_marg_RD$estimate),
  marg_SE = c(sqrt(res0_marg_RD$variance), sqrt(res1_marg_RD$variance), 
              sqrt(res2_marg_RD$variance)),
  marg_pval = c(as.data.frame(res0_marg_RD$contrast_mat)["1 v.s. 0", "Pr(>|z|)"],
                as.data.frame(res1_marg_RD$contrast_mat)["1 v.s. 0", "Pr(>|z|)"],
                as.data.frame(res2_marg_RD$contrast_mat)["1 v.s. 0", "Pr(>|z|)"]),
  cond_RD = c(as.numeric(res0_cond_RD[2, "condi_estimate"]), 
              as.numeric(res1_cond_RD[2, "condi_estimate"]),
              as.numeric(res2_cond_RD[2, "condi_estimate"])),
  cond_rb_SE = c(as.numeric(res0_cond_RD[2, "condi_rb_se"]), 
                 as.numeric(res1_cond_RD[2, "condi_rb_se"]),
                 as.numeric(res2_cond_RD[2, "condi_rb_se"])),
  cond_rb_pval = c(as.numeric(res0_cond_RD[2, "condi_rb_pval"]), 
                   as.numeric(res1_cond_RD[2, "condi_rb_pval"]),
                   as.numeric(res2_cond_RD[2, "condi_rb_pval"]))
)

RD_results_df


######################
##### Odds ratio #####
######################

### robin_glm marginal approach for odds ratio ###

# Marginal log odds ratio, no covariate adjustment
res0_marg_log_OR <- robin_glm(
  Y ~ trt,
  data = dat,
  treatment = trt ~ sp(1),
  contrast = "log_odds_ratio",
  family = binomial(link = "logit")
)

# Marginal log odds ratio, adjusting for X1
res1_marg_log_OR <- robin_glm(
  Y ~ trt + X1,
  data = dat,
  treatment = trt ~ sp(1),
  contrast = "log_odds_ratio",
  family = binomial(link = "logit")
)

# Marginal log odds ratio, adjusting for X1 and X2
res2_marg_log_OR <- robin_glm(
  Y ~ trt + X1 + X2,
  data = dat,
  treatment = trt ~ sp(1),
  contrast = "log_odds_ratio",
  family = binomial(link = "logit")
)

# Convert marginal log odds ratios to marginal odds ratios
res0_marg_OR <- exp(res0_marg_log_OR$estimate)
res1_marg_OR <- exp(res1_marg_log_OR$estimate)
res2_marg_OR <- exp(res2_marg_log_OR$estimate)

# Convert standard error of marginal log ORs to standard error of marginal ORs
res0_marg_OR_SE <- res0_marg_OR * sqrt(res0_marg_log_OR$variance)
res1_marg_OR_SE <- res1_marg_OR * sqrt(res1_marg_log_OR$variance)
res2_marg_OR_SE <- res2_marg_OR * sqrt(res2_marg_log_OR$variance)

### Robust SE approach for conditional odds ratio ###

# Conditional log odds ratio, no covariate adjustment
res0_cond_OR <- conditional_summary(data = dat, formula = "Y ~ trt", 
                                    type = "OR")

# Conditional log odds ratio, adjusting for X1
res1_cond_OR <- conditional_summary(data = dat, formula = "Y ~ trt + X1", 
                                    type = "OR")

# Conditional log odds ratio, adjusting for X1 and X2
res2_cond_OR <- conditional_summary(data = dat, formula = "Y ~ trt + X1 + X2", 
                                    type = "OR")

# Combine results. Now conditional OR estimates change more than marginal OR 
# estimates (not targeting same estimand). 
# Conditional and marginal p-values are still similar.
OR_results_df <- data.frame(
  formula = c("Y ~ trt",  "Y ~ trt + X1", "Y ~ trt + X1 + X2"),
  marg_OR = c(res0_marg_OR, res1_marg_OR, res2_marg_OR),
  marg_SE = c(res0_marg_OR_SE, res1_marg_OR_SE, res2_marg_OR_SE),
  marg_pval = c(as.data.frame(res0_marg_log_OR$contrast_mat)["1 v.s. 0", "Pr(>|z|)"],
                as.data.frame(res1_marg_log_OR$contrast_mat)["1 v.s. 0", "Pr(>|z|)"],
                as.data.frame(res2_marg_log_OR$contrast_mat)["1 v.s. 0", "Pr(>|z|)"]),
  cond_OR = c(as.numeric(res0_cond_OR[2, "condi_estimate"]), 
              as.numeric(res1_cond_OR[2, "condi_estimate"]),
              as.numeric(res2_cond_OR[2, "condi_estimate"])),
  cond_rb_SE = c(as.numeric(res0_cond_OR[2, "condi_rb_se"]), 
                 as.numeric(res1_cond_OR[2, "condi_rb_se"]),
                 as.numeric(res2_cond_OR[2, "condi_rb_se"])),
  cond_rb_pval = c(as.numeric(res0_cond_OR[2, "condi_rb_pval"]), 
                   as.numeric(res1_cond_OR[2, "condi_rb_pval"]),
                   as.numeric(res2_cond_OR[2, "condi_rb_pval"]))
)

OR_results_df
