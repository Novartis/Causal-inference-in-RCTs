# Estimate hypothetical estimand (effect of treatment in a hypothetical trial 
# where switching to rescue medication is not possible) using G-computation.
# This code simulates data, computes the point estimate, and constructs
# a bootstrapped CI.

library(tidyverse)
library(data.table)
library(progress)

#########################
# Simulate or load data 
#########################

# option 1: simulated data
#simulate covariates
set.seed(2023)
n <- 500 # number of subjects in each arm
X1 <- rnorm(n*2, 0, 1)
X2 <- rnorm(n*2, 0, 1)
X3 <- rnorm(n*2, 0, 1)
X4_W26 <- rnorm(n*2, 0, 1)
trt <- rbinom(n*2, 1, 0.5)

analysis_data <- data.frame(usubjid = seq(1,n), trt = trt, switch = NA, 
                            X1 = X1, X2 = X2, X3 = X3, X4_W26 = X4_W26, 
                            Y = NA)

#model for treatment switching
predictors <- cbind(X1, X2, X3, X4_W26)
pre_trt <- 1 + predictors %*% c(0.01, -0.01, -0.01, -0.02)
pre_con <- 0.5 + predictors %*% c(0.02, -0.02, -0.03, -0.04)
ps_trt <- exp(pre_trt) / (1 + exp(pre_trt))
ps_con <- exp(pre_con) / (1 + exp(pre_con))
prob_switch <- analysis_data$trt*ps_trt + (1 - analysis_data$trt)*ps_con
switch <- 1 - rbinom(nrow(analysis_data), 1, prob_switch)

analysis_data$switch <- switch

#outcome model

link <- 0.2 + 0.5*trt + 0.1*X1 -0.1* X2 + 0.05*X3 -0.05 * X4_W26 + 
  0.1* trt * X1 + (-0.3)*trt*X2 + 0.2*trt*X3 + (-0.2) *trt * X4_W26
response_probs <- exp(link)/(1 + exp(link))
Y_nonswitch <- rbinom(nrow(analysis_data[switch==0,]), 1, response_probs[switch == 0])
analysis_data$Y[switch == 0] <- Y_nonswitch

#the analysis_data is ready for performing the analysis


# option 2: user loads their own data, should keep the same data structure
# analysis_data <- readRDS(file = "") #the path to the dataset

###########################################
# Get point estimates using g-computation
###########################################

# Fit outcome model only using patients who do not switch treatment
fit_logistic_noswitch <- glm(data = analysis_data %>% filter(switch == 0), 
                             formula = Y ~ trt + 
                               trt*X1 + 
                               trt*X2 +
                               trt*X3 + 
                               trt*X4_W26, 
                             family = binomial(link = "logit"))

summary(fit_logistic_noswitch)

# Prediction based on model that only included patients without treatment switch
# Prediction based only on patients from within the same group
analysis_data$predict_hypo <- predict(object = fit_logistic_noswitch, 
                                      newdata  = analysis_data %>% 
                                        select(trt, X1, X2, X3, X4_W26),
                                      type = "response") 

# Calculate marginal mean under hypothetical scenario
analysis_data %>% 
  group_by(trt) %>% 
  summarize(mean_resp = mean(predict_hypo), .groups = "keep")


# Remove predict_hypo column
analysis_data$predict_hypo <- NULL

#################################################
# Calculate confidence intervals using bootstrap
#################################################

# Set up bootstrap parameters
B <- 15000
boot_list <- vector(mode = "list", length = B)

# Define progress bar
pb <- progress_bar$new(
  format = "  Iteration :current of :total [:bar] :percent eta: :eta",
  total = B, clear = FALSE, width = 60)

# Set seed
set.seed(950)

for(b in 1:B) {
  
  # Update progress bar
  pb$tick()
  
  # Resample within treatment group
  analysis_data_boot <- analysis_data %>% 
    group_by(trt) %>% 
    sample_n(size = n(), replace = TRUE) %>% 
    ungroup
  
  # Fit outcome model only using patients who do not switch treatment
  fit_logistic_switch_boot <- 
    glm(data = analysis_data_boot %>% filter(switch == 0), 
        formula = Y ~ trt + 
          trt*X1 + 
          trt*X2 +
          trt*X3 + 
          trt*X4_W26, 
        family = binomial(link = "logit"))
  
  # Prediction based on model that only included patients without trt switching
  # Prediction based only on patients from within the same group
  analysis_data_boot$predict_hypo <- 
    predict(object = fit_logistic_switch_boot, 
            newdata  = analysis_data_boot %>% 
              dplyr::select(trt, X1, X2, X3, X4_W26),
            type = "response") 
  
  # Calculate marginal mean under hypothetical scenario
  boot_list[[b]] <- analysis_data_boot %>% 
    group_by(trt) %>% 
    summarize(mean_resp = mean(predict_hypo), .groups = "keep") %>% 
    mutate(run = b) %>% 
    ungroup
  
  # Remove analysis_data_boot from memory
  rm("analysis_data_boot")
  
}

# E[Y(Z_0 = 0, Z_1 = 0)] and E[Y(Z_0 = 1, Z_1 = 0)] estimate summaries
bind_rows(boot_list) %>% 
  group_by(trt) %>% 
  summarize(q025 = quantile(mean_resp, 0.025),
            q975 = quantile(mean_resp, 0.975),
            mean = mean(mean_resp),
            median = median(mean_resp),
            sd = sd(mean_resp), 
            .groups = "keep")

# E[Y(Z_0 = 0, Z_1 = 0)] - E[Y(Z_0 = 1, Z_1 = 0)] estimate summaries
bind_rows(boot_list) %>% 
  pivot_wider(values_from = "mean_resp", names_from = "trt",
              names_prefix = "trt_") %>% 
  mutate(trt1_vs_trt0 = trt_1 - trt_0) %>% 
  summarize(effect = mean(trt1_vs_trt0),
            lcl = quantile(trt1_vs_trt0, probs = 0.025),
            ucl = quantile(trt1_vs_trt0, probs = 0.975))


# Save bootstrapped data
#saveRDS(boot_list, file = "") #should have the correct path and file name to store the bootstrap results
