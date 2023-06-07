# Estimate hypothetical estimand (effect of treatment in a hypothetical trial 
# where switching to rescue medication is not possible) using IPW.
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

# option 2: user load their own data, should keep the same data structure
#analysis_data <-  readRDS(file = "") #the path to the dataset

################################
# Get point estimates using IPW
################################

# Treatment arms
arms <- unique(analysis_data$trt)

# Empty data frame to store weights
wts <- data.frame(usubjid = NULL, trt = NULL, switch = NULL, wt = NULL, 
                  Y = NULL)

# Calculate weights for each arm separately
for (i in 1:length(arms)){
  
  # Pull data of people in arm i
  ipwdata <- analysis_data %>% dplyr::filter(trt == arms[i])
  
  # Fit model for the probability of no trt switching
  ipw <- glm(1 - switch ~ X1 + X2 +
               X3 + X4_W26, 
             data = ipwdata, family = "binomial")
  
  # Compute weights of 1 / P(Z_1 = 0 | Z_0, X_0, X_1) 
  ipwdata$wt <- 1 / ipw$fitted.values
  
  # Store weights for Z_0 = arms[i]
  wts <- rbind(wts, ipwdata %>% select(usubjid, trt, switch, wt, Y))
}

# Only keep patients with Z_1 = 0 (no switching)
wts_sub <- wts %>% filter(switch == 0)

# Get weighted outcomes of non-switchers in treatment arm
trt_data <- wts_sub %>% filter(trt == 1)

trt_out <- sum(trt_data$wt * trt_data$Y) / sum(trt_data$wt)

# Get weighted outcomes of non-switchers in placebo arm
pbo_data <- wts_sub %>% filter(trt == 0)

pbo_out <- sum(pbo_data$wt * pbo_data$Y) / sum(pbo_data$wt)


#################################################
# Calculate confidence intervals using bootstrap
#################################################

# Set seed
set.seed(1234)

# Number of bootstrapped samples
B <- 15000

# Dataframe to store bootstrapped point estimates
boot_res <- data.frame(bootid = 1:B, trt = NA, pbo = NA)

# Define progress bar
pb <- progress_bar$new(
  format = "  Iteration :current of :total [:bar] :percent eta: :eta",
  total = B, clear = FALSE, width = 60)

for(b in 1:B){ 
  
  # Update progress bar
  pb$tick()
  
  # Resample within treatment group
  analysis_data_boot <- analysis_data %>% 
    group_by(trt) %>% 
    sample_n(size = n(), replace = TRUE) %>% 
    ungroup
  
  # Treatment arms
  arms <- unique(analysis_data_boot$trt)
  
  # Empty data frame to store weights
  wts_boot <- data.frame(usubjid = NULL, trt = NULL, switch = NULL, 
                         wt = NULL, Y = NULL)
  
  # Calculate weights for each arm separately
  for (i in 1:length(arms)){
    
    # Pull data of people in arm i
    ipwdata <- analysis_data_boot %>% dplyr::filter(trt == arms[i])
    
    # Fit model for the probability of no trt switching
    ipw <- glm(1 - switch ~ X1 + X2 +
                 X3 + X4_W26, 
               data = ipwdata, family = "binomial")
    
    # Compute weights of 1 / P(Z_1 = 0 | Z_0, X_0, X_1) 
    ipwdata$wt <- 1 / ipw$fitted.values
    
    # Store weights for Z_0 = arms[i]
    wts_boot <- 
      rbind(wts_boot, 
            ipwdata %>% select(usubjid, trt, switch, wt, Y))
  }
  
  # Only keep patients with Z_1 = 0 (no switching)
  wts_bootsub <- wts_boot %>% filter(switch == 0)
  
  # Get weighted outcomes of non-switchers in treatment arm
  trt_data <- wts_bootsub %>% filter(trt == 1)
  
  trt_out <- sum(trt_data$wt * trt_data$Y) / sum(trt_data$wt)
  
  # Get weighted outcomes of non-switchers in placebo arm
  pbo_data <- wts_bootsub %>% filter(trt == 0)
  
  pbo_out <- sum(pbo_data$wt * pbo_data$Y) / sum(pbo_data$wt)
  
  # comparisons
  boot_res[b, 2:3] <- c(trt_out, pbo_out)
  
  # Remove analysis_data_boot from memory
  rm("analysis_data_boot")
}

# E[Y(Z_0 = 1, Z_1 = 0)] quantiles
round(quantile(boot_res$trt, c(0.025, 0.975)), 2)

# E[Y(Z_0 = 0, Z_1 = 0)] quantiles
round(quantile(boot_res$pbo, c(0.025, 0.975)), 2)

# E[Y(Z_0 = 1, Z_1 = 0)] - E[Y(Z_0 = 0, Z_1 = 0)] quantiles
round(quantile(boot_res$trt - boot_res$pbo, c(0.025, 0.975)), 2)

# Save bootstrapped data
#saveRDS(boot_res, file = "") #should have the correct path and file name to store the bootstrap results
