#########################
##### Simulate data #####
#########################

#* Here we are going to generate sim data using benchtm package.
#* The output data is stored in conditional_marginal/data/toy_data.rds.
#* The formula for the data is as below:
#* binary outcome
#* logit(p) = 1*(X1=='Y') + 0.3*X2 + 0.3*trt
#* the prognostic factors are X1 and X2
#* the rest of the covariates are not prognostic
#* patients were randomised to pbo (0) or trt (1) group in 1:1 ratio


#* load libs
# devtools::install_github("Sophie-Sun/benchtm") # install benchtm if not installed
library(benchtm)
library(tidyverse)

#* prognostic factors
prog <- "1*(X1=='Y') + 0.3*X2"

#* generate baseline covariate data
set.seed(2023)

X <- generate_X_dist(n = 500, p = 10, rho = 0.5)

#* generate treatment data
set.seed(51)

trt <- generate_trt(n = nrow(X), p_trt = 0.5)

#* generate outcome data
dat <- generate_y(X, trt,
  prog = prog,
  pred = "0", 
  b0 = 0.3, # treatment effect
  b1 = 0, # no predictive component (no interaction of covars and trt)
  type = "binary"
)

dat$trt <- factor(dat$trt)

#* save data
saveRDS(dat, file = file.path("conditional_marginal/data/toy_data.rds"))
