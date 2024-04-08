#########################
##### Simulate data #####
#########################

#* here we are going to generate sim data using benchtm package
#* the formula for the data is as below:
#* binary data
#* logit(p) = 0.5*(X1=='Y') + 1*X3 + 0.3*trt
#* the prognostic factor would be X1, X2, X3 with different effect size
#* the rest of the covariates are not prognostic
#* patients were randomised to pbo (0) or trt (1) group in 1:1 ratio

#* install benchtm package for generating simulated data
devtools::install_github("Sophie-Sun/benchtm")

#* load libs
library(benchtm)
library(tidyverse)

set.seed(2023)
#* prognostic factors
prog <- "0.5*(X1=='Y') + 1*X3"
#* generating data
X <- generate_X_dist(n = 500, p = 10, rho = 0.5)
# observed data set
trt <- generate_trt(n = nrow(X), p_trt = 0.5)
dat <- generate_y(X, trt,
  prog = prog,
  pred = "X3>0", b0 = 0.3, b1 = 0,
  type = "binary"
)
dat$trt <- factor(dat$trt)
saveRDS(dat, file = file.path("conditional_marginal/data/toy_data.rds"))
