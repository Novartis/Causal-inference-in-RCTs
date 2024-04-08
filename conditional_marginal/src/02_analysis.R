#* load libs and funs
library(tidyverse)
source(file.path("conditional_marginal/funs/funs.R"))
set.seed(2023)
dat <- readRDS(file = file.path("conditional_marginal/data/toy_data.rds"))

#* unadjusted
formula0 <- as.formula("Y ~ trt")
#* adjust for one weak/moderate prognostic factor X1
formula1 <- as.formula("Y ~ trt + X1")
#* adjust for one weak/moderate prognostic factors X1 and a strong prognostic factor X3
formula2 <- as.formula("Y ~ trt + X1 + X3")
#* adjusted with both prognostic factors and an un-prognostic factor X8
formula3 <- as.formula("Y ~ trt + X1 + X3 + X8")

#* num for bootstrap replicates
nsim <- 5000

for (type in c("OR", "RD")) {
  res0 <- summary_Estimate(data = dat, formula0, nsim = nsim, trt.var = "trt", type = type)
  res1 <- summary_Estimate(data = dat, formula1, nsim = nsim, trt.var = "trt", type = type)
  res2 <- summary_Estimate(data = dat, formula2, nsim = nsim, trt.var = "trt", type = type)
  res3 <- summary_Estimate(data = dat, formula3, nsim = nsim, trt.var = "trt", type = type)

  rbind(res0, res1, res2, res3) %>%
    write.csv(file.path("conditional_marginal/output", paste0("Result_condi_margin_", type, ".csv")))
}
