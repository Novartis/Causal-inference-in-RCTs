# Compute IPW estimate for heart transplant example.

library(data.table)

# Create long version of heart transplant example data
heart_df <- data.table(X = c(rep(1, 12), rep(0, 8)),
                       Y = c(rep(1, 5), rep(0, 7), rep(1, 3), rep(0, 5)),
                       Z = c(rep(1, 3), rep(0, 2), rep(1, 6), rep(0, 1),
                             rep(1, 1), rep(0, 2), rep(1, 3), rep(0, 2)))

# Check counts of each (X, Y, Z) combination
stopifnot(nrow(heart_df[X == 1 & Y == 1 & Z == 1]) == 3)
stopifnot(nrow(heart_df[X == 1 & Y == 1 & Z == 0]) == 2)
stopifnot(nrow(heart_df[X == 1 & Y == 0 & Z == 1]) == 6)
stopifnot(nrow(heart_df[X == 1 & Y == 0 & Z == 0]) == 1)
stopifnot(nrow(heart_df[X == 0 & Y == 1 & Z == 1]) == 1)
stopifnot(nrow(heart_df[X == 0 & Y == 1 & Z == 0]) == 2)
stopifnot(nrow(heart_df[X == 0 & Y == 0 & Z == 1]) == 3)
stopifnot(nrow(heart_df[X == 0 & Y == 0 & Z == 0]) == 2)

##########################
##### Point estimate #####
##########################

# FIT propensity score model (conditional mean of Z at each value of X)
prop_score_mod <- glm(Z ~ X, data = heart_df, family = "binomial")

# CHECK propensity score predictions at all levels of X

## Dataset with both levels of X
example_df <- data.frame(X = c(1, 0))

rownames(example_df) <- c("X == 1", "X == 0")

## Predictions at all X
predict(prop_score_mod, type = "response",
        newdata = example_df)

## Check that predictions at all X have correct values
all.equal(unname(predict(prop_score_mod, type = "response",
                         newdata = example_df)), 
          c(3/4, 1/2))

# COMPUTE IPW estimate

## Get propensity score P(Z = 1 | X) for each observation
heart_df$pi_hat <- predict(prop_score_mod, newdata = data.frame(X = heart_df$X),
                           type = "response")

## Estimate E[Y(1)] and E[Y(0)]
heart_df[, EY1_comp := Z * Y / pi_hat]
heart_df[, EY0_comp := (1 - Z) * Y / (1 - pi_hat)]

EY1_est <- mean(heart_df$EY1_comp)
EY0_est <- mean(heart_df$EY0_comp)

## IPW estimate of ATE
est_ipw <- EY1_est - EY0_est

## Check that IPW estimate of ATE equals -0.3
all.equal(round(est_ipw, 2), -0.3)

###########################
##### Bootstrapped CI #####
###########################

# Set seed
set.seed(20230113)

# Number of bootstrapped samples
n_boot <- 1000

# Vector to store bootstrapped point estimates
boot_vec <- rep(NA, n_boot)

for(b in 1:n_boot) {
  
  # Randomly sample indices
  boot_ind <- sample(1:nrow(heart_df), size = nrow(heart_df), replace = TRUE)

  # Get bootstrapped sample
  boot_df <- heart_df[boot_ind]

  # Fit propensity score model (conditional mean of Z at each value of X)
  boot_prop_score_mod <- glm(Z ~ X, data = boot_df, family = "binomial")

  # Get propensity score P(Z = 1 | X) for each observation
  boot_df$pi_hat <- predict(boot_prop_score_mod, 
                            newdata = data.frame(X = boot_df$X),
                            type = "response")

  # Estimate E[Y(1)] and E[Y(0)]
  boot_df[, EY1_comp := Z * Y / pi_hat]
  boot_df[, EY0_comp := (1 - Z) * Y / (1 - pi_hat)]

  EY1_est <- mean(boot_df$EY1_comp)
  EY0_est <- mean(boot_df$EY0_comp)

  # Compute IPW estimate of ATE
  boot_vec[b] <- EY1_est - EY0_est

}

# Percentile bootstrapped CI
round(quantile(boot_vec, c(0.025, 0.975)), 2)
