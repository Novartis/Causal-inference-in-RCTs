# Compute G-computation estimate for heart transplant example.

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

# FIT outcome model (conditional mean of Y at each (X, Z) combination)
outcome_mod <- glm(Y ~ X*Z, data = heart_df, family = "binomial")

# CHECK outcome model predictions on all combinations of (X, Z)

## Dataset with all combinations of (X, Z)
example_df <- data.frame(Z = c(1, 1, 0, 0), X = c(1, 0, 1, 0))

rownames(example_df) <- 
  paste0("(Z == ", example_df$Z, ", X == ", example_df$X, ")")

## Predictions at all combinations of (X, Z)
predict(outcome_mod, type = "response", newdata = example_df)

## Check that predictions at all combinations of (X, Z) have correct values
all.equal(unname(predict(outcome_mod, type = "response",
                         newdata = example_df)), 
          c(3/9, 1/4, 2/3, 2/4))

# COMPUTE G-computation estimate

## Predictions under Z = 1
heart_df$Q1 <- predict(outcome_mod, type = "response", 
                       newdata = data.frame(Z = 1, X = heart_df$X))

## Predictions under Z = 0
heart_df$Q0 <- predict(outcome_mod, type = "response",
                       newdata = data.frame(Z = 0, X = heart_df$X))

## Compute G-computation estimate of ATE
est_gcomp <- mean(heart_df$Q1) - mean(heart_df$Q0)

## Check that G-computation estimate of ATE equals -0.3
all.equal(round(est_gcomp, 2), -0.3)

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

  # Fit outcome model (conditional mean of Y at each (X, Z) combination)
  boot_outcome_mod <- glm(Y ~ X*Z, data = boot_df, family = "binomial")

  # Predictions under Z = 1
  boot_df$Q1 <- predict(boot_outcome_mod, type = "response", 
                        newdata = data.frame(Z = 1, X = boot_df$X))

  # Predictions under Z = 0
  boot_df$Q0 <- predict(boot_outcome_mod, type = "response",
                        newdata = data.frame(Z = 0, X = boot_df$X))

  # Compute G-computation estimate of ATE
  boot_vec[b] <- mean(boot_df$Q1) - mean(boot_df$Q0)

}

# Percentile bootstrapped CI
round(quantile(boot_vec, c(0.025, 0.975)), 2)
