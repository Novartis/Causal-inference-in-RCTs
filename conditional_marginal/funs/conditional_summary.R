#' Calculate conditional treatment effect estimates with SEs and p-vals
#' 
#' @description Calculate conditional risk difference, conditional risk ratio,
#' or conditional odds ratio.
#' 
#' @param data input data
#' @param formula formula including response, treatment, and adjustment 
#' covariate(s), e.g., y ~ trt + X1
#' @param type type of estimate, "OR", "RD", "RR"
#'
#' @return data.frame containing the following estimates for each term in the 
#' model: conditional effect estimate, SE of conditional effect estimate,
#' p-value of conditional effect estimate, 
#' robust SE of conditional effect estimate, and
#' robust p-value of conditional effect estimate.
#' @export
#'
#' @examples conditional_summary(data = dat, formula = "y ~ trt + X1", type = "OR")
conditional_summary <- function(data, formula, type = "OR") {

  # Fit model to target conditional risk difference, risk ratio, or odds ratio
  if (type == "RD") {
    m0 <- lm(formula = formula, data = data)
  } else if (type == "RR") {
    m0 <- try(glm(formula = formula, data = data, family = binomial("log")))
  } else if (type == "OR") {
    m0 <- glm(formula = formula, family = binomial, data = data)
  }

  if ("try-error" %in% class(m0)) {
    stop("cannot get conditional estimate, please use the marginal estimate\n")
  }

  # Robust sandwich estimator of variance-covariance matrix
  sandwich <- sandwich(m0)

  # Get robust p-values for coefficients, based on robust sandwich estimator
  # of standard errors
  condi_robust <- coeftest(m0, sandwich)

  if(type %in% c("RR", "OR")) {
    condi_robust <- as.data.frame(condi_robust[, c(1, 2, 4)]) %>%
      tibble::rownames_to_column("term") %>%
      dplyr::select(term, condi_rb_pval = `Pr(>|z|)`)
  } else if(type == "RD") {
    condi_robust <- as.data.frame(condi_robust[, c(1, 2, 4)]) %>%
      tibble::rownames_to_column("term") %>%
      dplyr::select(term, condi_rb_pval = `Pr(>|t|)`)
  }
  
  # Convert model output to data frame
  model.df <- broom::tidy(m0)

  # Extra point estimate, SE of point estimate, p-value of point estimate,
  # robust SE of point estimate, and robust p-value of point estimate for
  # each term
  if (type == "OR") {

    # Convert treatment effect estimates from log OR scale to OR scale.
    # Get summary statistics on OR scale.
    condi_results <- model.df %>%
      mutate(
        or = exp(estimate), # Odds ratio/gradient
        var.diag = diag(vcov(m0)), # Variance of OR of each coefficient
        conditional_or.se = sqrt(or^2 * var.diag), # SE of OR of each coefficient
        var.diag_rb = diag(sandwich), # Robust variance of OR of each coefficient
        condi_rb_or.se = sqrt(or^2 * var.diag_rb) # Robust SE of OR of each coefficient
      ) %>%
      left_join(condi_robust, by = "term") %>%
      dplyr::select(term, # Term (name of variable)
        condi_estimate = or, # Conditional OR estimate
        conditional_se = conditional_or.se, # SE of conditional OR estimate
        conditional_pval = p.value, # p-value of conditional OR estimate
        condi_rb_se = condi_rb_or.se, # Robust SE of conditional OR estimate
        condi_rb_pval) # Robust p-value of conditional OR estimate

  } else {

    # Get summary statistics of treatment effects
    condi_results <- model.df %>%
      dplyr::select(term, # Term (name of variable)
                    condi_estimate = estimate, # Conditional RD or RR estimate
                    conditional_se = std.error, # SE of conditional estimate
                    conditional_pval = p.value) %>% # p-value of conditional estimate
      mutate(condi_rb_se = sqrt(diag(sandwich))) %>% # Robust SE of conditional estimate
      left_join(condi_robust, by = "term") # Join to get robust p-value of conditional estimate

  }

  return(condi_results)
} 