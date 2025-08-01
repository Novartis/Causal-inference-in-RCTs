#* load used libs
library(mgcv)
require(lmtest)
require(sandwich)
library(future.apply)
# uncommented the next line for parallel
# plan(multisession, workers = 2) ## Run in parallel


#' Calculate conditional and marginal estimates with SE and p-vals
#' 
#' This is the wrapper function for calculating both conditional and marginal estimates.
#' It utilized the separate functions including
#' * marginal.OR to get marginal estimate
#' * marginal_delta for delta method SE 
#' * marginal.Stats for parametric bootstrap SE
#' * mar_est and marginal.Stats_npboot for non-parametric bootstrap SE
#' 
#' 
#' @param data input data
#' @param formula formula including response and adjusted covariate, e.g. y ~ trt + X1
#' @param trt.var name of treatment variable in the formula
#' @param type type of estimate, "OR", "RD", "RR"
#' @param nsim bootstrap replication number
#'
#' @return data.frame,   including used formula, conditional estimates
#'  with (robust) SE and p-val, marginal estimates with SEs and pvals
#'  from different methods (parametric bootstrap, non-parametric bootstrap and delta).
#'  For the two bootstrap methods, p-values were calculated using both z-score and tail-pvalue.
#' @export
#'
#' @examples summary_Estimate(data = dat, formula = "y ~ trt + X1", nsim = 100, trt.var = "trt", type = "OR")
summary_Estimate <- function(data, formula, trt.var = "trt", type = "OR", nsim = 5000) {
  trt_n <- levels(as.data.frame(data)[, trt.var])
  cat("calculating conditional estimate\n")
  m0_error <- FALSE
  if (type == "RD") {
    m0 <- try(glm(formula = formula, data = data, family = binomial("identity")))
  } else if (type == "RR") {
    m0 <- try(glm(formula = formula, data = data, family = binomial("log")))
  } else {
    m0 <- glm(formula = formula, family = binomial, data = data)
  }

  if ("try-error" %in% class(m0)) {
    m0 <- glm(formula = formula, family = binomial, data = data)
    m0_error <- TRUE
  }
  #* adjusted
  sandwich <- sandwich(m0)
  condi_robust <- coeftest(m0, sandwich)

  condi_robust <- as.data.frame(condi_robust[, c(1, 2, 4)]) %>%
    tibble::rownames_to_column("term") %>%
    dplyr::select(term, condi_rb_pval = `Pr(>|z|)`)

  model.df <- broom::tidy(m0)

  if (type == "OR") {
    condi_res <- model.df %>%
      mutate(
        or = exp(estimate), # Odds ratio/gradient
        var.diag = diag(vcov(m0)), # Variance of each coefficient
        conditional_or.se = sqrt(or^2 * var.diag),
        var.diag_rb = diag(sandwich), # Variance of each coefficient
        condi_rb_or.se = sqrt(or^2 * var.diag_rb)
      ) %>%
      left_join(condi_robust, by = "term") %>%
      dplyr::select(term,
        condi_estimate = or, conditional_se = conditional_or.se,
        conditional_pval = p.value,
        condi_rb_se = condi_rb_or.se, condi_rb_pval)
  } else {
    condi_res <- model.df %>%
      dplyr::select(term, condi_estimate = estimate,
                    conditional_se = std.error,
                    conditional_pval = p.value) %>%
      mutate(condi_rb_se = sqrt(diag(sandwich))) %>%
      left_join(condi_robust, by = "term")
  }

  if (m0_error) {
    condi_res[, c(2, 3, 5)] <- NA
    cat("cannot get conditional estimate, please use the marginal estimate\n")
  }

  #* marginal estimate - standardized
  cat("calculating marginal estimate\n")
  res_0 <- marginal_delta(data, formula = formula, trt.pbo = trt_n[1], trt = trt.var)

  #* delta method
  cat("calculating delta-method SE\n")
  if (type == "OR") {
    delta_res <- res_0 %>% dplyr::select(term,
      delta_Estimate = odds_ratio,
      delta_SE = se_delta_or, delta_pval = pval
    )
  } else if (type == "RD") {
    delta_res <- res_0 %>% dplyr::select(term,
      delta_Estimate = risk_difference,
      delta_SE = se_delta_rd, delta_pval = pval_rd
    )
  } else if (type == "RR") {
    delta_res <- res_0 %>% dplyr::select(term,
      delta_Estimate = risk_ratio,
      delta_SE = se_delta_rr, delta_pval = pval
    )
  }

  #* parametric bootstrap
  cat("calculating parametric bootstrap SE\n")
  m0 <- glm(formula = formula, family = binomial, data = data)
  m0marg <- marginal.Stats(nsim = nsim, mod = m0, data = data, trt.var = trt.var, type = type, robust = TRUE)
  parametric_res <- m0marg[-4]
  colnames(parametric_res) <- c("term", paste0("parametric_", colnames(parametric_res)[-1]))

  #* non parametric bootstrap
  cat("calculating non-parametric bootstrap SE\n")
  nonpar_res <- marginal.Stats_npboot(
    nsim = nsim, data = data,
    trt.var = trt.var, trt.ref = trt_n[1], trt.trt = trt_n[-1], type = type,
    formula = formula
  )

  nonpar_res <- nonpar_res[-4]
  colnames(nonpar_res) <- c("term", paste0("nonparametric_", colnames(nonpar_res)[2:5]))

  #* combine all result
  tmp <- condi_res %>%
    filter(stringr::str_detect(term, trt.var)) %>%
    mutate(term = gsub(trt.var, "", term)) %>%
    left_join(parametric_res, by = "term") %>%
    left_join(nonpar_res, by = "term") %>%
    left_join(delta_res, by = "term") %>%
    mutate(formula_cov = as.character(Reduce(paste, deparse(formula)))) %>%
    dplyr::select(formula = formula_cov,
                  treatment.var = term,
                  conditional.estimate = condi_estimate, conditional.se = conditional_se,
                  conditional.pval = conditional_pval,
                  conditional.se.robust = condi_rb_se, conditional.pval.robust = condi_rb_pval,
                  marginal.estimate = parametric_Estimate,
                  parametric.se = parametric_se, parametric.pval = parametric_Pval,
                  parametric.tail.pval = parametric_tailP,
                  nonparametric.se = nonparametric_se,nonparametric.pval = nonparametric_Pval,
                  nonparametric.tail.pval = nonparametric_tailP,
                  delta.se = delta_SE, delta.pval = delta_pval)
  return(tmp)
}


#' marginal_delta function for binary data using delta method for SE
#'
#' @param data input data
#' @param formula formula including response and adjusted covariate, e.g. y ~ trt + X1
#' @param trt name of treatment variable in the formula
#' @param trt.pbo ref level of treatment variable
#' @return data.frame including OR, RD, RR and corresponding SE, p-val for different treatment levels
#' @export
#'
#' @examples marginal_delta(data, formula = "y ~ trt", trt = "trt", trt.pbo = "Placebo")
marginal_delta <- function(data,
                            formula,
                            trt = "trt",
                            trt.pbo = "Placebo") {
  glmfit <- glm(formula, data = data, family = binomial(link = "logit"))
  
  df <- glmfit$model
  # if(!all(df[,2]%in%c(0,1)))stop("The treatment variable is not valid")
  
  ### estimate the proportion respectively
  
  coef_hat <- coef(glmfit)
  
  trt_n <- levels(as.data.frame(data)[, trt])
  df_list <- list()
  for (trt.l in trt_n) {
    df[, 2] <- trt.l
    df_list <- c(df_list, list(df))
  }
  df_list <- do.call("rbind", df_list)
  df_list[, 2] <- factor(df_list[, 2], levels = trt_n)
  mat_TC <- model.matrix(glmfit$formula, data = df_list)
  
  ### assume all subjects are in control group
  index_pbo <- which(df_list[, 2] == trt.pbo)
  mat_C <- mat_TC[index_pbo, ]
  Prob_control <- plogis(mat_C %*% coef_hat)
  Prob_c_mean <- mean(Prob_control)
  
  trt.list <- list()
  
  for (trt.trt in setdiff(trt_n, trt.pbo)) {
    ### assume all subjects are in treatment group
    index_trt <- which(df_list[, 2] == trt.trt)
    mat_T <- mat_TC[index_trt, ]
    Prob_treat <- plogis(mat_T %*% coef_hat)
    Prob_t_mean <- mean(Prob_treat)
    
    ### calculate target estimands
    odds_ratio <- (Prob_t_mean / (1 - Prob_t_mean)) / (Prob_c_mean / (1 -
                                                                        Prob_c_mean))
    risk_difference <- Prob_t_mean - Prob_c_mean
    risk_ratio <- Prob_t_mean / Prob_c_mean
    
    coef_cov <- vcov(glmfit)
    
    n <- nrow(mat_T)
    A <- Prob_t_mean / (1 - Prob_t_mean)
    B <- (1 - Prob_c_mean) / Prob_c_mean
    derive_or <- (B / (1 - Prob_t_mean)^2 * (t(Prob_treat * (1 - Prob_treat)) %*%
                                               mat_T) / n -
                    A / (Prob_c_mean)^2 * (t(Prob_control * (1 - Prob_control)) %*%
                                             mat_C) / n)
    
    var_delta_or <- derive_or %*% coef_cov %*% t(derive_or)
    se_delta_or <- sqrt(var_delta_or)
    
    ss <- Prob_treat / (Prob_t_mean * (1 - Prob_t_mean)) - Prob_control /
      (Prob_c_mean * (1 - Prob_c_mean))
    se_delta_or_eq <- sqrt(var_delta_or + var(ss) / n)
    
    derive_rd <- t(Prob_treat * (1 - Prob_treat)) %*% mat_T / n - t(Prob_control *
                                                                      (1 - Prob_control)) %*% mat_C / n
    se_delta_rd <- sqrt(derive_rd %*% coef_cov %*% t(derive_rd))
    
    derive_rr <- (
      Prob_c_mean * t(Prob_treat * (1 - Prob_treat)) %*% mat_T / n - Prob_t_mean *
        t(Prob_control * (1 - Prob_control)) %*% mat_C / n
    ) / Prob_c_mean^2
    se_delta_rr <- sqrt(derive_rr %*% coef_cov %*% t(derive_rr))
    
    
    pval <- 2 * min(pnorm(log(odds_ratio) / (se_delta_or / odds_ratio)), pnorm(-log(odds_ratio) /
                                                                                 (se_delta_or / odds_ratio)))
    pval_rd <- 2 * min(
      pnorm(risk_difference / se_delta_rd),
      pnorm(-risk_difference / se_delta_rd)
    )
    
    
    tmp <- data.frame(
      term = trt.trt,
      odds_ratio = odds_ratio,
      risk_difference = risk_difference,
      risk_ratio = risk_ratio,
      se_delta_or = se_delta_or,
      se_delta_rd = se_delta_rd,
      se_delta_rr = se_delta_rr,
      pval = pval,
      pval_rd = pval_rd
    )
    trt.list <- c(trt.list, list(tmp))
  }
  trt.list <- do.call("rbind", trt.list)
  return(trt.list)
}

#' marginal estimate
#'
#' @param mod fitted model
#' @param data input data frame, including used covariates
#' @param trt.var name of treatment variable
#' @param trt.ref reference level of trt
#' @param trt.trt active treatment level
#' @param type type of estimator, OR, RD, RR, LOR
#'
#' @return scalar, marginal estimate
#' @export
#'
#' @examples
marginal.OR <- function(mod, data, trt.var = "trt", trt.ref = 0, trt.trt = 1, type = "LOR") {
  data <- as.data.frame(data)
  pred.data <- data[data[, trt.var] %in% c(trt.ref, trt.trt), ]
  pred.data <- data
  pred.data[, trt.var] <- trt.trt
  preds1 <- predict(mod, newdata = pred.data, type = "response")

  pred.data[, trt.var] <- trt.ref
  preds2 <- predict(mod, newdata = pred.data, type = "response")

  (pr1 <- mean(preds1, na.rm = T))
  (pr2 <- mean(preds2, na.rm = T))

  if (type == "RD") res <- pr1 - pr2
  if (type == "RR") res <- pr1 / pr2
  if (type == "OR") res <- (pr1 / (1 - pr1)) / (pr2 / (1 - pr2))
  if (type == "LOR") res <- log((pr1 / (1 - pr1)) / (pr2 / (1 - pr2)))

  return(res)
}


#' Parametric bootstrap for marginal SE
#'
#' @param nsim replicates of bootstrap
#' @param mod fitted model
#' @param data input data frame including used covariates
#' @param trt.var name of treatment variable
#' @param type type of estimator, OR, RD, RR, LOR
#' @param robust if robust estimator is used during generating parametric bootstrap sample
#'
#' @return data.frame, including treatment levels and marginal estimate, parametric bootstrap SE,
#' z-statistic and associated p-val, bootstrap tail p-value
#' @export
#'
#' @examples
marginal.Stats <- function(nsim = 10000, mod, data, trt.var = "trt", type = "LOR", robust = FALSE) {

  # use robust variance
  if (robust) {
    v_mat <- sandwich(mod)
  } else {
    v_mat <- vcov(mod)
  }
  # set.seed(2022)
  simmat <- mgcv::rmvn(nsim, coef(mod), v_mat)

  cond_p <- colMeans(simmat > 0)
  cond_p <- unlist(lapply(cond_p, function(x) 2 * min(x, 1 - x)))
  names(cond_p) <- gsub(trt.var, "", names(cond_p))

  trt_n <- levels(as.data.frame(data)[, trt.var])

  re <- list()
  for (trt_v in trt_n[-1]) {
    Thetaout <- future_apply(simmat, 1, function(x) {
      tmpmod <- mod
      tmpmod$coefficients <- x
      marginal.OR(mod = tmpmod, data = data, trt.var = trt.var, trt.ref = trt_n[1], trt.trt = trt_v, type = type)
    })

    Est <- marginal.OR(mod = mod, data = data, trt.var = trt.var, trt.ref = trt_n[1], trt.trt = trt_v, type = type)
    SEout <- sd(Thetaout)
    if (type %in% c("LOR", "RD")) {
      Zscore <- Est / sd(Thetaout)
    } else {
      Zscore <- log(Est) / sd(log(Thetaout))
    }
    Pval <- 2 * min(pnorm(Zscore), pnorm(-Zscore))

    # add monte-carlo tail p-value
    if (type == "RD") tailP1 <- mean(Thetaout > 0)
    if (type == "RR") tailP1 <- mean(Thetaout > 1)
    if (type == "OR") tailP1 <- mean(Thetaout > 1)
    if (type == "LOR") tailP1 <- mean(Thetaout > 0)

    tailP2 <- 2 * min(tailP1, 1 - tailP1)

    tmp <- data.frame(term = trt_v, Estimate = Est, se = SEout, Z = Zscore, Pval = Pval, tailP = tailP2)
    re <- c(re, list(tmp))
  }

  do.call(rbind, re) %>% mutate(cond_bt_p = cond_p[trt_n[-1]])
}


#' Marginal estimate function for bootstrap
#'
#' @param formula y ~ trt+X1
#' @param data input data frame, including treatment and used covariates
#' @param indices
#' @param trt.var name of treatment variable in the formula
#' @param trt.ref ref level of treatment variable
#' @param trt.trt active level of treatment variable
#' @param type type of estimator, OR, RD, RR, LOR
#'
#' @return scalar, marginal estimate
#' @export
#'
#' @examples
mar_est <- function(formula, data, indices, trt.var = "trt", trt.ref = 0, trt.trt = 1, type = "LOR") {
  data_bt <- data[indices, ] # selecting sample with boot
  fit <- glm(formula, data = data_bt, family = binomial)
  re0 <- sapply(trt.trt, function(x) {
    marginal.OR(fit,
      data = data_bt, trt.var = trt.var,
      trt.ref = trt.ref, trt.trt = x, type = type
    )
  })
  return(re0)
}


#' Non-parametric bootstrap for marginal SE
#'
#' @param nsim replicates of bootstrap
#' @param data input data frame including used covariates
#' @param trt.var name of treatment variable in the formula
#' @param trt.ref ref level of treatment variable
#' @param trt.trt active level of treatment variable
#' @param type type of estimator, OR, RD, RR, LOR
#' @param formula y ~ trt+X1
#'
#' @return data.frame, including treatment levels and marginal estimate, non-parametric bootstrap SE,
#' z-statistic and associated p-val, bootstrap tail p-value
#' @examples
marginal.Stats_npboot <- function(nsim = 10000, data, trt.var = "trt", trt.ref = 0, trt.trt = 1, type = "LOR", formula) {
  res <- boot::boot(
    data = data, statistic = mar_est, R = nsim,
    formula = formula, trt.var = trt.var, trt.ref = trt.ref,
    trt.trt = trt.trt, type = type, parallel = "multicore"
  )

  Est <- res$t0
  Thetaout <- res$t
  SEout <- apply(Thetaout, 2, sd)
  if (type %in% c("RD", "LOR")) {
    Zscore <- Est / SEout
  } else {
    Zscore <- log(Est) / apply(log(Thetaout), 2, sd)
  }

  Pval <- sapply(Zscore, function(x) {
    (2 * min(pnorm(x), pnorm(-x)))
  })


  # add monte-carlo tail p-value
  if (type %in% c("RD", "LOR")) {
    mu <- 0
  } else if (type %in% c("RR", "OR")) {
    mu <- 1
  }

  tailP <- apply(Thetaout, 2, function(x) {
    2 * min(mean(x > mu), 1 - mean(x > mu))
  })
  data.frame(term = names(Est), Estimate = Est, se = SEout, Z = Zscore, Pval = Pval, tailP = tailP)
}
