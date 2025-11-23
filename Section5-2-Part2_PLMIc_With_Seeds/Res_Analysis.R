library(ggplot2)
library(tidyr)
library(patchwork)
library(purrr)
library(tidyverse)
library(testit)
source(
  "~/BayesianRecordLinkage/helper_fns/helper_fns_R/Rubin_interv.R"
)
source(
  "~/BayesianRecordLinkage/helper_fns/helper_fns_R/getHpi.R"
)

# looping variable --------------------------------------------------------
# looping_grid <- expand.grid(0:99, c(1,3), c(0.5), c(0.3,0.9), 1:5, c(0.1, 0.05, 0.01))
looping_grid <- rbind(
  expand.grid(0:99, c(1), c(0.5), c(0.3), 1:5, c(0.1, 0.05, 0.01)),
  expand.grid(0:99, c(3), c(0.5), c(0.9), 1:5, c(0.1, 0.05, 0.01))
)
MCMC_length <- 900

# Confidence Alpha --------------------------------------------------------
alpha.conf <- 0.1
# Inference of models -----------------------------------------------------
inference_df <- map_dfr(
  1:nrow(looping_grid),
  function(loop_var) {
    dataset <- looping_grid[loop_var, 1]
    n_error <- looping_grid[loop_var, 2]
    overlap_proportion <- looping_grid[loop_var, 3]
    R_sq <- looping_grid[loop_var, 4]
    sub_id <- looping_grid[loop_var, 5]
    seed_proportion <- looping_grid[loop_var, 6]

    true_beta_0 <- 3
    true_beta_1 <- 3

    # EM and TS.OLS ----------------------------------------------------
    result_file_name <- paste0(
      "RegResult_dataset",
      dataset,
      "(",
      sub_id,
      ")",
      "_nerrors",
      n_error,
      "_percdups",
      overlap_proportion * 100,
      "_percRsq",
      R_sq * 100,
      "_seedproportion",
      seed_proportion * 100,
      "_modelEM_TSols_Perfect.RData"
    )

    load(
      paste0(
        "~/BayesianRecordLinkage/Section5-2-Part2_PLMIc_With_Seeds/Simu_Result/",
        result_file_name
      ),
      temp_env <- new.env()
    )
    result_data <- as.list(temp_env)

    # TS.OLS regression
    beta_0.TS.OLS.estimates <- sapply(1:MCMC_length, function(i) {
      result_data[["res"]][["EM_and_TSols_result"]][[i]][[
        "TS.OLS_coefficients"
      ]][1]
    })
    beta_0.TS.OLS.var <- sapply(1:MCMC_length, function(i) {
      result_data[["res"]][["EM_and_TSols_result"]][[i]][["TS.OLS_cov"]][1, 1]
    })
    beta_0.TS.OLS_rubin_interv <- Rubin_interv(
      beta_0.TS.OLS.estimates,
      beta_0.TS.OLS.var,
      alpha = alpha.conf
    )
    beta_0.TS.OLS.lb <- beta_0.TS.OLS_rubin_interv[1]
    beta_0.TS.OLS.ub <- beta_0.TS.OLS_rubin_interv[2]
    beta_0.TS.OLS.pt <- mean(beta_0.TS.OLS.estimates)
    beta_0.TS.OLS.length <- beta_0.TS.OLS.ub - beta_0.TS.OLS.lb
    beta_0.TS.OLS.cover <- ((beta_0.TS.OLS.ub >= true_beta_0) &
      (beta_0.TS.OLS.lb <= true_beta_0))
    beta_0.TS.OLS.SE <- (beta_0.TS.OLS.pt - true_beta_0)^2
    beta_0.TS.OLS.bias <- abs(beta_0.TS.OLS.pt - true_beta_0)

    beta_1.TS.OLS.estimates <- sapply(1:MCMC_length, function(i) {
      result_data[["res"]][["EM_and_TSols_result"]][[i]][[
        "TS.OLS_coefficients"
      ]][2]
    })
    beta_1.TS.OLS.var <- sapply(1:MCMC_length, function(i) {
      result_data[["res"]][["EM_and_TSols_result"]][[i]][["TS.OLS_cov"]][2, 2]
    })
    beta_1.TS.OLS_rubin_interv <- Rubin_interv(
      beta_1.TS.OLS.estimates,
      beta_1.TS.OLS.var,
      alpha = alpha.conf
    )
    beta_1.TS.OLS.lb <- beta_1.TS.OLS_rubin_interv[1]
    beta_1.TS.OLS.ub <- beta_1.TS.OLS_rubin_interv[2]
    beta_1.TS.OLS.pt <- mean(beta_1.TS.OLS.estimates)
    beta_1.TS.OLS.length <- beta_1.TS.OLS.ub - beta_1.TS.OLS.lb
    beta_1.TS.OLS.cover <- ((beta_1.TS.OLS.ub >= true_beta_1) &
      (beta_1.TS.OLS.lb <= true_beta_1))
    beta_1.TS.OLS.SE <- (beta_1.TS.OLS.pt - true_beta_1)^2
    beta_1.TS.OLS.bias <- abs(beta_1.TS.OLS.pt - true_beta_1)

    # Perfect OLS results
    perfect.coeffs <- result_data[["res"]][["Perfect.OLS_result"]][[
      "Perfect.OLS_coefficients"
    ]]
    perfect.cov <- result_data[["res"]][["Perfect.OLS_result"]][[
      "Perfect.OLS_cov"
    ]]

    beta_0.perfect.pt <- perfect.coeffs[1]
    beta_0.perfect.lb <- beta_0.perfect.pt -
      qt(1 - alpha.conf / 2, df = 500 * overlap_proportion - 2) *
        sqrt(perfect.cov[1, 1])
    beta_0.perfect.ub <- beta_0.perfect.pt +
      qt(1 - alpha.conf / 2, df = 500 * overlap_proportion - 2) *
        sqrt(perfect.cov[1, 1])
    beta_0.perfect.length <- beta_0.perfect.ub - beta_0.perfect.lb
    beta_0.perfect.SE <- (beta_0.perfect.pt - true_beta_0)^2
    beta_0.perfect.bias <- abs(beta_0.perfect.pt - true_beta_0)
    beta_0.perfect.cover <- (beta_0.perfect.ub >= true_beta_0) &
      (beta_0.perfect.lb <= true_beta_0)

    beta_1.perfect.pt <- perfect.coeffs[2]
    beta_1.perfect.lb <- beta_1.perfect.pt -
      qt(1 - alpha.conf / 2, df = 500 * overlap_proportion - 2) *
        sqrt(perfect.cov[2, 2])
    beta_1.perfect.ub <- beta_1.perfect.pt +
      qt(1 - alpha.conf / 2, df = 500 * overlap_proportion - 2) *
        sqrt(perfect.cov[2, 2])
    beta_1.perfect.length <- beta_1.perfect.ub - beta_1.perfect.lb
    beta_1.perfect.SE <- (beta_1.perfect.pt - true_beta_1)^2
    beta_1.perfect.bias <- abs(beta_1.perfect.pt - true_beta_1)
    beta_1.perfect.cover <- (beta_1.perfect.ub >= true_beta_1) &
      (beta_1.perfect.lb <= true_beta_1)

    # EM results
    beta_0.EM.estimates <- sapply(1:MCMC_length, function(i) {
      result_data[["res"]][["EM_and_TSols_result"]][[i]][["EM_result"]][[
        "par"
      ]][1]
    })
    beta_0.EM.var <- sapply(1:MCMC_length, function(i) {
      result_data[["res"]][["EM_and_TSols_result"]][[i]][["EM_result"]][[
        "laplace_variance"
      ]][1, 1]
    })
    beta_0_rubin_interv <- Rubin_interv(
      beta_0.EM.estimates,
      beta_0.EM.var,
      alpha = alpha.conf
    )
    beta_0.EM.lb <- beta_0_rubin_interv[1]
    beta_0.EM.ub <- beta_0_rubin_interv[2]
    beta_0.EM.pt <- mean(beta_0.EM.estimates)
    beta_0.EM.length <- beta_0.EM.ub - beta_0.EM.lb
    beta_0.EM.cover <- ((beta_0.EM.ub >= true_beta_0) &
      (beta_0.EM.lb <= true_beta_0))
    beta_0.EM.SE <- (beta_0.EM.pt - true_beta_0)^2
    beta_0.EM.bias <- abs(beta_0.EM.pt - true_beta_0)

    beta_1.EM.estimates <- sapply(1:MCMC_length, function(i) {
      result_data[["res"]][["EM_and_TSols_result"]][[i]][["EM_result"]][[
        "par"
      ]][2]
    })
    beta_1.EM.var <- sapply(1:MCMC_length, function(i) {
      result_data[["res"]][["EM_and_TSols_result"]][[i]][["EM_result"]][[
        "laplace_variance"
      ]][2, 2]
    })
    beta_1_rubin_interv <- Rubin_interv(
      beta_1.EM.estimates,
      beta_1.EM.var,
      alpha = alpha.conf
    )
    beta_1.EM.lb <- beta_1_rubin_interv[1]
    beta_1.EM.ub <- beta_1_rubin_interv[2]
    beta_1.EM.pt <- mean(beta_1.EM.estimates)
    beta_1.EM.length <- beta_1.EM.ub - beta_1.EM.lb
    beta_1.EM.cover <- ((beta_1.EM.ub >= true_beta_1) &
      (beta_1.EM.lb <= true_beta_1))
    beta_1.EM.SE <- (beta_1.EM.pt - true_beta_1)^2
    beta_1.EM.bias <- abs(beta_1.EM.pt - true_beta_1)

    return(data.frame(
      dataset,
      n_error,
      overlap_proportion,
      R_sq,
      sub_id,
      seed_proportion,
      true_beta_0,
      true_beta_1,
      beta_0.TS.OLS.pt,
      beta_0.TS.OLS.length,
      beta_0.TS.OLS.lb,
      beta_0.TS.OLS.ub,
      beta_0.TS.OLS.cover,
      beta_0.TS.OLS.SE,
      beta_0.TS.OLS.bias,
      beta_1.TS.OLS.pt,
      beta_1.TS.OLS.length,
      beta_1.TS.OLS.lb,
      beta_1.TS.OLS.ub,
      beta_1.TS.OLS.cover,
      beta_1.TS.OLS.SE,
      beta_1.TS.OLS.bias,
      beta_0.EM.pt,
      beta_0.EM.length,
      beta_0.EM.lb,
      beta_0.EM.ub,
      beta_0.EM.cover,
      beta_0.EM.SE,
      beta_0.EM.bias,
      beta_1.EM.pt,
      beta_1.EM.length,
      beta_1.EM.lb,
      beta_1.EM.ub,
      beta_1.EM.cover,
      beta_1.EM.SE,
      beta_1.EM.bias,
      beta_0.perfect.pt,
      beta_0.perfect.length,
      beta_0.perfect.lb,
      beta_0.perfect.ub,
      beta_0.perfect.cover,
      beta_0.perfect.SE,
      beta_0.perfect.bias,
      beta_1.perfect.pt,
      beta_1.perfect.length,
      beta_1.perfect.lb,
      beta_1.perfect.ub,
      beta_1.perfect.cover,
      beta_1.perfect.SE,
      beta_1.perfect.bias
    ))
  },
  .progress = TRUE
)


result_beta1 <- inference_df %>%
  group_by(seed_proportion, n_error, R_sq) %>%
  summarise(
    mean_EM_cover = mean(beta_1.EM.cover, na.rm = TRUE),
    mean_TS.OLS_cover = mean(beta_1.TS.OLS.cover, na.rm = TRUE),
    mean_perfect_cover = mean(beta_1.perfect.cover, na.rm = TRUE),
    mean_EM_length = median(beta_1.EM.length, na.rm = TRUE),
    mean_TS.OLS_length = median(beta_1.TS.OLS.length, na.rm = TRUE),
    mean_perfect_length = median(beta_1.perfect.length, na.rm = TRUE),
    mean_EM_MSE = mean(beta_1.EM.SE, na.rm = TRUE),
    mean_TS.OLS_MSE = mean(beta_1.TS.OLS.SE, na.rm = TRUE),
    mean_perfect_MSE = mean(beta_1.perfect.SE, na.rm = TRUE),
    mean_EM_bias = median(beta_1.EM.bias, na.rm = TRUE),
    mean_TS.OLS_bias = median(beta_1.TS.OLS.bias, na.rm = TRUE),
    mean_perfect_bias = median(beta_1.perfect.bias, na.rm = TRUE),
    .groups = "keep"
  )

result_beta0 <- inference_df %>%
  group_by(seed_proportion, n_error, R_sq) %>%
  summarise(
    mean_EM_cover = mean(beta_0.EM.cover, na.rm = TRUE),
    mean_TS.OLS_cover = mean(beta_0.TS.OLS.cover, na.rm = TRUE),
    mean_perfect_cover = mean(beta_0.perfect.cover, na.rm = TRUE),
    mean_EM_length = median(beta_0.EM.length, na.rm = TRUE),
    mean_TS.OLS_length = median(beta_0.TS.OLS.length, na.rm = TRUE),
    mean_perfect_length = median(beta_0.perfect.length, na.rm = TRUE),
    mean_EM_MSE = mean(beta_0.EM.SE, na.rm = TRUE),
    mean_TS.OLS_MSE = mean(beta_0.TS.OLS.SE, na.rm = TRUE),
    mean_perfect_MSE = mean(beta_0.perfect.SE, na.rm = TRUE),
    mean_EM_bias = median(beta_0.EM.bias, na.rm = TRUE),
    mean_TS.OLS_bias = median(beta_0.TS.OLS.bias, na.rm = TRUE),
    mean_perfect_bias = median(beta_0.perfect.bias, na.rm = TRUE),
    .groups = "keep"
  )
