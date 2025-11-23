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
looping_grid <- expand.grid(
  0:99, # dataset
  c(1), # n_error
  c(0.5), # overlap_proportion
  c(0.05), # seed_proportion
  c(0.5), # x_2.prob
  c(1:5), # dataset_var
  c(1, 10) # beta_2
)

MCMC_length <- 900

# Confidence Alpha --------------------------------------------------------
alpha.conf <- 0.1
# Inference of models -----------------------------------------------------
inference_df <- map_dfr(
  1:nrow(looping_grid),
  function(loop_var) {
    # browser()
    dataset <- looping_grid[loop_var, 1]
    n_error <- looping_grid[loop_var, 2]
    overlap_proportion <- looping_grid[loop_var, 3]
    seed_proportion <- looping_grid[loop_var, 4]
    x_2.prob <- looping_grid[loop_var, 5]
    dataset_var <- looping_grid[loop_var, 6]
    beta_2 <- looping_grid[loop_var, 7]
    
    true_beta_0 <- 3
    true_beta_1 <- 0.5
    true_beta_2 <- beta_2
    
    result_file_name <- paste(
      "RegResult",
      "_dataset",
      dataset,
      "(",
      dataset_var,
      ")",
      "_nerrors",
      n_error,
      "_percdups",
      overlap_proportion * 100,
      '_x2prob',
      x_2.prob,
      '_beta2',
      beta_2,
      '_model',
      'EM_TSols_Perfect',
      '.RData',
      sep = ""
    )
    
    load(
      paste0(
        "~/BayesianRecordLinkage/",
        "Section5-6_Multiple_Regression/Simu_Result/",
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
    
    beta_2.TS.OLS.estimates <- sapply(1:MCMC_length, function(i) {
      result_data[["res"]][["EM_and_TSols_result"]][[i]][[
        "TS.OLS_coefficients"
      ]][3]
    })
    beta_2.TS.OLS.var <- sapply(1:MCMC_length, function(i) {
      result_data[["res"]][["EM_and_TSols_result"]][[i]][["TS.OLS_cov"]][3, 3]
    })
    beta_2.TS.OLS_rubin_interv <- Rubin_interv(
      beta_2.TS.OLS.estimates,
      beta_2.TS.OLS.var,
      alpha = alpha.conf
    )
    beta_2.TS.OLS.lb <- beta_2.TS.OLS_rubin_interv[1]
    beta_2.TS.OLS.ub <- beta_2.TS.OLS_rubin_interv[2]
    beta_2.TS.OLS.pt <- mean(beta_2.TS.OLS.estimates)
    beta_2.TS.OLS.length <- beta_2.TS.OLS.ub - beta_2.TS.OLS.lb
    beta_2.TS.OLS.cover <- ((beta_2.TS.OLS.ub >= true_beta_2) &
                              (beta_2.TS.OLS.lb <= true_beta_2))
    beta_2.TS.OLS.SE <- (beta_2.TS.OLS.pt - true_beta_2)^2
    beta_2.TS.OLS.bias <- abs(beta_2.TS.OLS.pt - true_beta_2)
    
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
    
    beta_2.perfect.pt <- perfect.coeffs[3]
    beta_2.perfect.lb <- beta_2.perfect.pt -
      qt(1 - alpha.conf / 2, df = 500 * overlap_proportion - 2) *
      sqrt(perfect.cov[3, 3])
    beta_2.perfect.ub <- beta_2.perfect.pt +
      qt(1 - alpha.conf / 2, df = 500 * overlap_proportion - 2) *
      sqrt(perfect.cov[3, 3])
    beta_2.perfect.length <- beta_2.perfect.ub - beta_2.perfect.lb
    beta_2.perfect.SE <- (beta_2.perfect.pt - true_beta_2)^2
    beta_2.perfect.bias <- abs(beta_2.perfect.pt - true_beta_2)
    beta_2.perfect.cover <- (beta_2.perfect.ub >= true_beta_2) &
      (beta_2.perfect.lb <= true_beta_2)
    
    # PLMIc results
    beta_0.PLMIc.estimates <- sapply(1:MCMC_length, function(i) {
      result_data[["res"]][["EM_and_TSols_result"]][[i]][["PLMIc_result"]][[
        "par"
      ]][1]
    })
    beta_0.PLMIc.var <- sapply(1:MCMC_length, function(i) {
      result_data[["res"]][["EM_and_TSols_result"]][[i]][["PLMIc_result"]][[
        "laplace_variance"
      ]][1, 1]
    })
    beta_0_rubin_interv <- Rubin_interv(
      beta_0.PLMIc.estimates,
      beta_0.PLMIc.var,
      alpha = alpha.conf
    )
    beta_0.PLMIc.lb <- beta_0_rubin_interv[1]
    beta_0.PLMIc.ub <- beta_0_rubin_interv[2]
    beta_0.PLMIc.pt <- mean(beta_0.PLMIc.estimates)
    beta_0.PLMIc.length <- beta_0.PLMIc.ub - beta_0.PLMIc.lb
    beta_0.PLMIc.cover <- ((beta_0.PLMIc.ub >= true_beta_0) &
                          (beta_0.PLMIc.lb <= true_beta_0))
    beta_0.PLMIc.SE <- (beta_0.PLMIc.pt - true_beta_0)^2
    beta_0.PLMIc.bias <- abs(beta_0.PLMIc.pt - true_beta_0)
    
    beta_1.PLMIc.estimates <- sapply(1:MCMC_length, function(i) {
        return(result_data[["res"]][["EM_and_TSols_result"]][[i]][["PLMIc_result"]][["par"]][2])
    })
    beta_1.PLMIc.var <- sapply(1:MCMC_length, function(i) {
      result_data[["res"]][["EM_and_TSols_result"]][[i]][["PLMIc_result"]][[
        "laplace_variance"
      ]][2, 2]
    })
    beta_1_rubin_interv <- Rubin_interv(
      beta_1.PLMIc.estimates,
      beta_1.PLMIc.var,
      alpha = alpha.conf
    )
    beta_1.PLMIc.lb <- beta_1_rubin_interv[1]
    beta_1.PLMIc.ub <- beta_1_rubin_interv[2]
    beta_1.PLMIc.pt <- mean(beta_1.PLMIc.estimates)
    beta_1.PLMIc.length <- beta_1.PLMIc.ub - beta_1.PLMIc.lb
    beta_1.PLMIc.cover <- ((beta_1.PLMIc.ub >= true_beta_1) &
                          (beta_1.PLMIc.lb <= true_beta_1))
    beta_1.PLMIc.SE <- (beta_1.PLMIc.pt - true_beta_1)^2
    beta_1.PLMIc.bias <- abs(beta_1.PLMIc.pt - true_beta_1)
    
    beta_2.PLMIc.estimates <- sapply(1:MCMC_length, function(i) {
      return(result_data[["res"]][["EM_and_TSols_result"]][[i]][["PLMIc_result"]][["par"]][3])
    })
    beta_2.PLMIc.var <- sapply(1:MCMC_length, function(i) {
      result_data[["res"]][["EM_and_TSols_result"]][[i]][["PLMIc_result"]][[
        "laplace_variance"
      ]][3, 3]
    })
    beta_2_rubin_interv <- Rubin_interv(
      beta_2.PLMIc.estimates,
      beta_2.PLMIc.var,
      alpha = alpha.conf
    )
    beta_2.PLMIc.lb <- beta_2_rubin_interv[1]
    beta_2.PLMIc.ub <- beta_2_rubin_interv[2]
    beta_2.PLMIc.pt <- mean(beta_2.PLMIc.estimates)
    beta_2.PLMIc.length <- beta_2.PLMIc.ub - beta_2.PLMIc.lb
    beta_2.PLMIc.cover <- ((beta_2.PLMIc.ub >= true_beta_2) &
                             (beta_2.PLMIc.lb <= true_beta_2))
    beta_2.PLMIc.SE <- (beta_2.PLMIc.pt - true_beta_2)^2
    beta_2.PLMIc.bias <- abs(beta_2.PLMIc.pt - true_beta_2)
    
    # PLMI results
    beta_0.PLMI.estimates <- sapply(1:MCMC_length, function(i) {
      result_data[["res"]][["EM_and_TSols_result"]][[i]][["PLMI_result"]][[
        "param"
      ]][1]
    })
    beta_0.PLMI.var <- sapply(1:MCMC_length, function(i) {
      result_data[["res"]][["EM_and_TSols_result"]][[i]][["PLMI_result"]][[
        "laplace_variance"
      ]][1, 1]
    })
    beta_0_rubin_interv <- Rubin_interv(
      beta_0.PLMI.estimates,
      beta_0.PLMI.var,
      alpha = alpha.conf
    )
    beta_0.PLMI.lb <- beta_0_rubin_interv[1]
    beta_0.PLMI.ub <- beta_0_rubin_interv[2]
    beta_0.PLMI.pt <- mean(beta_0.PLMI.estimates)
    beta_0.PLMI.length <- beta_0.PLMI.ub - beta_0.PLMI.lb
    beta_0.PLMI.cover <- ((beta_0.PLMI.ub >= true_beta_0) &
                             (beta_0.PLMI.lb <= true_beta_0))
    beta_0.PLMI.SE <- (beta_0.PLMI.pt - true_beta_0)^2
    beta_0.PLMI.bias <- abs(beta_0.PLMI.pt - true_beta_0)
    
    beta_1.PLMI.estimates <- sapply(1:MCMC_length, function(i) {
      return(result_data[["res"]][["EM_and_TSols_result"]][[i]][["PLMI_result"]][["param"]][2])
    })
    beta_1.PLMI.var <- sapply(1:MCMC_length, function(i) {
      result_data[["res"]][["EM_and_TSols_result"]][[i]][["PLMI_result"]][[
        "laplace_variance"
      ]][2, 2]
    })
    beta_1_rubin_interv <- Rubin_interv(
      beta_1.PLMI.estimates,
      beta_1.PLMI.var,
      alpha = alpha.conf
    )
    beta_1.PLMI.lb <- beta_1_rubin_interv[1]
    beta_1.PLMI.ub <- beta_1_rubin_interv[2]
    beta_1.PLMI.pt <- mean(beta_1.PLMI.estimates)
    beta_1.PLMI.length <- beta_1.PLMI.ub - beta_1.PLMI.lb
    beta_1.PLMI.cover <- ((beta_1.PLMI.ub >= true_beta_1) &
                             (beta_1.PLMI.lb <= true_beta_1))
    beta_1.PLMI.SE <- (beta_1.PLMI.pt - true_beta_1)^2
    beta_1.PLMI.bias <- abs(beta_1.PLMI.pt - true_beta_1)
    
    beta_2.PLMI.estimates <- sapply(1:MCMC_length, function(i) {
      return(result_data[["res"]][["EM_and_TSols_result"]][[i]][["PLMI_result"]][["param"]][3])
    })
    beta_2.PLMI.var <- sapply(1:MCMC_length, function(i) {
      result_data[["res"]][["EM_and_TSols_result"]][[i]][["PLMI_result"]][[
        "laplace_variance"
      ]][3, 3]
    })
    beta_2_rubin_interv <- Rubin_interv(
      beta_2.PLMI.estimates,
      beta_2.PLMI.var,
      alpha = alpha.conf
    )
    beta_2.PLMI.lb <- beta_2_rubin_interv[1]
    beta_2.PLMI.ub <- beta_2_rubin_interv[2]
    beta_2.PLMI.pt <- mean(beta_2.PLMI.estimates)
    beta_2.PLMI.length <- beta_2.PLMI.ub - beta_2.PLMI.lb
    beta_2.PLMI.cover <- ((beta_2.PLMI.ub >= true_beta_2) &
                            (beta_2.PLMI.lb <= true_beta_2))
    beta_2.PLMI.SE <- (beta_2.PLMI.pt - true_beta_2)^2
    beta_2.PLMI.bias <- abs(beta_2.PLMI.pt - true_beta_2)
    
    
    return(data.frame(
      dataset,
      n_error,
      overlap_proportion,
      dataset_var,
      seed_proportion,
      x_2.prob,
      true_beta_0,
      true_beta_1,
      true_beta_2,
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
      beta_2.TS.OLS.pt,
      beta_2.TS.OLS.length,
      beta_2.TS.OLS.lb,
      beta_2.TS.OLS.ub,
      beta_2.TS.OLS.cover,
      beta_2.TS.OLS.SE,
      beta_2.TS.OLS.bias,
      beta_0.PLMIc.pt,
      beta_0.PLMIc.length,
      beta_0.PLMIc.lb,
      beta_0.PLMIc.ub,
      beta_0.PLMIc.cover,
      beta_0.PLMIc.SE,
      beta_0.PLMIc.bias,
      beta_1.PLMIc.pt,
      beta_1.PLMIc.length,
      beta_1.PLMIc.lb,
      beta_1.PLMIc.ub,
      beta_1.PLMIc.cover,
      beta_1.PLMIc.SE,
      beta_1.PLMIc.bias,
      beta_2.PLMIc.pt,
      beta_2.PLMIc.length,
      beta_2.PLMIc.lb,
      beta_2.PLMIc.ub,
      beta_2.PLMIc.cover,
      beta_2.PLMIc.SE,
      beta_2.PLMIc.bias,
      beta_0.PLMI.pt, # PLMI
      beta_0.PLMI.length,
      beta_0.PLMI.lb,
      beta_0.PLMI.ub,
      beta_0.PLMI.cover,
      beta_0.PLMI.SE,
      beta_0.PLMI.bias,
      beta_1.PLMI.pt,
      beta_1.PLMI.length,
      beta_1.PLMI.lb,
      beta_1.PLMI.ub,
      beta_1.PLMI.cover,
      beta_1.PLMI.SE,
      beta_1.PLMI.bias,
      beta_2.PLMI.pt,
      beta_2.PLMI.length,
      beta_2.PLMI.lb,
      beta_2.PLMI.ub,
      beta_2.PLMI.cover,
      beta_2.PLMI.SE,
      beta_2.PLMI.bias,
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
      beta_1.perfect.bias,
      beta_2.perfect.pt,
      beta_2.perfect.length,
      beta_2.perfect.lb,
      beta_2.perfect.ub,
      beta_2.perfect.cover,
      beta_2.perfect.SE,
      beta_2.perfect.bias
    ))
  },
  .progress = TRUE
)


result_beta1 <- inference_df %>%
  group_by(n_error, overlap_proportion, seed_proportion, x_2.prob, true_beta_2) %>%
  summarise(
    mean_PLMIc_cover = mean(beta_1.PLMIc.cover, na.rm = TRUE),
    mean_PLMI_cover = mean(beta_1.PLMI.cover, na.rm = TRUE),
    mean_TS.OLS_cover = mean(beta_1.TS.OLS.cover, na.rm = TRUE),
    mean_perfect_cover = mean(beta_1.perfect.cover, na.rm = TRUE),
    mean_PLMIc_length = median(beta_1.PLMIc.length, na.rm = TRUE),
    mean_PLMI_length = median(beta_1.PLMI.length, na.rm = TRUE),
    mean_TS.OLS_length = median(beta_1.TS.OLS.length, na.rm = TRUE),
    mean_perfect_length = median(beta_1.perfect.length, na.rm = TRUE),
    mean_PLMIc_MSE = mean(beta_1.PLMIc.SE, na.rm = TRUE),
    mean_PLMI_MSE = mean(beta_1.PLMI.SE, na.rm = TRUE),
    mean_TS.OLS_MSE = mean(beta_1.TS.OLS.SE, na.rm = TRUE),
    mean_perfect_MSE = mean(beta_1.perfect.SE, na.rm = TRUE),
    mean_PLMIc_bias = median(beta_1.PLMIc.bias, na.rm = TRUE),
    mean_PLMI_bias = median(beta_1.PLMI.bias, na.rm = TRUE),
    mean_TS.OLS_bias = median(beta_1.TS.OLS.bias, na.rm = TRUE),
    mean_perfect_bias = median(beta_1.perfect.bias, na.rm = TRUE),
    .groups = "keep"
  )

result_beta2 <- inference_df %>%
  group_by(n_error, overlap_proportion, seed_proportion, x_2.prob, true_beta_2) %>%
  summarise(
    mean_PLMIc_cover = mean(beta_2.PLMIc.cover, na.rm = TRUE),
    mean_PLMI_cover = mean(beta_2.PLMI.cover, na.rm = TRUE),
    mean_TS.OLS_cover = mean(beta_2.TS.OLS.cover, na.rm = TRUE),
    mean_perfect_cover = mean(beta_2.perfect.cover, na.rm = TRUE),
    mean_PLMIc_length = median(beta_2.PLMIc.length, na.rm = TRUE),
    mean_PLMI_length = median(beta_2.PLMI.length, na.rm = TRUE),
    mean_TS.OLS_length = median(beta_2.TS.OLS.length, na.rm = TRUE),
    mean_perfect_length = median(beta_2.perfect.length, na.rm = TRUE),
    mean_PLMIc_MSE = mean(beta_2.PLMIc.SE, na.rm = TRUE),
    mean_PLMI_MSE = mean(beta_2.PLMI.SE, na.rm = TRUE),
    mean_TS.OLS_MSE = mean(beta_2.TS.OLS.SE, na.rm = TRUE),
    mean_perfect_MSE = mean(beta_2.perfect.SE, na.rm = TRUE),
    mean_PLMIc_bias = median(beta_2.PLMIc.bias, na.rm = TRUE),
    mean_PLMI_bias = median(beta_2.PLMI.bias, na.rm = TRUE),
    mean_TS.OLS_bias = median(beta_2.TS.OLS.bias, na.rm = TRUE),
    mean_perfect_bias = median(beta_2.perfect.bias, na.rm = TRUE),
    .groups = "keep"
  )

result_beta0 <- inference_df %>%
  group_by(n_error, overlap_proportion, seed_proportion, x_2.prob, true_beta_2) %>%
  summarise(
    mean_PLMIc_cover = mean(beta_0.PLMIc.cover, na.rm = TRUE),
    mean_PLMI_cover = mean(beta_0.PLMI.cover, na.rm = TRUE),
    mean_TS.OLS_cover = mean(beta_0.TS.OLS.cover, na.rm = TRUE),
    mean_perfect_cover = mean(beta_0.perfect.cover, na.rm = TRUE),
    mean_PLMIc_length = mean(beta_0.PLMIc.length, na.rm = TRUE),
    mean_PLMI_length = mean(beta_0.PLMI.length, na.rm = TRUE),
    mean_TS.OLS_length = mean(beta_0.TS.OLS.length, na.rm = TRUE),
    mean_perfect_length = mean(beta_0.perfect.length, na.rm = TRUE),
    mean_PLMIc_MSE = mean(beta_0.PLMIc.SE, na.rm = TRUE),
    mean_PLMI_MSE = mean(beta_0.PLMI.SE, na.rm = TRUE),
    mean_TS.OLS_MSE = mean(beta_0.TS.OLS.SE, na.rm = TRUE),
    mean_perfect_MSE = mean(beta_0.perfect.SE, na.rm = TRUE),
    mean_PLMIc_bias = mean(beta_0.PLMIc.bias, na.rm = TRUE),
    mean_PLMI_bias = mean(beta_0.PLMI.bias, na.rm = TRUE),
    mean_TS.OLS_bias = mean(beta_0.TS.OLS.bias, na.rm = TRUE),
    mean_perfect_bias = mean(beta_0.perfect.bias, na.rm = TRUE),
    .groups = "keep"
  )