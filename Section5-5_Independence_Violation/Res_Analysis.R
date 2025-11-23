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
  c(0.3), # R_sq,
  c(0.05), # seed_proportion
  c(1:5), # sub_id
  c(1, 3, 5) # mean_x_nonlink
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
    R_sq <- looping_grid[loop_var, 4]
    seed_proportion <- looping_grid[loop_var, 5]
    dataset_var <- looping_grid[loop_var, 6]
    mean_x_nonlink <- looping_grid[loop_var, 7]
    
    true_beta_0 <- 3
    true_beta_1 <- 3
    
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
      "_percRsq",
      R_sq * 100,
      '_meanXNonLink',
      mean_x_nonlink,
      '_model',
      'EM_TSols_Perfect',
      '.RData',
      sep = ""
    )
    
    load(
      paste0(
        "~/BayesianRecordLinkage/",
        "Section5-5_Independence_Violation/Simu_Result/",
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
    
    return(data.frame(
      dataset,
      n_error,
      overlap_proportion,
      R_sq,
      dataset_var,
      seed_proportion,
      mean_x_nonlink,
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
  group_by(n_error, overlap_proportion, R_sq, seed_proportion, mean_x_nonlink) %>%
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

result_beta0 <- inference_df %>%
  group_by(n_error, overlap_proportion, R_sq, seed_proportion, mean_x_nonlink) %>%
  summarise(
    mean_PLMIc_cover = mean(beta_0.PLMIc.cover, na.rm = TRUE),
    mean_PLMI_cover = mean(beta_0.PLMI.cover, na.rm = TRUE),
    mean_TS.OLS_cover = mean(beta_0.TS.OLS.cover, na.rm = TRUE),
    mean_perfect_cover = mean(beta_0.perfect.cover, na.rm = TRUE),
    mean_PLMIc_length = median(beta_0.PLMIc.length, na.rm = TRUE),
    mean_PLMI_length = median(beta_0.PLMI.length, na.rm = TRUE),
    mean_TS.OLS_length = median(beta_0.TS.OLS.length, na.rm = TRUE),
    mean_perfect_length = median(beta_0.perfect.length, na.rm = TRUE),
    mean_PLMIc_MSE = mean(beta_0.PLMIc.SE, na.rm = TRUE),
    mean_PLMI_MSE = mean(beta_0.PLMI.SE, na.rm = TRUE),
    mean_TS.OLS_MSE = mean(beta_0.TS.OLS.SE, na.rm = TRUE),
    mean_perfect_MSE = mean(beta_0.perfect.SE, na.rm = TRUE),
    mean_PLMIc_bias = median(beta_0.PLMIc.bias, na.rm = TRUE),
    mean_PLMI_bias = median(beta_0.PLMI.bias, na.rm = TRUE),
    mean_TS.OLS_bias = median(beta_0.TS.OLS.bias, na.rm = TRUE),
    mean_perfect_bias = median(beta_0.perfect.bias, na.rm = TRUE),
    .groups = "keep"
  )