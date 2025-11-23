library(purrr)
library(tidyverse)
library(testit)
library(ggplot2)
source(
  "~/BayesianRecordLinkage/helper_fns/helper_fns_R/Rubin_interv.R"
)
source(
  "~/BayesianRecordLinkage/helper_fns/helper_fns_R/getHpi.R"
)

MCMC_length <- 500
alpha.conf <- 0.1

# Load result data
load(
  paste0(
    "~/BayesianRecordLinkage/",
    "Section6_Genuine_Data_Application/res_RealData_withName.RData"
  ),
  temp_env <- new.env()
)
result_data <- as.list(temp_env)

load(paste0(
  "~/BayesianRecordLinkage/",
  "Section6_Genuine_Data_Application/preprocessed_data_withName.RData"
))

linked_data.perfect <- merge(
  data$file_1 |> mutate(id = data$file_1.id),
  data$file_2 |> mutate(id = data$file_2.id),
  by = 'id'
)

# Number of coefficients (intercept + covariates)
n_coef <- length(result_data$res$Perfect.OLS_result$Perfect.OLS_coefficients)

# Function to get TS.OLS results
get_TSOLS <- function(j) {
  est <- sapply(1:MCMC_length, function(i) {
    result_data$res$EM_and_TSols_result[[i]]$TS.OLS_coefficients[j]
  })
  var <- sapply(1:MCMC_length, function(i) {
    result_data$res$EM_and_TSols_result[[i]]$TS.OLS_cov[j, j]
  })
  ci <- Rubin_interv(est, var, alpha = alpha.conf)
  tibble(
    Method = "TS.OLS",
    Parameter = paste0("beta_", j - 1),
    Point_Estimate = median(est),
    Lower_Bound = ci[1],
    Upper_Bound = ci[2]
  )
}

# Function to get Perfect OLS results
get_Perfect <- function(j) {
  coeffs <- result_data$res$Perfect.OLS_result$Perfect.OLS_coefficients
  covmat <- result_data$res$Perfect.OLS_result$Perfect.OLS_cov
  pt <- coeffs[j]
  se <- sqrt(covmat[j, j])
  lb <- pt - qt(1 - alpha.conf/2, df = nrow(linked_data.perfect) - n_coef) * se
  ub <- pt + qt(1 - alpha.conf/2, df = nrow(linked_data.perfect) - n_coef) * se
  tibble(
    Method = "Perfect",
    Parameter = paste0("beta_", j - 1),
    Point_Estimate = pt,
    Lower_Bound = lb,
    Upper_Bound = ub
  )
}

# Generic function for PLMI and PLMIc
get_PLMI_like <- function(j, type = c("PLMIc", "PLMI")) {
  type <- match.arg(type)
  # if (type=='PLMI' & j==2) {browser()}
  est <- sapply(1:MCMC_length, function(i) {
    if (type == "PLMIc") {
      result_data$res$EM_and_TSols_result[[i]][[paste0(type, "_result")]]$par[j]
    } else {
      result_data$res$EM_and_TSols_result[[i]][[paste0(type, "_result")]]$param[j]
    }
  })
  iter <- sapply(1:MCMC_length, function(i) {
    result_data$res$EM_and_TSols_result[[i]][[paste0(type, "_result")]]$iter
  })
  var <- sapply(1:MCMC_length, function(i) {
    result_data$res$EM_and_TSols_result[[i]][[paste0(type, "_result")]]$laplace_variance[j, j]
  })
  est <- est[iter < 2000]
  var <- var[iter < 2000]
  ci <- Rubin_interv(est, var, alpha = alpha.conf)
  tibble(
    Method = type,
    Parameter = paste0("beta_", j - 1),
    Point_Estimate = median(est),
    Lower_Bound = ci[1],
    Upper_Bound = ci[2]
  )
}

# Collect results for all coefficients & methods
inference_results <- map_dfr(1:n_coef, function(j) {
  bind_rows(
    get_TSOLS(j),
    get_Perfect(j),
    get_PLMI_like(j, "PLMIc"),
    get_PLMI_like(j, "PLMI")
  )
})

# Add CI width
inference_results <- inference_results %>%
  mutate(CI_Width = Upper_Bound - Lower_Bound)

# print(inference_results)

# Wide format for easier comparison
inference_wide <- inference_results %>%
  pivot_wider(
    id_cols = Parameter,
    names_from = Method,
    values_from = c(Point_Estimate, Lower_Bound, Upper_Bound, CI_Width),
    names_sep = "_"
  )

# Plot linkage accuracy ---------------------------------------------------
Z <- result_data$res$Z
n1 <- nrow(result_data$res$file_1); n2 <- nrow(result_data$res$file_2)

Z_true <- rep(n1 + (1:n2), n2)  # 默认为 non-link
for (f2_pos in 1:n2) {
  id_in_file2 <- result_data$res$file_2.id[f2_pos]
  matching_positions <- which(result_data$res$file_1.id == id_in_file2)

  if (length(matching_positions) > 0) {
    Z_true[f2_pos] <- matching_positions[1]
  }
}

n_links_claimed <- apply(Z, 1, function(z) sum(z <= n1))
n_correct_links <- apply(Z, 1, function(z) sum(z <= n1 & z == Z_true))
true_n_links <- sum(Z_true <= n1)
cat(
  'Average linking accuracy:',
  mean(n_correct_links/n_links_claimed),
  '\n Average links claimed:',
  mean(n_links_claimed),
  '\n Average true links:',
  mean(n_correct_links)
)