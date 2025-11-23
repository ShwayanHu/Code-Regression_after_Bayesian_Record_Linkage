packages <- c(
  'RecordLinkage',
  'tidyverse',
  'purrr',
  'gtools',
  'mvtnorm',
  'pracma',
  'coda',
  'testit',
  'VaRES',
  'parallel',
  'Rcpp',
  'BRL',
  'psych',
  'tictoc'
)
suppressPackageStartupMessages(lapply(packages, library, character.only = TRUE))
rm(packages)

HOME_DIR <- "~/BayesianRecordLinkage/Section6_Genuine_Data_Application/"
PROJECT_PATH <- "~/BayesianRecordLinkage/"
select <- dplyr::select
purrr::walk(
  list.files(
    paste0(PROJECT_PATH, "helper_fns/helper_fns_R"),
    pattern = "*.R$",
    full.names = TRUE
  ),
  source,
  .GlobalEnv
)

set.seed(20011231)


# Control Variables ------------- ---------------------------------------
Iter_start <- 501
nIter <- 1000
num_seeds <- 50

load(paste0(HOME_DIR, 'preprocessed_data_withName.RData'))
file_1 <- data$file_1; file_2 <- data$file_2
file_1.id <- data$file_1.id; file_2.id <- data$file_2.id
overlap_ids <- data$overlap_ids
known_links <- data$known_links
Z_known <- data$Z_known

# Perfect Regression ------------------------------------------------------
linked_data.perfect <- merge(
  file_1 %>% mutate(NQUEST = file_1.id),
  file_2 %>% mutate(NQUEST = file_2.id),
  by = 'NQUEST',
  all = FALSE
)

fit <- lm(Y ~ X1 + X2 + X3 + X4, linked_data.perfect)
print(summary(fit))
pairs(linked_data.perfect |> select(Y, X1, X2, X3, X4))

n1 <- nrow(file_1); n2 <- nrow(file_2)
cellinds <- expand.grid(1:n1, 1:n2)
records <- rbind(file_1 |> select(starts_with("f")), file_2 |> select(starts_with("f")))

gamma_in_c <- compareRecords(
  file_1 |> select(starts_with("f")), 
  file_2 |> select(starts_with("f")), 
  flds = paste0("f", 1:7),
  types = c('bi', 'bi', 'bi', 'bi','bi', 'lv', 'lv'),       
  breaks = c(0, .001, .25, .5)
)

gamma <- data.frame(
  record_1 = rep(1:n1, n2),
  record_2 = rep(1:n2, each = n1),
  f1 = as.numeric(!gamma_in_c$comparisons[,1]),
  f2 = as.numeric(!gamma_in_c$comparisons[,3]),
  f3 = as.numeric(!gamma_in_c$comparisons[,5]),
  f4 = as.numeric(!gamma_in_c$comparisons[,7]),
  f5 = as.numeric(!gamma_in_c$comparisons[,9]),
  f6 = sapply(1:nrow(gamma_in_c$comparisons), function(i){which(gamma_in_c$comparisons[i,11:15])-1}),
  f7 = sapply(1:nrow(gamma_in_c$comparisons), function(i){which(gamma_in_c$comparisons[i,16:20])-1})
)

# known links - CORRECTED VERSION
direct_matches <- intersect(file_1.id, file_2.id)

Z_true <- rep(n1 + (1:n2), n2)  # 默认为 non-link
for (f2_pos in 1:n2) {
  id_in_file2 <- file_2.id[f2_pos]
  matching_positions <- which(file_1.id == id_in_file2)
  
  if (length(matching_positions) > 0) {
    Z_true[f2_pos] <- matching_positions[1]
  }
}

tic("MCMC Time")
mcmc <- bipartiteGibbs_withSeed(gamma_in_c, Z_known, nIter=nIter, a=1, b=1, aBM=1, bBM=1, seed=20011231)
toc()

# organize result
Z = t(mcmc$Z)
# plot(1:nrow(Z), apply(Z, 1, function(z){sum(z <= n1)}), 'l', main = 'Number of links claimed')
# plot(1:nrow(Z), apply(Z, 1, function(z){sum(z <= n1 & z==Z_true)}), 'l', col='blue')

m <- lapply(
  1:nIter, 
  function(s) {
    res <- list(
      f1 = mcmc$m[1:2,s], 
      f2 = mcmc$m[3:4,s],
      f3 = mcmc$m[5:6,s],
      f4 = mcmc$m[7:8,s],
      f5 = mcmc$m[9:10,s],
      f6 = mcmc$m[11:15,s], 
      f7 = mcmc$m[16:20,s]
    )
  }
)
u <- lapply(
  1:nIter, 
  function(s) {
    res <- list(
      f1 = mcmc$u[1:2,s],
      f2 = mcmc$u[3:4,s],
      f3 = mcmc$u[5:6,s],
      f4 = mcmc$u[7:8,s],
      f5 = mcmc$u[9:10,s],
      f6 = mcmc$u[11:15,s], 
      f7 = mcmc$u[16:20,s]
    )
  }
)

# Posterior Inference ------------------------------------------------------
if ("Y" %in% colnames(file_1)) {
  p.y <- approxfun(density(file_1$Y))
} else {
  p.y <- approxfun(density(file_2$Y))
}


EM_and_TSols_history <- lapply(
  Iter_start:nIter,
  function(s) {
    tic(msg = paste0("Iteration ", s, " Time: "))
    Z_tmp <- Z[s, ]
    m_tmp <- m[[s]]
    u_tmp <- u[[s]]
    linked_data_tmp <- z2regdata_with_goodness_and_seed_flelxible(Z_tmp, file_1, file_2, m_tmp, u_tmp, gamma, Z_known)
    PLMIc_result <- EM_estimation_Seed_ConfidenceMeasure_4covariates(
      X1 = as.numeric(linked_data_tmp$X$X1),
      X2 = as.numeric(linked_data_tmp$X$X2),
      X3 = as.numeric(linked_data_tmp$X$X3),
      X4 = as.numeric(linked_data_tmp$X$X4),
      Y = as.numeric(linked_data_tmp$Y$Y),
      G = linked_data_tmp$G,
      p.y = p.y,
      Seed.vec = as.numeric(linked_data_tmp$Seed.vec)
    )
    # PLMI
    PLMI_result <- EM_estimation_withSeed_flexible(
      X = linked_data_tmp$X %>% mutate(
        X2 = as.numeric(X2),
        X3 = as.numeric(X3)
      ),
      Y = as.numeric(linked_data_tmp$Y$Y),
      Seed.vec = as.numeric(linked_data_tmp$Seed.vec),
      p.y = p.y
    )
    
    # run loop to calculate and save TS.OLS result
    X_mat_tmp <- cbind(
      rep(1, length(linked_data_tmp$Y$Y)),
      linked_data_tmp$X$X1,
      linked_data_tmp$X$X2,
      linked_data_tmp$X$X3,
      linked_data_tmp$X$X4
    )
    Y_mat_tmp <- matrix(
      as.numeric(linked_data_tmp$Y$Y),
      ncol = 1
    )
    TS.OLS_coefficients <- solve(t(X_mat_tmp) %*% X_mat_tmp) %*%
      t(X_mat_tmp) %*%
      Y_mat_tmp
    TS.OLS_cov <- sum((Y_mat_tmp - X_mat_tmp %*% TS.OLS_coefficients)^2) /
      (nrow(X_mat_tmp) - ncol(X_mat_tmp)) *
      solve(t(X_mat_tmp) %*% X_mat_tmp)
    res <- list(
      TS.OLS_coefficients = TS.OLS_coefficients,
      TS.OLS_cov = TS.OLS_cov
    )
    
    results <- list(
      PLMIc_result = PLMIc_result,
      PLMI_result = PLMI_result,
      TS.OLS_coefficients = TS.OLS_coefficients,
      TS.OLS_cov = TS.OLS_cov
    )
    toc()
    return(results)
  }
)

# Perfect Regression ------------------------------------------------------
linked_data.perfect <- merge(
  file_1 %>% mutate(NQUEST = file_1.id),
  file_2 %>% mutate(NQUEST = file_2.id),
  by = 'NQUEST',
  all = FALSE
)

fit <- lm(Y ~ X1 + X2 + X3 + X4, linked_data.perfect)
print(summary(fit))

Y_tmp.perfect <- matrix(linked_data.perfect$Y, ncol = 1)
X_tmp.perfect <- cbind(
  rep(1, nrow(linked_data.perfect)),
  linked_data.perfect$X1,
  linked_data.perfect$X2,
  linked_data.perfect$X3,
  linked_data.perfect$X4
)
Perfect.OLS_coefficients <- solve(t(X_tmp.perfect) %*% X_tmp.perfect) %*%
  t(X_tmp.perfect) %*%
  Y_tmp.perfect
Perfect.OLS_cov <- sum(
  (Y_tmp.perfect - X_tmp.perfect %*% Perfect.OLS_coefficients)^2
) /
  (nrow(X_tmp.perfect) - ncol(X_tmp.perfect)) *
  solve(t(X_tmp.perfect) %*% X_tmp.perfect)
Perfect.OLS_result <- list(
  Perfect.OLS_coefficients = Perfect.OLS_coefficients,
  Perfect.OLS_cov = Perfect.OLS_cov
)


# Save result -------------------------------------------------------------
res <- list(
  EM_and_TSols_result = EM_and_TSols_history,
  Perfect.OLS_result = Perfect.OLS_result,
  Z = Z,
  m = m,
  u = u,
  file_1 = file_1,
  file_2 = file_2,
  file_1.id = file_1.id,
  file_2.id = file_2.id,
  gamma = gamma,
  linked_data.perfect = linked_data.perfect
)

save(res, file = paste0(HOME_DIR, "res_RealData_withName.RData"))
