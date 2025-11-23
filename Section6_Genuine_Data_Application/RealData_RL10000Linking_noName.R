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
purrr::walk(
  list.files(
    paste0(PROJECT_PATH, "helper_fns/helper_fns_R"),
    pattern = "*.R$",
    full.names = TRUE
  ),
  source,
  .GlobalEnv
)

set.seed(5250)


# Control Variables ------------- ---------------------------------------
Iter_start <- 501
nIter <- 1000
num_seeds <- 50

load(paste0(HOME_DIR, "preprocessed_data_withName.RData"))
file_1 <- data$file_1 |> select(-f6, -f7)
file_2 <- data$file_2 |> select(-f6, -f7)
file_1.id <- data$file_1.id
file_2.id <- data$file_2.id
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

gamma <- data.frame(
  record_1 = cellinds$Var1,
  record_2 = cellinds$Var2,
  f1 = AgrLevBinComp(records$f1, cellinds),
  f2 = AgrLevBinComp(records$f2, cellinds),
  f3 = AgrLevBinComp(records$f3, cellinds),
  f4 = AgrLevBinComp(records$f4, cellinds),
  f5 = AgrLevBinComp(records$f5, cellinds)
)

gamma_in_c <- compareRecords(
  file_1 |> select(starts_with("f")), 
  file_2 |> select(starts_with("f")), 
  flds = paste0("f", 1:5),
  types = c('bi', 'bi', 'bi', 'bi','bi')
)

direct_matches <- intersect(file_1.id, file_2.id)
cat("Number of true links:", length(direct_matches), "\n")
cat("Number of records in file 1:", n1, "\n")
cat("Number of records in file 2:", n2, "\n")
cat("Number of seeds:", num_seeds, "\n")


Z_known <- rep(-1, n2)
for(seed_id in known_links) {
  f1_pos <- which(file_1.id == seed_id)
  f2_pos <- which(file_2.id == seed_id)
  Z_known[f2_pos] <- f1_pos
}

Z_true <- rep(n1 + (1:n2), n2)
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


Z = t(mcmc$Z)

m <- lapply(
  1:nIter, 
  function(s) {
    res <- list(
      f1 = mcmc$m[1:2,s], 
      f2 = mcmc$m[3:4,s],
      f3 = mcmc$m[5:6,s],
      f4 = mcmc$m[7:8,s],
      f5 = mcmc$m[9:10,s]
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
      f5 = mcmc$u[9:10,s]
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
      # X2 = as.numeric(linked_data_tmp$X$X2),
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

save(res, file = paste0(HOME_DIR, "res_RealData_noName.RData"))