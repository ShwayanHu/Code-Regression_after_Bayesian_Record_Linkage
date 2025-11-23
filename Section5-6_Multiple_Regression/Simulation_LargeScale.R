# Packages, seeds, and Helper Functions -----------------------------------
PROJECT_PATH <- "/BayesianRecordLinkage/"

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
  'BRL'
)
suppressPackageStartupMessages(lapply(packages, library, character.only = TRUE))
rm(packages)

purrr::walk(
  list.files(
    paste0(PROJECT_PATH, "helper_fns/helper_fns_R"),
    pattern = "*.R$",
    full.names = TRUE
  ),
  source,
  .GlobalEnv
)

# Loops -------------------------------------------------------------------
# looping variables
looping_grid <- expand.grid(
  0:99, # dataset
  c(1), # n_error
  c(0.5), # overlap_proportion
  c(0.05), # seed_proportion
  c(0.5), # x_2.prob
  c(10, 1) # beta_2
  )

lapply(
  1:nrow(looping_grid),
  # 1,
  function(loop_var) {
    # control random seed explicitly
    set.seed(1231 + loop_var)
    
    # loop control variable
    dataset <- looping_grid[loop_var, 1]
    n_error <- looping_grid[loop_var, 2]
    overlap_proportion <- looping_grid[loop_var, 3]
    seed_proportion <- looping_grid[loop_var, 4]
    x_2.prob <- looping_grid[loop_var, 5]
    beta_2 <- looping_grid[loop_var, 6]
    
    # load data
    records <- read.csv(
      paste(
        PROJECT_PATH,
        "Simulation Datafiles/",
        "percdups50_nerrors",
        n_error,
        "_dataset",
        dataset,
        "_4fields_twofiles_BRLpaper.csv",
        sep = ""
      )
    )
    records$file <- rep(2:1, length.out = dim(records)[1])
    records <- records[records$file %in% c(1, 2), ] # pick original record and first dup
    
    # Preprocess data (create 5 sub-datasets) -----------------------------
    preprocessed_data <- PreProcess_SimulationData_withSeed_MLRandNonNormal(
      records = records,
      num_simu = 5,
      n1 = 500,
      n2 = 500,
      overlap_proportion = overlap_proportion,
      seed_proportion = seed_proportion,
      x_2.prob = x_2.prob,
      beta_0=3,
      beta_1 = 0.5,
      beta_2 = beta_2,
      sigma = 1 # standard deviation of regression error term
    )
    
    # Simulation for each sub-dataset in `preprocessed_data`---------------
    for (dataset_var in 1:5) {
      # Result File Name --------------------------------------------------------
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
      if (
        result_file_name %in%
        list.files(paste0(PROJECT_PATH, "Section5-6_Multiple_Regression/Simu_Result"))
      ) {
        next
      }
      
      file_1 <- preprocessed_data[[dataset_var]]$file_1
      file_2 <- preprocessed_data[[dataset_var]]$file_2
      gamma <- preprocessed_data[[dataset_var]]$gamma
      Z_known <- preprocessed_data[[dataset_var]]$Z_known
      
      # Run MCMC with Sadinle's code (in C language) ------------------------
      gamma_in_c <- compareRecords(file_1, file_2, flds=c("gname", "fname", "age", "occup"), types=c("lv","lv","bi","bi"), breaks=c(0,.25,.5))
      mcmc <- bipartiteGibbs_withSeed(gamma_in_c, Z_known, nIter=1000, a=1, b=1, aBM=1, bBM=1, seed=loop_var)
      Z = t(mcmc$Z)
      m <- lapply(1:1000, function(s) {res <- list(f1 = mcmc$m[1:4,s], f2 = mcmc$m[5:8,s], f3 = mcmc$m[9:10,s], f4 = mcmc$m[11:12,s])})
      u <- lapply(1:1000, function(s) {res <- list(f1 = mcmc$u[1:4,s], f2 = mcmc$u[5:8,s], f3 = mcmc$u[9:10,s], f4 = mcmc$u[11:12,s])})
      
      # Fit EM Model and TS.OLS ------------------------------------------------------------
      # estimate p.y
      p.y <- approxfun(density(file_2$y))
      
      EM_and_TSols_history <- lapply(
        101:1000, # first 100 iters are burn-in
        function(s) {
          print(s)
          # run loop to calculate and save EM result
          Z_tmp <- Z[s, ]
          m_tmp <- m[[s]]
          u_tmp <- u[[s]]
          linked_data_tmp <- z2regdata_with_goodness_and_seed_MLRandNonNormal(Z_tmp, file_1, file_2, m_tmp, u_tmp, gamma, Z_known)
          # PLMIc
          PLMIc_result <- EM_estimation_Seed_ConfidenceMeasure_MLRandNonNormal(
            X1 = linked_data_tmp$X1,
            X2 = linked_data_tmp$X2,
            Y = as.numeric(linked_data_tmp$Y),
            G = linked_data_tmp$G,
            p.y = p.y,
            Seed.vec = as.numeric(linked_data_tmp$Seed.vec)
          )
          # PLMI
          PLMI_result <- EM_estimation_withSeed_MLRandNonNormal(
            X1.vec = linked_data_tmp$X1,
            X2.vec = linked_data_tmp$X2,
            Y.vec = as.numeric(linked_data_tmp$Y),
            Seed.vec = as.numeric(linked_data_tmp$Seed.vec),
            p.y = p.y
          )
          
          # run loop to calculate and save TS.OLS result
          X_mat_tmp <- cbind(
            rep(1, length(linked_data_tmp$Y)),
            linked_data_tmp$X1,
            linked_data_tmp$X2
          )
          Y_mat_tmp <- matrix(
            as.numeric(linked_data_tmp$Y),
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
          return(results)
        }
      )
      
      # Perfect Regression ------------------------------------------------------
      linked_data.perfect <- merge(
        file_1,
        file_2,
        by = 'rec.id',
        all = FALSE
      ) %>%
        select(x_1, x_2, y)
      Y_tmp.perfect <- matrix(linked_data.perfect$y, ncol = 1)
      X_tmp.perfect <- cbind(
        rep(1, nrow(linked_data.perfect)),
        linked_data.perfect$x_1,
        linked_data.perfect$x_2
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
      
      res <- list(
        EM_and_TSols_result = EM_and_TSols_history,
        Perfect.OLS_result = Perfect.OLS_result
      )
      
      save(
        res,
        file = paste(
          PROJECT_PATH,
          "Section5-6_Multiple_Regression/Simu_Result/",
          result_file_name,
          sep = ''
        )
      )
    } # end of `dataset_var` (sub-dataset)
  } # end of `loop_var` (original dataset)
  # mc.cores = 36
)
