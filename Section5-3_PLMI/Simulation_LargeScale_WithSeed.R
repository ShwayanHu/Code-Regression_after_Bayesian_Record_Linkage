options(warn = 2)
# Packages, seeds, and Helper Functions -----------------------------------
PROJECT_PATH <- "~/BayesianRecordLinkage/"

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

# looping variables
looping_grid <- expand.grid(0:99, c(1), c(0.5), c(0.3), c(0.01,0.05, 0.1))

lapply(
  1:nrow(looping_grid),
  function(loop_var) {
    print(loop_var)
    # control random seed explicitly
    set.seed(1231 + loop_var)
    
    # loop control variable
    dataset <- looping_grid[loop_var, 1]
    n_error <- looping_grid[loop_var, 2]
    overlap_proportion <- looping_grid[loop_var, 3]
    R_sq <- looping_grid[loop_var, 4]
    seed_proportion <- looping_grid[loop_var, 5]
    if (seed_proportion != 0.05) {return(NULL)}
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
    preprocessed_data <- PreProcess_SimulationData_withSeed(
      records = records,
      num_simu = 5,
      n1 = 500,
      n2 = 500,
      overlap_proportion = overlap_proportion,
      R_sq = R_sq,
      seed_proportion = seed_proportion,
      beta_0 = 3,
      beta_1= 3
    )
    
    # Simulation for each sub-dataset in `preprocessed_data`---------------
    for (dataset_var in 1:5) {
      file_1 <- preprocessed_data[[dataset_var]]$file_1
      file_2 <- preprocessed_data[[dataset_var]]$file_2
      gamma <- preprocessed_data[[dataset_var]]$gamma
      Z_known <- preprocessed_data[[dataset_var]]$Z_known
      
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
        "_percSeed",
        seed_proportion * 100,
        '_model',
        'TSols_Perfect',
        '.RData',
        sep = ""
      )
      
      # Run MCMC with Sadinle's code (in C language) ------------------------
      gamma_in_c <- compareRecords(file_1, file_2, flds=c("gname", "fname", "age", "occup"), types=c("lv","lv","bi","bi"), breaks=c(0,.25,.5))
      mcmc <- bipartiteGibbs_withSeed(gamma_in_c, Z_known, nIter=1000, a=1, b=1, aBM=1, bBM=1, seed=loop_var)
      Z = t(mcmc$Z)
      m <- lapply(1:1000, function(s) {res <- list(f1 = mcmc$m[1:4,s], f2 = mcmc$m[5:8,s], f3 = mcmc$m[9:10,s], f4 = mcmc$m[11:12,s])})
      u <- lapply(1:1000, function(s) {res <- list(f1 = mcmc$u[1:4,s], f2 = mcmc$u[5:8,s], f3 = mcmc$u[9:10,s], f4 = mcmc$u[11:12,s])})
      # Fit EM Model and TS.OLS ------------------------------------------------------------
      # estimate p.y
      mu.y <- mean(file_2$y)
      sigma.y <- sqrt(sum((file_2$y-mu.y)^2)/nrow(file_2))
      
      EM_and_TSols_history <- lapply(
        101:1000, # first 100 iters are burn-in
        function(s){
          # run loop to calculate and save EM result
          Z_tmp <- Z[s,]
          m_tmp <- m[[s]]
          u_tmp <- u[[s]]
          linked_data_tmp <- z2regdata_withSeed(Z_tmp, file_1, file_2, Z_known)
          X_tmp <- linked_data_tmp$X
          Y_tmp <- linked_data_tmp$Y
          Seed.vec_tmp <- linked_data_tmp$Seed.vec
          EM_result <- EM_estimation_withSeed(
            X = as.numeric(X_tmp[,2]),
            Y = as.numeric(Y_tmp),
            Seed.vec = as.numeric(Seed.vec_tmp),
            p.y = function(y){dnorm(y, mean = mu.y, sd = sigma.y)}
          )
          # run loop to calculate and save TS.OLS result
          TS.OLS_coefficients <- solve(t(X_tmp) %*% X_tmp) %*% t(X_tmp) %*% Y_tmp
          TS.OLS_cov <- sum((Y_tmp - X_tmp %*% TS.OLS_coefficients)^2)/(nrow(X_tmp)-ncol(X_tmp)) * solve(t(X_tmp) %*% X_tmp)
          results <- list(
            EM_result = EM_result,
            TS.OLS_coefficients = TS.OLS_coefficients,
            TS.OLS_cov = TS.OLS_cov
          )
          return(results)
        }
      )

      # Perfect Regression ------------------------------------------------------
      linked_data.perfect <- merge(file_1, file_2, by = 'rec.id', all=FALSE) %>% select(x,y)
      Y_tmp.perfect <- matrix(linked_data.perfect$y, ncol=1)
      X_tmp.perfect <- matrix(
        c(rep(1, length(Y_tmp.perfect)), linked_data.perfect$x),
        ncol=2,
        byrow = FALSE
      )
      Perfect.OLS_coefficients <- solve(t(X_tmp.perfect) %*% X_tmp.perfect) %*% t(X_tmp.perfect) %*% Y_tmp.perfect
      Perfect.OLS_cov <- sum((Y_tmp.perfect - X_tmp.perfect %*% Perfect.OLS_coefficients)^2)/(nrow(X_tmp.perfect)-ncol(X_tmp.perfect)) * solve(t(X_tmp.perfect) %*% X_tmp.perfect)
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
        file=paste(
          PROJECT_PATH,
          "Section5-3_PLMI/Simu_Result_NOEM/",
          result_file_name,
          sep = ""
        )
      )
    } # end of `dataset_var` (sub-dataset)
  } # end of `loop_var` (original dataset)
)




















































