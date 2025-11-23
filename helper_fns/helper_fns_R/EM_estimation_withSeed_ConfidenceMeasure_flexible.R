EM_estimation_Seed_ConfidenceMeasure_flexible <- function(X, Y, G, p.y, Seed.vec) {
  # Detect number of covariates
  n_covariates <- ncol(X) - 1  # Subtract 1 for intercept column
  n_beta <- ncol(X)
  n_params <- n_covariates + 1 + 1 + 2  # betas + sigma + 2 etas
  
  # Helper functions
  # Safe h with clipping
  h <- function(g, eta0, eta1) {
    p <- 1 / (1 + exp(-(eta0 + eta1 * g)))
    p <- pmin(pmax(p, 1e-8), 1 - 1e-8)
    return(p)
  }
  
  # Safe phi (density of N(X %*% beta, sigma^2))
  phi <- function(y, X_row, beta_vec, sigma) {
    mu <- as.numeric(X_row %*% beta_vec)
    dnorm(y, mean = mu, sd = sigma)
  }
  
  # Observed data log-likelihood for a single observation
  log_likelihood.1 <- function(param, p.y, X_row, y, g) {
    n_betas <- ncol(X)
    beta_vec <- param[1:n_betas]
    sigma <- param[n_betas + 1]
    eta_0 <- param[n_betas + 2]
    eta_1 <- param[n_betas + 3]
    
    mu <- as.numeric(X_row %*% beta_vec)
    term.1 <- (2*pi*sigma^2)^(-1/2) * exp(-1/(2*sigma^2) * (y-mu)^2) * h(g, eta_0, eta_1)
    term.2 <- p.y(y) * (1 - h(g, eta_0, eta_1))
    return(log(term.1 + term.2))
  }
  
  # E-step function
  prob_t_is_1 <- function(X, y, g, par, p.y, Seed.vec = NULL) {
    n_betas <- ncol(X)
    beta_vec <- par[1:n_betas]
    sigma <- max(par[n_betas + 1], 1e-4)
    eta_0 <- par[n_betas + 2]
    eta_1 <- par[n_betas + 3]
    
    h_val <- h(g, eta_0, eta_1)
    phi_val <- sapply(1:length(y), function(i) {
      pmax(phi(y[i], X[i, , drop=FALSE], beta_vec, sigma), 1e-8)
    })
    p_y_safe <- pmax(p.y(y), 1e-8)
    
    numer <- h_val * phi_val
    denom <- numer + (1 - h_val) * p_y_safe
    
    t_hat <- numer / denom
    if (!is.null(Seed.vec)) t_hat[Seed.vec==1] <- 1
    return(t_hat)
  }
  
  # M-step objective function
  Q <- function(X, Y, G, T_Estep, par) {
    n_betas <- ncol(X)
    beta_vec <- par[1:n_betas]
    sigma <- max(par[n_betas + 1], 1e-4)  # enforce positivity
    eta_0 <- par[n_betas + 2]
    eta_1 <- par[n_betas + 3]
    
    h_val <- h(G, eta_0, eta_1)
    phi_val <- sapply(1:length(Y), function(i) {
      pmax(phi(Y[i], X[i, , drop=FALSE], beta_vec, sigma), 1e-8)
    })
    
    loglik <- T_Estep * (log(h_val + 1e-8) + log(phi_val)) +
      (1 - T_Estep) * log(1 - h_val + 1e-8)
    return(-sum(loglik))
  }
  
  Q_regression <- function(par_beta_sigma, X, Y, T_Estep) {
    n_betas <- ncol(X)
    beta_vec <- par_beta_sigma[1:n_betas]
    sigma <- max(par_beta_sigma[n_betas + 1], 1e-4)
    mu <- as.numeric(X %*% beta_vec)
    -sum(T_Estep * dnorm(Y, mu, sigma, log=TRUE))  # regression part only
  }
  
  Q_logistic <- function(par_eta, G, T_Estep) {
    eta0 <- par_eta[1]
    eta1 <- par_eta[2]
    h_val <- h(G, eta0, eta1)
    -sum(T_Estep * log(h_val) + (1 - T_Estep) * log(1 - h_val))  # logistic part only
  }
  
  # model: Y = X %*% beta + epsilon. epsilon ~ N(0,sigma^2)
  # Pr(t = 1|g) = 1/(1 + exp(-(eta_0 + eta_1 * g)))
  
  # Initialization
  fit_lm <- lm.fit(x = as.matrix(X[Seed.vec==1, , drop=FALSE]), y = Y[Seed.vec==1])
  
  par_iter <- c(fit_lm$coefficients, 
                sqrt(sum(fit_lm$residuals^2)/fit_lm$df.residual),
                0, 0)
  par_prev <- rep(0, n_params)
  iter <- 1
  
  # EM iterations
  while (sqrt(sum((par_iter - par_prev)^2)) > 1e-7 & iter <= 2e3) {
    # E-step
    t <- prob_t_is_1(X, Y, G, par_iter, p.y); t[Seed.vec==1] <- 1 # Force the latent class probability of seed to be 1
    
    # M-step
    n_betas <- ncol(X)
    fit1 <- optim(par = par_iter[1:(n_betas + 1)], fn = Q_regression, method = "BFGS", 
                  X = X, Y = Y, T_Estep = t)  # M1-step
    fit2 <- optim(par = par_iter[(n_betas + 2):(n_betas + 3)], fn = Q_logistic, method = "BFGS", 
                  G = G, T_Estep = t)  # M2-step
    
    ## Update parameters
    par_prev <- par_iter
    par_iter <- c(fit1$par, fit2$par)
    
    iter <- iter + 1
  }
  
  # Variance estimation (moved outside convergence check)
  variance_estimates <- NA
  laplace_variance <- NA
  
  # Complete data score function
  complete_score <- function(par, X, Y, G) {
    n_betas <- ncol(X)
    beta_vec <- par[1:n_betas]
    sigma <- max(par[n_betas + 1], 1e-4)
    eta0  <- par[n_betas + 2]
    eta1  <- par[n_betas + 3]
    
    mu <- as.numeric(X %*% beta_vec)
    h_val <- h(G, eta0, eta1)
    phi_val <- dnorm(Y, mu, sigma)
    p_y_val <- p.y(Y)
    
    # E-step latent weights
    t <- (h_val * phi_val) / (h_val * phi_val + (1 - h_val) * p_y_val)
    
    # Initialize score matrix
    score_matrix <- matrix(0, nrow = length(Y), ncol = n_params)
    
    # Regression score contributions (betas)
    for (j in 1:n_betas) {
      score_matrix[, j] <- t * X[, j] * (Y - mu) / sigma^2
    }
    
    # Sigma score
    score_matrix[, n_betas + 1] <- t * (-1/sigma + (Y - mu)^2 / sigma^3)
    
    # Logistic part
    score_matrix[, n_betas + 2] <- (t - h_val)  # d_eta0
    score_matrix[, n_betas + 3] <- (t - h_val) * G  # d_eta1
    
    return(score_matrix)
  }
  # browser()
  
  # Complete data information - CORRECTED VERSION
  complete_info <- function(par, X, Y, G, Seed.vec = NULL) {
    n_betas <- ncol(X)
    beta_vec <- par[1:n_betas]
    sigma <- max(par[n_betas + 1], 1e-4)
    eta0  <- par[n_betas + 2]
    eta1  <- par[n_betas + 3]
    
    mu <- as.numeric(X %*% beta_vec)
    h_val <- h(G, eta0, eta1)
    t <- prob_t_is_1(X, Y, G, par, p.y, Seed.vec)  # 强制 seed
    W <- diag(as.numeric(t))
    I_beta <- t(X) %*% W %*% X / sigma^2
    
    I_sigma <- sum(t * 2 / sigma^2)
    
    w <- as.numeric(h_val * (1 - h_val))
    I_eta <- matrix(0, 2, 2)
    I_eta[1,1] <- sum(w)
    I_eta[1,2] <- I_eta[2,1] <- sum(w * G)
    I_eta[2,2] <- sum(w * G^2)
    
    info <- matrix(0, n_betas + 3, n_betas + 3)
    info[1:n_betas, 1:n_betas] <- I_beta
    info[n_betas + 1, n_betas + 1] <- I_sigma
    info[(n_betas + 2):(n_betas + 3), (n_betas + 2):(n_betas + 3)] <- I_eta
    
    # small ridge for numerical stability
    eps <- 1e-6
    info <- info + diag(eps, nrow(info))
    
    return(info)
  }
  
  scores <- complete_score(par_iter, X, Y, G)
  B <- t(scores) %*% scores
  C <- complete_info(par_iter, X, Y, G, Seed.vec)
  eig <- eigen(C, symmetric=TRUE)$values
  cond <- max(abs(eig)) / min(abs(eig))
  if (is.nan(cond) || cond > 1e12) {
    warning("C is ill-conditioned: adding ridge")
    C <- C + diag(1e-6, nrow(C))
  }
  laplace_variance <- tryCatch(solve(C), error = function(e) MASS::ginv(C))

  # Return results
  return(list(
    par = par_iter,
    iter = iter,
    # converged = (iter < 100),
    laplace_variance = MASS::ginv(C),
    n_covariates = n_covariates
  ))
}