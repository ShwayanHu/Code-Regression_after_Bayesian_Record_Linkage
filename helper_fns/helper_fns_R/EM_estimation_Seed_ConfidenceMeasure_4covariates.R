EM_estimation_Seed_ConfidenceMeasure_4covariates <- function(X1, X2, X3, X4, Y, G, p.y, Seed.vec) {
  # Helper functions
  # Safe h with clipping
  h <- function(g, eta0, eta1) {
    p <- 1 / (1 + exp(-(eta0 + eta1 * g)))
    p <- pmin(pmax(p, 1e-8), 1 - 1e-8)
    return(p)
  }
  
  # Safe phi (density of N(beta_0 + beta_1*x, sigma^2))
  phi <- function(y, x1, x2, x3, x4, beta0, beta1, beta2, beta3, beta4, sigma) {
    mu <- beta0 + beta1 * x1 + beta2 * x2 + beta3 * x3 + beta4 * x4
    dnorm(y, mean = mu, sd = sigma)
  }
  
  # Observed data log-likelihood for a single observation
  log_likelihood.1 <- function(param, p.y) {
    beta_0 <- param[1]
    beta_1 <- param[2]
    beta_2 <- param[3]
    beta_3 <- param[4]
    beta_4 <- param[5]
    sigma <- param[6]
    eta_0 <- param[7]
    eta_1 <- param[8]
    x1 <- param[9]
    x2 <- param[10]
    x3 <- param[11]
    x4 <- param[12]
    y <- param[13]
    g <- param[14]
    
    term.1 <- (2*pi*sigma^2)^(-1/2) * exp(-1/(2*sigma^2) * (y-beta_0-beta_1*x1-beta_2*x2-beta_3*x3-beta_4*x4)^2) * h(g, eta_0, eta_1)
    term.2 <- p.y(y) * (1 - h(g, eta_0, eta_1))
    return(log(term.1 + term.2))
  }
  
  # E-step function
  prob_t_is_1 <- function(x1, x2, x3, x4, y, g, par, p.y) {
    beta_0 <- par[1]
    beta_1 <- par[2]
    beta_2 <- par[3]
    beta_3 <- par[4]
    beta_4 <- par[5]
    sigma <- max(par[6], 1e-4)
    eta_0 <- par[7]
    eta_1 <- par[8]
    # browser()
    h_val <- h(g, eta_0, eta_1)
    phi_val <- pmax(phi(y, x1, x2, x3, x4, beta_0, beta_1, beta_2, beta_3, beta_4, sigma), 1e-8)
    p_y_safe <- pmax(p.y(y), 1e-8)
    
    numer <- h_val * phi_val
    denom <- numer + (1 - h_val) * p_y_safe
    
    return(numer / denom)
  }
  
  # M-step objective function
  Q <- function(X1, X2, X3, X4, Y, G, T_Estep, par) {
    theta_0 <- par[1]
    theta_1 <- par[2]
    theta_2 <- par[3]
    theta_3 <- par[4]
    theta_4 <- par[5]
    sigma <- max(par[6], 1e-4)  # enforce positivity
    eta_0 <- par[7]
    eta_1 <- par[8]
    
    h_val <- h(G, eta_0, eta_1)
    phi_val <- pmax(phi(Y, X1, X2, X3, X4, theta_0, theta_1, theta_2, theta_3, theta_4, sigma), 1e-8)  # avoid log(0)
    
    loglik <- T_Estep * (log(h_val + 1e-8) + log(phi_val)) +
      (1 - T_Estep) * log(1 - h_val + 1e-8)
    return(-sum(loglik))
  }
  
  Q_regression <- function(par_beta_sigma, X1, X2, X3, X4, Y, T_Estep) {
    beta0 <- par_beta_sigma[1]
    beta1 <- par_beta_sigma[2]
    beta2 <- par_beta_sigma[3]
    beta3 <- par_beta_sigma[4]
    beta4 <- par_beta_sigma[5]
    sigma <- max(par_beta_sigma[6], 1e-4)
    mu <- beta0 + beta1 * X1 + beta2 * X2 + beta3 * X3 + beta4 * X4
    -sum(T_Estep * dnorm(Y, mu, sigma, log=TRUE))  # 仅回归部分
  }
  
  Q_logistic <- function(par_eta, G, T_Estep) {
    eta0 <- par_eta[1]
    eta1 <- par_eta[2]
    h_val <- h(G, eta0, eta1)
    -sum(T_Estep * log(h_val) + (1 - T_Estep) * log(1 - h_val))  # 仅Logistic部分
  }
  # model: Y = beta_0 + beta_1 * X + epsilon. epsilon ~ N(0,sigma^2)
  # Pr(t = 1|g) = 1/(1 + exp(-(eta_0 + eta_1 * g)))
  
  # Initialization
  # browser()
  fit_lm <- lm(Y[Seed.vec==1] ~ X1[Seed.vec==1] + X2[Seed.vec==1] + X3[Seed.vec==1] + X4[Seed.vec==1] + 0)
  # fit_lm <- lm(Y ~ X1 + X2 + X3 + X4)
  par_iter <- as.numeric(
    c(fit_lm$coefficients, sqrt(sum(fit_lm$residuals^2)/fit_lm$df.residual), 0, 0, 0)
  )
  par_iter[is.na(par_iter)] <- 1
  
  par_prev <- rep(0, length(par_iter))
  iter <- 1
  
  # EM iterations
  while (any(abs(par_iter - par_prev) > 1e-3) & iter <= 2000) {
    # E-step
    t <- prob_t_is_1(X1, X2, X3, X4, Y, G, par_iter, p.y); t[Seed.vec==1] <- 1 # Force the latent class probability of seed to be 1
    
    # Extract current eta values from par_iter
    eta0 <- par_iter[7]
    eta1 <- par_iter[8]
    
    # M-step
    # browser()
    fit1 <- optim(par = par_iter[1:6], fn = Q_regression, method = "BFGS", X1 = X1, X2 = X2, X3 = X3, X4 = X4, Y=Y, T_Estep=t)  # M1-step
    fit2 <- optim(par = par_iter[7:8], fn = Q_logistic, method = "BFGS", G=G, T_Estep=t)  # M2-step
    
    ## Update parameters
    par_prev <- par_iter
    par_iter <- c(fit1$par, fit2$par)
    
    iter <- iter + 1
  }
  
  # Variance estimation (moved outside convergence check)
  variance_estimates <- NA
  laplace_variance <- NA
  
  # Complete data score function
  complete_score <- function(par, X1, X2, X3, X4, Y, G) {
    beta0 <- par[1]
    beta1 <- par[2]
    beta2 <- par[3]
    beta3 <- par[4]
    beta4 <- par[5]
    sigma <- max(par[6], 1e-4)
    eta0  <- par[7]
    eta1  <- par[8]
    
    mu <- beta0 + beta1 * X1 + beta2 * X2 + beta3 * X3 + beta4 * X4
    h_val <- h(G, eta0, eta1)
    phi_val <- dnorm(Y, mu, sigma)
    p_y_val <- p.y(Y)
    
    # E-step latent weights
    t <- (h_val * phi_val) / (h_val * phi_val + (1 - h_val) * p_y_val)
    
    # Regression score contributions
    d_beta0 <- t * (Y - mu) / sigma^2
    d_beta1 <- t * X1 * (Y - mu) / sigma^2
    d_beta2 <- t * X2 * (Y - mu) / sigma^2
    d_beta3 <- t * X3 * (Y - mu) / sigma^2
    d_beta4 <- t * X4 * (Y - mu) / sigma^2
    d_sigma <- t * (-1/sigma + (Y - mu)^2 / sigma^3)
    
    # Logistic part
    d_eta0 <- (t - h_val)
    d_eta1 <- (t - h_val) * G
    
    cbind(d_beta0, d_beta1, d_beta2, d_beta3, d_beta4, d_sigma, d_eta0, d_eta1)
  }
  
  # Complete data information
  complete_info <- function(par, X1, X2, X3, X4, Y, G) {
    beta0 <- par[1]
    beta1 <- par[2]
    beta2 <- par[3]
    beta3 <- par[4]
    beta4 <- par[5]
    sigma <- max(par[6], 1e-4)
    eta0  <- par[7]
    eta1  <- par[8]
    
    mu <- beta0 + beta1 * X1 + beta2 * X2 + beta3 * X3 + beta4 * X4
    h_val <- h(G, eta0, eta1)
    t <- prob_t_is_1(X1, X2, X3, X4, Y, G, par, p.y)
    
    ## --- Regression block (β0, β1, β2) ---
    Xmat <- cbind(1, X1, X2, X3, X4)   # design matrix
    I_beta <- matrix(0, 5, 5)
    for (i in 1:length(Y)) {
      xi <- Xmat[i, , drop=FALSE]
      I_beta <- I_beta + t[i] * (t(xi) %*% xi) / sigma^2
    }
    
    ## --- Sigma block ---
    I_sigma <- sum(t * ( -1/sigma^2 + 2*(Y - mu)^2/sigma^4 ))   # expected info
    
    ## --- Logistic block (eta0, eta1) ---
    w <- h_val * (1 - h_val)  # logistic variance weights
    I_eta <- matrix(0, 2, 2)
    I_eta[1,1] <- sum(w)
    I_eta[1,2] <- I_eta[2,1] <- sum(w * G)
    I_eta[2,2] <- sum(w * G^2)
    
    ## --- Combine ---
    info <- matrix(0, 8, 8)
    info[1:5, 1:5] <- I_beta
    info[6, 6]     <- I_sigma
    info[7:8, 7:8] <- I_eta
    
    info
  }
  
  scores <- complete_score(par_iter, X1, X2, X3, X4, Y, G)
  outer_prod <- apply(scores, 1, function(x) x %*% t(x))
  B <- matrix(rowSums(outer_prod), 8, 8)
  C <- complete_info(par_iter, X1, X2, X3, X4, Y, G)
  
  # Return results
  return(list(
    par = par_iter,
    iter = iter,
    # converged = (iter < 100),
    laplace_variance = MASS::ginv(C)
  ))
}