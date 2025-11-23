# EM Algorithm Implementation with Variance Estimation
# Main EM estimation function
EM_estimation_Seed_ConfidenceMeasure_MLRandNonNormal <- function(X1, X2, Y, G, p.y, Seed.vec) {
  # Helper functions
  # Safe h with clipping
  h <- function(g, eta0, eta1) {
    p <- 1 / (1 + exp(-(eta0 + eta1 * g)))
    p <- pmin(pmax(p, 1e-8), 1 - 1e-8)
    return(p)
  }
  
  # Safe phi (density of N(beta_0 + beta_1*x, sigma^2))
  phi <- function(y, x1, x2, beta0, beta1, beta2, sigma) {
    mu <- beta0 + beta1 * x1 + beta2 * x2
    dnorm(y, mean = mu, sd = sigma)
  }
  
  # Observed data log-likelihood for a single observation
  log_likelihood.1 <- function(param, p.y) {
    beta_0 <- param[1]
    beta_1 <- param[2]
    beta_2 <- param[3]
    sigma <- param[4]
    eta_0 <- param[5]
    eta_1 <- param[6]
    x1 <- param[7]
    x2 <- param[8]
    y <- param[9]
    g <- param[10]
    
    term.1 <- (2*pi*sigma^2)^(-1/2) * exp(-1/(2*sigma^2) * (y-beta_0-beta_1*x1-beta_2*x2)^2) * h(g, eta_0, eta_1)
    term.2 <- p.y(y) * (1 - h(g, eta_0, eta_1))
    return(log(term.1 + term.2))
  }
  
  # E-step function
  prob_t_is_1 <- function(x1, x2, y, g, par, p.y) {
    beta_0 <- par[1]
    beta_1 <- par[2]
    beta_2 <- par[3]
    sigma <- max(par[4], 1e-4)
    eta_0 <- par[5]
    eta_1 <- par[6]
    
    h_val <- h(g, eta_0, eta_1)
    phi_val <- pmax(phi(y, x1, x2, beta_0, beta_1, beta_2, sigma), 1e-8)
    p_y_safe <- pmax(p.y(y), 1e-8)
    
    numer <- h_val * phi_val
    denom <- numer + (1 - h_val) * p_y_safe
    
    return(numer / denom)
  }
  
  # M-step objective function
  Q <- function(X1, X2, Y, G, T_Estep, par) {
    theta_0 <- par[1]
    theta_1 <- par[2]
    theta_2 <- par[3]
    sigma <- max(par[4], 1e-4)  # enforce positivity
    eta_0 <- par[5]
    eta_1 <- par[6]
    
    h_val <- h(G, eta_0, eta_1)
    phi_val <- pmax(phi(Y, X1, X2, theta_0, theta_1, theta_2, sigma), 1e-8)  # avoid log(0)
    
    loglik <- T_Estep * (log(h_val + 1e-8) + log(phi_val)) +
      (1 - T_Estep) * log(1 - h_val + 1e-8)
    return(-sum(loglik))
  }
  
  Q_regression <- function(par_beta_sigma, X1, X2, Y, T_Estep) {
    beta0 <- par_beta_sigma[1]
    beta1 <- par_beta_sigma[2]
    beta2 <- par_beta_sigma[3]
    sigma <- max(par_beta_sigma[4], 1e-4)
    mu <- beta0 + beta1 * X1 + beta2 * X2
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
  if (length(unique(X2[Seed.vec==1]))==1) { # if seeds' X2 are equal
    fit_lm <- lm(Y[Seed.vec==1] ~ X1[Seed.vec==1])
    par_iter <- as.numeric(
      c(fit_lm$coefficients, 1, sqrt(sum(fit_lm$residuals^2)/fit_lm$df.residual),
        0, 0)
    )
  } else {
    fit_lm <- lm(Y[Seed.vec==1] ~ X1[Seed.vec==1] + X2[Seed.vec==1])
    par_iter <- as.numeric(
      c(fit_lm$coefficients, 
        sqrt(sum(fit_lm$residuals^2)/fit_lm$df.residual),
        0,0)
    )
  }
  par_prev <- rep(0, 6)
  iter <- 1
  
  # EM iterations
  while (sqrt(sum((par_iter - par_prev)^2)) > 0.01 & iter < 100) {
    # E-step
    t <- prob_t_is_1(X1, X2, Y, G, par_iter, p.y); t[Seed.vec==1] <- 1 # Force the latent class probability of seed to be 1
    
    # Extract current eta values from par_iter
    eta0 <- par_iter[4]
    eta1 <- par_iter[5]
    
    # M-step
    fit1 <- optim(par = par_iter[1:4], fn = Q_regression, method = "BFGS", X1 = X1, X2 = X2, Y=Y, T_Estep=t)  # M1-step
    fit2 <- optim(par = par_iter[5:6], fn = Q_logistic, method = "BFGS", G=G, T_Estep=t)  # M2-step
    # ## M1-step: Update theta parameters (beta0, beta1, sigma)
    # Q1 <- function(par_beta_sigma) {
    #   Q(X, Y, G, t, c(par_beta_sigma, eta0, eta1))
    # }
    # fit1 <- optim(par = par_iter[1:3], fn = Q1, method = "BFGS")
    # 
    # ## M2-step: Update eta parameters (eta0, eta1)
    # Q2 <- function(par_eta) {
    #   Q(X, Y, G, t, c(fit1$par, par_eta))
    # }
    # fit2 <- optim(par = par_iter[4:5], fn = Q2, method = "BFGS")
    
    ## Update parameters
    par_prev <- par_iter
    par_iter <- c(fit1$par, fit2$par)
    
    iter <- iter + 1
  }
  
  # Variance estimation (moved outside convergence check)
  variance_estimates <- NA
  laplace_variance <- NA
  
  # Complete data score function
  complete_score <- function(par, X1, X2, Y, G) {
    beta0 <- par[1]
    beta1 <- par[2]
    beta2 <- par[3]
    sigma <- max(par[4], 1e-4)
    eta0  <- par[5]
    eta1  <- par[6]
    
    mu <- beta0 + beta1 * X1 + beta2 * X2
    h_val <- h(G, eta0, eta1)
    phi_val <- dnorm(Y, mu, sigma)
    p_y_val <- p.y(Y)
    
    # E-step latent weights
    t <- (h_val * phi_val) / (h_val * phi_val + (1 - h_val) * p_y_val)
    
    # Regression score contributions
    d_beta0 <- t * (Y - mu) / sigma^2
    d_beta1 <- t * X1 * (Y - mu) / sigma^2
    d_beta2 <- t * X2 * (Y - mu) / sigma^2
    d_sigma <- t * (-1/sigma + (Y - mu)^2 / sigma^3)
    
    # Logistic part
    d_eta0 <- (t - h_val)
    d_eta1 <- (t - h_val) * G
    
    cbind(d_beta0, d_beta1, d_beta2, d_sigma, d_eta0, d_eta1)
  }
  
  # Complete data information
  complete_info <- function(par, X1, X2, Y, G) {
    beta0 <- par[1]
    beta1 <- par[2]
    beta2 <- par[3]
    sigma <- max(par[4], 1e-4)
    eta0  <- par[5]
    eta1  <- par[6]
    
    mu <- beta0 + beta1 * X1 + beta2 * X2
    h_val <- h(G, eta0, eta1)
    t <- prob_t_is_1(X1, X2, Y, G, par, p.y)
    
    ## --- Regression block (β0, β1, β2) ---
    Xmat <- cbind(1, X1, X2)   # design matrix
    I_beta <- matrix(0, 3, 3)
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
    info <- matrix(0, 6, 6)
    info[1:3, 1:3] <- I_beta
    info[4, 4]     <- I_sigma
    info[5:6, 5:6] <- I_eta
    
    info
  }
  
  scores <- complete_score(par_iter, X1, X2, Y, G)
  outer_prod <- apply(scores, 1, function(x) x %*% t(x))
  B <- matrix(rowSums(outer_prod), 6, 6)
  C <- complete_info(par_iter, X1, X2, Y, G)
  
  # Return results
  return(list(
    par = par_iter,
    em_iter = iter,
    converged = (iter < 100),
    laplace_variance = MASS::ginv(C)
  ))
}


# Calculate Fisher information --------------------------------------------

calculate_observed_fisher <- function(X, Y, G, p.y, mle_params) {
  # Define the negative log-likelihood function for all observations
  neg_loglik <- function(params) {
    beta_0 <- params[1]
    beta_1 <- params[2]
    sigma <- max(params[3], 1e-4)  # Ensure sigma > 0
    eta_0 <- params[4]
    eta_1 <- params[5]
    
    # Calculate h(g) values
    h_vals <- h(G, eta_0, eta_1)
    
    # Calculate phi(y|x) values
    phi_vals <- pmax(phi(Y, X, beta_0, beta_1, sigma), 1e-8)
    
    # Calculate p.y(y) values
    p_y_vals <- pmax(p.y(Y), 1e-8)
    
    # Calculate the complete data negative log-likelihood
    term1 <- h_vals * phi_vals
    term2 <- (1 - h_vals) * p_y_vals
    loglik <- log(term1 + term2)
    
    return(-sum(loglik))
  }
  
  # Compute Hessian numerically
  require(numDeriv)
  hessian <- hessian(neg_loglik, mle_params)
  
  # Observed Fisher information is the Hessian of the negative log-likelihood
  observed_fisher <- hessian
  
  # Inverse of Fisher information is the variance-covariance matrix
  cov_matrix <- solve(observed_fisher)
  
  return(list(observed_fisher = observed_fisher, 
              covariance_matrix = cov_matrix))
}