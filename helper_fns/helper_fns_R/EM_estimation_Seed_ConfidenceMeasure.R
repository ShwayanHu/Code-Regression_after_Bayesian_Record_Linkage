# EM Algorithm Implementation with Variance Estimation

# Safe h with clipping
h <- function(g, eta0, eta1) {
  p <- 1 / (1 + exp(-(eta0 + eta1 * g)))
  p <- pmin(pmax(p, 1e-8), 1 - 1e-8)
  return(p)
}

# Safe phi (density of N(beta_0 + beta_1*x, sigma^2))
phi <- function(y, x, beta0, beta1, sigma) {
  mu <- beta0 + beta1 * x
  dnorm(y, mean = mu, sd = sigma)
}

# Observed data log-likelihood for a single observation
log_likelihood.1 <- function(param, p.y) {
  beta_0 <- param[1]
  beta_1 <- param[2]
  sigma <- param[3]
  eta_0 <- param[4]
  eta_1 <- param[5]
  x <- param[6]
  y <- param[7]
  g <- param[8]
  
  term.1 <- (2*pi*sigma^2)^(-1/2) * exp(-1/(2*sigma^2) * (y-beta_0-beta_1*x)^2) * h(g, eta_0, eta_1)
  term.2 <- p.y(y) * (1 - h(g, eta_0, eta_1))
  return(log(term.1 + term.2))
}

# E-step function
prob_t_is_1 <- function(x, y, g, par, p.y) {
  beta_0 <- par[1]
  beta_1 <- par[2]
  sigma <- max(par[3], 1e-4)
  eta_0 <- par[4]
  eta_1 <- par[5]
  
  h_val <- h(g, eta_0, eta_1)
  phi_val <- pmax(phi(y, x, beta_0, beta_1, sigma), 1e-8)
  p_y_safe <- pmax(p.y(y), 1e-8)
  
  numer <- h_val * phi_val
  denom <- numer + (1 - h_val) * p_y_safe
  
  return(numer / denom)
}

# M-step objective function
Q <- function(X, Y, G, T_Estep, par) {
  theta_0 <- par[1]
  theta_1 <- par[2]
  sigma <- max(par[3], 1e-4)  # enforce positivity
  eta_0 <- par[4]
  eta_1 <- par[5]
  
  h_val <- h(G, eta_0, eta_1)
  phi_val <- pmax(phi(Y, X, theta_0, theta_1, sigma), 1e-8)  # avoid log(0)
  
  loglik <- T_Estep * (log(h_val + 1e-8) + log(phi_val)) +
    (1 - T_Estep) * log(1 - h_val + 1e-8)
  return(-sum(loglik))
}

Q_regression <- function(par_beta_sigma, X, Y, T_Estep) {
  beta0 <- par_beta_sigma[1]
  beta1 <- par_beta_sigma[2]
  sigma <- max(par_beta_sigma[3], 1e-4)
  mu <- beta0 + beta1 * X
  -sum(T_Estep * dnorm(Y, mu, sigma, log=TRUE))  # 仅回归部分
}

Q_logistic <- function(par_eta, G, T_Estep) {
  eta0 <- par_eta[1]
  eta1 <- par_eta[2]
  h_val <- h(G, eta0, eta1)
  -sum(T_Estep * log(h_val) + (1 - T_Estep) * log(1 - h_val))  # 仅Logistic部分
}

# Main EM estimation function
EM_estimation_Seed_ConfidenceMeasure <- function(X, Y, G, p.y, Seed.vec) {
  # model: Y = beta_0 + beta_1 * X + epsilon. epsilon ~ N(0,sigma^2)
  # Pr(t = 1|g) = 1/(1 + exp(-(eta_0 + eta_1 * g)))
  
  # Initialization
  fit_lm <- lm(Y[Seed.vec==1] ~ X[Seed.vec==1])
  
  par_iter <- as.numeric(
    c(fit_lm$coefficients, 
      sqrt(sum(fit_lm$residuals^2)/fit_lm$df.residual),
      0,0)
  )
  par_prev <- rep(0, 5)
  iter <- 1
  
  # EM iterations
  while (sqrt(sum((par_iter - par_prev)^2)) > 0.01 & iter < 100) {
    # E-step
    t <- prob_t_is_1(X, Y, G, par_iter, p.y); t[Seed.vec==1] <- 1 # Force the latent class probability of seed to be 1
    
    # Extract current eta values from par_iter
    eta0 <- par_iter[4]
    eta1 <- par_iter[5]
    
    # M-step
    fit1 <- optim(par = par_iter[1:3], fn = Q_regression, method = "BFGS", X = X, Y=Y, T_Estep=t)  # M1-step
    fit2 <- optim(par = par_iter[4:5], fn = Q_logistic, method = "BFGS", G=G, T_Estep=t)  # M2-step
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
  complete_score <- function(par, X, Y, G) {
    beta0 <- par[1]
    beta1 <- par[2]
    sigma <- max(par[3], 1e-4)
    eta0 <- par[4]
    eta1 <- par[5]
    
    # Calculate all needed terms
    mu <- beta0 + beta1 * X
    h_val <- h(G, eta0, eta1)
    phi_val <- dnorm(Y, mu, sigma)
    p_y_val <- p.y(Y)
    
    # t_i estimates
    t <- (h_val * phi_val) / (h_val * phi_val + (1 - h_val) * p_y_val)
    
    # Derivatives
    # For normal model
    d_beta0 <- t * (Y - mu)/sigma^2
    d_beta1 <- t * X * (Y - mu)/sigma^2
    d_sigma <- t * (-1/sigma + (Y - mu)^2/sigma^3)
    
    # For logistic model
    dh <- h_val * (1 - h_val)
    d_eta0 <- (t - h_val) * dh / (h_val * (1 - h_val))
    d_eta1 <- (t - h_val) * G * dh / (h_val * (1 - h_val))
    
    cbind(d_beta0, d_beta1, d_sigma, d_eta0, d_eta1)
  }
  
  # Complete data information
  complete_info <- function(par, X, Y, G) {
    beta0 <- par[1]
    beta1 <- par[2]
    sigma <- max(par[3], 1e-4)
    eta0 <- par[4]
    eta1 <- par[5]
    
    n <- length(Y)
    mu <- beta0 + beta1 * X
    h_val <- h(G, eta0, eta1)
    t <- prob_t_is_1(X, Y, G, par, p.y)
    
    # Information matrix components
    I_beta <- matrix(0, 2, 2)
    I_beta[1,1] <- sum(t)/sigma^2
    I_beta[1,2] <- I_beta[2,1] <- sum(t * X)/sigma^2
    I_beta[2,2] <- sum(t * X^2)/sigma^2
    
    # Corrected sigma information (for sigma^2)
    I_sigma <- sum(t * (1/sigma^2 - 2*(Y - mu)^2/sigma^4 + (Y - mu)^2/sigma^6))
    
    I_eta <- matrix(0, 2, 2)
    w <- h_val * (1 - h_val)  # weights for logistic model
    I_eta[1,1] <- sum(w)
    I_eta[1,2] <- I_eta[2,1] <- sum(w * G)
    I_eta[2,2] <- sum(w * G^2)
    
    # Combine into full matrix
    info <- matrix(0, 5, 5)
    info[1:2,1:2] <- I_beta
    info[3,3] <- I_sigma
    info[4:5,4:5] <- I_eta
    
    info
  }
  
  scores <- complete_score(par_iter, X, Y, G)
  outer_prod <- apply(scores, 1, function(x) x %*% t(x))
  B <- matrix(rowSums(outer_prod), 5, 5)
  C <- complete_info(par_iter, X, Y, G)
  
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