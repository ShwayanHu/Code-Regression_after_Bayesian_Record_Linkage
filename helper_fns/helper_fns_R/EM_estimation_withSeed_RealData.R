EM_estimation_withSeed_flexible <- function(X, Y, Seed.vec, p.y) {
  # X is a n by p matrix, with first column 1's
  
  # Helper functions
  # phi <- function(y, x1, x2, x3, beta_0, beta_1, beta_2, beta_3, sigma) {
  #   return(dnorm(y, mean=beta_0+beta_1*x1 + beta_2*x2 + beta_3*x3, sd=sigma))
  # }
  # 
  # # Flexible version that works with matrix X
  # phi_flexible <- function(y, x_row, beta_vec, sigma) {
  #   return(dnorm(y, mean=sum(x_row * beta_vec), sd=sigma))
  # }
  
  logit <- function(p) log(p / (1 - p))
  
  log_likelihood_withSeed <- function(param, X, Y, t.vec, p.y) { 
    # param = (beta_0, beta_1, ..., beta_p, sigma, prob). 
    # sigma and prob are log-transformed and logit-transformed respectively
    p <- ncol(X)
    beta_vec <- param[1:p]
    sigma <- exp(param[p+1]) # sigma is log-transformed
    prob <- plogis(param[p+2]) # prob is logit-transformed
    
    m <- sum(t.vec) # number of links
    X.link <- X[t.vec, , drop=FALSE]
    Y.link <- Y[t.vec] # linked regression data
    mu.link <- X.link %*% beta_vec
    term.1 <- m*log(prob) + sum(dnorm(Y.link, mean=mu.link, sd=sigma, log=TRUE))
    
    n <- sum(!t.vec) # number of non-links
    y.nonlink <- Y[!t.vec] # non-linked regression data
    term.2 <- n*log(1-prob) + sum(log(p.y(y.nonlink)))
    return(term.1 + term.2)
  }
  
  prob_t_is_1_withSeed <- function(x_row, y, param, p.y) {
    p <- ncol(X)
    beta_vec <- param[1:p]
    sigma <- exp(param[p+1]) # sigma is log-transformed
    prob <- pmin(pmax(plogis(param[p+2]), 1e-8), 1-1e-8) # prob is logit-transformed
    
    mu <- sum(x_row * beta_vec)
    numerator <- prob * dnorm(y, mean=mu, sd=sigma)
    denomenator <- prob * dnorm(y, mean=mu, sd=sigma) + (1-prob)*p.y(y)
    return(numerator/denomenator)
  }
  
  # M-step objective function
  Q_withSeed <- function(X, Y, T_Estep.vec, param, p.y) {
    # browser()
    p <- ncol(X)
    beta_vec <- param[1:p]
    sigma <- exp(param[p+1]) # sigma is log-transformed
    prob <- pmin(pmax(plogis(param[p+2]), 1e-8), 1-1e-8) # prob is logit-transformed
    
    mu <- as.matrix(X) %*% as.matrix(beta_vec,ncol=1)
    if (any(is.nan(dnorm(Y, mean = mu, sd = sigma, log = TRUE)))){browser()}
    return(sum(
      T_Estep.vec*(log(prob) + dnorm(Y, mean=mu, sd=sigma, log=TRUE)) +
        (1-T_Estep.vec) * (log(p.y(Y)) + log(1-prob))
    ))
  }
  
  # calculate_observed_fisher_withSeed <- function(X, Y, Seed.vec, p.y, mle_params) {
  #   # Compute the Hessian of the negative log-likelihood at MLE
  #   # Using numerical differentiation
  # 
  #   # Define the negative log-likelihood function
  #   neg_loglik <- function(params) {
  #     p <- ncol(X)
  #     beta_vec <- params[1:p]
  #     sigma <- params[p+1]
  #     prob <- pmin(pmax(plogis(params[p+2]), 1e-8), 1-1e-8) # prob is logit-transformed
  # 
  #     # Calculate t_i probabilities at current parameters
  #     t_probs <- sapply(1:nrow(X), function(i) {
  #       prob_t_is_1_withSeed(x_row=X[i,], y=Y[i], param = params, p.y = p.y)
  #     })
  #     t_probs[Seed.vec==1] <- 1 # For seeds, t_i is always 1
  # 
  #     # Calculate the complete data negative log-likelihood
  #     mu <- as.matrix(X) %*% as.matrix(beta_vec, col=1)
  #     term1 <- sum(t_probs * dnorm(Y, mean = mu, sd = sigma, log = TRUE))
  #     term2 <- sum((1-t_probs) * log(p.y(Y)))
  #     term3 <- sum(t_probs * log(prob) + (1-t_probs) * log(1-prob))
  # 
  #     return(-(term1 + term2 + term3))
  #   }
  # 
  #   # Compute Hessian numerically
  #   hessian <- numDeriv::hessian(neg_loglik, mle_params)
  # 
  #   # Observed Fisher information is the Hessian of the negative log-likelihood
  #   observed_fisher <- hessian
  # 
  #   # Inverse of Fisher information is the variance-covariance matrix
  #   cov_matrix <- MASS::ginv(observed_fisher)
  # 
  #   return(list(observed_fisher = observed_fisher,
  #               covariance_matrix = cov_matrix))
  # }
  
  calculate_observed_fisher_withSeed <- function(X, Y, Seed.vec, p.y, mle_params) {
    # Compute the Hessian of the negative log-likelihood at MLE
    # Using numerical differentiation
    
    # Define the negative log-likelihood function
    neg_loglik <- function(params) {
      p <- ncol(X)
      beta_vec <- params[1:p]
      sigma <- exp(params[p+1])  # ✓ 修复：添加exp变换
      prob <- pmin(pmax(plogis(params[p+2]), 1e-8), 1-1e-8) # prob is logit-transformed
      
      # Calculate t_i probabilities at current parameters
      t_probs <- sapply(1:nrow(X), function(i) {
        prob_t_is_1_withSeed(x_row=X[i,], y=Y[i], param = params, p.y = p.y)
      })
      t_probs[Seed.vec==1] <- 1 # For seeds, t_i is always 1
      
      # Calculate the complete data negative log-likelihood
      mu <- as.matrix(X) %*% as.matrix(beta_vec, col=1)
      term1 <- sum(t_probs * dnorm(Y, mean = mu, sd = sigma, log = TRUE))
      term2 <- sum((1-t_probs) * log(pmax(p.y(Y), 1e-10)))  # ✓ 修复：添加pmax避免log(0)
      term3 <- sum(t_probs * log(prob) + (1-t_probs) * log(1-prob))
      
      return(-(term1 + term2 + term3))
    }
    
    # ✓ 修复：需要将mle_params转换为变换空间
    p <- ncol(X)
    params_transformed <- mle_params
    params_transformed[p+1] <- log(mle_params[p+1])  # sigma -> log(sigma)
    params_transformed[p+2] <- qlogis(mle_params[p+2])  # prob -> logit(prob)
    
    # Compute Hessian numerically
    hessian <- numDeriv::hessian(neg_loglik, params_transformed)  # ✓ 修复：使用变换后的参数
    
    # ✓ 修复：检查并修正非正定Hessian
    eigen_vals <- eigen(hessian, symmetric = TRUE, only.values = TRUE)$values
    if (any(eigen_vals <= 1e-8)) {
      ridge <- abs(min(eigen_vals)) + 0.1
      hessian <- hessian + diag(nrow(hessian)) * ridge
    }
    
    # Observed Fisher information is the Hessian of the negative log-likelihood
    observed_fisher <- hessian
    
    # Inverse of Fisher information is the variance-covariance matrix
    cov_matrix <- tryCatch({
      solve(observed_fisher)
    }, error = function(e) {
      MASS::ginv(observed_fisher)
    })
    
    # ✓ 修复：Delta方法转换回原始参数空间
    J <- diag(p+2)
    J[p+1, p+1] <- mle_params[p+1]  # d(exp(log_sigma))/d(log_sigma) = sigma
    J[p+2, p+2] <- mle_params[p+2] * (1 - mle_params[p+2])  # d(logit^-1)/d(logit)
    cov_matrix <- J %*% cov_matrix %*% t(J)
    
    return(list(observed_fisher = observed_fisher,
                covariance_matrix = cov_matrix))
  }
  
  # model: model: Y = beta_0 + beta_1 * X1 + ... + beta_p * Xp + epsilon. epsilon ~ N(0,sigma^2)
  # Seed.vec: Bool vector indicating seeds
  # Initialization
  # browser()
  seed_data <- data.frame(
    Y = Y[Seed.vec==1],
    X[Seed.vec==1, -1]
  )
  fit <- lm(Y ~ ., data = seed_data)
  p <- ncol(X)
  par_iter <- as.numeric(c(fit$coefficients, log(sqrt(sum(fit$residuals^2)/fit$df.residual)), logit(0.5))) # beta_0, ..., beta_p, log(sigma), logit(prob)
  par_iter[is.na(par_iter)] <- 1
  par_prev <- rep(0, p+2)
  iter <- 1
  
  while (any(abs(par_iter - par_prev) > 1e-3) & iter < 2000) {
    # E-step
    ## Calculate probability that t==1
    T_Estep.vec <- sapply(1:nrow(X), function(i) {
      prob_t_is_1_withSeed(x_row=X[i,], y=Y[i], param = par_iter, p.y = p.y)
    })
    T_Estep.vec[Seed.vec==1] <- 1
    
    optimization_target <- function(param) {
      return(Q_withSeed(X, Y, T_Estep.vec, param, p.y))
    }
    fit <- optim(par=par_iter, fn=optimization_target, method="BFGS", control=list(fnscale=-1))
    par_prev <- par_iter; par_iter <- fit$par; iter <- iter + 1
  }
  
  par_iter[p+1] <- exp(par_iter[p+1]) # back-transform sigma
  par_iter[p+2] <- plogis(par_iter[p+2]) # back-transform prob
  laplace_variance <- calculate_observed_fisher_withSeed(X, Y, Seed.vec, p.y, par_iter)$covariance_matrix
  
  result <- list(
    param = par_iter,
    iter = iter,
    laplace_variance = laplace_variance
  )
  return(result)
}