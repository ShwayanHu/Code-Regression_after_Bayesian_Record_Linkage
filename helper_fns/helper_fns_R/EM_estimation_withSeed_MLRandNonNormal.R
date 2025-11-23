# Introduction ------------------------------------------------------------
# Implementation of EM algorithm without confidence measure, but with seeds.

EM_estimation_withSeed_MLRandNonNormal <- function(X1.vec, X2.vec, Y.vec, Seed.vec, p.y) {
  # Helper functions
  phi <- function(y, x1, x2, beta_0, beta_1, beta_2 ,sigma) {
    return(dnorm(y, mean=beta_0+beta_1*x1 + beta_2*x2, sd=sigma))
  }
  
  logit <- function(p) log(p / (1 - p))
  
  log_likelihood_withSeed <- function(param, X1.vec, X2.vec, Y.vec, t.vec, p.y) { 
    # param = (beta_0, beta_1, sigma, prob). 
    #    `prob` models the latent class variable via I(t_i=1)~Bernoulli(prob)
    # X.vec, Y.vec = regression column vectors
    # t.vec = vector indicating latent class for data. For seed pairs, the t.vec==1 and seed.vec==TRUE.
    # phi = distribution of y for pairs
    # p.y = distribution of y for non-pairs
    beta_0 <- param[1]
    beta_1 <- param[2]
    beta_2 <- param[3]
    sigma <- exp(param[4]) # sigma is log-transformed
    prob <- plogis(param[5]) # prob is logit-transformed
    
    m <- sum(t.vec) # number of links
    x1.link <- X1.vec[t.vec]; x2.link <- X2.vec[t.vec]; y.link <- Y.vec[t.vec] # linked regression data
    term.1 <- m*log(prob) + sum(dnorm(y.link, mean=beta_0+beta_1*x1.link + beta_2*x2.link, sd=sigma, log=TRUE))
    n <- sum(!t.vec) # number of non-links
    y.nonlink <- Y.vec[!t.vec] # non-linked regression data
    term.2 <- n*log(1-prob) + sum(log(p.y(y.nonlink)))
    return(term.1 + term.2)
  }
  
  prob_t_is_1_withSeed <- function(x1, x2, y, param, p.y) {
    beta_0 <- param[1]
    beta_1 <- param[2]
    beta_2 <- param[3]
    sigma <- exp(param[4]) # sigma is log-transformed
    prob <- plogis(param[5]) # prob is logit-transformed
    
    numerator <- prob * dnorm(y, mean=beta_0 + beta_1*x1 + beta_2*x2, sd=sigma)
    denomenator <- prob * dnorm(y, mean=beta_0 + beta_1*x1 + beta_2*x2, sd=sigma) + (1-prob)*p.y(y)
    return(numerator/denomenator)
  }
  
  # M-step objective function
  Q_withSeed <- function(X1.vec, X2.vec, Y.vec, T_Estep.vec, param, p.y) {
    beta_0 <- param[1]
    beta_1 <- param[2]
    beta_2 <- param[3]
    sigma <- exp(param[4]) # sigma is log-transformed
    prob <- plogis(param[5]) # prob is logit-transformed
    if (any(is.nan(dnorm(Y.vec, mean = beta_0 + beta_1 * X1.vec + beta_2 * X2.vec, sd = sigma, log = TRUE)))){browser()}
    return(sum(
      T_Estep.vec*(log(prob) + dnorm(Y.vec, mean=beta_0+beta_1*X1.vec + beta_2*X2.vec, sd=sigma, log=TRUE)) +
        (1-T_Estep.vec) * (log(p.y(Y.vec)) + log(1-prob))
    ))
  }
  
  calculate_observed_fisher_withSeed <- function(X1.vec, X2.vec, Y.vec, Seed.vec, p.y, mle_params) {
    # Compute the Hessian of the negative log-likelihood at MLE
    # Using numerical differentiation
    
    # Define the negative log-likelihood function
    neg_loglik <- function(params) {
      beta_0 <- params[1]
      beta_1 <- params[2]
      beta_2 <- params[3]
      sigma <- exp(params[4]) # sigma is log-transformed
      prob <- plogis(params[5]) # prob is logit-transformed
      
      # Calculate t_i probabilities at current parameters
      t_probs <- sapply(1:length(X1.vec), function(i) {prob_t_is_1_withSeed(x1=X1.vec[i], x2=X2.vec[i], y=Y.vec[i], param = params, p.y = p.y)})
      t_probs[Seed.vec==1] <- 1 # For seeds, t_i is always 1
      
      # Calculate the complete data negative log-likelihood
      term1 <- sum(t_probs * dnorm(Y.vec, mean = beta_0 + beta_1*X1.vec + beta_2*X2.vec, sd = sigma, log = TRUE))
      term2 <- sum((1-t_probs) * log(p.y(Y.vec)))
      term3 <- sum(t_probs * log(prob) + (1-t_probs) * log(1-prob))
      
      return(-(term1 + term2 + term3))
    }
    
    # Compute Hessian numerically
    hessian <- numDeriv::hessian(neg_loglik, mle_params)
    
    # Observed Fisher information is the Hessian of the negative log-likelihood
    observed_fisher <- hessian
    
    # Inverse of Fisher information is the variance-covariance matrix
    cov_matrix <- solve(observed_fisher)
    
    return(list(observed_fisher = observed_fisher, 
                covariance_matrix = cov_matrix))
  }
  
  # model: model: Y = beta_0 + beta_1 * X + epsilon. epsilon ~ N(0,sigma^2)
  # Seed.vec: Bool vector indicating seeds
  # Initialization
  fit <- lm(Y.vec[Seed.vec==1] ~ X1.vec[Seed.vec==1] + X2.vec[Seed.vec==1])
  par_iter <- as.numeric(c(fit$coefficients, log(sqrt(sum(fit$residuals^2)/fit$df.residual)), logit(0.5))) # beta_0, beta_1, log(sigma), logit(prob)
  par_prev <- rep(0,5)
  iter <- 1
  while (sqrt(sum((par_iter-par_prev)^2)) > 0.01 & iter < 100) {
    # E-step
    ## Calculate probability that t==1
    T_Estep.vec <- sapply(1:length(X1.vec), function(i) {prob_t_is_1_withSeed(x1=X1.vec[i], x2=X2.vec[i], y=Y.vec[i], param = par_iter, p.y = p.y)})
    T_Estep.vec[Seed.vec==1] <- 1
    
    optimization_target <- function(param) {
      return(-Q_withSeed(X1.vec, X2.vec, Y.vec, T_Estep.vec, param, p.y))
    }
    fit <- optim(par=par_iter, fn=optimization_target, method="BFGS")
    par_prev <- par_iter; par_iter <- fit$par; iter <- iter + 1
  }
  
  laplace_variance <- calculate_observed_fisher_withSeed(X1.vec, X2.vec, Y.vec, Seed.vec, p.y, par_iter)$covariance_matrix
  
  par_iter[4] <- exp(par_iter[4]) # back-transform sigma
  par_iter[5] <- plogis(par_iter[5]) # back-transform prob
  result <- list(
    param = par_iter,
    iter = iter,
    laplace_variance = laplace_variance
  )
  return(result)
}