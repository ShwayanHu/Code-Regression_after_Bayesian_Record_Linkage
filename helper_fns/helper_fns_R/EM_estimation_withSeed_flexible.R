EM_estimation_withSeed_flexible <- function(X, Y, Seed.vec, p.y) {
  phi <- function(y, x, param) {
    return(dnorm(y, mean=sum(x * param[1:(length(param)-2)]), sd=param[length(param)-1]))
  }
  
  logit <- function(p) log(p / (1 - p))
  
  log_likelihood_withSeed <- function(param, X, Y, t.vec, p.y) {
    m <- sum(t.vec)
    X.link <- X[t.vec == 1, , drop=FALSE]
    Y.link <- Y[t.vec == 1]
    term.1 <- m*log(param[length(param)]) + 
      sum(
        dnorm(
          Y.link, 
          mean = X.link %*% matrix(param[1:(length(param)-2)], ncol=1),
          sd = param[length(param)-1],
          log = TRUE
        )
      )
    n <- sum(!t.vec)
    Y.nonlink <- Y[!t.vec]
    term.2 <- n * log(1-param[length(param)]) +
      sum(log(p.y(Y.nonlink)))
    return(term.1 + term.2)
  }
  
  prob_t_is_1_withSeed <- function(X, Y, param, p.y, Seed.vec) {
    numerator <- param[length(param)] *
      dnorm(Y, mean = X %*% matrix(param[1:(length(param)-2)], ncol=1), sd = param[length(param)-1])
    denominator <- numerator + (1 - param[length(param)]) * p.y(Y)
    res <- numerator/denominator
    res[Seed.vec==1] <- 1
    return(res)
  }
  
  Q_withSeed <- function(X, Y, T_Estep.vec, param, p.y) {
    return(
      sum(
        T_Estep.vec * (
          log(param[length(param)]) +
            dnorm(Y, mean = X %*% matrix(param[1:(length(param)-2)], ncol=1), sd = param[length(param)-1], log = TRUE)
        ) +
          (1 - T_Estep.vec) * (
            log(p.y(Y)) + log(1 - param[length(param)])
          )
      )
    )
  }
  
  calculate_observed_fisher_withSeed <- function(X, Y, Seed.vec, p.y, mle_params) {
    neg_loglik <- function(params) {
      return(-log_likelihood_withSeed(params, X, Y, Seed.vec, p.y))
    }
    
    hessian <- numDeriv::hessian(neg_loglik, mle_params)
    observed_fisher <- hessian
    cov_matrix <- tryCatch({MASS::ginv(observed_fisher)}, error = function(e) {matrix(NA, nrow=length(mle_params), ncol=length(mle_params))})
    return(list(
      observed_fisher = observed_fisher,
      cov_matrix = cov_matrix
    ))
  }
  
  # EM algorithm
  ## Initialization
  fit <- lm(Y[Seed.vec==1] ~ X[Seed.vec==1, -1])
  par_iter <- as.numeric(c(fit$coefficients, summary(fit)$sigma, 0.5))
  par_prev <- rep(0, length(par_iter))
  iter <- 1
  
  while (sqrt(sum((par_iter - par_prev)^2)) > 1e-4 & iter <= 2e3) {
    # E-step
    T_Estep.vec <- sapply(
      1:nrow(X),
      function(i) {
        if (Seed.vec[i] == 1) {
          return(1)
        } else {
          return(prob_t_is_1_withSeed(X[i, , drop=FALSE], Y[i], par_iter, p.y))
        }
      }
    )
    
    # M-step
    optimization_target <- function(param) {
      return(-Q_withSeed(X, Y, T_Estep.vec, param, p.y))
    }
    fit <- optim(par = par_iter, fn = optimization_target, method = "L-BFGS-B", 
                 lower = c(rep(-Inf, length(par_iter)-2), 0.001, 0.001), 
                 upper = c(rep(Inf, length(par_iter)-2), Inf, 0.999))
    par_prev <- par_iter
    par_iter <- fit$par
    iter <- iter + 1
  }
  
  laplace_variance <- calculate_observed_fisher_withSeed(X, Y, Seed.vec, p.y, par_iter)$cov_matrix
  result <- list(
    param = par_iter,
    iter = iter,
    laplace_variance = laplace_variance
  )
  return(result)
}