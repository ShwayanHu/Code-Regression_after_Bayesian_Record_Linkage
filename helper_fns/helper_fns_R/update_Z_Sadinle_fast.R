update_Z_Sadinle_fast <- function(Z_tmp, n1, n2, gamma, m_tmp, u_tmp, beta_pi, alpha_pi) {
  gamma_array <- array(NA, dim = c(n1, n2, 4))
  
  for (f in 1:4) {
    field <- paste0("f", f)
    for (i in 1:nrow(gamma)) {
      r1 <- gamma$record_1[i]
      r2 <- gamma$record_2[i]
      gamma_array[r1, r2, f] <- gamma[[field]][i]
    }
  }
  
  for (j in 1:n2) {
    Z_minus_j <- Z_tmp[-j]
    n12 <- sum(Z_minus_j <= n1)  # fast n_12 function
    
    possible_qs <- c(1:n1, n1 + j)
    w_j <- numeric(length(possible_qs))
    
    linked_set <- Z_minus_j[Z_minus_j <= n1]
    
    for (idx in seq_along(possible_qs)) {
      q <- possible_qs[idx]
      if (q <= n1) {
        if (q %in% linked_set) {
          w_j[idx] <- 0
        } else {
          g <- gamma_array[q, j, ] + 1  # add 1 for indexing
          log_ratio <- sum(log(mapply(function(m, u, l) m[l]/u[l], 
                                      m_tmp, u_tmp, g)))
          w_j[idx] <- exp(log_ratio)
        }
      } else {
        prior_term <- (n1 - n12) * (n2 - n12 - 1 + beta_pi) / (n12 + alpha_pi)
        w_j[idx] <- prior_term
      }
    }
    
    w_j <- w_j / sum(w_j)
    if (j == 1) {
      print(w_j[1:min(10, length(w_j))])
    }
    Z_tmp[j] <- sample(possible_qs, 1, prob = w_j)
  }
  
  return(Z_tmp)
}