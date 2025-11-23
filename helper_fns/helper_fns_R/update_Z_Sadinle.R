update_Z_Sadinle <- function(Z_tmp, n1, n2, gamma, m_tmp, u_tmp, beta_pi, alpha_pi) {
  for (j in 1:n2) {
    Z_minus_j <- Z_tmp[-j]
    n12 <- n_12(Z_minus_j, n1)
    # Define possible matches: 1:n1 (links) or n1 + j (non-link)
    possible_qs <- c(1:n1, n1 + j)
    w_j <- sapply(
      possible_qs,
      function(q) {
        if (q <= n1) {  # Link case (q is in file 1)
          if (q %in% Z_minus_j) {
            return(0)  # Ensure one-to-one: q already linked elsewhere
          } else {
            # Calculate likelihood ratio for gamma fields
            case <- gamma %>% dplyr::filter(record_1 == q & record_2 == j)
            log_lik_ratio <- sum(
              sapply(1:4, function(f) {
                l <- as.numeric(pull(case, num2field(f)))
                log(m_tmp[[num2field(f)]][l + 1]) - log(u_tmp[[num2field(f)]][l + 1])
              })
            )
            return(exp(log_lik_ratio))
          }
        } else {  # Non-link case (q = n1 + j)
          # Prior term with correction factor (see notes)
          prior_term <- (n1 - n12) * 
            (n2 - n12 - 1 + beta_pi) / (n12 + alpha_pi)
          # (n12 + 1) / (n2 - n12)
          return(prior_term)
        }
      }
    )
    # Normalize probabilities to avoid underflow/overflow
    w_j <- w_j / sum(w_j)
    
    # Draw new Z_j from possible_qs with weights w_j
    Z_j_tmp <- sample(possible_qs, size = 1, prob = w_j)
    Z_tmp[j] <- Z_j_tmp  # Update Z_tmp in place
  }
  return(Z_tmp)
}
