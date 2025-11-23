update_m_fast <- function(m, Z, gamma, L, alpha) {
  for (f in 1:4) {
    m[[num2field(f)]] <- as.numeric(rdirichlet(
      n = 1,
      alpha = sapply(0:L[f], function(l) {
        gamma_tmp <- gamma[as.numeric(pull(gamma, num2field(f))) == l, ]
        index_i <- gamma_tmp$record_1
        index_j <- gamma_tmp$record_2
        if (length(index_i) == 0) {
          return(0)
        } else {
          return(sum(index_i == Z[index_j]))
        }
      }) + alpha[[num2field(f)]]
    ))
  }
  return(m)
}