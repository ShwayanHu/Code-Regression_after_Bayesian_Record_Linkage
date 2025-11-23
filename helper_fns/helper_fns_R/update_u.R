update_u <- function(u, Z, gamma, L, beta) {
  for (f in 1:4) {
    u[[num2field(f)]] = as.numeric(rdirichlet(
      n = 1, 
      alpha = sapply(
        0:L[f],
        function(l){
          slot <- as.numeric(pull(gamma, num2field(f)))==l
          gamma_tmp <- gamma[slot,]
          index_i <- as.numeric(gamma_tmp$record_1)
          index_j <- as.numeric(gamma_tmp$record_2)
          num_cases <- nrow(gamma_tmp)
          if (num_cases==0) {
            return(0)
          } else {
            return(sum(sapply(1:num_cases, function(k) {sum(!(index_i[k] %in% Z[[index_j[k]]]))})))
          }
        }
      ) + beta[[num2field(f)]]
    ))
  }
  return(u)
} 