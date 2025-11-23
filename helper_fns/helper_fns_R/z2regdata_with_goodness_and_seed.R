z2regdata_with_goodness_and_seed <- function(Z, file_1, file_2, m, u, gamma, Z_known, Z_nonlink){
  n1 <- nrow(file_1)
  n2 <- nrow(file_2)
  X <- Y <- G <- c()
  for (j in 1:length(Z)){
    if (Z[j]==n1+j) {next}
    slot_1 <- Z[j]
    X <- c(X, file_1$x[slot_1])
    Y <- c(Y, file_2$y[j])
    case <- gamma %>% filter(record_1==slot_1 & record_2==j)
    log_lik_ratio <- sum(
      sapply(1:4, function(f) {
        l <- as.numeric(pull(case, num2field(f)))
        log(m[[num2field(f)]][l+1]) - log(u[[num2field(f)]][l+1])
      })
    )
    G <- c(G, log_lik_ratio)
  }
  if (length(X)==0) {return(NULL)}
  slot_2 <- (1:n2)[Z <= n1]
  slot_1 <- Z[slot_2]
  slot_known <- Z_known[slot_2]
  Seed.vec <- slot_known!=-1
  
  if (!is.null(Z_nonlink)){
    # calculate the G for non-links
    X.nonlink <- file_1$x[Z_nonlink$from_file_1]
    Y.nonlink <- file_2$y[Z_nonlink$from_file_2]
    G.nonlink <- sapply(
      1:length(Z_nonlink$from_file_1),
      function(i) {
        slot_1 <- Z_nonlink$from_file_1[i]
        slot_2 <- Z_nonlink$from_file_2[i]
        case <- gamma %>% filter(record_1==slot_1 & record_2==slot_2)
        log_lik_ratio <- sum(
          sapply(1:4, function(f) {
            l <- as.numeric(pull(case, num2field(f)))
            log(m[[num2field(f)]][l+1]) - log(u[[num2field(f)]][l+1])
          })
        )
        return(log_lik_ratio)
      }
    )
    
    non_link <- list(
      X = X.nonlink,
      Y = Y.nonlink,
      G = G.nonlink
    )
  } else {
    non_link <- NULL
  }
  
  res <- list(
    X = cbind(matrix(rep(1, length(X)), ncol = 1), matrix(X, ncol=1)),
    Y = matrix(Y, ncol = 1),
    G = G,
    Seed.vec = Seed.vec,
    non_link = non_link
  )
  return(res)
}
