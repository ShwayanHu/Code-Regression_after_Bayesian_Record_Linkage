z2regdata_with_goodness_and_seed_RealData <- function(Z, file_1, file_2, m, u, gamma, Z_known){
  n1 <- nrow(file_1)
  n2 <- nrow(file_2)
  X1 <- X2 <- X3 <- Y <- G <- c()
  for (j in 1:length(Z)){
    if (Z[j]==n1+j) {next}
    slot_1 <- Z[j]
    X1 <- c(X1, file_1$X1[slot_1])
    X2 <- c(X2, file_1$X2[slot_1])
    X3 <- c(X3, file_1$X3[slot_1])
    Y <- c(Y, file_2$Y[j])
    case <- gamma %>% filter(record_1==slot_1 & record_2==j)
    log_lik_ratio <- sum(
      sapply(1:3, function(f) {
        l <- as.numeric(pull(case, num2field(f)))
        log(m[[num2field(f)]][l+1]) - log(u[[num2field(f)]][l+1])
      })
    )
    G <- c(G, log_lik_ratio)
  }
  if (length(X1)==0 | length(X2)==0) {return(NULL)}
  slot_2 <- (1:n2)[Z <= n1]
  slot_1 <- Z[slot_2]
  Seed.vec <- Z_known[slot_2]!=-1
  
  res <- list(
    X1 = matrix(X1, ncol = 1),
    X2 = matrix(X2, ncol = 1),
    X3 = matrix(X3, ncol = 1),
    Y = matrix(Y, ncol = 1),
    G = G,
    Seed.vec = Seed.vec
  )
  return(res)
}