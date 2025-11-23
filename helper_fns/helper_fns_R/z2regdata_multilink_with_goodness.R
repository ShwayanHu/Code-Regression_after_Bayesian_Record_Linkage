z2regdata_multilink_with_goodness <- function(Z, file_1, file_2, m, u, gamma){
  n1 <- nrow(file_1)
  n2 <- nrow(file_2)
  X <- Y <- G <- c()
  for (j in 1:length(Z)) {
    if ((length(Z[[j]])==1) && (Z[[j]] == n1+j)) {next}
    slot_1 <- Z[[j]]
    x_tmp <- file_1$x[slot_1]
    y_tmp <- rep(file_2$y[j], length(slot_1))
    X <- c(X, x_tmp)
    Y <- c(Y, y_tmp)
    for (i in slot_1) {
      case <- gamma %>% filter(record_1 == i & record_2 == j)
      log_lik_ratio <- sum(
        sapply(1:4, function(f) {
          l <- as.numeric(pull(case, num2field(f)))
          log(m[[num2field(f)]][l + 1]) - log(u[[num2field(f)]][l + 1])
        })
      )
      G <- c(G, log_lik_ratio)
    }
  }
  if (length(X)==0) {return(NULL)}
  res <- list(
    X = cbind(matrix(rep(1, length(X)), ncol = 1), matrix(X, ncol=1)),
    Y = matrix(Y, ncol = 1),
    G = G
  )
  return(res)
}
