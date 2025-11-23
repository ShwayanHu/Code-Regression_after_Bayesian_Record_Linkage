z2regdata_multilink <- function(Z, file_1, file_2){
  n1 <- nrow(file_1)
  n2 <- nrow(file_2)
  X <- Y <- c()
  for (j in 1:length(Z)) {
    if ((length(Z[[j]])==1) && (Z[[j]] == n1+j)) {next}
    slot_1 <- Z[[j]]
    x_tmp <- file_1$x[slot_1]
    y_tmp <- rep(file_2$y[j], length(slot_1))
    X <- c(X, x_tmp)
    Y <- c(Y, y_tmp)
  }
  if (length(X)==0) {return(NULL)}
  res <- list(
    X = cbind(matrix(rep(1, length(X)), ncol = 1), matrix(X, ncol=1)),
    Y = matrix(Y, ncol = 1)
  )
  return(res)
}
