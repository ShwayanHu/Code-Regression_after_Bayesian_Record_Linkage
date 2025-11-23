z2regdata <- function(Z, file_1, file_2){
  n1 <- nrow(file_1)
  n2 <- nrow(file_2)
  slot_2 <- (1:n2)[Z <= n1]
  slot_1 <- Z[slot_2]
  Y <- matrix(file_2$y[slot_2], ncol = 1)
  X <- cbind(matrix(rep(1, length(slot_1)), ncol=1), matrix(file_1$x[slot_1], ncol = 1))
  return(list(X = X, Y = Y))
}
