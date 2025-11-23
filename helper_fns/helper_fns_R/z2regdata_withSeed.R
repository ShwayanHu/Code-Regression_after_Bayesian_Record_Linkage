z2regdata_withSeed <- function(Z, file_1, file_2, Z_known) {
  n1 <- nrow(file_1)
  n2 <- nrow(file_2)
  slot_2 <- (1:n2)[Z <= n1]
  slot_1 <- Z[slot_2]
  slot_known <- Z_known[slot_2]
  X <- matrix(
    c(rep(1, length(slot_1)), file_1$x[slot_1]),
    ncol = 2,
    byrow = FALSE
  )
  Y <- matrix(
    c(file_2$y[slot_2]),
    ncol = 1,
    byrow = FALSE
  )
  Seed.vec <- slot_known!=-1
  res <- list(
    X = X,
    Y = Y,
    Seed.vec = Seed.vec
  )
  return(res)
}