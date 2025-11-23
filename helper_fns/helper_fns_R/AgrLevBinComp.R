# function computes agreement levels for binary comparisons, used inside for
AgrLevBinComp <- function(x, cellinds) {
  same <- (x[cellinds[, 1]] == x[cellinds[, 2]])
  AgrLev <- 1 * same
  AgrLev[!same] <- 2
  AgrLev <- AgrLev - 1
  return(AgrLev)
}