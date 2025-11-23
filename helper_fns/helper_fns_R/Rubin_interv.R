Rubin_interv <- function(estimates, variances, alpha){
  valid_slot <- !is.na(estimates) & !is.na(variances) & variances>0
  estimates <- estimates[valid_slot]
  variances <- variances[valid_slot]
  stopifnot(length(estimates)==length(variances))
  m <- length(estimates)
  Q.bar <- mean(estimates)
  U.bar <- mean(variances)
  B.m <- sum((estimates-Q.bar)^2)/(m-1)
  T.m <- U.bar + (1+1/m)*B.m
  r.m <- (1+1/m)*B.m/U.bar
  v <- (m-1)*(1+1/r.m)^2
  upperbound <- Q.bar + qt(p=1-alpha/2, df = v)*sqrt(T.m)
  lowerbound <- Q.bar - qt(p=1-alpha/2, df = v)*sqrt(T.m)
  return(c(lowerbound, upperbound))
}
