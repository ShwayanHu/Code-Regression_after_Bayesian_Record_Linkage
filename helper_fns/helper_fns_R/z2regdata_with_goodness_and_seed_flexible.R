z2regdata_with_goodness_and_seed_flelxible <- function(Z, file_1, file_2, m, u, gamma, Z_known){
  n1 <- nrow(file_1)
  n2 <- nrow(file_2)
  file_1_rv <- file_1 %>% select(starts_with("Y") | starts_with("X"))
  file_2_rv <- file_2 %>% select(starts_with("Y") | starts_with("X"))
  
  G <- file_1_pair <- file_2_pair <- NULL
  
  for (j in 1:length(Z)){
    if (Z[j]==n1+j) {next}
    slot_1 <- Z[j]
    file_1_pair <- rbind(file_1_pair, file_1_rv[slot_1,])
    file_2_pair <- rbind(file_2_pair, file_2_rv[j,])
    case <- gamma %>% filter(record_1==slot_1 & record_2==j)
    if (nrow(case)==0){browser()}
    log_lik_ratio <- sum(
      sapply(1:(file_1 |> select(starts_with("f")) |> ncol()), function(f) {
        l <- as.numeric(pull(case, num2field(f)))
        log(m[[num2field(f)]][l+1]) - log(u[[num2field(f)]][l+1])
      })
    )
    G <- rbind(G, log_lik_ratio)
  }
  # stopifnot(nrow(X)==nrow(Y) & nrow(X)==nrow(G))
  slot_2 <- (1:n2)[Z <= n1]
  slot_1 <- Z[slot_2]
  Seed.vec <- Z_known[slot_2]!=-1
  X_and_Y <- cbind(
    X0 = rep(1,nrow(file_1_pair)),
    file_1_pair %>%
      as.data.frame() %>%
      `colnames<-`(colnames(file_1_rv)), 
    file_2_pair %>%
      as.data.frame() %>%
      `colnames<-`(colnames(file_2_rv))
  )
  # browser()
  res <- list(
    X = X_and_Y %>% select(starts_with("X")),
    Y = X_and_Y %>% select(starts_with("Y")),
    G = G,
    Seed.vec = Seed.vec
  )
  return(res)
}

