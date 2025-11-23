suppressPackageStartupMessages(library(RecordLinkage))
suppressPackageStartupMessages(library(tidyverse))

simulate_data <- function(num_true_pair, num_individuals, beta_0 = 3, beta_1 = 0.4, sigma = 1) {
  # browser()
  # Generate simulation data =====
  data(RLdata10000)
  RLdata10000_with_id <- RLdata10000 %>% mutate(id = identity.RLdata10000)
  duplicates <- RLdata10000_with_id %>%
    subset(id %in% identity.RLdata10000[duplicated(identity.RLdata10000)])
  non_duplicates <- RLdata10000_with_id %>%
    subset(!(id %in% identity.RLdata10000[duplicated(identity.RLdata10000)]))

  ## sample true pairs =====
  A_1_pair <- NULL
  A_2_pair <- NULL
  duplicates_id <- sample(unique(duplicates$id), num_true_pair, replace = FALSE)
  for (i in duplicates_id) {
    x_tmp <- rnorm(1, mean = 0, sd = 1)
    y_tmp <- rnorm(1, mean = beta_0 + beta_1*x_tmp, sd = sigma)
    A_1_pair <- rbind(A_1_pair, (duplicates %>% filter(id == i))[1, ] %>% mutate(x = x_tmp))
    A_2_pair <- rbind(A_2_pair, (duplicates %>% filter(id == i))[2, ] %>% mutate(y = y_tmp))
  }
  reg_with_only_true_link <- lm(A_2_pair$y ~ A_1_pair$x)

  ## sample non-duplicates =====
  A_1_nonpair <- NULL
  A_2_nonpair <- NULL
  non_duplicates_id <- sample(non_duplicates$id, 2*num_individuals, replace = FALSE)
  for (i in 1:num_individuals) {
    x_tmp <- rnorm(1, mean = 0, sd = 1)
    A_1_nonpair <- rbind(A_1_nonpair, non_duplicates %>% filter(id == non_duplicates_id[i]) %>% mutate(x = x_tmp))
    
    x_tmp <- rnorm(1, mean = 0, sd = 1)
    y_tmp <- rnorm(1, mean = beta_0 + beta_1 * x_tmp, sd = sigma)
    A_2_nonpair <- rbind(A_2_nonpair, non_duplicates %>% filter(id == non_duplicates_id[i + num_individuals]) %>% mutate(y = y_tmp))
  }
  reg_with_non_link <- lm(A_2_nonpair$y ~ A_1_nonpair$x)
  
  ## combine =====
  A_1 <- rbind(A_1_pair, A_1_nonpair)
  A_2 <- rbind(A_2_pair, A_2_nonpair)
  
  ## shuffle =====
  A_1 <- A_1[sample(nrow(A_1), size = nrow(A_1), replace=FALSE), ]
  A_2 <- A_2[sample(nrow(A_2), size = nrow(A_2), replace=FALSE), ]

  ## data =====
  A_1_id <- A_1$id
  A_2_id <- A_2$id
  A_1 <- A_1 %>% dplyr::select(-id)
  A_2 <- A_2 %>% dplyr::select(-id)

  return(list(
    file_1 = A_1,
    file_2 = A_2,
    id_1 = A_1_id,
    id_2 = A_2_id,
    reg_with_only_true_link = reg_with_only_true_link,
    reg_with_non_link = reg_with_non_link
  ))
}
