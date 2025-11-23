# Assumptions:
# - n_1 ≤ n_2, i.e. datafile D_1 has less datafiles than D_2
# - Datafile D_2 is already very exhausted, if a record in D_1 is not "similar" to any record in D_2,
# then it it is highly possible that it doesn't link to anyone in D_2.
# - Resulting from previous assumption, every record in D_2 is linked to some record, but may not be
# contained in D_1. Therefore, every y in D_2 follows same marginal distribution.

PreProcess_SimulationData_withSeed_IndepViolation <- function(
  records,
  num_simu,
  n1,
  n2,
  overlap_proportion,
  R_sq,
  seed_proportion,
  mean_x_nonlink,
  mean_x_link = 0,
  beta_0 = 3,
  beta_1 = 3
) {
  # get entity ID
  records$rec.id <- gsub("^rec-", "", records$rec.id)
  records$rec.id <- as.numeric(gsub("-[[:alnum:]]+", "", records$rec.id))

  organized_data <- lapply(
    1:num_simu,
    function(simu_var) {
      if (seed_proportion > overlap_proportion) {
        stop(
          "The `seed_proportion` must be less than or equal to overlap_proportion."
        )
      }
      if (n1 * seed_proportion != round(n1 * seed_proportion)) {
        stop("The number of seed must be an integer.")
      }
      if (n2 * overlap_proportion != round(n2 * overlap_proportion)) {
        stop("The number of overlap records must be an integer.")
      }

      # number of true links already known
      num_seed <- n1 * seed_proportion

      # identify which records will actually belong to each sample
      sample_1 <- sample(0:999, size = n1, replace = FALSE)
      records$samp1 <- records$rec.id %in% sample_1
      common_sample <- sample(
        sample_1,
        size = n1 * overlap_proportion,
        replace = FALSE
      )
      sample_2 <- sample(
        setdiff(0:999, sample_1),
        size = n2 - n1 * overlap_proportion,
        replace = FALSE
      )
      records$samp2 <- records$rec.id %in% c(sample_2, common_sample)

      # drop extra records
      records <- records[records$samp1 | records$file == 2, ]
      records <- records[records$samp2 | records$file == 1, ]

      # don't need these two columns any more
      records$samp1 <- NULL
      records$samp2 <- NULL

      records <- records[order(records$file), ]

      # generate datafile 1 and 2 (including regresion data: x and y)
      file_1 <- records |> filter(file == 1)
      file_2 <- records |> filter(file == 2)
      sigma <- sqrt(beta_1^2 / R_sq - beta_1^2)
      file_1 <- file_1 |>
        mutate(
          x = sapply(
            1:n1,
            function(row_id) {
              if (file_1$rec.id[row_id] %in% file_2$rec.id) {
                return(rnorm(1, mean = 0, sd = 1))
              } else {
                return(rnorm(1, mean = 5, sd = 1))
              }
            }
          )
        )
      file_2 <- file_2 |>
        mutate(
          y = sapply(
            1:n2,
            function(row_id) {
              if (file_2$rec.id[row_id] %in% file_1$rec.id) {
                return(rnorm(
                  1,
                  mean = beta_0 + beta_1 * file_1[file_1$rec.id == file_2$rec.id[row_id], 'x'],
                  sd = sigma
                ))
              } else {
                return(rnorm(
                  1,
                  mean = beta_0,
                  sd = sqrt(sigma^2 + beta_1^2)
                ))
              }
            }
          )
        )

      file_1.id <- file_1$rec.id
      file_2.id <- file_2$rec.id
      file_1 <- file_1 %>% select(rec.id, gname, fname, age, occup, x)
      file_2 <- file_2 %>% select(rec.id, gname, fname, age, occup, y)

      # known links
      known_links <- sample(common_sample, size = num_seed, replace = FALSE)
      Z_known <- rep(-1, n2)
      Z_known[match(known_links, file_2.id)] <- match(known_links, file_1.id) # if element == -1, then it is not a link; otherwise it is a known link

      n <- dim(records)[1]
      cellinds <- expand.grid(1:n1, 1:n2)

      # computing agreement levels for binary comparisons
      AgrLevAge <- AgrLevBinComp(records$age, cellinds)
      AgrLevOccup <- AgrLevBinComp(records$occup, cellinds)

      # computing agreement levels for comparisons based on Levenshtein distance
      AgrLevGname <- AgrLevLevenshtein(records$gname, cellinds)
      AgrLevFname <- AgrLevLevenshtein(records$fname, cellinds)

      gamma <- data.frame(
        record_1 = cellinds[, 1],
        record_2 = cellinds[, 2],
        f1 = AgrLevGname,
        f2 = AgrLevFname,
        f3 = AgrLevAge,
        f4 = AgrLevOccup
      )

      res <- list(
        file_1 = file_1,
        file_2 = file_2,
        gamma = gamma,
        Z_known = Z_known
      )
      return(res)
    }
  )
}
