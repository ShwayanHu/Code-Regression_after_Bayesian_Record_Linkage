# Preprocessing Sadinle's data
PreProcess_SimulationData <- function(records, num_simu, n1, n2, overlap_proportion, R_sq, beta_0 = 3, beta_1 = 3) {
  # get entity ID
  records$rec.id <- gsub("^rec-", "", records$rec.id)
  records$rec.id <- as.numeric(gsub("-[[:alnum:]]+", "", records$rec.id))
  
  organized_data <- lapply(
    1:num_simu,
    function(simu_var) {
      # identify which records will actually belong to each sample
      sample_1 <- sample(0:999, size = n1, replace=FALSE)
      records$samp1 <- records$rec.id %in% sample_1
      common_sample <- sample(sample_1, size = n1 * overlap_proportion, replace=FALSE)
      sample_2 <- sample(setdiff(0:999, sample_1), size = n2 - n1*overlap_proportion, replace=FALSE)
      records$samp2 <- records$rec.id %in% c(sample_2, common_sample)
      
      # drop extra records
      records <- records[records$samp1 | records$file == 2, ]
      records <- records[records$samp2 | records$file == 1, ]
      
      # don't need these two columns any more
      records$samp1 <- NULL
      records$samp2 <- NULL
      
      records <- records[order(records$file), ]
      
      # add regression data
      sigma <- sqrt(beta_1^2/R_sq - beta_1^2)
      file_1 <- records %>% filter(file == 1) %>% mutate(x = rnorm(n1, mean = 0, sd = 1))
      file_2 <- records %>% filter(file == 2)
      file_2 <- file_2 %>% mutate(y = sapply(1:n2, function(row_id) {
        id <- file_2[row_id, 'rec.id']
        if (id %in% file_1$rec.id) {
          return(rnorm(1, mean = beta_0 + beta_1 * file_1[file_1$rec.id == id, 'x'], sd = sigma))
        } else{
          return(rnorm(1, mean = beta_0, sd = sqrt(sigma^2 + beta_1 ^ 2)))
        }
      }))
      
      file_1.id <- file_1$rec.id
      file_2.id <- file_2$rec.id
      file_1 <- file_1 %>% select(rec.id, gname, fname, age, occup, x)
      file_2 <- file_2 %>% select(rec.id, gname, fname, age, occup, y)
      
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
        gamma = gamma
      )
      return(res)
    }
  )
}
