packages <- c(
  'RecordLinkage',
  'tidyverse',
  'purrr',
  'gtools',
  'mvtnorm',
  'pracma',
  'coda',
  'testit',
  'VaRES',
  'parallel',
  'Rcpp',
  'BRL',
  'psych',
  'tictoc'
)
suppressPackageStartupMessages(lapply(packages, library, character.only = TRUE))
rm(packages)

HOME_DIR <- "~/BayesianRecordLinkage/Section6_Genuine_Data_Application/"
PROJECT_PATH <- "~/BayesianRecordLinkage/"
select <- dplyr::select
purrr::walk(
  list.files(
    paste0(PROJECT_PATH, "helper_fns/helper_fns_R"),
    pattern = "*.R$",
    full.names = TRUE
  ),
  source,
  .GlobalEnv
)

set.seed(20011231)


# Control Variables ------------- ---------------------------------------
Iter_start <- 501
nIter <- 1000
num_seeds <- 50

# Preprocess --------------------------------------------------------------
file_95 <- list.files(paste0(HOME_DIR,"ind95_ascii"), full.names = TRUE)
for (i in seq_along(file_95)) {
  if (i == 1) {
    data.95 <- read.csv(file_95[i]) |> 
      distinct(NQUEST, NORD, .keep_all = TRUE)
    next
  }
  temp <- read.csv(file_95[i])
  if ("NORD" %in% colnames(temp)){
    data.95 <- left_join(data.95, temp, by=c("NQUEST","NORD")) |> 
      distinct(NQUEST, NORD, .keep_all = TRUE)
  } else {
    data.95 <- left_join(data.95, temp, by = c("NQUEST")) |> 
      distinct(NQUEST, .keep_all = TRUE)
  }
}

data.95 <- data.95 %>% 
  filter(NORD == 1)

file_98 <- list.files(paste0(HOME_DIR,"ind98_ascii"), full.names = TRUE)
for (i in seq_along(file_98)) {
  if (i == 1) {
    data.98 <- read.csv(file_98[i]) |> 
      distinct(NQUEST, NORD, .keep_all = TRUE)
    next
  }
  temp <- read.csv(file_98[i])
  if ("NORD" %in% colnames(temp)){
    data.98 <- left_join(data.98, temp, by=c("NQUEST","NORD")) |> 
      distinct(NQUEST, NORD, .keep_all = TRUE)
  } else {
    data.98 <- left_join(data.98, temp, by = c("NQUEST")) |> 
      distinct(NQUEST, .keep_all = TRUE)
  }
}

data.98 <- data.98 %>% 
  filter(NORD == 1)

rm(temp)

# Matching with RLdata10000
RLdata <- RLdata10000 |>
  mutate(id = identity.RLdata10000) |>
  select(id, fname_c1, lname_c1, by, bm, bd) |>
  drop_na()

id_counts <- RLdata %>%
  count(id, name = "n")

unique_ids <- id_counts %>% filter(n == 1) %>% pull(id)
duplicate_ids <- id_counts %>% filter(n > 1) %>% pull(id)

RL_unique <- RLdata %>%
  filter(id %in% unique_ids) %>% 
  arrange(id)

RL_dup <- RLdata %>%
  filter(id %in% duplicate_ids) %>%
  group_by(id) %>%
  mutate(row_num = row_number()) %>%
  ungroup()

RL_dup1 <- RL_dup %>% filter(row_num %% 2 == 1) %>% select(-row_num) %>% arrange(id)
RL_dup2 <- RL_dup %>% filter(row_num %% 2 == 0) %>% select(-row_num) %>% arrange(id)

common_NQUEST <- intersect(data.95$NQUEST, data.98$NQUEST)
overlap_NQUEST <- sample(common_NQUEST, min(1000, length(common_NQUEST)))
n_overlap <- length(overlap_NQUEST)

data.98_NQUEST <- c(
  sample(setdiff(data.98$NQUEST, common_NQUEST), n_overlap),
  overlap_NQUEST
)
data.95_NQUEST <- c(
  sample(setdiff(data.95$NQUEST, common_NQUEST), n_overlap),
  overlap_NQUEST
)

cat("Need:", n_overlap, "ID's\n")
cat("Num. of available ID's in RL_dup:", length(unique(RL_dup1$id)), "\n")
cat("Num. of available ID's in RL_unique:", nrow(RL_unique), "\n")

overlap_ids <- sample(unique(RL_dup1$id), n_overlap)

overlap_rl_for_95 <- overlap_ids  
overlap_rl_for_98 <- overlap_ids  

unique_pool <- sample(RL_unique$id, 2 * n_overlap)
data.95_rl <- c(unique_pool[1:(n_overlap)], overlap_rl_for_95)
data.98_rl <- c(unique_pool[(n_overlap + 1):(2*n_overlap)], overlap_rl_for_98)


cat("\n=== Check ID allocation ===\n")
cat("data.95_rl length:", length(data.95_rl), "\n")
cat("data.98_rl length:", length(data.98_rl), "\n")
cat("Expected num. of overlap ID's:", n_overlap, "\n")
cat("Actural num. of overlap ID's:", length(intersect(data.95_rl, data.98_rl)), "\n")

linkingvar_95 <- data.frame(
  NQUEST = data.95_NQUEST,
  rl_id = data.95_rl
) %>%
  left_join(
    rbind(
      RL_dup1 %>% select(id, fname_c1, lname_c1, by, bm, bd),
      RL_unique %>% select(id, fname_c1, lname_c1, by, bm, bd)
    ),
    by = c('rl_id' = 'id')
  )

linkingvar_98 <- data.frame(
  NQUEST = data.98_NQUEST,
  rl_id = data.98_rl
) %>%
  left_join(
    rbind(
      RL_dup2 %>% select(id, fname_c1, lname_c1, by, bm, bd),
      RL_unique %>% select(id, fname_c1, lname_c1, by, bm, bd)
    ),
    by = c('rl_id' = 'id')
  )


data.95 <- linkingvar_95 %>%
  left_join(
    data.95 %>% 
      mutate(EDU = pmax(STUPCF, STUMCF, na.rm = TRUE)) %>% 
      select(NQUEST, Y, EDU, SEX, SCORTA),
    by = 'NQUEST'
  ) %>%
  drop_na()

data.98 <- linkingvar_98 %>%
  left_join(
    data.98 %>% 
      mutate(EDU = pmax(STUPCF, STUMCF, na.rm = TRUE)) %>% 
      select(NQUEST, Y, SEX, EDU),
    by = 'NQUEST'
  ) %>%
  drop_na()

Y_upper <- quantile(data.98$Y, 0.95); Y_lower <- quantile(data.98$Y, 0.05)
data.98 <- data.98 |> filter(Y<Y_upper & Y>Y_lower)
AR_upper <- quantile(data.95$Y, 0.95); AR_lower <- quantile(data.95$Y, 0.05)
data.95 <- data.95 |> filter(Y<AR_upper & Y>AR_lower)

file_1 <- data.98 %>% 
  mutate(
    Y = scale(Y)[,1],
    f1 = by,
    f2 = bm,
    f3 = bd,
    f4 = as.numeric(SEX == 1),
    f5 = as.numeric(EDU >= 3),
    f6 = fname_c1,
    f7 = lname_c1,
  ) %>%
  select(f1:f7, Y, NQUEST) %>% 
  drop_na()

file_1.id <- file_1$NQUEST
file_1 <- file_1 %>% select(-NQUEST)

file_2 <- data.95 %>% 
  mutate(
    X1 = scale(Y)[,1],
    X2 = as.numeric(EDU >= 3),
    X3 = as.numeric(SEX == 1),
    X4 = scale(SCORTA)[,1],
    f1 = by,
    f2 = bm,
    f3 = bd,
    f4 = as.numeric(SEX == 1),
    f5 = as.numeric(EDU >= 3),
    f6 = fname_c1,
    f7 = lname_c1
  ) %>% 
  select(NQUEST, f1:f7, X1:X4) %>% 
  drop_na()

file_2.id <- file_2$NQUEST
file_2 <- file_2 %>% select(-NQUEST)

if (nrow(file_1)<nrow(file_2)){
  tmp <- file_1; file_1 <- file_2; file_2 <- tmp
  tmp <- file_1.id; file_1.id <- file_2.id; file_2.id <- tmp
}

n1 <- nrow(file_1); n2 <- nrow(file_2)
overlap_ids <- intersect(file_1.id, file_2.id)
known_links <- sample(overlap_ids, num_seeds, replace = FALSE)

Z_known <- rep(-1, n2)
for(seed_id in known_links) {
  f1_pos <- which(file_1.id == seed_id)
  f2_pos <- which(file_2.id == seed_id)
  Z_known[f2_pos] <- f1_pos
}

data <- list(
  file_1 = file_1,
  file_1.id = file_1.id,
  file_2 = file_2,
  file_2.id = file_2.id,
  overlap_ids = overlap_ids,
  known_links = known_links,
  Z_known = Z_known
)
save(data, file = paste0(HOME_DIR, "preprocessed_data_withName.RData"))