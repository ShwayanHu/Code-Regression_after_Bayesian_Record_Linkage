suppressPackageStartupMessages(library(stringdist))
ld2category <- function(ld) {
  if (ld == 0) {
    return(3)
  } else if (ld > 0 & ld <= .25) {
    return(2)
  } else if (ld > .25 & ld <= .5) {
    return(1)
  } else if (ld > .5) {
    return(0)
  }
}

compare_field <- function(record_1, record_2) {
  i <- record_1$idx
  j <- record_2$idx
  dist_fname <- ld2category(stringdist(record_1$fname_c1, record_2$fname_c1, method = "lv") / max(length(record_1$fname_c1), length(record_2$fname_c1)))
  dist_lname <- ld2category(stringdist(record_1$lname_c1, record_2$lname_c1, method = "lv") / max(length(record_1$lname_c1), length(record_2$lname_c1)))
  year <- (record_1$by_proc == record_2$by_proc) + 0
  day <- (record_1$bd_proc == record_2$bd_proc) + 0
  return(c(i, j, dist_fname, dist_lname, year, day))
}
