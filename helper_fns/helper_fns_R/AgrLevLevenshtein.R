# function computes agreement levels for Levenshtein comparisons, used inside for
AgrLevLevenshtein <- function(
  x,
  cellinds,
  breaks = c(-Inf, .001, .25, .5, Inf)
) {
  x <- as.character(x)
  LevenshteinSim <- 1 - levenshteinSim(x[cellinds[, 1]], x[cellinds[, 2]])
  AgrLev <- cut(
    LevenshteinSim,
    breaks = breaks,
    labels = seq_len(length(breaks) - 1)
  )
  AgrLev <- as.numeric(AgrLev) - 1
  return(AgrLev)
}

Levenshtein_similarity <- function(x, cellinds) {
  x <- as.character(x)
  LevenshteinSim <- levenshteinSim(x[cellinds[, 1]], x[cellinds[, 2]])
  return(LevenshteinSim)
}
