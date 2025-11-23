getHpi <- function(density, prob) {
  df <- data.frame(x = density$x, y = density$y)
  df_sorted <- df[order(-df$y), ]
  dx <- diff(density$x)[1]
  df_sorted$cum_area <- cumsum(df_sorted$y) * dx
  threshold_index <- which(df_sorted$cum_area >= prob)[1]
  y_cutoff <- df_sorted$y[threshold_index]
  hpi_region <- df[df$y >= y_cutoff, ]
  hpi_lower <- min(hpi_region$x)
  hpi_upper <- max(hpi_region$x)
  return(
    c(hpi_lower, hpi_upper)
  )
}