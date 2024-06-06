#' Wrapper for Kolmogorov-Smirnov fitting
#'
#' @description
#' Wrapper function around the base R function `ks.test()`
#'
#' @param x data
#' @param distribution name of the cumulative distribution function
#' @param ... additional arguments to forward into the cumulative distribution
#'   function
#' @return list with fields `ks_distance` for the KS-distance, and `ks_p` for
#'   the corresponding p-value
#' @export
#' @examples
#' x <- stats::rnorm(100, mean = 10, sd = 2)
#' kolmogorov_smirnov(x, "pnorm", mean = 10, sd = 3)
#'
kolmogorov_smirnov <- function(x, distribution, ...) {
  ft <- ks.test(x, distribution, ...)
  list(
    ks_distance = as.vector(ft$statistic),
    ks_p = as.vector(ft$p.value)
  )
}


#' Prepare dataframes for Kolmogorov-Smirnov plot
#'
#' @description
#' Prepare three dataframes for Kolmogorov-Smirnov plotting.
#'
#' @md
#' @inheritParams kolmogorov_smirnov
#' @param xp x-positions at which to plot the fit line
#' @return list with three dataframes containing `x` and `y` data
#'   * `data`: cumulative measured data
#'   * `fit`: cumulative fit
#'   * `ks`: line showing location and size of largest KS distance
#' @export
#' @examples
#' x <- stats::rnorm(50, mean = 10, sd = 2)
#' l <- kolmogorov_smirnov_plot(x, "pnorm", mean = 10, sd = 2)
#'
#' plot(l$data$x, l$data$y, "l")
#' lines(l$fit$x, l$fit$y, col = "red")
#' lines(l$ks$x, l$ks$y, col = "blue")
#'
kolmogorov_smirnov_plot <- function(
    x,
    distribution,
    xp = seq(min(x), max(x), l = 251),
    ...
) {
  # number of points
  n <- length(x)
  # cumulative experimental
  df_measured <- data.frame(
    x = rep(sort(x), each = 2),
    y = rep(seq(0, 1, l = n + 1), each = 2)[2:(2*n + 1)]
  )
  # cumulative fit
  fun <- get(distribution)
  df_fit <- data.frame(
    x = xp,
    y = fun(xp, ...)
  )
  # dataframe with ks line
  df_measured$yp <- fun(df_measured$x, ...)
  i <- which.max(abs(df_measured$y - df_measured$yp))
  df_ks <- data.frame(
    x = rep(df_measured$x[i], 2),
    y = c(df_measured$yp[i], df_measured$y[i])
  )
  # return list of dataframes
  list(
    data = df_measured[, c("x", "y")],
    fit = df_fit,
    ks = df_ks
  )
}
