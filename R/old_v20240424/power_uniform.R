# Fit power-law with uniform distributions of residuals
# 26/03/2024 - G. J. Meijer


#' Fit power-law assuming uniform distribution of residuals
#'
#' @description
#' Fit the best power-law curve through a series of `x`,`y` data, where
#'
#'   y_predict = multiplier*x^exponent
#'
#' The ratio y/y_predict is assumed to be uniformly distributed with a width
#' (distance between min and max) `width`.
#'
#' @inheritParams power_weibull_fit
#' @return list with fields for the (weighted) loglikelihood of (y/y_fit,
#'   field `loglikelihood`), the power-law parameters (`multiplier` and
#'   `exponent`) and the width of the distribution of y/y_fit (`width`).
#' @export
#' @examples
#' # generate data
#' x <- seq(1, 10, l = 51)
#' multiplier <- 20
#' exponent <- -0.5
#' width <- 0.5
#' y <- multiplier*x^exponent
#' y <- y*stats::runif(length(x), 1 - 0.5*width, 1 + 0.5*width)
#'
#' # fit
#' ft <- power_uniform_fit(x, y)
#' xp <- seq(min(x), max(x), l = 251)
#' yp <- ft$multiplier*xp^ft$exponent
#'
#' # plot
#' plot(x, y)
#' lines(xp, yp, col = "red")
#' lines(xp, yp*(1 - 0.5*ft$width), col = "red")
#' lines(xp, yp*(1 + 0.5*ft$width), col = "red")
#' lines(x[ft$index], y[ft$index], "p", col = "red")
#'
#' # log plot
#' plot(log(x), log(y))
#' lines(log(xp), log(yp), col = "red")
#' lines(log(xp), log(yp*(1 - 0.5*ft$width)), col = "red")
#' lines(log(xp), log(yp*(1 + 0.5*ft$width)), col = "red")
#' lines(log(x[ft$index]), log(y[ft$index]), "p", col = "red")
#'
power_uniform_fit <- function(x, y, weights = rep(1, length(x))) {
  # create a convex hull i log-space (clockwise order)
  lx <- log(x)
  ly <- log(y)
  i <- grDevices::chull(lx, ly)
  # differences in log(x) and log(y) for points on hull
  dlx <- diff(lx[c(i, i[1])])
  dly <- diff(ly[c(i, i[1])])
  # best fit must be parallel to one of the hull vertices. For each vertex,
  # calculate the intercept of the vertex line with the y-axis, and also the
  # intercepts of lines with the same gradient through all points on the hull.
  # The max difference in gradient is positively related to the likelihood.
  dist <- apply(
    cbind(dlx, dly, deparse.level = 0),
    1,
    function(z) {
      tmp <- ly[i] - lx[i]*z[2]/z[1]
      max(tmp, na.rm = TRUE) - min(tmp, na.rm = TRUE)
    }
  )
  # minimum distance between intercepts = best fit (thinnest band)
  imin <- which.min(dist)
  # calculate fit parameters
  tmp <- ly[i] - lx[i]*dly[imin]/dlx[imin]
  ly_min <- min(tmp)
  ly_max <- max(tmp)
  exponent <- dly[imin]/dlx[imin]
  width <- 2*(exp(ly_max) - exp(ly_min))/(exp(ly_max) + exp(ly_min))
  # indices of points on edges
  index1 <- c(i, i[1])[imin + c(0, 1)]
  yrel <- y/(x^exponent)
  if (yrel[index1[1]] > mean(yrel)) {
    index2 <- which.min(yrel)
  } else {
    index2 <- which.max(yrel)
  }
  # return
  list(
    loglikelihood = -sum(weights)*log(width),
    multiplier = 0.5*exp(ly_min) + 0.5*exp(ly_max),
    exponent = exponent,
    width = width
  )
}


#' Covariance matrix of fitting parameters for `power_normal_fit()`
#'
#' @description
#' Estimate the variance-covariance matrix for power law fitting parameters, based
#' on Fisher information matrix, for normally distributed residuals
#'
#' @inheritParams power_uniform_fit
#' @param multiplier,exponent power-law multiplier and power coefficient for the
#'   mean
#' @param width width of the (uniform) range of y/y_predict
#' @return 2*2 matrix (for parameters `multiplier` and `exponent`)
#' @examples
#' # generate data
#' x <- seq(1, 10, l = 51)
#' multiplier <- 20
#' exponent <- -0.5
#' width <- 0.5
#' y <- multiplier*x^exponent
#' y <- y*stats::runif(length(x), 1 - 0.5*width, 1 + 0.5*width)
#'
#' # fit
#' ft <- power_uniform_fit(x, y)
#'
#' # covariance matrix
#' power_uniform_covariancematrix(
#'   x, y,
#'   ft$multiplier, ft$exponent, ft$width, ft$index
#' )
#'
#' ft <- power_weibull_fit(x, y)
#' power_weibull_covariancematrix(x, y, ft$multiplier, ft$exponent, ft$shape)
#'
power_uniform_covariancematrix <- function(
    x,
    y,
    multiplier,
    exponent,
    width,
    index,
    weights = rep(1, length(x))
) {
  ## WRONG !!!
  # find pivot points
  xp <- x[index]
  yp <- y[index]
  # width of log-transformed band: z
  z <- max(log(yp) - exponent*log(xp)) - min(log(yp) - exponent*log(xp))
  mx <- max(yp/xp^exponent)
  dz_da <- 2/mx
  dz_db <- max(abs(log(xp[1:2]) - log(xp[3])))
  d2z_da2 <- -4/mx^2
  d2z_dadb <- 0
  d2z_db2 <- 0
  # derivatives of loglikelihood with respect to log-transformed with
  dL_dz <- 2*exp(z)/(1 - exp(2*z))*sum(weights)
  d2L_dz2 <- 2*(1 + exp(2*z))*exp(z)/((1 - exp(2*z))^2)*sum(weights)
  # combined derivatives
  d2L_da2 <- d2L_dz2*(dz_da)^2 + dL_dz*d2z_da2
  d2L_dadb <- d2L_dz2*dz_da*dz_db + dL_dz*d2z_dadb
  d2L_db2 <- d2L_dz2*(dz_db)^2 + dL_dz*d2z_db2
  # Fisher information matrix - multiplier & power
  fisher <- -matrix(c(d2L_da2, d2L_dadb, d2L_dadb, d2L_db2), nrow = 2)
  # return matrix
  solve(fisher)
}


#' Generate prediction interval for uniform fit
#'
#' @description
#' Generate prediction interval with specified level of confidence, based on
#' power law fit with uniform residuals
#'
#' @inheritParams power_uniform_covariance
#' @param x root diameters at which to predict interval
#' @param level confidence level (fraction)
#' @return dataframe with field for diameter (`x`), average power law strength
#'   (`y`) and the lower and upper bound of the prediction interval (`ymin`,
#'   `ymax`)
#' @examples
#' x <- seq(1, 10, l = 10)
#' power_uniform_predictioninterval(x, 20, -0.5, 0.3)
#'
power_uniform_predictioninterval <- function(
    x,
    multiplier,
    exponent,
    width,
    level = 0.95
) {
  y <- multiplier*x^exponent
  data.frame(
    x = x,
    y = y,
    ymin = y*(1 - level*0.5*width),
    ymax = y*(1 + level*0.5*width)
  )
}
