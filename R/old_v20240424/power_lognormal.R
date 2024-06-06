# Fit power-law with log-normal distributions of residuals
# 26/03/2024 - G. J. Meijer


#' Fit power-law function using linear fitting on the log-log-transform
#'
#' @description
#' Fit a power law function between `x` and `y` in the form y = a*x^b.
#'
#' The fit is generated using the `stats::lm()` linear least-squares
#' regression function in R, between the log-transformed x and y data.
#'
#' The power-law function describing the 'mean' data is given by
#'
#'   y = exp(0.5*a*sdlog^2)*x^b
#'
#' where a is the intercept, b is the gradient of the linear fit on the
#' log-log-transform. sdlog is the standard deviation of the residuals of the
#' log-log fit. The correction exp(0.5*sdlog^2) is required since the mean of
#' the log(y) is not equal to the mean of y.
#'
#' @inheritParams power_weibull_fit
#' @param correction if `TRUE`, the required means correction is used, correction
#'   for the fact that the mean of log-transformed value is not the same as the
#'   log-transformed mean. If `FALSE`, this correction is not made and values are
#'   returned as implemented in the existing literature (incorrect approach).
#' @return list with the power-law multiplier (`multiplier`) and exponent
#'   (`exponent`), the standard deviation of the log-transformed residuals
#'   (`sdlog`) and the loglikelihood score (`loglikelihood`)
#' @export
#' @examples
#' # input
#' x <- seq(2, 10, l = 25)
#' y <- 10*x^-0.5
#' y <- y*stats::rweibull(length(x), shape = 4, scale = 1/gamma(1 + 1/4))
#'
#' # fit
#' ft <- power_lognormal_fit(x, y)
#' xp <- seq(min(x), max(x), l = 101)
#' yp <- ft$multiplier*xp^ft$exponent
#'
#' # plot
#' plot(x, y)
#' lines(xp, yp, col = "red")
#'
power_lognormal_fit <- function(
    x,
    y,
    weights = rep(1, length(x)),
	  correction = TRUE
) {
  # log-log fitting
  ft <- stats::lm(log(y) ~ log(x), weights = weights)
  par <- as.vector(ft$coefficients)
  # residuals and variance
  logyp <- stats::predict(ft)
  res <- log(y) - logyp
  sdlog <- sqrt(sum(weights*res^2)/sum(weights))
  # calculate multiplier
  if (correction == TRUE) {
	  multiplier <- exp(par[1] + 0.5*sdlog^2)
  } else {
	  multiplier <- exp(par[1])
  }
  # log-probability (assumes normal distribution of residuals)
  logpi <- stats::dnorm(ft$residuals, mean = 0, sd = sdlog, log = TRUE)
  # return
  list(
    loglikelihood = sum(weights*logpi),
    multiplier = multiplier,
    exponent = par[2],
    sdlog = sdlog
  )
}


#' Covariance matrix of fitting parameters for `power_lognormal_fit()`
#'
#' @description
#' Estimate the variance-covariance matrix for power law fitting parameters, based
#' on Fisher information matrix, for log-normally distributed residuals
#'
#' @inheritParams power_lognormal_fit
#' @param multiplier,exponent power-law multiplier and exponent for the
#'   mean
#' @param sdlog standard deviation for log-transformed residuals
#' @return 2*2 Fisher information matrix (for parameter `multiplier` and
#'   `exponent`)
#' @examples
#' # input
#' x <- seq(2, 10, l = 25)
#' y <- 10*x^-0.5
#' y <- y*stats::rweibull(length(x), shape = 4, scale = 1/gamma(1 + 1/4))
#'
#' # fit
#' ft <- power_lognormal_fit(x, y)
#'
#' # covariance
#' power_lognormal_covariancematrix(x, y, ft$multiplier, ft$exponent, ft$sdlog)
#'
power_lognormal_covariancematrix <- function(
    x,
    y,
    multiplier,
    exponent,
    sdlog,
    weights = rep(1, length(x)),
    correction = TRUE
) {
  # calculate log-normal means: mu
  if (correction == TRUE) {
    mu <- log(multiplier) + exponent*log(x) - 0.5*sdlog^2
  } else {
    mu <- log(multiplier) + exponent*log(x)
  }
  # derivatives of mu with respect to multiplier and exponent
  dmu_da <- 1/multiplier
  dmu_db <- log(x)
  d2mu_da2 <- -1/multiplier^2
  d2mu_dadb <- 0
  d2mu_db2 <- 0
  # derivatives of log-probability with respect to log-mean mu
  dlogp_dmu <- (log(y) - mu)/sdlog^2
  d2logp_dmu2 <- -1/sdlog^2
  # derivatives of log-probability
  dlogp_da2 <- d2logp_dmu2*dmu_da^2 + dlogp_dmu*d2mu_da2
  dlogp_dadb <- d2logp_dmu2*dmu_da*dmu_db + dlogp_dmu*d2mu_dadb
  dlogp_db2 <- d2logp_dmu2*dmu_db^2 + dlogp_dmu*d2mu_db2
  # derivatives of log-likelihood
  d2L_da2 <- sum(weights*dlogp_da2)
  d2L_dadb <- sum(weights*dlogp_dadb)
  d2L_db2 <- sum(weights*dlogp_db2)
  # Fisher information matrix - multiplier & exponent
  fisher <- -matrix(c(d2L_da2, d2L_dadb, d2L_dadb, d2L_db2), nrow = 2)
  # return matrix
  solve(fisher)
}


#' Generate prediction interval for lognormal fit
#'
#' @description
#' Generate prediction interval with specified level of confidence, based on
#' power law fit with lognormally distributed residuals
#'
#' @inheritParams power_lognormal_covariancematrix
#' @param x root diameters at which to predict interval
#' @param level confidence level (fraction)
#' @return dataframe with field for diameter (`x`), average power law strength
#'   (`y`) and the lower and upper bound of the prediction interval (`ymin`,
#'   `ymax`)
#' @examples
#' x <- seq(2, 10, l = 25)
#' power_lognormal_predictioninterval(x, 20, -0.5, 0.2)
#'
power_lognormal_predictioninterval <- function(
    x,
    multiplier,
    exponent,
    sdlog,
    level = 0.95,
    correction = TRUE
) {
  if (correction == TRUE) {
    mu <- -0.5*sdlog^2
  } else {
    mu <- 0
  }
  y <- multiplier*x^exponent
  data.frame(
    x = x,
    y = y,
    ymin = y*stats::qlnorm(0.5 - 0.5*level, mu, sdlog),
    ymax = y*stats::qlnorm(0.5 + 0.5*level, mu, sdlog)
  )
}
