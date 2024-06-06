#' Fit power-law distribution to data
#'
#' @md
#' @description
#' Fit a power probability distribution to a series of observations. The
#' probability density function is defined as:
#'
#'   p(x) = a^x^b
#'
#' The function is only defined between lower bound x = xmin and upper bound
#' x = xmax. The multiplier `a` is defined such so that the total probability
#' across the interval equals 1.
#'
#' The fit is conducted by finding the value of power coefficient `b` that sets
#' the derivative of the loglikelihood function equal to zero. This is done
#' using the bisection solver inbuilt in R (`stats::uniroot()`).
#'
#' @param x vector with observations (must all be finite and >0)
#' @param xmin,xmax min and max limits of the density function
#' @param weights weighting for each observation in `x`
#' @param eps small value, to catch case where b == -1
#' @param range initial search range around best initial guess, used for the
#'   bisection algorithm
#' @export
#' @return a list with fields:
#'   * loglikelihood: loglikelihood of the best fit
#'   * xmin: lower bound of fit domain
#'   * xmax: upper bound of fit domain
#'   * multiplier: fitted power-law multiplier
#'   * exponent: fitted power-law exponent
#' @examples
#' b <- 1.5
#' xmin <- 2
#' xmax <- 7
#' y <- stats::runif(100, 0, 1)
#' x <- (y*(xmax^(b + 1) - xmin^(b + 1)) + xmin^(b + 1))^(1/(b + 1))
#' power_fit(x)
#'
power_fit <- function(
    x,
    xmin = min(x),
    xmax = max(x),
    weights = rep(1, length(x)),
    eps = .Machine$double.eps^0.5,
    range = c(-1, 1)
) {
  # initial guess
  b0 <- 0
  # fit
  b <- stats::uniroot(
    power_root,
    interval = b0 + range,
    x = x,
    xmin = xmin,
    xmax = xmax,
    weights = weights,
    eps = eps,
    extendInt = "downX"
  )$root
  # calculate multiplier
  if (abs(b + 1) < eps) {
    a <- 1/log(xmax/xmin)
  } else {
    a <- (b + 1)/(xmax^(b + 1) - xmin^(b + 1))
  }
  # return list
  list(
    loglikelihood = power_loglikelihood(
      b, x, xmin = xmin, xmax = xmax,
      weights = weights, eps = eps
    ),
    xmin = xmin,
    xmax = xmax,
    multiplier = a,
    exponent = b
  )
}


#' Loglikelihood function for power law probability density
#'
#' @description
#' Returns the loglikelihood score for a power law probability density
#' function, given the power coefficient `b`.
#'
#' @inheritParams power_fit
#' @param b power coefficient
#' @return loglikelihood score
#' @keywords internal
#'
power_loglikelihood <- function(
    b,
    x,
    xmin = min(x),
    xmax = max(x),
    weights = rep(1, length(x)),
    eps = .Machine$double.eps^0.5
) {
  # calculate multiplier
  if (abs(b + 1) < eps) {
    a <- 1/log(xmax/xmin)
  } else {
    a <- (b + 1)/(xmax^(b + 1) - xmin^(b + 1))
  }
  # calculate log-probability
  logpi <- log(a) + b*log(x)
  # return (weighted) likelihood
  sum(weights*logpi)
}


#' Derivative of loglikelihood score for power law probability density
#'
#' @description
#' Returns the derivative of the loglikelihood function for a power law
#' probability distribution function. By solving this function == 0, the
#' power coefficient at which the loglikelihood score is maximised can be
#' found.
#'
#' @inheritParams power_fit
#' @param b value for exponent
#' @return derivative of loglikelihood score
#' @keywords internal
#'
power_root <- function(
    b,
    x,
    xmin = min(x),
    xmax = max(x),
    weights = rep(1, length(x)),
    eps = .Machine$double.eps^0.5
) {
  # calculate multiplier, and derivatives
  if (abs(b + 1) < eps) {
    a <- 1/log(xmax/xmin)
    da_db <- 0.5 - log(xmax)/log(xmax/xmin)
  } else {
    a <- (b + 1)/(xmax^(b + 1) - xmin^(b + 1))
    da_db <- 1/(xmax^(b + 1) - xmin^(b + 1)) - (b + 1)*(log(xmax)*xmax^(b + 1) - log(xmin)*xmin^(b + 1))/(xmax^(b + 1) - xmin^(b + 1))^2
  }
  # calculate log-probabilities
  dlogpi_db <- da_db/a + log(x)
  # return derivative of loglikelihoood
  sum(weights*dlogpi_db)
}
