#' Fit power-law assuming lognormally distributed residuals
#'
#' @md
#' @description
#' Solve power-law fit with log-normally distributed residuals to a series
#' of (x, y) data, using the loglikelihood method.
#'
#' The power-law is described by:
#'
#'   y0*x^beta
#'
#' where y0 is the power law multiplier and beta the power law exponent.
#'
#' The lognormal distribution is characterised by log-mean parameter muL and
#' log-standard deviation sdlog.
#'
#' The mean (at each value of x) must satisfy:
#'
#'   exp(muL + 0.5*sdlog^2) = y0*x^beta
#'
#' And therefore muL can be described as:
#'
#'   muL = log(y0) + beta*log(x) - 0.5*sdlog^2
#'
#' @inheritParams powerlaw_gamma_fit
#' @return a list containing the fields
#'   * `loglikelihood`: the log-likelihood score of the fit
#'   * `multiplier`: fitted power-law multiplier
#'   * `sdlog`: log-standard deviation
#' @export
#' @examples
#' # generate data
#' y0 <- 20
#' beta <- -0.5
#' sdlog <- 0.3
#' muL <- log(y0) - 0.5*sdlog^2
#' x <- seq(1, 8, l = 101)
#' y <- x^beta*rlnorm(length(x), muL, sdlog)
#'
#' # fit
#' ft <- powerlaw_lognormal_fit(x, y)
#' ft
#' xp <- seq(min(x), max(x), l = 101)
#' yp <- ft$multiplier*xp^ft$exponent
#'
#' # plot
#' plot(x, y)
#' lines(xp, yp, col = "red")
#'
powerlaw_lognormal_fit <- function(
    x,
    y,
    weights = rep(1, length(x))
) {
  # coefficients
  c1 <- sum(weights)
  c2 <- sum(weights*log(x))
  c3 <- sum(weights*log(y))
  c4 <- sum(weights*log(x)^2)
  c5 <- sum(weights*log(x)*log(y))
  c6 <- sum(weights*log(y)^2)
  # solution
  beta <- -(c2*c3 - c1*c5)/(- c2^2 + c1*c4)
  sdlog <- sqrt(c6/c1 - (c1*c5^2 - 2*c2*c3*c5 + c3^2*c4)/(c1*(c1*c4 - c2^2)))
  y0 <- exp((c3 - beta*c2)/c1 + sdlog^2/2)
  # loglikelihood
  logL <- powerlaw_lognormal_loglikelihood(c(y0, beta, sdlog), x, y, weights = weights)
  # return list
  list(
    loglikelihood = logL,
    multiplier = y0,
    exponent = beta,
    sdlog = sdlog
  )
}


#' Calculate power-law lognormal loglikelihood
#'
#' @description
#' Calculate the weighted loglikelihood for a power-law fit with lognormal
#' residuals. Can also be used to calculate the first or second partial
#' derivative with respect to fitting parameters.
#'
#' @inheritParams powerlaw_lognormal_fit
#' @param par vector with fitting parameters (power-law multiplier, power-law
#'   exponent, and log-standard deviation)
#' @param deriv order of partial derivative requested
#' @return loglikelihood, or its partial derivatives to order `deriv`
#' @keywords internal
#'
powerlaw_lognormal_loglikelihood <- function(
    par,
    x,
    y,
    weights = rep(1, length(x)),
    deriv = 0
) {
  # split input
  y0 <- par[1]
  beta <- par[2]
  sdlog <- par[3]
  # coefficients
  c1 <- sum(weights)
  c2 <- sum(weights*log(x))
  c3 <- sum(weights*log(y))
  c4 <- sum(weights*log(x)^2)
  c5 <- sum(weights*log(x)*log(y))
  c6 <- sum(weights*log(y)^2)
  # return result
  if (deriv == 0) {
    (log(y0)/sdlog^2 - 1)*(c3 - beta*c2 - c1*log(y0)/2) -
      c1*(log(sdlog) + log(2*pi)/2 + sdlog^2/8) -
      (beta*c2 + c3)/2 -
      (beta^2*c4 - 2*beta*c5 + c6)/(2*sdlog^2)
  } else if (deriv == 1) {
    c(
      (c3 - beta*c2  - c1*log(y0))/(y0*sdlog^2) + c1/(2*y0),
      0.5*c2 - (beta*c4 - c5 + c2*log(y0))/(sdlog^2),
      (c1*log(y0)^2 - 2*log(y0)*(c3 - beta*c2) + beta^2*c4 - 2*beta*c5 + c6)/sdlog^3 - c1/sdlog - c1*sdlog/4
    )
  } else if (deriv == 2) {
    d2logL_dy02 <- -((c1*sdlog^2)/2 + c1 + c3 - beta*c2 - c1*log(y0))/(sdlog^2*y0^2)
    d2logL_dy0dbeta <- -c2/(sdlog^2*y0)
    d2logL_dy0dsdlog <- (2*(beta*c2 - c3 + c1*log(y0)))/(sdlog^3*y0)
    d2logL_dbeta2 <- -c4/sdlog^2
    d2logL_dbetadsdlog <- (2*(beta*c4 - c5 + c2*log(y0)))/sdlog^3
    d2logL_dsdlog2 <- -3*(c1*log(y0)^2 - 2*log(y0)*(c3 - beta*c2) + beta^2*c4 - 2*beta*c5 + c6)/sdlog^4 + c1/sdlog^2 - c1/4
    matrix(
      c(
        d2logL_dy02, d2logL_dy0dbeta, d2logL_dy0dsdlog,
        d2logL_dy0dbeta, d2logL_dbeta2, d2logL_dbetadsdlog,
        d2logL_dy0dsdlog, d2logL_dbetadsdlog, d2logL_dsdlog2
      ),
      nrow = 3
    )
  }
}


