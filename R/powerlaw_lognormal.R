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
#' @examples
#' # generate data
#' y0 <- 20
#' beta <- -0.5
#' sdlog <- 0.3
#' muL <- log(y0) - 0.5*sdlog^2
#' x <- seq(1, 8, l = 1001)
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
#' # compare fit to linear regression (including means correction)
#' ftl <- stats::lm(log(y) ~ log(x))
#' betaL <- ftl$coef[2]
#' alphaL <- ftl$coef[1]
#' sdL <- sqrt(1/length(x)*sum((log(y) - alphaL - betaL*log(x))^2))
#' c(exp(alphaL + sdL^2/2), betaL)
#' c(ft$multiplier, ft$exponent)
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
  # kolmogorov-smirnov distance
  ks <- powerlaw_lognormal_ks(x, y, y0, beta, sdlog)
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
#' @examples
#' # generate some data
#' y0 <- 20
#' beta <- -0.5
#' sdlog <- 0.2
#' x <- seq(1, 7, l = 51)
#' y <- y0*x^beta*rlnorm(length(x), -0.5*sdlog^2, sdlog)
#' w <- runif(length(x), 0.8, 1.2)
#'
#' # check likelihood calculation
#' par <- c(y0, beta, sdlog)
#' powerlaw_lognormal_loglikelihood(par, x, y, weights = w)
#' sum(w*dlnorm(y, log(y0) + beta*log(x) - 0.5*sdlog^2, sdlog, log = TRUE))
#'
#' # check first derivative
#' eps <- 1e-6
#' f0 <- powerlaw_lognormal_loglikelihood(par, x, y, weights = w, deriv = 0)
#' f1 <- powerlaw_lognormal_loglikelihood(par + c(eps, 0, 0), x, y, weights = w, deriv = 0)
#' f2 <- powerlaw_lognormal_loglikelihood(par + c(0, eps, 0), x, y, weights = w, deriv = 0)
#' f3 <- powerlaw_lognormal_loglikelihood(par + c(0, 0, eps), x, y, weights = w, deriv = 0)
#' (c(f1, f2, f3) - f0)/eps
#' powerlaw_lognormal_loglikelihood(par, x, y, weights = w, deriv = 1)
#'
#' # check second derivative
#' eps <- 1e-6
#' f0 <- powerlaw_lognormal_loglikelihood(par, x, y, weights = w, deriv = 1)
#' f1 <- powerlaw_lognormal_loglikelihood(par + c(eps, 0, 0), x, y, weights = w, deriv = 1)
#' f2 <- powerlaw_lognormal_loglikelihood(par + c(0, eps, 0), x, y, weights = w, deriv = 1)
#' f3 <- powerlaw_lognormal_loglikelihood(par + c(0, 0, eps), x, y, weights = w, deriv = 1)
#' (cbind(f1, f2, f3) - f0)/eps
#' powerlaw_lognormal_loglikelihood(par, x, y, weights = w, deriv = 2)
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


#' Calculate Kolmogorov-Smirnov parameters for power-law fit + lognormal
#'
#' @description
#' Calculate Kolmogorov-Smirnov parameter for power-law fit with lognormal
#' distributed residuals
#'
#' @md
#' @inheritParams powerlaw_normal_fit
#' @param multiplier,exponent multiplier and exponent for power-law fit
#'   describing the mean
#' @param sdlog standard deviation for for power-law fit
#'   describing the log-standard deviation of residuals
#' @return list with fields
#'   * `ks_distance`: Kolmogorov-Smirnov distance
#' @examples
#' y0 <- 20
#' beta <- -0.5
#' sdlog <- 0.3
#' muL <- log(y0) - 0.5*sdlog^2
#' x <- seq(1, 8, l = 101)
#' y <- x^beta*rlnorm(length(x), muL, sdlog)
#'
#' ft <- powerlaw_lognormal_fit(x, y)
#'
#' powerlaw_lognormal_ks(x, y, ft$multiplier, ft$exponent, ft$sdlog)
#'
powerlaw_lognormal_ks <- function(
    x,
    y,
    multiplier,
    exponent,
    sdlog
) {
  # prediction
  C2 <- sort(stats::plnorm(
    y,
    log(multiplier) + exponent*log(x) - 0.5*sdlog^2,
    sdlog
  ))
  # real cumulative
  C1_lower <- seq(length(x))/(1 + length(x))
  C1_upper <- C1_lower + 1/(1 + length(x))
  # distance
  ks_distance <- max(abs(c(C1_lower - C2, C1_upper - C2)))
  # return list (in case of future expansion)
  list(
    ks_distance = ks_distance
  )
}
