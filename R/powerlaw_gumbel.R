#' Fit power-law assuming Gumbel distributed residuals
#'
#' @description
#' Fit a power-law regression to a series of (x, y) data using weighted
#' loglikelihood optimalisation.
#'
#' The power law curve describes the mean of the gumbel distribution. The shape
#' parameter of the gumbel distribution is assumed to scale with the power-law
#' mean.
#'
#' @inheritParams powerlaw_fit
#' @param gamma Euler-Mascheroni constant
#' @param scale_min minimum value for the scale parameters during making an
#'   initial guess. Value is needed to ensure algorithm doesn't crash
#' @param start (optimal) initial guess for power-law exponent and gumbel shape
#'   parameter. If not defined, an educated guess is made using the function
#'   `powerlaw_gumbel_initialguess()`.
#' @return a list containing the fields
#'   * `loglikelihood`: the log-likelihood score of the fit
#'   * `multiplier`: fitted power-law multiplier
#'   * `exponent`: fitted power-law exponent
#'   * `scale`: fitted gumbel scale parameter, at x = 1
#' @export
#' @examples
#' # generate some data
#' y0 <- 20
#' beta <- -0.5
#' theta0 <- 5
#' gamma <- 0.57721566490153286060651209008240243104215933593992
#' x <- seq(1, 6, l = 51)
#' y <- y0*x^beta*rgumbel(
#'   length(x),
#'   1 - gamma*theta0/y0,
#'   theta0/y0
#' )
#'
#' # fit
#' ft <- powerlaw_gumbel_fit(x, y)
#' ft
#' xp <- seq(min(x), max(x), l = 101)
#' yp <- ft$multiplier*xp^ft$exponent
#'
#' # plot
#' plot(x, y)
#' lines(xp, yp, col = "red")
#'
powerlaw_gumbel_fit <- function(
    x,
    y,
    weights = rep(1, length(x)),
    start = NULL,
    gamma = 0.57721566490153286060651209008240243104215933593992,
    method = "bisection",
    range = c(-1, 1),
    scale_min = 0.1
) {
  # initial guess
  if (is.null(start)) {
    start <- powerlaw_gumbel_initialguess(
      x,
      y,
      weights = weights,
      method = method,
      range = range,
      scale_min = scale_min
    )
  }
  # root solve
  sol <- rootSolve::multiroot(
    powerlaw_gumbel_root,
    start,
    jacfunc = powerlaw_gumbel_root_jacobian,
    x = x,
    y = y,
    weights = weights
  )
  beta <- sol$root[1]
  theta0 <- sol$root[2]
  # get multiplier
  c1 <- sum(weights)
  c6 <- sum(weights*exp(-y/(theta0*x^beta)))
  y0 <- theta0*(log(c1) - log(c6) + gamma)
  # calculate loglikelihood
  par <- c(y0, beta, theta0)
  logL <- powerlaw_gumbel_loglikelihood(par, x, y, weights = weights)
  # return list
  list(
    loglikelihood = logL,
    multiplier = y0,
    exponent = beta,
    scale = theta0
  )
}


#' Calculate power-law gumbel loglikelihood
#'
#' @description
#' Calculate the weighted loglikelihood for a power-law fit with gumbel
#' residuals. Can also be used to calculate the first or second partial
#' derivative with respect to fitting parameters.
#'
#' @inheritParams powerlaw_gumbel_fit
#' @param par vector with fitting parameters (power-law multiplier, power-law
#'   exponent, and gumbel shape parameter at x = 1)
#' @param deriv order of partial derivative requested
#' @return loglikelihood, or its partial derivatives to order `deriv`
#' @keywords internal
#'
powerlaw_gumbel_loglikelihood <- function(
    par,
    x,
    y,
    weights = rep(1, length(x)),
    deriv = 0,
    gamma = 0.57721566490153286060651209008240243104215933593992
) {
  # unpack parameters
  y0 <- par[1]
  beta <- par[2]
  theta0 <- par[3]
  # coefficients
  c1 <- sum(weights)
  c2 <- sum(weights*log(x))
  c3 <- sum(weights*y/(x^beta))
  c4 <- sum(weights*y*log(x)/(x^beta))
  c5 <- sum(weights*y*log(x)^2/(x^beta))
  c6 <- sum(weights*exp(-y/(theta0*x^beta)))
  c7 <- sum(weights*y/(x^beta)*exp(-y/(theta0*x^beta)))
  c8 <- sum(weights*y*log(x)/(x^beta)*exp(-y/(theta0*x^beta)))
  c9 <- sum(weights*y*log(x)^2/(x^beta)*exp(-y/(theta0*x^beta)))
  c10 <- sum(weights*y^2/(x^(2*beta))*exp(-y/(theta0*x^beta)))
  c11 <- sum(weights*y^2*log(x)/(x^(2*beta))*exp(-y/(theta0*x^beta)))
  c12 <- sum(weights*y^2*log(x)^2/(x^(2*beta))*exp(-y/(theta0*x^beta)))
  # loglikelihood
  if (deriv == 0) {
    c1*(y0/theta0 - log(theta0) - gamma) - beta*c2 -
      c3/theta0 - exp(y0/theta0 - gamma)*c6
  } else if (deriv == 1) {
    dlogL_dy0 <- c1/theta0 - c6/theta0*exp(y0/theta0 - gamma)
    dlogL_dbeta <- c4/theta0 - c2 - c8/theta0*exp(y0/theta0 - gamma)
    dlogL_dtheta0 <- (c3 - c1*y0)/theta0^2 -
      c1/theta0 +
      (c6*y0 - c7)/theta0^2*exp(y0/theta0 - gamma)
    c(dlogL_dy0, dlogL_dbeta, dlogL_dtheta0)
  } else if (deriv == 2) {
    d2logL_dy02 <- -c6/theta0^2*exp(y0/theta0 - gamma)
    d2logL_dy0dbeta <- -c8/theta0^2*exp(y0/theta0 - gamma)
    d2logL_dy0dtheta0 <- 1/theta0^2*(c6*(1 + y0/theta0) - c7/theta0)*exp(y0/theta0 - gamma) -
      c1/theta0^2
    d2logL_dbeta2 <- 1/theta0*(c9 - c12/theta0)*exp(y0/theta0 - gamma) -
      c5/theta0
    d2logL_dbetadtheta0 <- -c4/theta0^2 +
      1/theta0^2*(c8*(1 + y0/theta0) - c11/theta0)*exp(y0/theta0 - gamma)
    d2logL_dtheta02 <- c1/theta0^2 -
      2*(c3 - c1*y0)/theta0^3 -
      2*(c6*y0 - c7)/theta0^3*exp(y0/theta0 - gamma) +
      (2*c7*y0 - c10 - c6*y0^2)/theta0^4*exp(y0/theta0 - gamma)
    matrix(
      c(
        d2logL_dy02, d2logL_dy0dbeta, d2logL_dy0dtheta0,
        d2logL_dy0dbeta, d2logL_dbeta2, d2logL_dbetadtheta0,
        d2logL_dy0dtheta0, d2logL_dbetadtheta0, d2logL_dtheta02
      ),
      nrow = 3
    )
  }
}


#' Generate initial guess for `powerlaw_gumbel_fit()`
#'
#' @description
#' Generate an initial guess for power-law exponent and Gumbel scale parameter
#' at x = 1, to be used in function `powerlaw_gumbel_fit()`.
#'
#' The power-law exponent is estimated using linear regression on
#' log-transformed x and y data. The Gumbel shape parameter is subsequently
#' estimated from Gumbel fitting of the scaled data (y/x^beta)
#'
#' @inheritParams powerlaw_gumbel_fit
#' @return estimate for power-law exponent and Gumbel scale parameter at x = 1
#' @keywords internal
#'
powerlaw_gumbel_initialguess <- function(
    x,
    y,
    weights = rep(1, length(x)),
    method = "bisection",
    range = c(-1, 1),
    scale_min = 0.01
) {
  # initial guess for beta - linear regression on log data
  ft0 <- stats::lm(log(y) ~ log(x), weights = weights)
  beta <- as.numeric(ft0$coef[2])
  # get scale parameter based gumbel fitting of y/x^beta
  ft1 <- gumbel_fit(
    y/x^beta,
    method = method,
    range = range,
    scale_min = scale_min,
    weights = weights
  )
  # return
  c(beta, ft1$scale)
}


#' Root to solve for power-law Gumbel fitting
#'
#' @description
#' Root equation to solve for power-law Gumbel fitting
#'
#' @inheritParams powerlaw_gumbel_fit
#' @param par fitting parameter (power-law exponent, and gumbel shape parameter
#'   at x = 1)
#' @return two-parameter vector with roots
#' @keywords internal
#'
powerlaw_gumbel_root <- function(
    par,
    x,
    y,
    weights = rep(1, length(x))
) {
  # unpack parameters
  beta <- par[1]
  theta0 <- par[2]
  # coefficients
  c1 <- sum(weights)
  c2 <- sum(weights*log(x))
  c3 <- sum(weights*y/(x^beta))
  c4 <- sum(weights*y*log(x)/(x^beta))
  c6 <- sum(weights*exp(-y/(theta0*x^beta)))
  c7 <- sum(weights*y/(x^beta)*exp(-y/(theta0*x^beta)))
  c8 <- sum(weights*y*log(x)/(x^beta)*exp(-y/(theta0*x^beta)))
  # roots
  dlogL_dbeta <- c4/theta0 - c1*c8/(theta0*c6) - c2
  dlogL_dtheta0 <- c3/theta0^2 - c1*c7/(theta0^2*c6) - c1/theta0
  # return
  c(dlogL_dbeta, dlogL_dtheta0)
}


#' Jacobian of root to solve for power-law Gumbel fitting
#'
#' @description
#' Jacobian of root equation to solve for power-law Gumbel fitting.
#'
#' @inheritParams powerlaw_gumbel_root
#' @return derivative of function `powerlaw_gumbel_root()` with respect to input
#'   argument `par`
#' @keywords internal
#'
powerlaw_gumbel_root_jacobian <- function(
    par,
    x,
    y,
    weights = rep(1, length(x))
) {
  # unpack parameters
  beta <- par[1]
  theta0 <- par[2]
  # coefficients
  c1 <- sum(weights)
  c3 <- sum(weights*y/(x^beta))
  c4 <- sum(weights*y*log(x)/(x^beta))
  c5 <- sum(weights*y*log(x)^2/(x^beta))
  c6 <- sum(weights*exp(-y/(theta0*x^beta)))
  c7 <- sum(weights*y/(x^beta)*exp(-y/(theta0*x^beta)))
  c8 <- sum(weights*y*log(x)/(x^beta)*exp(-y/(theta0*x^beta)))
  c9 <- sum(weights*y*log(x)^2/(x^beta)*exp(-y/(theta0*x^beta)))
  c10 <- sum(weights*y^2/(x^(2*beta))*exp(-y/(theta0*x^beta)))
  c11 <- sum(weights*y^2*log(x)/(x^(2*beta))*exp(-y/(theta0*x^beta)))
  c12 <- sum(weights*y^2*log(x)^2/(x^(2*beta))*exp(-y/(theta0*x^beta)))
  # roots
  d2logL_dbeta2 <- -c5/theta0 - c1*c12/(theta0^2*c6) + c1*c9/(theta0*c6) +
    c1*c8^2/(theta0^2*c6^2)
  d2logL_dbetadtheta0 <- -c1*c11/(theta0^3*c6) + c1*c7*c8/(theta0^3*c6^2) -
    c4/theta0^2 + c1*c8/(theta0^2*c6)
  d2logL_dtheta02 <- -2*c3/theta0^3 + 2*c1*c7/(theta0^3*c6) -
    c1*c10/(theta0^4*c6) + c1*c7^2/(theta0^4*c6^2) + c1/theta0^2
  # return
  matrix(
    c(d2logL_dbeta2, d2logL_dbetadtheta0, d2logL_dbetadtheta0, d2logL_dtheta02),
    nrow = 2
  )
}

