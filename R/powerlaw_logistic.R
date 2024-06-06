#' Fit power law + gamma using loglikelihood fitting
#'
#' @description
#' Fit a power law distribution plus a logistic gamma distribution to a series
#' of test data (x, y), using the log-likelihood method.
#'
#' The mean logistic distribution is assumed to be equal
#' to the power-law describing the mean. The shape parameter is assumed
#' to scale with the mean. The magnitude of the shape is defined at x = 1.
#'
#' @md
#' @inheritParams powerlaw_fit
#' @param start vector with (optimal) initial guess for power-law multiplier,
#'   power-law exponent and logistic shape parameter (at x = 1).
#'   If not defined, an educated guess is made using the function
#'   `powerlaw_logistic_initialguess()`.
#' @return a list containing the fields
#'   * `loglikelihood`: the log-likelihood score of the fit
#'   * `multiplier`: fitted power-law multiplier
#'   * `exponent`: fitted power coefficient
#'   * `scale`: fitted logistic scale parameter, defined at x = 1
#' @export
#' @examples
#' y0 <- 20
#' beta <- -0.5
#' s0 <- 4
#' x <- seq(1, 6, l = 101)
#' y <- abs(rlogis(length(x), y0*x^beta, s0*x^beta))
#' w <- runif(length(x), 0.8, 1.2)
#'
#' ft <- powerlaw_logistic_fit(x, y, weights = w)
#' ft
#' xp <- seq(min(x), max(x), l = 101)
#' yp <- ft$multiplier*xp^ft$exponent
#'
#' plot(x, y)
#' lines(xp, yp, col = "red")
#'
powerlaw_logistic_fit <- function(
    x,
    y,
    weights = rep(1, length(x)),
    start = NULL,
    method = "bisection",
    range = c(-1, 1)
) {
  # initial guess
  if (is.null(start)) {
    start <- powerlaw_logistic_initialguess(
      x, y, weights = weights,
      method = method, range = range
    )
  }
  # solve
  sol <- rootSolve::multiroot(
    function(par) {
      powerlaw_logistic_loglikelihood(par, x, y, weights = weights, deriv = 1)
    },
    start,
    jacfunc = function(par) {
      powerlaw_logistic_loglikelihood(par, x, y, weights = weights, deriv = 2)
    }
  )
  y0 <- sol$root[1]
  beta <- sol$root[2]
  s0 <- sol$root[3]
  # loglikelihood
  logL <- powerlaw_logistic_loglikelihood(sol$root, x, y, weights = weights)
  # return
  list(
    loglikelihood = logL,
    multiplier = y0,
    exponent = beta,
    scale = s0
  )
}


#' Generate initial guess for fitting in `powerlaw_logistic_fit()`
#'
#' @description
#' Generate an initial guess for logistic fitting parameter (power-law
#' multiplier, power-law exponent, and logistic shape parameter (at x = 1).
#'
#' The power-law function describing the mean is estimated using a power-law
#' regression with normally distributed standard deviations that scale with
#' the mean.
#'
#' The scale parameter is then estimated from the standard deviation using
#' the method of moments, i.e. scale = sd*sqrt(3)/pi
#'
#' @inheritParams powerlaw_logistic_fit
#' @return vector with initial guess for power-law multiplier, power-law
#'   exponent, and logistic shape parameter (at x = 1).
#'
powerlaw_logistic_initialguess <- function(
    x,
    y,
    weights = rep(1, length(x)),
    method = "bisection",
    range = c(-1, 1)
) {
  # guess mean from fit with normal regression
  ft0 <- powerlaw_normal_fit(
    x, y, weights = weights,
    sd_exponent = "scaled",
    method = method, range = range
  )
  # return
  c(
    ft0$multiplier,
    ft0$exponent,
    ft0$sd_multiplier*sqrt(3)/pi
  )
}


#' Calculate power-law logistic loglikelihood
#'
#' @description
#' Calculate the weighted loglikelihood for a power-law fit with logistic
#' residuals. Can also be used to calculate the first or second partial
#' derivative with respect to fitting parameters.
#'
#' @inheritParams powerlaw_logistic_fit
#' @param par vector with fitting parameters (power-law multiplier, power-law
#'   exponent, and logistic shape parameter)
#' @param deriv order of partial derivative requested
#' @return loglikelihood, or its partial derivatives to order `deriv`
#' @examples
#' # parameters and data
#' y0 <- 20
#' beta <- -0.5
#' s0 <- 4
#' x <- seq(1, 8, l = 100)
#' y <- stats::rlogis(length(x), y0*x^beta, s0*x^beta)
#' w <- runif(length(x), 0.8, 1.2)
#'
#' # check loglikelihood function
#' par <- c(y0, beta, s0)
#' powerlaw_logistic_loglikelihood(par, x, y, weights = w, deriv = 0)
#' sum(w*stats::dlogis(y, y0*x^beta, s0*x^beta, log = TRUE))
#'
#' # test first derivative
#' eps <- 1e-6
#' f0 <- powerlaw_logistic_loglikelihood(par, x, y, weights = w, deriv = 0)
#' f1 <- powerlaw_logistic_loglikelihood(par + c(eps, 0, 0), x, y, weights = w, deriv = 0)
#' f2 <- powerlaw_logistic_loglikelihood(par + c(0, eps, 0), x, y, weights = w, deriv = 0)
#' f3 <- powerlaw_logistic_loglikelihood(par + c(0, 0, eps), x, y, weights = w, deriv = 0)
#' (c(f1, f2, f3) - f0)/eps
#' powerlaw_logistic_loglikelihood(par, x, y, weights = w, deriv = 1)
#'
#' # test second derivative
#' f0 <- powerlaw_logistic_loglikelihood(par, x, y, weights = w, deriv = 1)
#' f1 <- powerlaw_logistic_loglikelihood(par + c(eps, 0, 0), x, y, weights = w, deriv = 1)
#' f2 <- powerlaw_logistic_loglikelihood(par + c(0, eps, 0), x, y, weights = w, deriv = 1)
#' f3 <- powerlaw_logistic_loglikelihood(par + c(0, 0, eps), x, y, weights = w, deriv = 1)
#' (cbind(f1, f2, f3) - f0)/eps
#' powerlaw_logistic_loglikelihood(par, x, y, weights = w, deriv = 2)
#'
powerlaw_logistic_loglikelihood <- function(
    par,
    x,
    y,
    weights = rep(1, length(x)),
    deriv = 0
) {
  # split input
  y0 <- par[1]
  beta <- par[2]
  s0 <- par[3]
  # inermediate parameters
  eta <- 1/cosh((y/x^beta - y0)/(2*s0))
  zeta <- tanh((y/x^beta - y0)/(2*s0))
  # weighted loglikelihood
  if (deriv == 0) {
    logp <- -log(4*s0) - beta*log(x) + 2*log(eta)
    sum(weights*logp)
  } else if (deriv == 1) {
    # derivatives of loglikelihood
    dlogpL_dy0 <- zeta/s0
    dlogp_dbeta <- y*log(x)*zeta/(s0*x^beta) - log(x)
    dlogp_ds0 <- (y/x^beta - y0)*zeta/s0^2 - 1/s0
    # return
    c(
      sum(weights*dlogpL_dy0),
      sum(weights*dlogp_dbeta),
      sum(weights*dlogp_ds0)
    )
  } else if (deriv == 2) {
    # derivatives of zeta
    dzeta_dy0 <- -eta^2/(2*s0)
    dzeta_dbeta <- -y*log(x)*eta^2/(2*s0*x^beta)
    dzeta_ds0 <- -(y/x^beta - y0)*eta^2/(2*s0^2)
    # derivatives of probability
    d2p_dy02 <- dzeta_dy0/s0
    d2p_dy0dbeta <- dzeta_dbeta/s0
    d2p_dy0ds0 <- (dzeta_ds0 - zeta/s0)/s0
    d2p_dbeta2 <- y*log(x)*(dzeta_dbeta - log(x)*zeta)/(s0*x^beta)
    d2p_dbetads0 <- y*log(x)*(dzeta_ds0 - zeta/s0)/(s0*x^beta)
    d2p_ds02 <- (y/x^beta - y0)*(dzeta_ds0 - 2*zeta/s0)/s0^2 + 1/s0^2
    # derivatives of loglikelihood
    d2logL_dy02 <- sum(weights*d2p_dy02)
    d2logL_dy0dbeta <- sum(weights*d2p_dy0dbeta)
    d2logL_dy0ds0 <- sum(weights*d2p_dy0ds0)
    d2logL_dbeta2 <- sum(weights*d2p_dbeta2)
    d2logL_dbetads0 <- sum(weights*d2p_dbetads0)
    d2logL_ds02 <- sum(weights*d2p_ds02)
    # return
    matrix(
      c(
        d2logL_dy02, d2logL_dy0dbeta, d2logL_dy0ds0,
        d2logL_dy0dbeta, d2logL_dbeta2, d2logL_dbetads0,
        d2logL_dy0ds0, d2logL_dbetads0, d2logL_ds02
      ),
      nrow = 3
    )
  }
}


#' Calculate Kolmogorov-Smirnov parameters for power-law fit + logistic
#'
#' @description
#' Calculate Kolmogorov-Smirnov parameter for power-law fit with logistic
#' distributed residuals
#'
#' @md
#' @inheritParams powerlaw_normal_fit
#' @param multiplier,exponent multiplier and exponent for power-law fit
#'   describing the mean
#' @param scale scale parameter for standard deviation of residuals, at x = 1
#' @return list with fields
#'   * `ks_distance`: Kolmogorov-Smirnov distance
#' @examples
#' y0 <- 20
#' beta <- -0.5
#' s0 <- 4
#' x <- seq(1, 6, l = 51)
#' y <- abs(rlogis(length(x), y0*x^beta, s0*x^beta))
#' w <- runif(length(x), 0.8, 1.2)
#'
#' ft <- powerlaw_logistic_fit(x, y, weights = w)
#'
#' powerlaw_logistic_ks(x, y, ft$multiplier, ft$exponent, ft$scale)
#'
powerlaw_logistic_ks <- function(
    x,
    y,
    multiplier,
    exponent,
    scale
) {
  # prediction
  C2 <- sort(stats::plogis(
    y,
    multiplier*x^exponent,
    scale*x^exponent
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
