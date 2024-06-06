#' Fit power law + gamma using loglikelihood fitting
#'
#' @description
#' Fit a power law distribution plus a gamma distribution to a series
#' of test data (x, y), using the log-likelihood method.
#'
#' The scale parameter of the gamma distribution is assumed to be equal
#' to the power-law describing the mean. The shape parameter is assumed
#' constant with x
#'
#' @md
#' @inheritParams powerlaw_weibull_fit
#' @param start initial guess for power-log exponent. If not defined, an
#'   initial guess is made using linear regression on log-transformed x and y
#'   values
#' @param method choose `newton` for gradient descent solving, using the
#'   `rootSolve::multiroot()` function, or `bisection` for bisection root
#'   solving algorithm using the `stats::uniroot()` function.
#' @param range two-value array to add to best guess, to define initial
#'   interval for bisection algorithm
#' @return a list containing the fields:
#'   * `loglikelihood`: the log-likelihood score of the fit
#'   * `multiplier`: fitted power-law multiplier
#'   * `exponent`: fitted power coefficient
#'   * `shape`: fitted gamma shape parameter
#' @export
#' @examples
#' # parameters
#' x <- seq(1, 7, l = 101)
#' y0 <- 25
#' beta <- -0.5
#' k <- 10
#' theta <- y0*x^beta/k
#' y <- stats::rgamma(length(x), shape = k, scale = theta)
#'
#' # fit
#' ft <- powerlaw_gamma_fit(x, y)
#' ft
#' xp <- seq(min(x), max(x), l = 251)
#' yp <- ft$multiplier*xp^ft$exponent
#'
#' # plot data and best fit
#' plot(x, y)
#' lines(xp, yp, col = "red")
#'
powerlaw_gamma_fit <- function(
    x,
    y,
    weights = rep(1, length(x)),
    start = NULL,
    method = "bisection",
    range = c(-1, 1)
) {
  # initial guess for beta
  if (is.null(start)) {
    ft1 <- stats::lm(log(y) ~ log(x), weights = weights)
    start <- as.numeric(ft1$coef[2])
  }
  # find beta, using root solving: c1*c5 = c2*c4
  if (method == "newton") {
    beta <- rootSolve::multiroot(
      powerlaw_gamma_root_beta,
      start,
      jacfunc = powerlaw_gamma_root_beta_jacobian,
      x = x,
      y = y,
      weights = weights
    )$root
  } else if (method == "bisection") {
    beta <- stats::uniroot(
      powerlaw_gamma_root_beta,
      start + range,
      extendInt = "downX",
      x = x,
      y = y,
      weights = weights
    )$root
  } else {
    stop("`method` not recognised")
  }
  # initial guess for k
  c1 <- sum(weights)
  c2 <- sum(weights*log(x))
  c3 <- sum(weights*log(y))
  c4 <- sum(weights*y*x^(-beta))
  k0 <- 0.5*c1/(beta*c2 - c3 + c1*log(c4) - c1*log(c1))
  # find shape parameter k, using root solving
  if (method == "newton") {
    k <- rootSolve::multiroot(
      powerlaw_gamma_root_k,
      k0,
      jacfunc = powerlaw_gamma_root_k_jacobian,
      beta = beta,
      x = x,
      y = y,
      weights = weights
    )$root
  } else if (method == "bisection") {
    k0_range <- pmax(1e-2, k0 + range)  # to ensure k > 0
    k <- stats::uniroot(
      powerlaw_gamma_root_k,
      k0_range,
      extendInt = "downX",
      beta = beta,
      x = x,
      y = y,
      weights = weights
    )$root
  } else {
    stop("`method` not recognised")
  }
  # multiplier <y0>
  c1 <- sum(weights)
  c4 <- sum(weights*y*x^(-beta))
  y0 <- c4/c1
  # likelihood calculation
  logL <- powerlaw_gamma_loglikelihood(c(y0, beta, k), x, y, weights = weights)
  # return dataframe
  list(
    loglikelihood = logL,
    multiplier = y0,
    exponent = beta,
    shape = k
  )
}


#' Calculate power-law gamma loglikelihood
#'
#' @description
#' Calculate the weighted loglikelihood for a power-law fit with gamma
#' residuals. Can also be used to calculate the first or second partial
#' derivative with respect to fitting parameters.
#'
#' @inheritParams powerlaw_gamma_fit
#' @param par vector with fitting parameters (power-law multiplier, power-law
#'   exponent, and gamma shape parameter)
#' @param deriv order of partial derivative requested
#' @return loglikelihood, or its partial derivatives to order `deriv`
#' @keywords internal
#'
powerlaw_gamma_loglikelihood <- function(
    par,
    x,
    y,
    weights = rep(1, length(x)),
    deriv = 0
) {
  # split input
  y0 <- par[1]
  beta <- par[2]
  k <- par[3]
  # coefficients
  c1 <- sum(weights)
  c2 <- sum(weights*log(x))
  c3 <- sum(weights*log(y))
  c4 <- sum(weights*x^(-beta)*y)
  c5 <- sum(weights*x^(-beta)*y*log(x))
  c6 <- sum(weights*x^(-beta)*y*log(x)^2)
  # gamma functions
  logg <- lgamma(k)
  p <- digamma(k)
  q <- psigamma(k, deriv = 1)
  # loglikelihood
  if (deriv == 0) {
    c1*(k*log(k) - logg - k*log(y0)) - beta*k*c2 + (k - 1)*c3 - k*c4/y0
  } else if (deriv == 1) {
    c(
      k/y0*(c4/y0 - c1),
      k*(c5/y0 - c2),
      c1*(1 + log(k) - p - log(y0)) - beta*c2 + c3 - c4/y0
    )
  } else if (deriv == 2) {
    d2logL_dy02 <- k/y0^2*(c1 - 2*c4/y0)
    d2logL_dy0dbeta <- -k*c5/y0^2
    d2logL_dy0dk <- 1/y0*(c4/y0 - c1)
    d2logL_dbeta2 <- -k*c6/y0
    d2logL_dbetadk <- c5/y0 - c2
    d2logL_dk2 <- c1*(1/k - q)
    matrix(
      c(
        d2logL_dy02, d2logL_dy0dbeta, d2logL_dy0dk,
        d2logL_dy0dbeta, d2logL_dbeta2, d2logL_dbetadk,
        d2logL_dy0dk, d2logL_dbetadk, d2logL_dk2
      ),
      nrow = 3
    )
  }
}


#' Root solve exponent in gamma power-law
#'
#' @description
#' Root solve equation to obtain the exponent in power-law fitting
#' with gamma-distribution
#'
#' @inheritParams powerlaw_gamma_fit
#' @param beta power-law exponent to solve for
#' @return value of root function to solve
#' @keywords internal
#'
powerlaw_gamma_root_beta <- function(
    beta,
    x,
    y,
    weights = rep(1, length(x))
) {
  # coefficients
  c1 <- sum(weights)
  c2 <- sum(weights*log(x))
  c4 <- sum(weights*y*x^(-beta))
  c5 <- sum(weights*y*x^(-beta)*log(x))
  # derivatives of root
  c1*c5/c4 - c2
}


#' Jacobian of function `powerlaw_gamma_root_beta()`
#'
#' @description
#' Returns the derivative of the function `powerlaw_gamma_root_beta()`
#' with respect to input argument `beta`
#'
#' @inheritParams powerlaw_gamma_root_beta
#' @return derivative of `powerlaw_gamma_root_beta()` with respect to `beta`
#' @keywords internal
#'
powerlaw_gamma_root_beta_jacobian <- function(
    beta,
    x,
    y,
    weights = rep(1, length(x))
) {
  # coefficients
  c1 <- sum(weights)
  c4 <- sum(weights*y*x^(-beta))
  c5 <- sum(weights*y*x^(-beta)*log(x))
  c6 <- sum(weights*y*x^(-beta)*log(x)^2)
  # derivative of root
  c1/c4*(c5^2/c4 - c6)
}


#' Root solve shape in gamma power-law
#'
#' @description
#' Root solve equation to obtain the Gamma distribution shape parameter
#' in power-law fitting with gamma-distribution
#'
#' @inheritParams powerlaw_gamma_fit
#' @param k unknown gamma distribution shape parameter
#' @param beta known power-law exponen
#' @return value of root function to solve
#' @keywords internal
#'
powerlaw_gamma_root_k <- function(
    k,
    beta,
    x,
    y,
    weights = rep(1, length(x))
) {
  # coefficients
  c1 <- sum(weights)
  c2 <- sum(weights*log(x))
  c3 <- sum(weights*log(y))
  c4 <- sum(weights*y*x^(-beta))
  # root
  c1*(log(k) - digamma(k) - log(c4) + log(c1)) - beta*c2 + c3
}


#' Jacobian of function `powerlaw_gamma_root_k()`
#'
#' @description
#' Returns the derivative of the function `powerlaw_gamma_root_k()`
#' with respect to input argument `k`
#'
#' @inheritParams powerlaw_gamma_root_k
#' @return derivative of `powerlaw_gamma_root_beta()` with respect to `k`
#' @keywords internal
#'
powerlaw_gamma_root_k_jacobian <- function(
    k,
    beta,
    x,
    y,
    weights = rep(1, length(x))
) {
  # coefficients
  c1 <- sum(weights)
  # derivative of root
  c1*(1/k - psigamma(k, deriv = 1))
}

