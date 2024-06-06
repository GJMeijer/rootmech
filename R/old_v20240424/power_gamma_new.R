#' Fit power law + gamma using loglikelihood fitting
#'
#' @description
#' Fit a power law distribution plus a gamma distribution to a series
#' of test data (x, y), using the log-likelihood method. The power law
#' describing the mean value of y for every value of x is given by:
#'
#'   y_fit = multiplier*x^exponent
#'
#' and the ratio y/y_fit is assumed to follow a gamma distribution with
#' a mean value of 1.
#'
#' @param x measured x-values (e.g. root diameters). These values are assumed
#'   to be already normalised by a reference value to ensure a unitless
#'   parameter.
#' @param y measured y-values (e.g. root tensile strength)
#' @param weights weighting for each measurement. Default = 1, but a strong
#'   case can be made for weighting with `weights = x^2` because of the
#'   large effect of thick roots on root-reinforcement
#' @param method choose `newton` for gradient descent solving, using the
#'   `rootSolve::multiroot()` function, or `bisection` for bisection root
#'   solving algorithm using the `stats::uniroot()` function.
#' @param range two-value array to add to best guess, to define initial
#'   interval for bisection algorithm
#' @return a list containing the fields
#'   * `loglikelihood`: the log-likelihood score of the fit
#'   * `multiplier`: fit power-law multiplier
#'   * `exponent`: fit power coefficient
#'   * `shape`: fit Gamma shape parameter
#' @export
#' @examples
#' # parameters
#' x <- seq(1, 7, l = 101)
#' y0 <- 25
#' beta <- -0.5
#' k <- 4
#' y <- y0*x^beta*stats::rgamma(length(x), shape = k, scale = 1/k)
#'
#' # fit
#' ft <- power_gamma_fit(x, y)
#' xp <- seq(min(x), max(x), l = 251)
#' yp <- ft$multiplier*xp^ft$exponent
#'
#' # plot data and best fit
#' plot(x, y)
#' lines(xp, yp, col = "red")
#'
power_gamma_fit <- function(
    x,
    y,
    weights = rep(1, length(x)),
    method = "newton",
    range = c(-1, 1)
) {
  # initial guess for beta
  ft1 <- stats::lm(log(y) ~ log(x), weights = weights)
  beta0 <- as.numeric(ft1$coef[2])
  # find beta, using root solving: c1*c5 = c2*c4
  if (method == "newton") {
    beta <- rootSolve::multiroot(
      power_gamma_root_beta,
      beta0,
      jacfunc = power_gamma_root_beta_jacobian,
      x = x,
      y = y,
      weights = weights
    )$root
  } else if (method == "bisection") {
    beta <- stats::uniroot(
      power_gamma_root_beta,
      beta0 + range,
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
      power_gamma_root_k,
      k0,
      jacfunc = power_gamma_root_k_jacobian,
      beta = beta,
      x = x,
      y = y,
      weights = weights
    )$root
  } else if (method == "bisection") {
    k <- stats::uniroot(
      power_gamma_root_k,
      k0 + range,
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
  yrel <- y/(y0*x^beta)
  logp <- stats::dgamma(yrel, shape = k, scale = 1/k, log = TRUE)
  # return dataframe
  list(
    loglikelihood = sum(weights*logp),
    multiplier = y0,
    exponent = beta,
    shape = k
  )
}


#' @examples
#' # parameters
#' x <- seq(1, 7, l = 101)
#' y0 <- 25
#' beta <- -0.5
#' k <- 4
#' y <- y0*x^beta*stats::rgamma(length(x), shape = k, scale = 1/k)
#' w <- stats::runif(length(x), 0.8, 1.2)
#'
#' # Test first derivative
#' eps <- 1e-6
#' par <- c(y0, beta, k)
#' f0 <- power_gamma_loglikelihood(par, x, y, weights = w)
#' f1 <- power_gamma_loglikelihood(par + c(eps, 0, 0), x, y, weights = w)
#' f2 <- power_gamma_loglikelihood(par + c(0, eps, 0), x, y, weights = w)
#' f3 <- power_gamma_loglikelihood(par + c(0, 0, eps), x, y, weights = w)
#' (c(f1, f2, f3) - f0)/eps
#' power_gamma_loglikelihood(par, x, y, weights = w, deriv = 1)
#'
#' # test second derivative
#' f0 <- power_gamma_loglikelihood(par, x, y, weights = w, deriv = 1)
#' f1 <- power_gamma_loglikelihood(par + c(eps, 0, 0), x, y, weights = w, deriv = 1)
#' f2 <- power_gamma_loglikelihood(par + c(0, eps, 0), x, y, weights = w, deriv = 1)
#' f3 <- power_gamma_loglikelihood(par + c(0, 0, eps), x, y, weights = w, deriv = 1)
#' (cbind(f1, f2, f3) - f0)/eps
#' power_gamma_loglikelihood(par, x, y, weights = w, deriv = 2)
#'
power_gamma_loglikelihood <- function(
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
  if (deriv == 0) {
    c1*(k*log(k) - log(gamma(k)) - k*log(y0)) + k*c3 - beta*k*c2 - k*c4/y0
  } else if (deriv == 1) {
    c(
      -c1*k/y0 + k*c4/y0^2,
      -k*c2 + k*c5/y0,
      c1*(1 + log(k) - log(y0) - digamma(k)) + c3 - beta*c2 - c4/y0
    )
  } else if (deriv == 2) {
    d2logL_dy02 <- c1*k/y0^2 - 2*k*c4/y0^3
    d2logL_dy0dbeta <- -k*c5/y0^2
    d2logL_dy0dk <- -c1/y0 + c4/y0^2
    d2logL_dbeta2 <- -k*c6/y0
    d2logL_dbetadk <- -c2 + c5/y0
    d2logL_dk2 <- c1*(1/k - psigamma(k, deriv = 1))
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
#' @inheritParams power_gamma_fit
#' @param beta power-law exponent to solve for
#' @return value of root function to solve
#'
power_gamma_root_beta <- function(
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


#' Jacobian of function `power_gamma_root_beta()`
#'
#' @description
#' Returns the derivative of the function `power_gamma_root_beta()`
#' with respect to input argument `beta`
#'
#' @inheritParams power_gamma_root_beta
#'
power_gamma_root_beta_jacobian <- function(
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
#' @inheritParams power_gamma_fit
#' @param k unknown gamma distribution shape parameter
#' @param beta known power-law exponen
#' @return value of root function to solve
#'
power_gamma_root_k <- function(
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
  c1*(log(k) + log(c1) - log(c4) - digamma(k)) + c3 - beta*c2
}


#' Jacobian of function `power_gamma_root_k()`
#'
#' @description
#' Returns the derivative of the function `power_gamma_root_k()`
#' with respect to input argument `k`
#'
#' @inheritParams power_gamma_root_beta
#'
power_gamma_root_k_jacobian <- function(
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


#' Covariance matrix of fitting parameters for `power_gamma_fit()`
#'
#' @description
#' Estimate the variance-covariance matrix for power law fitting parameters, based
#' on Fisher information matrix, for gamma-distributed residuals
#'
#' @inheritParams power_gamma_fit
#' @param multiplier,exponent power-law multiplier and power coefficient for the
#'   mean
#' @param shape Gamma shape parameter
#' @return 2*2 matrix (for parameters `multiplier` and `exponent`)
#' @examples
#' # parameters
#' multiplier <- 50
#' exponent <- -0.3
#' shape <- 4
#'
#' # generate data
#' n <- 51
#' x <- seq(2, 10, l = n)
#' y <- multiplier*x^exponent
#' y <- y*stats::rgamma(length(x), shape = shape, scale = 1/shape)
#'
#' # covariance matrix
#' power_gamma_covariancematrix(
#'   multiplier, exponent, shape,
#'   x, y
#' )
#'
power_gamma_covariancematrix <- function(
    y0,
    beta,
    k,
    x,
    y,
    weights = rep(1, length(x))
) {
  # 2nd derivative of loglikelihoodcoefficients
  J <- power_gamma_loglikelihood(
    y0, beta, k,
    x, y, weights = weights,
    deriv = 2
  )
  # return matrix
  solve(-J)
}
