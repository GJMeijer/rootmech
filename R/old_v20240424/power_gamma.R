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
  # multiplier <a> for scale parameter (theta = a*x^beta)
  c1 <- sum(weights)
  c4 <- sum(weights*y*x^(-beta))
  a <- c4/(c1*k)
  # likelihood calculation
  theta <- a*x^beta
  logp <- stats::dgamma(y, shape = k, scale = theta, log = TRUE)
  # return dataframe
  list(
    loglikelihood = sum(weights*logp),
    multiplier = k*a,
    exponent = beta,
    shape = k
  )
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
  -c1*c6/c4 + c1*c5^2/c4^2
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
  c3 + c1*(log(c1) + log(k) - log(c4) - digamma(k)) - beta*c2
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
  c1/k - c1*psigamma(k, 1)
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
#'   x, y,
#'   multiplier, exponent,
#'   shape
#' )
#'
power_gamma_covariancematrix <- function(
    x,
    y,
    multiplier,
    exponent,
    shape,
    weights = rep(1, length(x))
) {
  # coefficients
  c1 <- sum(weights)
  c4 <- sum(weights*y*x^(-exponent))
  c5 <- sum(weights*y*x^(-exponent)*log(x))
  c6 <- sum(weights*y*x^(-exponent)*log(x)^2)
  # second derivatives of log-transformed probabilities
  d2L_da2 <- c1*shape/multiplier^2 - 2*c4/multiplier^3
  d2L_dadb <- -c5/multiplier^2
  d2L_db2 <- -c6/multiplier
  # Fisher information matrix - multiplier & power
  fisher <- -matrix(c(dL2_da2, dL2_dadb, dL2_dadb, dL2_db2), nrow = 2)
  # return matrix
  solve(fisher)
}
