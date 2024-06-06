#' Power-law fit with normally distributed residuals - independent power law
#'
#' @description
#' Power-law fit of (x, y) data. For each value of x, the mean
#' of y is described by
#'
#'   mu = y0*x^beta
#'
#' where y0 is the power law multiplier, and beta the power law exponent.
#'
#' The standard deviations are normally distributed, where the standard
#' deviation scales with the mean, i.e.
#'
#'   sigma = sigma0*x^delta
#'
#' where sigma0 is the power-law multiplier for the standard deviation, and
#' delta the power-law exponent.
#'
#' The optimal fitting parameters y0, beta, sigma0 and delta are found by
#' maximising the (weighted) loglikelihood.
#'
#' @inheritParams power_gamma_fit
#' @param method Only used for making an initial guess for the standard
#'    deviation power-law exponent.
#'    Choose `newton` for gradient descent solving, using the
#'    `rootSolve::multiroot()` function, or `bisection` for bisection root
#'    solving algorithm using the `stats::uniroot()` function.
#' @param range two-value array to add to best guess, to define initial
#'   interval for bisection algorithm
#' @return a list containing the fields
#'   * `loglikelihood`: the log-likelihood score of the fit
#'   * `multiplier`: fitted power-law multiplier
#'   * `exponent`: fitted power coefficient
#'   * `sd_multiplier`: fitted standard deviation power-law multiplier
#'   * `sd_exponent`: fitted standard deviation power-law exponent
#'   * `ks`: komogorov-smirnov distance of best fit
#' @export
#' @examples
#' # parameters
#' y0 <- 20
#' beta <- -0.25
#' sigma0 <- 2.5
#' delta <- -1
#' x <- seq(1, 8, l = 10)
#' mu <- y0*x^beta
#' sigma <- sigma0*x^delta
#' y <- abs(stats::rnorm(length(x), mu, sigma))
#' weights <- stats::runif(length(x), 0.9, 1.1)
#'
#' # fit
#' ft <- power_normal_sdpower_fit(x, y, weights = weights)
#' ft
#' xp <- seq(min(x), max(x), l = 251)
#' yp <- ft$multiplier*xp^ft$exponent
#'
#' # plot data and best fit
#' plot(x, y)
#' lines(xp, yp, col = "red")
#'
power_normal_sdpower_fit <- function(
    x,
    y,
    weights = rep(1, length(x)),
    method = "bisection",
    range = c(-1, 1),
    start = NULL
) {
  # initial guess
  if (is.null(start)) {
    start <- power_normal_sdpower_initialguess(
      x, y, weights = weights,
      method = method, range = range
    )
  }
  # solve
  par <- rootSolve::multiroot(
    power_normal_sdpower_root,
    start,
    jacfunc = power_normal_sdpower_root_jacobian,
    x = x,
    y = y,
    weights = weights
  )$root
  beta <- par[1]
  delta <- par[2]
  # multiplier <y0> and standard deviation <sigma>
  c1 <- sum(weights)
  c3 <- sum(weights*x^(2*beta - 2*delta))
  c6 <- sum(weights*x^(beta - 2*delta)*y)
  c9 <- sum(weights*x^(-2*delta)*y^2)
  y0 <- c6/c3
  sigma0 <- sqrt(c9/c1 - c6^2/(c1*c3))
  # loglikelihood
  logL <- power_normal_sdpower_loglikelihood(
    c(y0, beta, sigma0, delta),
    x, y, weights = weights
  )
  # kolmogorov-smirnov distance
  ks <- stats::ks.test((y/(x^beta) - y0)/(x^delta), "pnorm", 0, sigma0)
  # return list
  list(
    loglikelihood = logL,
    multiplier = y0,
    exponent = beta,
    sd_multiplier = sigma0,
    sd_exponent = delta,
    ks = as.numeric(ks$statistic)
  )
}


#' Generate initial guess for `power_normal_sdpower_fit()`
#'
#' @description
#' Generate an intiial guess for mean and standard deviation power-law
#' exponents, to be used in function `power_normal_sdpower_fit()`.
#'
#' The mean power-law exponent (beta) is estimated using linear regression on
#' log-transformed x and y data.
#'
#' The standard deviation power-law exponent (delta) is subsequently fitted
#' using logliklihood root solving, assuming <beta> as fixed.
#'
#' @inheritParams power_normal_sdpower_fit
#' @return estimate for power-law exponents for the mean and standard deviation
#'
power_normal_sdpower_initialguess <- function(
    x,
    y,
    weights = rep(1, length(x)),
    method = "bisection",
    range = c(-1, 1)
) {
  # initial guess for beta - linear regression on log-transformed data
  ft0 <- stats::lm(log(y) ~ log(x), weights = weights)
  beta <- as.numeric(ft0$coef[2])
  # initial guess for delta - use root solving, with assumed beta
  if (method == "bisection") {
    delta <- stats::uniroot(
      power_normal_sdpower_initialguess_root,
      beta + range,
      beta = beta,
      x = x,
      y = y,
      weights = weights,
      extendInt = "downX"
    )$root
  } else if (method == "newton") {
    delta <- rootSolve::multiroot(
      power_normal_sdpower_initialguess_root,
      beta,
      jacfunc = power_normal_sdpower_initialguess_root_jacobian,
      beta = beta,
      x = x,
      y = y,
      weights = weights,
    )$root
  }
  # return starting guess
  c(beta, delta)
}


#' Root solving equation for finding initial guess for sd power-law exponent
#'
#' @description
#' Returns the first partial derivative of the loglikelihood with respect to
#' the standard deviation power-law exponent (delta). This is used to make
#' an initial guess for the function `power_normal_sdpower_fit()`
#'
#' @inheritParams power_normal_sdpower_fit
#' @return first partial derivative of loglikelihood with respect to delta
#'
power_normal_sdpower_initialguess_root <- function(
    delta,
    beta,
    x,
    y,
    weights = rep(1, length(x))
) {
  # coefficients
  c1 <- sum(weights)
  c2 <- sum(weights*log(x))
  c3 <- sum(weights*x^(2*beta - 2*delta))
  c4 <- sum(weights*x^(2*beta - 2*delta)*log(x))
  c6 <- sum(weights*x^(beta - 2*delta)*y)
  c7 <- sum(weights*x^(beta - 2*delta)*y*log(x))
  c9 <- sum(weights*x^(-2*delta)*y^2)
  c10 <- sum(weights*x^(-2*delta)*y^2*log(x))
  # temporary variables
  zeta <- 1/(c3^2*c9 - c3*c6^2)
  # return derivative of loglikelihood with respect to delta
  c1*(c3^2*c10 - 2*c3*c6*c7 + c4*c6^2)*zeta - c2
}


#' Jacobian of root to solve initial guess for delta
#'
#' @description
#' Jacobian of function `power_normal_sdpower_initialguess_root()`
#'
#' @inheritParams power_normal_sdpower_initialguess_root
#' @return derivative of function `power_normal_sdpower_initialguess_root()`
#'   with respect to input argument `delta`
#'
power_normal_sdpower_initialguess_root_jacobian <- function(
    delta,
    beta,
    x,
    y,
    weights = rep(1, length(x))
) {
  # coefficients
  c1 <- sum(weights)
  c2 <- sum(weights*log(x))
  c3 <- sum(weights*x^(2*beta - 2*delta))
  c4 <- sum(weights*x^(2*beta - 2*delta)*log(x))
  c5 <- sum(weights*x^(2*beta - 2*delta)*log(x)^2)
  c6 <- sum(weights*x^(beta - 2*delta)*y)
  c7 <- sum(weights*x^(beta - 2*delta)*y*log(x))
  c8 <- sum(weights*x^(beta - 2*delta)*y*log(x)^2)
  c9 <- sum(weights*x^(-2*delta)*y^2)
  c10 <- sum(weights*x^(-2*delta)*y^2*log(x))
  c11 <- sum(weights*x^(-2*delta)*y^2*log(x)^2)
  # temporary variables
  zeta <- 1/(c3^2*c9 - c3*c6^2)
  dzeta_ddelta <- -(-4*c3*c4*c9 - 2*c3^2*c10 + 2*c4*c6^2 + 4*c3*c6*c7)*zeta^2
  # return 2nd derivative of loglikelihood with respect to delta
  c1*(c3^2*c10 - 2*c3*c6*c7 + c4*c6^2)*dzeta_ddelta +
    c1*(-4*c3*c4*c10 - 2*c3^2*c11 + 4*c3*c7^2 + 4*c3*c6*c8 - 2*c5*c6^2)*zeta
}


#' Calculate power-law normal loglikelihood - separate power law for sd
#'
#' @description
#' Calculate the weighted loglikelihood for a power-law fit with normal
#' residuals. Standard deviation follows its own power-law.
#' Can also be used to calculate the first or second partial
#' derivative with respect to fitting parameters.
#'
#' @inheritParams power_normal_sdpower_fit
#' @param par vector with fitting parameters (power-law multiplier, power-law
#'   exponent, standard deviation power-law multiplier, and standard deviation
#'   power-law exponent)
#' @param deriv order of partial derivative requested
#' @return loglikelihood, or its partial derivatives to order `deriv`
#' @examples
#' # define parameters
#' y0 <- 20
#' beta <- -0.5
#' sigma0 <- 2
#' delta <- -0.75
#' x <- seq(1, 8, l = 101)
#' mu <- y0*x^beta
#' sigma <- sigma0*x^delta
#' y <- rnorm(length(x), mu, sigma)
#' w <- stats::runif(length(x), 0.8, 1.2)
#' par <- c(y0, beta, sigma0, delta)
#'
#' # check loglikelihood
#' mu <- y0*x^beta
#' sigma <- sigma0*x^delta
#' sum(w*stats::dnorm(y, mu, sigma, log = TRUE))
#' power_normal_sdpower_loglikelihood(par, x, y, weights = w)
#'
#' # test first derivative
#' eps <- 1e-6
#' f0 <- power_normal_sdpower_loglikelihood(par, x, y, weights = w, deriv = 0)
#' f1 <- power_normal_sdpower_loglikelihood(par + c(eps, 0, 0, 0), x, y, weights = w, deriv = 0)
#' f2 <- power_normal_sdpower_loglikelihood(par + c(0, eps, 0, 0), x, y, weights = w, deriv = 0)
#' f3 <- power_normal_sdpower_loglikelihood(par + c(0, 0, eps, 0), x, y, weights = w, deriv = 0)
#' f4 <- power_normal_sdpower_loglikelihood(par + c(0, 0, 0, eps), x, y, weights = w, deriv = 0)
#' (c(f1, f2, f3, f4) - f0)/eps
#' power_normal_sdpower_loglikelihood(par, x, y, weights = w, deriv = 1)
#'
#' # test second derivative
#' f0 <- power_normal_sdpower_loglikelihood(par, x, y, weights = w, deriv = 1)
#' f1 <- power_normal_sdpower_loglikelihood(par + c(eps, 0, 0, 0), x, y, weights = w, deriv = 1)
#' f2 <- power_normal_sdpower_loglikelihood(par + c(0, eps, 0, 0), x, y, weights = w, deriv = 1)
#' f3 <- power_normal_sdpower_loglikelihood(par + c(0, 0, eps, 0), x, y, weights = w, deriv = 1)
#' f4 <- power_normal_sdpower_loglikelihood(par + c(0, 0, 0, eps), x, y, weights = w, deriv = 1)
#' (cbind(f1, f2, f3, f4) - f0)/eps
#' power_normal_sdpower_loglikelihood(par, x, y, weights = w, deriv = 2)
#'
power_normal_sdpower_loglikelihood <- function(
    par,
    x,
    y,
    weights = rep(1, length(x)),
    deriv = 0
) {
  # split input
  y0 <- par[1]
  beta <- par[2]
  sigma0 <- par[3]
  delta <- par[4]
  # coefficients
  c1 <- sum(weights)
  c2 <- sum(weights*log(x))
  c3 <- sum(weights*x^(2*beta - 2*delta))
  c4 <- sum(weights*x^(2*beta - 2*delta)*log(x))
  c5 <- sum(weights*x^(2*beta - 2*delta)*log(x)^2)
  c6 <- sum(weights*x^(beta - 2*delta)*y)
  c7 <- sum(weights*x^(beta - 2*delta)*y*log(x))
  c8 <- sum(weights*x^(beta - 2*delta)*y*log(x)^2)
  c9 <- sum(weights*x^(-2*delta)*y^2)
  c10 <- sum(weights*x^(-2*delta)*y^2*log(x))
  c11 <- sum(weights*x^(-2*delta)*y^2*log(x)^2)
  # loglikelihood
  if (deriv == 0) {
    -c1*(log(sigma0) + 0.5*log(2*pi)) - delta*c2 -
      (c9 - 2*y0*c6 + c3*y0^2)/(2*sigma0^2)
  } else if (deriv == 1) {
    c(
      (c6 - y0*c3)/sigma0^2,
      y0*(c7 - c4*y0)/sigma0^2,
      -c1/sigma0 + (c9 - 2*c6*y0 + c3*y0^2)/sigma0^3,
      (c10 - 2*c7*y0 + c4*y0^2)/sigma0^2 - c2
    )
  } else if (deriv == 2) {
    d2logL_dy02 <- -c3/sigma0^2
    d2logL_dy0dbeta <- (c7 - 2*c4*y0)/sigma0^2
    d2logL_dy0dsigma0 <- 2*(c3*y0 - c6)/sigma0^3
    d2logL_dy0ddelta <- 2*(c4*y0 - c7)/sigma0^2
    d2logL_dbeta2 <- y0*(c8 - 2*c5*y0)/sigma0^2
    d2logL_dbetadsigma0 <- 2*y0*(c4*y0 - c7)/sigma0^3
    d2logL_dbetaddelta <- 2*y0*(c5*y0 - c8)/sigma0^2
    d2logL_dsigma02 <- c1/sigma0^2 - 3*(c9 - 2*c6*y0 + c3*y0^2)/sigma0^4
    d2logL_dsigma0ddelta <- 2*(-c10 + 2*c7*y0 - c4*y0^2)/sigma0^3
    d2logL_ddelta2 <- (-2*c11 + 4*c8*y0 - 2*c5*y0^2)/sigma0^2
    matrix(
      c(
        d2logL_dy02, d2logL_dy0dbeta, d2logL_dy0dsigma0, d2logL_dy0ddelta,
        d2logL_dy0dbeta, d2logL_dbeta2, d2logL_dbetadsigma0, d2logL_dbetaddelta,
        d2logL_dy0dsigma0, d2logL_dbetadsigma0, d2logL_dsigma02, d2logL_dsigma0ddelta,
        d2logL_dy0ddelta, d2logL_dbetaddelta, d2logL_dsigma0ddelta, d2logL_ddelta2
      ),
      nrow = 4
    )
  }
}


#' Root solve exponent in normal power-law - seperate power law for sd
#'
#' @description
#' Root solve equation to obtain the exponent in power-law fitting
#' with normal distribution with separate power-law describing the standard
#' deviation
#'
#' @inheritParams power_normal_sdpower_fit
#' @param par vector with mean and standard deviation power-law exponents
#' @return value of root functions to solve
#'
power_normal_sdpower_root <- function(
    par,
    x,
    y,
    weights = rep(1, length(x))
) {
  # unpack parameters
  beta <- par[1]
  delta <- par[2]
  # coefficients
  c1 <- sum(weights)
  c2 <- sum(weights*log(x))
  c3 <- sum(weights*x^(2*beta - 2*delta))
  c4 <- sum(weights*x^(2*beta - 2*delta)*log(x))
  c6 <- sum(weights*x^(beta - 2*delta)*y)
  c7 <- sum(weights*x^(beta - 2*delta)*y*log(x))
  c9 <- sum(weights*x^(-2*delta)*y^2)
  c10 <- sum(weights*x^(-2*delta)*y^2*log(x))
  # temporary variables
  zeta <- 1/(c3^2*c9 - c3*c6^2)
  # roots
  dlogL_dbeta <- c1*c6*(c3*c7 - c4*c6)*zeta
  dlogL_ddelta <- c1*(c3^2*c10 - 2*c3*c6*c7 + c4*c6^2)*zeta - c2
  # return
  c(dlogL_dbeta, dlogL_ddelta)
}


#' Jacobian of function `power_normal_sdpower_root()`
#'
#' @description
#' Returns the derivative of the function `power_normal_sdpower_root()`
#' with respect to input argument `par`
#'
#' @inheritParams power_normal_sdpower_root
#' @return derivative of `power_normal_sdpower_root()` with respect to
#'   elements in `par`
#' @examples
#' # parameters
#' y0 <- 20
#' beta <- -0.25
#' sigma0 <- 2.5
#' delta <- -0.5
#' x <- seq(1, 8, l = 50001)
#' mu <- y0*x^beta
#' sigma <- sigma0*x^delta
#' y <- abs(stats::rnorm(length(x), mu, sigma))
#' w <- stats::runif(length(x), 0.9, 1.1)
#' par <- c(beta, delta)
#'
#' eps <- 1e-6
#' f0 <- power_normal_sdpower_root(par, x, y, weights = w)
#' f1 <- power_normal_sdpower_root(par + c(eps, 0), x, y, weights = w)
#' f2 <- power_normal_sdpower_root(par + c(0, eps), x, y, weights = w)
#' (cbind(f1, f2) - f0)/eps
#' power_normal_sdpower_root_jacobian(par, x, y, weights = w)
#'
power_normal_sdpower_root_jacobian <- function(
    par,
    x,
    y,
    weights = rep(1, length(x))
) {
  # unpack parameters
  beta <- par[1]
  delta <- par[2]
  # coefficients
  c1 <- sum(weights)
  c2 <- sum(weights*log(x))
  c3 <- sum(weights*x^(2*beta - 2*delta))
  c4 <- sum(weights*x^(2*beta - 2*delta)*log(x))
  c5 <- sum(weights*x^(2*beta - 2*delta)*log(x)^2)
  c6 <- sum(weights*x^(beta - 2*delta)*y)
  c7 <- sum(weights*x^(beta - 2*delta)*y*log(x))
  c8 <- sum(weights*x^(beta - 2*delta)*y*log(x)^2)
  c9 <- sum(weights*x^(-2*delta)*y^2)
  c10 <- sum(weights*x^(-2*delta)*y^2*log(x))
  c11 <- sum(weights*x^(-2*delta)*y^2*log(x)^2)
  # temporary variables
  zeta <- 1/(c3^2*c9 - c3*c6^2)
  dzeta_dbeta <- -(4*c3*c4*c9 - 2*c4*c6^2 - 2*c3*c6*c7)*zeta^2
  dzeta_ddelta <- -(-4*c3*c4*c9 - 2*c3^2*c10 + 2*c4*c6^2 + 4*c3*c6*c7)*zeta^2
  # roots
  d2logL_dbeta2 <- c1*(c3*c7 - c4*c6)*(c6*dzeta_dbeta + c7*zeta) +
    c1*c6*(c4*c7 + c3*c8 - 2*c5*c6)*zeta
  d2logL_dbetaddelta <-   c1*(c3*c7 - c4*c6)*(c6*dzeta_ddelta - 2*c7*zeta) +
    2*c1*c6*(c5*c6 - c3*c8)*zeta
  d2logL_ddelta2 <- c1*(c3^2*c10 - 2*c3*c6*c7 + c4*c6^2)*dzeta_ddelta +
    c1*(-4*c3*c4*c10 - 2*c3^2*c11 + 4*c3*c7^2 + 4*c3*c6*c8 - 2*c5*c6^2)*zeta
  # return matrix
  matrix(
    c(d2logL_dbeta2, d2logL_dbetaddelta, d2logL_dbetaddelta, d2logL_ddelta2),
    nrow = 2
  )
}
