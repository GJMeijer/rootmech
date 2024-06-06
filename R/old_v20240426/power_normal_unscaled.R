#' Power-law fit with homoscedatic normally distributed residuals
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
#' deviation is constant with x:
#'
#'   sigma = sigma
#'
#' The optimal fitting parameters y0, beta and sigma are found by maximising
#' the (weighted) loglikelihood.
#'
#' @inheritParams power_gamma_fit
#' @return a list containing the fields
#'   * `loglikelihood`: the log-likelihood score of the fit
#'   * `multiplier`: fitted power-law multiplier
#'   * `exponent`: fitted power coefficient
#'   * `sd_multiplier`: fitted standard deviation
#'   * `sd_exponent`: assumed standard deviation power-law exponent (beta = 0)
#'   * `ks`: komogorov-smirnov distance of best fit
#' @export
#' @examples
#' # parameters
#' y0 <- 100
#' beta <- -0.8
#' sigma <- 4
#' x <- seq(1, 8, l = 51)
#' mu <- y0*x^beta
#' y <- abs(stats::rnorm(length(x), mu, sigma))
#' weights <- stats::runif(length(x), 0.9, 1.1)
#'
#' # fit
#' ft <- power_normal_unscaled_fit(x, y, weights = weights)
#' ft
#' xp <- seq(min(x), max(x), l = 251)
#' yp <- ft$multiplier*xp^ft$exponent
#'
#' # plot data and best fit
#' plot(x, y)
#' lines(xp, yp, col = "red")
#'
power_normal_unscaled_fit <- function(
    x,
    y,
    weights = rep(1, length(x)),
    method = "bisection",
    range = c(-1, 1),
    start = NULL
) {
  # initial guess for exponent
  if (is.null(start)) {
    ft0 <- stats::lm(log(y) ~ log(x), weights = weights)
    start <- as.numeric(ft0$coef[2])
  }
  # solve
  if (method == "newton") {
    beta <- rootSolve::multiroot(
      power_normal_unscaled_root,
      start,
      jacfunc = power_normal_unscaled_root_jacobian,
      x = x,
      y = y,
      weights = weights
    )$root
  } else if (method == "bisection") {
    beta <- stats::uniroot(
      power_normal_unscaled_root,
      start + range,
      extendInt = "downX",
      x = x,
      y = y,
      weights = weights
    )$root
  } else {
    stop("`method` not recognised")
  }
  # multiplier <y0> and standard deviation <sigma>
  c1 <- sum(weights)
  c2 <- sum(weights*x^(2*beta))
  c5 <- sum(weights*x^beta*y)
  c8 <- sum(weights*y^2)
  y0 <- c5/c2
  sigma <- sqrt(c8/c1 - c5^2/(c1*c2))
  # loglikelihood
  logL <- power_normal_unscaled_loglikelihood(
    c(y0, beta, sigma),
    x, y, weights
  )
  # kolmogorov-smirnov distance
  ks <- stats::ks.test(y - y0*x^beta, "pnorm", 0, sigma)
  # return list
  list(
    loglikelihood = logL,
    multiplier = y0,
    exponent = beta,
    sd_multiplier = sigma,
    sd_exponent = 0,
    ks = as.numeric(ks$statistic)
  )
}


#' Calculate power-law normal (unscaled) loglikelihood
#'
#' @description
#' Calculate the weighted loglikelihood for a power-law fit with normal
#' residuals. Can also be used to calculate the first or second partial
#' derivative with respect to fitting parameters.
#'
#' @inheritParams power_normal_unscaled_fit
#' @param par vector with fitting parameters (power-law multiplier, power-law
#'   exponent, and standard deviation)
#' @param deriv order of partial derivative requested
#' @return loglikelihood, or its partial derivatives to order `deriv`
#' @examples
#' # define parameters
#' y0 <- 20
#' beta <- -0.5
#' x <- seq(1, 8, l = 10001)
#' mu <- y0*x^beta
#' sigma <- 2
#' y <- rnorm(length(x), mu, sigma)
#' w <- stats::runif(length(x), 0.8, 1.2)
#' par <- c(y0, beta, sigma)
#'
#' # check loglikelihood
#' mu <- y0*x^beta
#' sum(w*stats::dnorm(y, mu, sigma, log = TRUE))
#' power_normal_unscaled_loglikelihood(par, x, y, weights = w)
#'
#' # test first derivative
#' eps <- 1e-6
#' f0 <- power_normal_unscaled_loglikelihood(par, x, y, weights = w, deriv = 0)
#' f1 <- power_normal_unscaled_loglikelihood(par + c(eps, 0, 0), x, y, weights = w, deriv = 0)
#' f2 <- power_normal_unscaled_loglikelihood(par + c(0, eps, 0), x, y, weights = w, deriv = 0)
#' f3 <- power_normal_unscaled_loglikelihood(par + c(0, 0, eps), x, y, weights = w, deriv = 0)
#' (c(f1, f2, f3) - f0)/eps
#' power_normal_unscaled_loglikelihood(par, x, y, weights = w, deriv = 1)
#'
#' # test second derivative
#' f0 <- power_normal_scaled_loglikelihood(par, x, y, weights = w, deriv = 1)
#' f1 <- power_normal_scaled_loglikelihood(par + c(eps, 0, 0), x, y, weights = w, deriv = 1)
#' f2 <- power_normal_scaled_loglikelihood(par + c(0, eps, 0), x, y, weights = w, deriv = 1)
#' f3 <- power_normal_scaled_loglikelihood(par + c(0, 0, eps), x, y, weights = w, deriv = 1)
#' (cbind(f1, f2, f3) - f0)/eps
#' power_normal_scaled_loglikelihood(par, x, y, weights = w, deriv = 2)
#'
power_normal_unscaled_loglikelihood <- function(
    par,
    x,
    y,
    weights = rep(1, length(x)),
    deriv = 0
) {
  # split input
  y0 <- par[1]
  beta <- par[2]
  sigma <- par[3]
  # coefficients
  c1 <- sum(weights)
  c2 <- sum(weights*x^(2*beta))
  c3 <- sum(weights*x^(2*beta)*log(x))
  c4 <- sum(weights*x^(2*beta)*log(x)^2)
  c5 <- sum(weights*x^beta*y)
  c6 <- sum(weights*x^beta*y*log(x))
  c7 <- sum(weights*x^beta*y*log(x)^2)
  c8 <- sum(weights*y^2)
  # loglikelihood
  if (deriv == 0) {
    -c1*(log(sigma) + 0.5*log(2*pi)) - (c8 - 2*y0*c5 + c2*y0^2)/(2*sigma^2)
  } else if (deriv == 1) {
    c(
      (c5 - c2*y0)/sigma^2,
      y0*(c6 - c3*y0)/sigma^2,
      (c8 - 2*c5*y0 + c2*y0^2)/sigma^3 - c1/sigma
    )
  } else if (deriv == 2) {
    d2logL_dy02 <- -c2/sigma^2
    d2logL_dy0dbeta <- (c6 - 2*c3*y0)/sigma^2
    d2logL_dy0dsigma <- 2*(c2*y0 - c5)/sigma^3
    d2logL_dbeta2 <- y0*(c7 - 2*c4*y0)/sigma^2
    d2logL_dbetadsigma <- 2*y0*(c3*y0 - c6)/sigma^3
    d2logL_dsigma2 <- c1/sigma^4 - 3*(c8 - 2*c5*y0 + c2*y0^2)/sigma^4
    matrix(
      c(
        d2logL_dy02, d2logL_dy0dbeta, d2logL_dy0dsigma,
        d2logL_dy0dbeta, d2logL_dbeta2, d2logL_dbetadsigma,
        d2logL_dy0dsigma, d2logL_dbetadsigma, d2logL_dsigma2
      ),
      nrow = 3
    )
  }
}


#' Root solve exponent in normal (unscaled) power-law
#'
#' @description
#' Root solve equation to obtain the exponent in power-law fitting
#' with normal (unscaled, i.e. homoscedatic standar deviations)-distribution
#'
#' @inheritParams power_normal_unscaled_fit
#' @param beta power-law exponent to solve for
#' @return value of root function to solve
#'
power_normal_unscaled_root <- function(
    beta,
    x,
    y,
    weights = rep(1, length(x))
) {
  # calculate coefficients
  c1 <- sum(weights)
  c2 <- sum(weights*x^(2*beta))
  c3 <- sum(weights*x^(2*beta)*log(x))
  c4 <- sum(weights*x^(2*beta)*log(x)^2)
  c5 <- sum(weights*x^beta*y)
  c6 <- sum(weights*x^beta*y*log(x))
  c7 <- sum(weights*x^beta*y*log(x)^2)
  c8 <- sum(weights*y^2)
  # root
  c1*c5*(c2*c6 - c3*c5)/(c2*(c2*c8 - c5^2))
}


#' Jacobian of function `power_normal_unscaled_root()`
#'
#' @description
#' Returns the derivative of the function `power_normal_unscaled_root()`
#' with respect to input argument `beta`
#'
#' @inheritParams power_normal_unscaled_root
#' @return derivative of `power_normal_unscaled_root()` with respect to `beta`
#' @examples
#' # define parameters
#' y0 <- 20
#' beta <- -0.5
#' sigma <- 2
#' x <- seq(1, 8, l = 101)
#' mu <- y0*x^beta
#' y <- rnorm(length(x), mu, sigma)
#' w <- stats::runif(length(x), 0.8, 1.2)
#' par <- c(y0, beta, sigma0)
#'
#' eps <- 1e-6
#' f0 <- power_normal_unscaled_root(beta, x, y, weights = w)
#' f1 <- power_normal_unscaled_root(beta + eps, x, y, weights = w)
#' (f1 - f0)/eps
#' power_normal_unscaled_root_jacobian(beta + eps, x, y, weights = w)
#'
power_normal_unscaled_root_jacobian <- function(
    beta,
    x,
    y,
    weights = rep(1, length(x))
) {
  # calculate coefficients
  c1 <- sum(weights)
  c2 <- sum(weights*x^(2*beta))
  c3 <- sum(weights*x^(2*beta)*log(x))
  c4 <- sum(weights*x^(2*beta)*log(x)^2)
  c5 <- sum(weights*x^beta*y)
  c6 <- sum(weights*x^beta*y*log(x))
  c7 <- sum(weights*x^beta*y*log(x)^2)
  c8 <- sum(weights*y^2)
  # jacobian of root
  #(c1*c6*(c2*c6 - c3*c5) + c1*c5*(2*c3*c6 + c2*c7 - 2*c4*c5 - c3*c6))/(c2*(c2*c8 - c5^2)) -
  #  c1*c5*(c2*c6 - c3*c5)*(2*c3*(c2*c8 - c5^2) + c2*(2*c3*c8 - 2*c5*c6))/((c2*(c2*c8 - c5^2))^2)
  c1*(c2*c5*c7 - 2*c4*c5^2 + c2*c6^2)/(c2*(c2*c8 - c5^2)) +
    2*c1*c5*(c2*c6 - c3*c5)*(c3*c5^2 + c2*c5*c6 - 2*c2*c3*c8)/((c2*(c2*c8 - c5^2))^2)
}
