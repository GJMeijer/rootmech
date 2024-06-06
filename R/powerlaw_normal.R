#' Power-law fit with normally distributed residuals
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
#' @inheritParams powerlaw_gamma_fit
#' @param exponent will fix the value of the mean power-law exponent to a
#'   specific numeric value, if defined. If defined as a character string,
#'   beta = delta is assumed
#' @param sd_exponent will fix the value of the standard deviation power-law
#'   exponent to a specific numeric value, if defined. If defined as a
#'   character string, beta = delta is assumed
#' @param method Only used when finding a single root.
#'    Choose `newton` for gradient descent solving, using the
#'    `rootSolve::multiroot()` function, or `bisection` for bisection root
#'    solving algorithm using the `stats::uniroot()` function.
#'    When using a multivariate root, `rootSolve::multiroot()` will always be
#'    used.
#' @param range two-value array to add to best guess, to define initial
#'   interval for bisection algorithm
#' @param sd_exponent_range search range around beta, for initial guess of
#'   delta in the case where both beta and delta are to be fitted. See function
#'   `powerlaw_normal_freebetadelta_initialguess()` for more details.
#' @param n_range number of guesses for delta to be tried on range
#'   `sd_exponent_range`. See function
#'   `powerlaw_normal_freebetadelta_initialguess()` for more details.
#' @return a list containing the fields
#'   * `loglikelihood`: the log-likelihood score of the fit
#'   * `multiplier`: fitted power-law multiplier
#'   * `exponent`: fitted power coefficient
#'   * `sd_multiplier`: fitted standard deviation power-law multiplier
#'   * `sd_exponent`: fitted standard deviation power-law exponent
#' @export
#' @examples
#' # parameters
#' y0 <- 20
#' beta <- -0.25
#' sigma0 <- 2.5
#' delta <- -1
#' x <- seq(1, 8, l = 101)
#' mu <- y0*x^beta
#' sigma <- sigma0*x^delta
#' y <- abs(stats::rnorm(length(x), mu, sigma))
#' weights <- stats::runif(length(x), 0.9, 1.1)
#'
#' # fit
#' ft <- powerlaw_normal_fit(x, y, weights = weights)
#' ft <- powerlaw_normal_fit(x, y, weights = weights, sd_exponent = -2)
#' ft <- powerlaw_normal_fit(x, y, weights = weights, sd_exponent = NULL)
#' ft <- powerlaw_normal_fit(x, y, weights = weights, sd_exponent = "linked")
#' ft
#' xp <- seq(min(x), max(x), l = 251)
#' yp <- ft$multiplier*xp^ft$exponent
#'
#' # plot data and best fit
#' plot(x, y)
#' lines(xp, yp, col = "red")
#'
powerlaw_normal_fit <- function(
    x,
    y,
    weights = rep(1, length(x)),
    method = "bisection",
    range = c(-1, 1),
    start = NULL,
    exponent = NULL,
    sd_exponent = 0,
    sd_exponent_range = c(-2, 2),
    n_range = 9
) {
  if (is.character(exponent) | is.character(sd_exponent)) {
    # beta = delta
    if (is.null(start)) {
      start <- powerlaw_normal_linkedbetadelta_initialguess(x, y, weights = weights)
    }
    if (method == "newton") {
      beta <- rootSolve::multiroot(
        powerlaw_normal_linkedbetadelta_root,
        start,
        jacfunc = powerlaw_normal_linkedbetadelta_root_jacobian,
        x = x,
        y = y,
        weights = weights
      )$root
    } else if (method == "bisection") {
      beta <- stats::uniroot(
        powerlaw_normal_linkedbetadelta_root,
        start + range,
        extendInt = "downX",
        x = x,
        y = y,
        weights = weights
      )$root
    } else {
      stop("`method` not recognised")
    }
    delta <- beta
  } else if (is.numeric(exponent)) {
    if (is.numeric(sd_exponent)) {
      # beta known, delta known
      beta <- exponent
      detal <- sd_exponent
    } else {
      # beta known, delta unknown
      beta <- exponent
      if (is.null(start)) {
        start <- powerlaw_normal_freedelta_initialguess(
          x, y, weights = weights, beta = beta
        )
      }
      if (method == "newton") {
        delta <- rootSolve::multiroot(
          powerlaw_normal_freedelta_root,
          start,
          jacfunc = powerlaw_normal_freedelta_root_jacobian,
          x = x,
          y = y,
          beta = beta,
          weights = weights
        )$root
      } else if (method == "bisection") {
        delta <- stats::uniroot(
          powerlaw_normal_freedelta_root,
          start + range,
          extendInt = "downX",
          x = x,
          y = y,
          beta = beta,
          weights = weights
        )$root
      } else {
        stop("`method` not recognised")
      }
    }
  } else {
    if (is.numeric(sd_exponent)) {
      # beta unknown, delta known
      delta <- sd_exponent
      if (is.null(start)) {
        start <- powerlaw_normal_freebeta_initialguess(
          x, y, weights = weights, delta = delta
        )
      }
      if (method == "newton") {
        beta <- rootSolve::multiroot(
          powerlaw_normal_freebeta_root,
          start,
          jacfunc = powerlaw_normal_freebeta_root_jacobian,
          x = x,
          y = y,
          delta = delta,
          weights = weights
        )$root
      } else if (method == "bisection") {
        beta <- stats::uniroot(
          powerlaw_normal_freebeta_root,
          start + range,
          extendInt = "downX",
          x = x,
          y = y,
          delta = delta,
          weights = weights
        )$root
      } else {
        stop("`method` not recognised")
      }
    } else {
      # beta unknown, delta unknown
      if (is.null(start)) {
        start <- powerlaw_normal_freebetadelta_initialguess(
          x,
          y,
          weights = weights,
          method = method,
          sd_exponent_range = sd_exponent_range,
          n_range = n_range
        )
      }
      par <- rootSolve::multiroot(
        powerlaw_normal_freebetadelta_root,
        start,
        jacfunc = powerlaw_normal_freebetadelta_root_jacobian,
        x = x,
        y = y,
        weights = weights
      )$root
      beta <- par[1]
      delta <- par[2]
    }
  }
  # get multipliers
  c1 <- sum(weights)
  c3 <- sum(weights*x^(2*beta - 2*delta))
  c6 <- sum(weights*x^(beta - 2*delta)*y)
  c9 <- sum(weights*x^(-2*delta)*y^2)
  y0 <- c6/c3
  sigma0 <- sqrt(c9/c1 - c6^2/(c1*c3))
  # loglikelihood
  logL <- powerlaw_normal_freebetadelta_loglikelihood(
    c(y0, beta, sigma0, delta),
    x,
    y,
    weights = weights,
    deriv = 0
  )
  # return list
  list(
    loglikelihood = logL,
    multiplier = y0,
    exponent = beta,
    sd_multiplier = sigma0,
    sd_exponent = delta
  )
}


#' Generate initial guess for `powerlaw_normal_fit()` with two unknown exponents
#'
#' @description
#' Generate an initial guess for mean and standard deviation power-law
#' exponents, to be used in function `powerlaw_normal_fit()`.
#'
#' The mean power-law exponent (beta) is estimated using linear regression on
#' log-transformed x and y data.
#'
#' The standard deviation power-law exponent (delta) is subsequently fitted
#' using loglikelihood root solving, assuming <beta> as fixed.
#'
#' @inheritParams powerlaw_normal_fit
#' @return estimate for power-law exponents for the mean and standard deviation
#'
powerlaw_normal_freebetadelta_initialguess_old <- function(
    x,
    y,
    weights = rep(1, length(x)),
    method = "bisection",
    range = c(-1, 1)
) {
  # initial guess for beta - linear regression on log-transformed data
  beta <- powerlaw_normal_freebeta_initialguess(
    x, y, weights = weights
  )
  delta <- 0
  # refine guess for beta - normal fit with heteroscedastic sds
  beta <- powerlaw_normal_fit(
    x,
    y,
    weights = weights,
    sd_exponent = delta,
    start = beta,
    method = method,
    range = range
  )$exponent
  # make guess for delta - assume known value of beta
  delta <- powerlaw_normal_fit(
    x,
    y,
    weights = weights,
    exponent = beta,
    sd_exponent = NULL,
    start = beta,
    method = method,
    range = range
  )$sd_exponent
  # return starting guess
  c(beta, delta)
}


#' Generate initial guess for `powerlaw_normal_fit()` with two unknown exponents
#'
#' @description
#' Generate an initial guess for mean and standard deviation power-law
#' exponents, to be used in function `powerlaw_normal_fit()`.
#'
#' The mean power-law exponent (beta) is estimated using linear regression on
#' log-transformed x and y data.
#'
#' Subsequently, `n` equally spaced guess for the standard deviation
#' power-law exponents are tried, on the domain
#' (beta_lm + `sd_exponent_range`). The beta and delta value of the fit
#' with the lowest loglikelihood are taken as the initial guess
#'
#' @inheritParams powerlaw_normal_fit
#' @return estimate for power-law exponents for the mean and standard deviation
#'
powerlaw_normal_freebetadelta_initialguess <- function(
    x,
    y,
    weights = rep(1, length(x)),
    method = "bisection",
    range = c(-1, 1),
    sd_exponent_range = c(-2, 2),
    n_range = 9
) {
  # initial guess for beta - linear regression on log-transformed data
  beta <- powerlaw_normal_freebeta_initialguess(
    x, y, weights = weights
  )
  # try a number of guesses for delta
  delta_guess <- beta + seq(
    sd_exponent_range[1],
    sd_exponent_range[2],
    l = n_range
  )
  # fit using assumed delta
  logL <- sapply(
    delta_guess,
    function(d) {
      ft <- powerlaw_normal_fit(
        x,
        y,
        weights = weights,
        sd_exponent = d,
        method = method,
        range = range
      )
      c(ft$loglikelihood, ft$exponent)
    }
  )
  delta <- delta_guess[which.max(logL[1, ])]
  beta <- logL[2, which.max(logL[1, ])]
  # return starting guess
  c(beta, delta)
}


#' Calculate power-law normal loglikelihood - separate power law for sd
#'
#' @description
#' Calculate the weighted loglikelihood for a power-law fit with normal
#' residuals. Standard deviation follows its own power-law.
#' Can also be used to calculate the first or second partial
#' derivative with respect to fitting parameters.
#'
#' @inheritParams powerlaw_normal_fit
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
#' powerlaw_normal_freebetadelta_loglikelihood(par, x, y, weights = w)
#'
#' # test first derivative
#' eps <- 1e-6
#' f0 <- powerlaw_normal_freebetadelta_loglikelihood(
#'   par, x, y, weights = w, deriv = 0)
#' f1 <- powerlaw_normal_freebetadelta_loglikelihood(
#'   par + c(eps, 0, 0, 0), x, y, weights = w, deriv = 0
#' )
#' f2 <- powerlaw_normal_freebetadelta_loglikelihood(
#'   par + c(0, eps, 0, 0), x, y, weights = w, deriv = 0
#' )
#' f3 <- powerlaw_normal_freebetadelta_loglikelihood(
#'   par + c(0, 0, eps, 0), x, y, weights = w, deriv = 0
#' )
#' f4 <- powerlaw_normal_freebetadelta_loglikelihood(
#'   par + c(0, 0, 0, eps), x, y, weights = w, deriv = 0
#' )
#' (c(f1, f2, f3, f4) - f0)/eps
#' powerlaw_normal_freebetadelta_loglikelihood(
#'   par, x, y, weights = w, deriv = 1
#' )
#'
#' # test second derivative
#' f0 <- powerlaw_normal_freebetadelta_loglikelihood(
#'   par, x, y, weights = w, deriv = 1
#' )
#' f1 <- powerlaw_normal_freebetadelta_loglikelihood(
#'   par + c(eps, 0, 0, 0), x, y, weights = w, deriv = 1
#' )
#' f2 <- powerlaw_normal_freebetadelta_loglikelihood(
#'   par + c(0, eps, 0, 0), x, y, weights = w, deriv = 1
#' )
#' f3 <- powerlaw_normal_freebetadelta_loglikelihood(
#'   par + c(0, 0, eps, 0), x, y, weights = w, deriv = 1
#' )
#' f4 <- powerlaw_normal_freebetadelta_loglikelihood(
#'   par + c(0, 0, 0, eps), x, y, weights = w, deriv = 1
#' )
#' (cbind(f1, f2, f3, f4) - f0)/eps
#' powerlaw_normal_freebetadelta_loglikelihood(
#'   par, x, y, weights = w, deriv = 2
#' )
#'
powerlaw_normal_freebetadelta_loglikelihood <- function(
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


#' Root solve exponent in normal power-law - separate power law for sd
#'
#' @description
#' Root solve equation to obtain the exponent in power-law fitting
#' with normal distribution with separate power-law describing the standard
#' deviation
#'
#' @inheritParams powerlaw_normal_freebetadelta_loglikelihood
#' @param par vector with mean and standard deviation power-law exponents
#' @return value of root functions to solve
#'
powerlaw_normal_freebetadelta_root <- function(
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


#' Jacobian of function `powerlaw_normal_freebetadelta_root()`
#'
#' @description
#' Returns the derivative of the function `powerlaw_normal_freebetadelta_root()`
#' with respect to input argument `par`
#'
#' @inheritParams powerlaw_normal_freebetadelta_root
#' @return derivative of `powerlaw_normal_freebetadelta_root()` with respect to
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
#' f0 <- powerlaw_normal_freebetadelta_root(par, x, y, weights = w)
#' f1 <- powerlaw_normal_freebetadelta_root(par + c(eps, 0), x, y, weights = w)
#' f2 <- powerlaw_normal_freebetadelta_root(par + c(0, eps), x, y, weights = w)
#' (cbind(f1, f2) - f0)/eps
#' powerlaw_normal_freebetadelta_root_jacobian(par, x, y, weights = w)
#'
powerlaw_normal_freebetadelta_root_jacobian <- function(
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


#' Generate initial guess for `powerlaw_normal_fit()` with known delta
#'
#' @description
#' Generate an initial guess for mean power-law exponen, to be used in
#' function `powerlaw_normal_fit()`. The standard deviation power-law exponent
#' (delta) is assumed to be known.
#'
#' The mean power-law exponent (beta) is estimated using linear regression on
#' log-transformed x and y data.
#'
#' @inheritParams powerlaw_normal_fit
#' @param delta known value for standard deviation power law exponent
#' @return estimate for power-law exponents for the mean
#'
powerlaw_normal_freebeta_initialguess <- function(
    x,
    y,
    delta = 0,
    weights = rep(1, length(x))
) {
  # simple guess - linear regression of log transformed x and y data
  ft1 <- stats::lm(log(y) ~ log(x), weights = weights)
  as.numeric(ft1$coef[2])
}


#' Calculate power-law normal loglikelihood - fixed delta
#'
#' @description
#' Calculate the weighted loglikelihood for a power-law fit with normal
#' residuals. The power-law exponent for the standard deviation (delta) is
#' assumed to be known.
#'
#' Can also be used to calculate the first or second partial
#' derivative with respect to fitting parameters.
#'
#' @inheritParams powerlaw_normal_fit
#' @param par vector with fitting parameters (power-law multiplier, power-law
#'   exponent, standard deviation power-law multiplier)
#' @param delta known standard deviation power-law exponent
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
#' par <- c(y0, beta, sigma0)
#'
#' # check loglikelihood
#' mu <- y0*x^beta
#' sigma <- sigma0*x^delta
#' sum(w*stats::dnorm(y, mu, sigma, log = TRUE))
#' powerlaw_normal_freebeta_loglikelihood(par, x, y, delta = delta, weights = w)
#'
#' # test first derivative
#' eps <- 1e-6
#' f0 <- powerlaw_normal_freebeta_loglikelihood(
#'   par, x, y, delta = delta, weights = w, deriv = 0
#' )
#' f1 <- powerlaw_normal_freebeta_loglikelihood(
#'   par + c(eps, 0, 0), x, y, delta = delta, weights = w, deriv = 0
#' )
#' f2 <- powerlaw_normal_freebeta_loglikelihood(
#'   par + c(0, eps, 0), x, y, delta = delta, weights = w, deriv = 0
#' )
#' f3 <- powerlaw_normal_freebeta_loglikelihood(
#'   par + c(0, 0, eps), x, y, delta = delta, weights = w, deriv = 0
#' )
#' (c(f1, f2, f3) - f0)/eps
#' powerlaw_normal_freebeta_loglikelihood(
#'   par, x, y, delta = delta, weights = w, deriv = 1
#' )
#'
#' # test second derivative
#' f0 <- powerlaw_normal_freebeta_loglikelihood(
#'   par, x, y, delta = delta, weights = w, deriv = 1
#' )
#' f1 <- powerlaw_normal_freebeta_loglikelihood(
#'   par + c(eps, 0, 0), x, y, delta = delta, weights = w, deriv = 1
#' )
#' f2 <- powerlaw_normal_freebeta_loglikelihood(
#'   par + c(0, eps, 0), x, y, delta = delta, weights = w, deriv = 1
#' )
#' f3 <- powerlaw_normal_freebeta_loglikelihood(
#'   par + c(0, 0, eps), x, y, delta = delta, weights = w, deriv = 1
#' )
#' (cbind(f1, f2, f3) - f0)/eps
#' powerlaw_normal_freebeta_loglikelihood(
#'   par, x, y, delta = delta, weights = w, deriv = 2
#' )
#'
powerlaw_normal_freebeta_loglikelihood <- function(
    par,
    x,
    y,
    delta = 0,
    weights = rep(1, length(x)),
    deriv = 0
) {
  # split input
  y0 <- par[1]
  beta <- par[2]
  sigma0 <- par[3]
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
    dlogL_dy0 <- 0
    dlogL_dbeta <- 0
    dlogL_dsigma0 <- 0
    dlogL_ddelta <- 0


      c(
      (c6 - y0*c3)/sigma0^2,
      y0*(c7 - c4*y0)/sigma0^2,
      -c1/sigma0 + (c9 - 2*c6*y0 + c3*y0^2)/sigma0^3
    )
  } else if (deriv == 2) {
    d2logL_dy02 <- -c3/sigma0^2
    d2logL_dy0dbeta <- (c7 - 2*c4*y0)/sigma0^2
    d2logL_dy0dsigma0 <- 2*(c3*y0 - c6)/sigma0^3
    d2logL_dbeta2 <- y0*(c8 - 2*c5*y0)/sigma0^2
    d2logL_dbetadsigma0 <- 2*y0*(c4*y0 - c7)/sigma0^3
    d2logL_dsigma02 <- c1/sigma0^2 - 3*(c9 - 2*c6*y0 + c3*y0^2)/sigma0^4
    d2logL_ddelta2 <- (-2*c11 + 4*c8*y0 - 2*c5*y0^2)/sigma0^2
    matrix(
      c(
        d2logL_dy02, d2logL_dy0dbeta, d2logL_dy0dsigma0,
        d2logL_dy0dbeta, d2logL_dbeta2, d2logL_dbetadsigma0,
        d2logL_dy0dsigma0, d2logL_dbetadsigma0, d2logL_dsigma02
      ),
      nrow = 3
    )
  }
}


#' Root solve exponent in normal power-law - known delta
#'
#' @description
#' Root solve equation to obtain the exponent in power-law fitting
#' with normal distribution. Power-law exponent for standard deviation (delta)
#' is assumed to be known
#'
#' @inheritParams powerlaw_normal_freebeta_loglikelihood
#' @param beta value for mean power-law exponent
#' @return value of root function to solve
#'
powerlaw_normal_freebeta_root <- function(
    beta,
    x,
    y,
    delta = 0,
    weights = rep(1, length(x))
) {
  # coefficients
  c1 <- sum(weights)
  c3 <- sum(weights*x^(2*beta - 2*delta))
  c4 <- sum(weights*x^(2*beta - 2*delta)*log(x))
  c6 <- sum(weights*x^(beta - 2*delta)*y)
  c7 <- sum(weights*x^(beta - 2*delta)*y*log(x))
  c9 <- sum(weights*x^(-2*delta)*y^2)
  # temporary variables
  zeta <- 1/(c3^2*c9 - c3*c6^2)
  # return root
  c1*c6*(c3*c7 - c4*c6)*zeta
}


#' Jacobian of function `powerlaw_normal_freebeta_root()`
#'
#' @description
#' Returns the derivative of the function `powerlaw_normal_freebeta_root()`
#' with respect to input argument `beta`
#'
#' @inheritParams powerlaw_normal_freebeta_root
#' @return derivative of `powerlaw_normal_freebeta_root()` with respect to
#'   `beta`#'
powerlaw_normal_freebeta_root_jacobian <- function(
    beta,
    x,
    y,
    delta = 0,
    weights = rep(1, length(x))
) {
  # coefficients
  c1 <- sum(weights)
  c3 <- sum(weights*x^(2*beta - 2*delta))
  c4 <- sum(weights*x^(2*beta - 2*delta)*log(x))
  c5 <- sum(weights*x^(2*beta - 2*delta)*log(x)^2)
  c6 <- sum(weights*x^(beta - 2*delta)*y)
  c7 <- sum(weights*x^(beta - 2*delta)*y*log(x))
  c8 <- sum(weights*x^(beta - 2*delta)*y*log(x)^2)
  c9 <- sum(weights*x^(-2*delta)*y^2)
  # temporary variables
  zeta <- 1/(c3^2*c9 - c3*c6^2)
  dzeta_dbeta <- -(4*c3*c4*c9 - 2*c4*c6^2 - 2*c3*c6*c7)*zeta^2
  # root
  c1*(c3*c7 - c4*c6)*(c6*dzeta_dbeta + c7*zeta) +
    c1*c6*(c4*c7 + c3*c8 - 2*c5*c6)*zeta
}


#' Generate initial guess for `powerlaw_normal_fit()` with known beta
#'
#' @description
#' Generate an initial guess for standard deviation mean power-law exponent,
#' to be used in function `powerlaw_normal_fit()`. The mean power-law exponent
#' (beta) is assumed to be known.
#'
#' Currently, a very simple guess delta = beta is made.
#'
#' @inheritParams powerlaw_normal_fit
#' @param beta known value for power law exponent
#' @return estimate for power-law exponents for the standard deviation
#'
powerlaw_normal_freedelta_initialguess <- function(
    x,
    y,
    beta = 0,
    weights = rep(1, length(x))
) {
  # simple guess - equal to mean power-law exponent
  beta
}


#' Calculate power-law normal loglikelihood - fixed beta
#'
#' @description
#' Calculate the weighted loglikelihood for a power-law fit with normal
#' residuals. The power-law exponent for the mean (beta) is
#' assumed to be known.
#'
#' Can also be used to calculate the first or second partial
#' derivative with respect to fitting parameters.
#'
#' @inheritParams powerlaw_normal_fit
#' @param par vector with fitting parameters (power-law multiplier,
#'   standard deviation power-law multiplier, standard deviation power-law
#'   exponent)
#' @param beta known power-law exponent for the mean
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
#' par <- c(y0, sigma0, delta)
#'
#' # check loglikelihood
#' mu <- y0*x^beta
#' sigma <- sigma0*x^delta
#' sum(w*stats::dnorm(y, mu, sigma, log = TRUE))
#' powerlaw_normal_freedelta_loglikelihood(par, x, y, beta = beta, weights = w)
#'
#' # test first derivative
#' eps <- 1e-6
#' f0 <- powerlaw_normal_freedelta_loglikelihood(
#'   par, x, y, beta = beta, weights = w, deriv = 0
#' )
#' f1 <- powerlaw_normal_freedelta_loglikelihood(
#'   par + c(eps, 0, 0), x, y, beta = beta, weights = w, deriv = 0
#' )
#' f2 <- powerlaw_normal_freedelta_loglikelihood(
#'   par + c(0, eps, 0), x, y, beta = beta, weights = w, deriv = 0
#' )
#' f3 <- powerlaw_normal_freedelta_loglikelihood(
#'   par + c(0, 0, eps), x, y, beta = beta, weights = w, deriv = 0
#' )
#' (c(f1, f2, f3) - f0)/eps
#' powerlaw_normal_freedelta_loglikelihood(
#'   par, x, y, beta = beta, weights = w, deriv = 1
#' )
#'
#' # test second derivative
#' f0 <- powerlaw_normal_freedelta_loglikelihood(
#'   par, x, y, beta = beta, weights = w, deriv = 1
#' )
#' f1 <- powerlaw_normal_freedelta_loglikelihood(
#'   par + c(eps, 0, 0), x, y, beta = beta, weights = w, deriv = 1
#' )
#' f2 <- powerlaw_normal_freedelta_loglikelihood(
#'   par + c(0, eps, 0), x, y, beta = beta, weights = w, deriv = 1
#' )
#' f3 <- powerlaw_normal_freedelta_loglikelihood(
#'   par + c(0, 0, eps), x, y, beta = beta, weights = w, deriv = 1
#' )
#' (cbind(f1, f2, f3) - f0)/eps
#' powerlaw_normal_freedelta_loglikelihood(
#'   par, x, y, beta = beta, weights = w, deriv = 2
#' )
#'
powerlaw_normal_freedelta_loglikelihood <- function(
    par,
    x,
    y,
    beta = 0,
    weights = rep(1, length(x)),
    deriv = 0
) {
  # split input
  y0 <- par[1]
  sigma0 <- par[2]
  delta <- par[3]
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
      -c1/sigma0 + (c9 - 2*c6*y0 + c3*y0^2)/sigma0^3,
      (c10 - 2*c7*y0 + c4*y0^2)/sigma0^2 - c2
    )
  } else if (deriv == 2) {
    d2logL_dy02 <- -c3/sigma0^2
    d2logL_dy0dsigma0 <- 2*(c3*y0 - c6)/sigma0^3
    d2logL_dy0ddelta <- 2*(c4*y0 - c7)/sigma0^2
    d2logL_dsigma02 <- c1/sigma0^2 - 3*(c9 - 2*c6*y0 + c3*y0^2)/sigma0^4
    d2logL_dsigma0ddelta <- 2*(-c10 + 2*c7*y0 - c4*y0^2)/sigma0^3
    d2logL_ddelta2 <- (-2*c11 + 4*c8*y0 - 2*c5*y0^2)/sigma0^2
    matrix(
      c(
        d2logL_dy02, d2logL_dy0dsigma0, d2logL_dy0ddelta,
        d2logL_dy0dsigma0, d2logL_dsigma02, d2logL_dsigma0ddelta,
        d2logL_dy0ddelta, d2logL_dsigma0ddelta, d2logL_ddelta2
      ),
      nrow = 3
    )
  }
}


#' Root solve exponent in normal power-law - known beta
#'
#' @description
#' Root solve equation to obtain the standard deviation exponent in power-law
#' fitting with normal distribution. Power-law exponent for mean (beta)
#' is assumed to be known.
#'
#' @inheritParams powerlaw_normal_freedelta_loglikelihood
#' @param delta value for standard deviation power-law exponent
#' @return value of root function to solve
#'
powerlaw_normal_freedelta_root <- function(
    delta,
    x,
    y,
    beta = 0,
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
  # root
  c1*(c3^2*c10 - 2*c3*c6*c7 + c4*c6^2)*zeta - c2
}


#' Jacobian of function `powerlaw_normal_freedelta_root()`
#'
#' @description
#' Returns the derivative of the function `powerlaw_normal_freedelta_root()`
#' with respect to input argument `delta`
#'
#' @inheritParams powerlaw_normal_freedelta_root
#' @return derivative of `powerlaw_normal_freedelta_root()` with respect to
#'   `delta`#'
powerlaw_normal_freedelta_root_jacobian <- function(
    delta,
    x,
    y,
    beta = 0,
    weights = rep(1, length(x))
) {
  # coefficients
  c1 <- sum(weights)
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
  # root
  c1*(c3^2*c10 - 2*c3*c6*c7 + c4*c6^2)*dzeta_ddelta +
    c1*(-4*c3*c4*c10 - 2*c3^2*c11 + 4*c3*c7^2 + 4*c3*c6*c8 - 2*c5*c6^2)*zeta
}


#' Generate initial guess for `powerlaw_normal_fit()` with beta=delta
#'
#' @description
#' Generate an initial guess for mean power-law exponent, to be used in
#' function `powerlaw_normal_fit()`. The standard deviation power-law exponent
#' (delta) is assumed to be equal to beta
#'
#' The mean power-law exponent (beta = delta) is estimated using linear
#' regression on log-transformed x and y data.
#'
#' @inheritParams powerlaw_normal_fit
#' @return estimate for power-law exponent for the mean
#'
powerlaw_normal_linkedbetadelta_initialguess <- function(
    x,
    y,
    weights = rep(1, length(x))
) {
  # simple guess - linear regression of log transformed x and y data
  ft1 <- stats::lm(log(y) ~ log(x), weights = weights)
  as.numeric(ft1$coef[2])
}


#' Calculate power-law normal loglikelihood - beta=delta
#'
#' @description
#' Calculate the weighted loglikelihood for a power-law fit with normal
#' residuals. The power-law exponent for the mean (beta) is assumed to be
#' equal to the power-law exponent for the mean (beta)
#'
#' Can also be used to calculate the first or second partial
#' derivative with respect to fitting parameters.
#'
#' @inheritParams powerlaw_normal_fit
#' @param par vector with fitting parameters (power-law multiplier, power-law
#'   exponent, standard deviation power-law multiplier)
#' @param deriv order of partial derivative requested
#' @return loglikelihood, or its partial derivatives to order `deriv`
#' @examples
#' # define parameters
#' y0 <- 20
#' beta <- -0.5
#' sigma0 <- 2
#' delta <- beta
#' x <- seq(1, 8, l = 101)
#' mu <- y0*x^beta
#' sigma <- sigma0*x^delta
#' y <- rnorm(length(x), mu, sigma)
#' w <- stats::runif(length(x), 0.8, 1.2)
#' par <- c(y0, beta, sigma0)
#'
#' # check loglikelihood
#' mu <- y0*x^beta
#' sigma <- sigma0*x^delta
#' sum(w*stats::dnorm(y, mu, sigma, log = TRUE))
#' powerlaw_normal_linkedbetadelta_loglikelihood(par, x, y, weights = w)
#'
#' # test first derivative
#' eps <- 1e-6
#' f0 <- powerlaw_normal_linkedbetadelta_loglikelihood(
#'   par, x, y, weights = w, deriv = 0
#' )
#' f1 <- powerlaw_normal_linkedbetadelta_loglikelihood(
#'   par + c(eps, 0, 0), x, y, weights = w, deriv = 0
#' )
#' f2 <- powerlaw_normal_linkedbetadelta_loglikelihood(
#'   par + c(0, eps, 0), x, y, weights = w, deriv = 0
#' )
#' f3 <- powerlaw_normal_linkedbetadelta_loglikelihood(
#'   par + c(0, 0, eps), x, y, weights = w, deriv = 0
#' )
#' (c(f1, f2, f3) - f0)/eps
#' powerlaw_normal_linkedbetadelta_loglikelihood(
#'   par, x, y, weights = w, deriv = 1
#' )
#'
#' # test second derivative
#' f0 <- powerlaw_normal_linkedbetadelta_loglikelihood(
#'   par, x, y, weights = w, deriv = 1
#' )
#' f1 <- powerlaw_normal_linkedbetadelta_loglikelihood(
#'   par + c(eps, 0, 0), x, y, weights = w, deriv = 1
#' )
#' f2 <- powerlaw_normal_linkedbetadelta_loglikelihood(
#'   par + c(0, eps, 0), x, y, weights = w, deriv = 1
#' )
#' f3 <- powerlaw_normal_linkedbetadelta_loglikelihood(
#'   par + c(0, 0, eps), x, y, weights = w, deriv = 1
#' )
#' (cbind(f1, f2, f3) - f0)/eps
#' powerlaw_normal_linkedbetadelta_loglikelihood(
#'   par, x, y, weights = w, deriv = 2
#' )
#'
powerlaw_normal_linkedbetadelta_loglikelihood <- function(
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
  delta <- beta
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
    dlogL_dy0 <- (c6 - y0*c3)/sigma0^2
    dlogL_dbeta <- y0*(c7 - c4*y0)/sigma0^2
    dlogL_dsigma0 <- -c1/sigma0 + (c9 - 2*c6*y0 + c3*y0^2)/sigma0^3
    dlogL_ddelta <- (c10 - 2*c7*y0 + c4*y0^2)/sigma0^2 - c2
    c(dlogL_dy0, dlogL_dbeta + dlogL_ddelta, dlogL_dsigma0)
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
        d2logL_dy02,
        d2logL_dy0dbeta + d2logL_dy0ddelta,
        d2logL_dy0dsigma0,
        d2logL_dy0dbeta + d2logL_dy0ddelta,
        d2logL_dbeta2 + 2*d2logL_dbetaddelta + d2logL_ddelta2,
        d2logL_dbetadsigma0 + d2logL_dsigma0ddelta,
        d2logL_dy0dsigma0,
        d2logL_dbetadsigma0 + d2logL_dsigma0ddelta,
        d2logL_dsigma02
      ),
      nrow = 3
    )
  }
}


#' Root solve exponent in normal power-law - beta=delta
#'
#' @description
#' Root solve equation to obtain the mean exponent in power-law
#' fitting with normal distribution (beta). Power-law exponent for the
#' standard deviation (delta) is assumed to be equal to beta.
#'
#' @inheritParams powerlaw_normal_linkedbetadelta_loglikelihood
#' @param beta value for mean power-law exponent
#' @return value of root function to solve
#'
powerlaw_normal_linkedbetadelta_root <- function(
    beta,
    x,
    y,
    weights = rep(1, length(x))
) {
  # unpack parameters
  delta <- beta
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
  dlogL_dbeta + dlogL_ddelta
}


#' Jacobian of function `powerlaw_normal_linkedbetadelta_root()`
#'
#' @description
#' Returns the derivative of the function `powerlaw_normal_linkedbetadelta_root()`
#' with respect to input argument `beta`
#'
#' @inheritParams powerlaw_normal_linkedbetadelta_root
#' @return derivative of `powerlaw_normal_linkedbetadelta_root()` with respect to
#'   `beta`#'
powerlaw_normal_linkedbetadelta_root_jacobian <- function(
    beta,
    x,
    y,
    weights = rep(1, length(x))
) {
  # unpack parameters
  delta <- beta
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
  # return
  d2logL_dbeta2 + 2*d2logL_dbetaddelta + d2logL_ddelta2
}


#' Calculate Kolmogorov-Smirnov parameters for power-law fit + normal
#'
#' @description
#' Calculate Kolmogorov-Smirnov parameter for power-law fit with normally
#' distributed residuals
#'
#' @md
#' @inheritParams powerlaw_normal_fit
#' @param multiplier,exponent multiplier and exponent for power-law fit
#'   describing the mean
#' @param sd_multiplier,sd_exponent multiplier and exponent for power-law fit
#'   describing the standard deviation of residuals
#' @return list with fields
#'   * `ks_distance`: Kolmogorov-Smirnov distance
#' @examples
#' x <- seq(1, 7, l = 251)
#' y <- stats::rnorm(length(x), 2, 0.5)
#' powerlaw_normal_ks(x, y, 2, 0, 0.5)
#'
powerlaw_normal_ks <- function(
    x,
    y,
    multiplier,
    exponent,
    sd_multiplier,
    sd_exponent = 0
) {
  # normalise
  yn <- sort((y - multiplier*x^exponent)/(sd_multiplier*x^sd_exponent))
  # prediction
  C2 <- rep(stats::pnorm(yn, 0, 1), each = 2)
  # real cumulative
  nx <- length(x)
  C1 <- rep(seq(0, 1, l = nx + 1), each = 2)[2:(2*nx + 1)]
  # distance
  distances <- C2 - C1
  i_max <- which.max(distances)
  # experimental cumulative trace
  df_exp <- data.frame(
    x = rep(yn, each = 2),
    y = C1
  )
  # fitted cumulative trace
  xp <- seq(min(yn, na.rm = TRUE), max(yn, na.rm = TRUE), l = n)
  df_fit <- data.frame(
    x = xp,
    y = stats::pnorm(xp, 0, 1)
  )
  # line for max
  df_con <- data.frame(
    x = rep(yn[floor((i_max + 1)/2)], 2),
    y = c(C1[i_max], C2[i_max])
  )
  # return list (in case of future expansion)
  list(
    ks_distance = distances[i_max],
    experimental = df_exp,
    fit = df_fit,
    difference = df_con
  )
}

