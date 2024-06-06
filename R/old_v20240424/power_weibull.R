# Fit power-law with weibull distributions of residuals
# 26/03/2024 - G. J. Meijer


#' Fit power law + weibull using loglikelihood fitting
#'
#' @description
#' Fit a power law distribution plus a Weibull shape parameter to a series
#' of test data (x, y), using the log-likelihood method.
#'
#' The variation in the data, defined as the ratio (y/y_fit) where y_fit is
#' predicted using a power law fit, is assumed to be Weibull distribution
#' with a mean of 1 and a to be fitted shape parameter.
#'
#' To find the best fitting parameters for the power law, the probability of
#' all measurements is maximised.
#'
#' Measurements may be weighted according to their x-values in the form:
#' probability_weighted = probability_{x,y}^weights, where `weights` is the
#' vector with weighting for each measured point. Given
#' that most biomechanical test relatively few thick roots, but thick roots
#' have a large effect on the calculated reinforcement, I suggest weighting
#' by the cross-sectional area, `weights = x^2`.
#'
#' The problem is solved by solving for the roots of the derivatives of the
#' loglikelihood function with respect to `exponent` and `shape`. The power-law
#' multiplier can be expressed in terms of those two parameters, thus
#' simplifying the solving method. Root solving is done using the
#' `rootSolve::multiroot()` function, with an analytically defined jacobian.
#'
#' Initial estimate for `exponent` is determined by linear fitting of the
#' log-log transformed data. Subsequently, a Weibull distribution is fitted
#' to the scaled data (scaled by `x^exponent`) using the function
#' `weibull_fit()`.
#'
#' @md
#' @param x measured x-values (e.g. root diameters). These values are assumed
#'   to be already normalised by a reference value to ensure a unitless
#'   parameter.
#' @param y measured y-values (e.g. root tensile strength)
#' @param weights weighting for each measurement. Default = 1, but a strong
#'   case can be made for weighting with `weights = x^2` because of the
#'   large effect of thick roots on root-reinforcement
#' @param start vector with user-defined initial guess for the exponent
#'   and Weibull shape parameter. If not defined, an initial
#'   guess is made based on linear fitting of log-transformed data (to find
#'   the exponent) and subsequent Weibull fitting of the residuals
#' @param factor reduction factor that is used when the algorithm gets stuck
#'   when very different weights are used (this makes it hard to get a good
#'   initial guess). The weighting is reduced by a exponent `factor` until
#'   a good solution is found. The weighting is then iteratively increased
#'   again (using solutions from previous steps) until a solution with the
#'   original weighting is obtained.
#' @return a list containing the fields
#'   * `loglikelihood`: the log-likelihood score of the fit (y/y_fit)
#'   * `multiplier`: fit power-law multiplier
#'   * `exponent`: fit power exponent
#'   * `shape`: fit Weibull shape parameter
#' @export
#' @examples
#' # input parameters
#' x <- seq(1, 10, l = 10001)
#' multiplier <- 5
#' exponent <- -0.4
#' shape <- 4
#' y <- multiplier*x^exponent
#' y <- y*stats::rweibull(length(x), shape, 1/gamma(1 + 1/shape))
#'
#' # fit
#' ft <- power_weibull_fit(x, y)
#' ft
#'
#' # plot
#' xp <- seq(min(x), max(x), l = 101)
#' yp <- ft$multiplier*xp^ft$exponent
#' plot(x, y)
#' lines(xp, yp, col = "red")
#'
power_weibull_fit <- function(
    x,
    y,
    weights = rep(1, length(x)),
    start = NULL,
    factor = 0.80,
    itermax = log(0.05)/log(factor)
) {
  # make initial guess for beta and kappa, if needed
  if (is.null(start)) {
    start <- power_weibull_initialguess(x, y, weights = weights)
  }
  # find beta and kappa using a root solve
  sol <- rootSolve::multiroot(
    power_weibull_root,
    start = start,
    jacfunc = power_weibull_root_jacobian,
    x = x,
    y = y,
    weights = weights
  )
  # if solution not satisfying, e.g. because of weighting, iterate weight scaling
  fac <- 1
  iter <- 1
  while (((sol$root[2] < 1) | is.na(sol$root[2])) & (iter <= itermax))  {
    fac <- fac*factor
    iter <- iter + 1
    sol <- rootSolve::multiroot(
      power_weibull_root,
      start = power_weibull_initialguess(x, y, weights = weights^fac),
      jacfunc = power_weibull_root_jacobian,
      x = x,
      y = y,
      weights = weights^fac
    )
  }
  while (fac < 1) {
    fac <- fac/factor
    sol <- rootSolve::multiroot(
      power_weibull_root,
      start = sol$root,
      jacfunc = power_weibull_root_jacobian,
      x = x,
      y = y,
      weights = weights^fac
    )
  }
  # assign
  exponent <- sol$root[1]
  shape <- sol$root[2]
  # get predictions
  alpha <- (sum(weights*y^shape*x^(-exponent*shape))/sum(weights))^(1/shape)
  multiplier = alpha*gamma(1 + 1/shape)
  yp <- multiplier*x^exponent
  # calculate loglikelihood
  logpi <- stats::dweibull(
    y/yp,
    shape = shape,
    scale = 1/gamma(1 + 1/shape),
    log = TRUE
  )
  # return dataframe
  list(
    loglikelihood = sum(weights*logpi),
    multiplier = multiplier,
    exponent = exponent,
    shape = shape
  )
}


#' Initial guess for `power_weibull_fit()`
#'
#' @description
#' Generate an initial guess for exponent and Weibull shape parameter, used
#' by the root solving methods in `power_weibull_fit()`.
#'
#' The exponent is found through linear fitting of log-transformed x and y
#' values. The Weibull shape parameter is found using seperate Weibull fitting
#' of the residuals.
#'
#' @inheritParams power_weibull_fit
#' @return two-parameter vector, with initial guesses for the exponent and
#'   Weibull shape parameter
#'
power_weibull_initialguess <- function(x, y, weights = rep(1, length(x))) {
  # fit exponent through linear regression on log-transformed data
  ft0 <- stats::lm(log(y) ~ log(x), weights = weights)
  power0 <- as.numeric(ft0$coefficients[2])
  # remove power-law effect from y-values
  yp <- y/(x^power0)
  # fit Weibull distribution to adjusted y-data
  shape0 <- weibull_fit(yp)$shape
  # return vector with initial guess
  c(power0, shape0)
}


#' Roots to solve in function `power_weibull_fit()`
#'
#' @description
#' Function that takes fitting parameters (`power` and `shape`), and returns
#' a two-parameter vector that must be zero for the most likely values of
#' fitting parameters
#'
#' @inheritParams power_weibull_fit
#' @param par two-parameter vector containing the values for `power` and `shape`
#' @return two-parameter vector, to be root-solved
#'
power_weibull_root <- function(par, x, y, weights = rep(1, length(x))) {
  # split parameters: beta = exponent, kappa = Weibull shape parameter
  beta <- par[1]
  kappa <- par[2]
  # calculate intermediate coefficients
  c1 <- sum(weights)
  c2 <- sum(weights*log(x))
  c3 <- sum(weights*log(y))
  c4 <- sum(weights*y^kappa*x^(-beta*kappa))
  c5 <- sum(weights*y^kappa*x^(-beta*kappa)*log(x))
  c6 <- sum(weights*y^kappa*x^(-beta*kappa)*log(y))
  # calculate roots
  r1 <- kappa*(c1*c5/c4 - c2)
  r2 <- c1/kappa + c1/c4*(beta*c5 - c6) - beta*c2 + c3
  # return roots
  c(r1, r2)
}


#' Jacobian of function `power_weibull_root()`
#'
#' @description
#' Returns the derivative of the function `power_weibull_root()` with respect
#' to its input arguments `par`
#'
#' @inheritParams power_weibull_root
#' @return a 2*2 matrix, returning the derivative of the root function (rows)
#'   with respect to the two arguments in `par` (columns)
#' @examples
#' # define some parameters
#' alpha <- 20
#' beta <- -2
#' kappa <- 6.5
#'
#' # generate some data
#' x <- seq(1, 10, l = 101)
#' y <- alpha*x^beta*stats::rweibull(length(x), kappa, 1/gamma(1 + 1/kappa))
#'
#' # check derivative - compare analytical and numerical solutions
#' par <- c(beta, kappa)
#' eps <- 1e-6
#' r0 <- power_weibull_root(par, x, y)
#' r1 <- power_weibull_root(par + c(eps, 0), x, y)
#' r2 <- power_weibull_root(par + c(0, eps), x, y)
#' (cbind(r1, r2) - r0)/eps
#' power_weibull_root_jacobian(par, x, y)
#'
power_weibull_root_jacobian <- function(par, x, y, weights = rep(1, length(x))) {
  # split parameters: beta = exponent, kappa = Weibull shape parameter
  beta <- par[1]
  kappa <- par[2]
  # calculate intermediate coefficients
  c1 <- sum(weights)
  c2 <- sum(weights*log(x))
  c3 <- sum(weights*log(y))
  c4 <- sum(weights*y^kappa*x^(-beta*kappa))
  c5 <- sum(weights*y^kappa*x^(-beta*kappa)*log(x))
  c6 <- sum(weights*y^kappa*x^(-beta*kappa)*log(y))
  c7 <- sum(weights*y^kappa*x^(-beta*kappa)*log(x)^2)
  c8 <- sum(weights*y^kappa*x^(-beta*kappa)*log(x)*log(y))
  c9 <- sum(weights*y^kappa*x^(-beta*kappa)*log(y)^2)
  # calculate derivatives of roots with respect to input arguments
  dr1_dbeta <- kappa^2*c1/c4*(c5^2/c4 - c7)
  dr1_dkappa <- c1/c4*(c5 + kappa*c5/c4*(beta*c5 - c6) + kappa*(c8 - beta*c7)) - c2
  dr2_dkappa <- c1*(c6 - beta*c5)^2/c4^2 - c1/c4*(c9 - 2*beta*c8 + beta^2*c7) - c1/kappa^2
  # return jacobian matrix
  matrix(
    c(dr1_dbeta, dr1_dkappa, dr1_dkappa, dr2_dkappa),
    nrow = 2,
    ncol = 2,
    byrow = TRUE
  )
}


#' Covariance matrix of fitting parameters for `power_weibull_fit()`
#'
#' @description
#' Estimate the variance-covariance matrix for power law fitting parameters, based
#' on Fisher information matrix, for Weibull fit
#'
#' @inheritParams power_weibull_fit
#' @param multiplier,exponent power-law multiplier and exponent for the mean
#' @param shape Weibull shape parameter
#' @return 2*2 matrix (for parameter `multiplier` and
#'   `exponent`)
#' @examples
#' #' input parameters
#' x <- seq(1, 10, l = 101)
#' multiplier <- 5
#' exponent <- -0.4
#' shape <- 4
#' y <- multiplier*x^exponent
#' y <- y*stats::rweibull(length(x), shape, 1/gamma(1 + 1/shape))
#'
#' # fit
#' ft <- power_weibull_fit(x, y)
#'
#' # covariance matrix
#' power_weibull_covariancematrix(x, y, ft$multiplier, ft$exponent, ft$shape)
#'
power_weibull_covariancematrix <- function(
    x,
    y,
    multiplier,
    exponent,
    shape,
    weights = rep(1, length(x))
) {
  # coefficients
  c1 <- sum(weights)
  c4 <- sum(weights*y^shape*x^(-exponent*shape))
  c5 <- sum(weights*y^shape*x^(-exponent*shape)*log(x))
  c7 <- sum(weights*y^shape*x^(-exponent*shape)*log(x)^2)
  # second derivatives of log-transformed probabilities
  dL2_da2 <- shape*c1/multiplier^2 - shape*(shape + 1)*c4*multiplier^(-shape - 2)
  dL2_dadb <- -shape^2*c5*multiplier^(-shape - 1)
  dL2_db2 <- -shape^2*multiplier^(-shape)*c7
  # Fisher information matrix - multiplier & exponent
  fisher <- -matrix(c(dL2_da2, dL2_dadb, dL2_dadb, dL2_db2), nrow = 2)
  # covariance matrix
  solve(fisher)
}


#' Generate prediction interval for Weibull fit
#'
#' @description
#' Generate prediction interval with specified level of confidence, based on
#' Weibull power law fit
#'
#' @inheritParams power_weibull_covariancematrix
#' @param x root diameters at which to predict interval
#' @param level confidence level (fraction)
#' @return dataframe with field for diameter (`x`), average power law strength
#'   (`y`) and the lower and upper bound of the prediction interval (`ymin`,
#'   `ymax`)
#' @examples
#' x <- seq(1, 10, l = 10)
#' multiplier <- 20
#' exponent <- -0.5
#' shape <- 4
#' power_weibull_predictioninterval(x, multiplier, exponent, shape)
#'
power_weibull_predictioninterval <- function(
    x,
    multiplier,
    exponent,
    shape,
    level = 0.95
) {
  y <- multiplier*x^exponent
  scale <- 1/gamma(1 + 1/shape)
  ymin <- y*stats::qweibull(0.5 - 0.5*level, shape, scale)
  ymax <- y*stats::qweibull(0.5 + 0.5*level, shape, scale)
  data.frame(x = x, y = y, ymin = ymin, ymax = ymax)
}

