# Fit Weibull distributions
# 25/04/2024 - G. J. Meijer


#' Weibull likelihood fitting
#'
#' @description
#' Determine best fitting Weibull distribution using weighted loglikelihood
#' fitting.
#'
#' @md
#' @param x array with observations
#' @param weights weighting factor for each observation in `x`
#' @param method choose `newton` for gradient descent solving, using the
#'   `rootSolve::multiroot()` function, or `bisection` for bisection root
#'   solving algorithm using the `stats::uniroot()` function.
#' @param range two-value array to add to best guess, to define initial
#'   interval for bisection algorithm
#' @param start initial guess for weibull shape parameter. If not defined, an
#'   estimate is made based on fitting of the linearised cumulative density
#'   distribution
#' @return a list containing the fields
#'   * `loglikelihood`: the log-likelihood score of the fit
#'   * `shape`: fitted weibull shape parameter
#'   * `scale`: fitted weibull scale parameter
#' @export
#' @examples
#' x <- stats::rweibull(1000, 4, 2)
#' weibull_fit(x)
#'
weibull_fit <- function(
    x,
    weights = rep(1, length(x)),
    method = "bisection",
    range = c(-1, 1),
    start = NULL
) {
  # initial guess, based on fitting of linearised survival function
  if (is.null(start)) {
    start <- weibull_initialguess(x, weights = weights)
  }
  # fit shape parameter
  if (method == "newton") {
    shape <- rootSolve::multiroot(
      weibull_root,
      start,
      jacfunc = weibull_root_jacobian,
      x = x,
      weights = weights
    )$root
  } else if (method == "bisection") {
    shape <- stats::uniroot(
      weibull_root,
      interval = start + range,
      extendInt = "downX",
      x = x,
      weights = weights
    )$root
  } else {
    stop("`method` not recognised")
  }
  # calculate scale parameter
  c1 <- sum(weights)
  c3 <- sum(weights*x^shape)
  scale <- (c3/c1)^(1/shape)
  # get likelihood
  logpi <- log(shape) - shape*log(scale) + (shape - 1)*log(x) - (x/scale)^shape
  # return
  list(
    loglikelihood = sum(weights*logpi),
    shape = shape,
    scale = scale
  )
}


#' Initial guess for weibull distribution fitting
#'
#' @description
#' Return an initial guess for the Weibull shape parameter, to be used as a
#' starting point in function `weibull_fit()`.
#'
#' The estimate is made using linear regression on the linearlised cumulative
#' density function, i.e by fitting
#'
#'   log(-log(1 - P)) = shape*log(x) - shape*log(scale)
#'
#' where P is the cumulative density for each value x
#'
#' @inheritParams weibull_fit
#' @return initial guess for weibull shape parameter
#'
weibull_initialguess <- function(x, weights = rep(1, length(x))) {
  x <- sort(x)
  n <- length(x)
  yp <- (2*seq(n, 1) - 1)/(2*n)
  as.numeric(
    stats::lm(log(-log(yp)) ~ log(x), weights = weights)$coef[2]
  )
}


#' Root to solve for Weibull fitting
#'
#' @description
#' Root equation to solve for Weibull fitting
#'
#' @inheritParams weibull_fit
#' @param shape shape parameter
#' @return root
#'
weibull_root <- function(shape, x, weights = rep(1, length(x))) {
  # coefficients
  c1 <- sum(weights)
  c2 <- sum(weights*log(x))
  c3 <- sum(weights*x^shape)
  c4 <- sum(weights*x^shape*log(x))
  # return root
  c1/shape - c1*c4/c3 + c2
}


#' Jacobian of root to solve for Weibull fitting
#'
#' @description
#' Jacobian of root equation to solve for Weibull fitting.
#'
#' @inheritParams weibull_root
#' @return derivative of function `weibull_root()` with respect to input
#'   argument `shape`
#' @examples
#' # parameters
#' shape <- 4
#'
#' # create some data
#' x <- stats::rweibull(50, shape, 2)
#' w <- stats::runif(length(x))
#'
#' # check derivative - compare analytical and numerical solutions
#' eps <- 1e-6
#' v0 <- weibull_root(shape, x, weights = w)
#' v1 <- weibull_root(shape + eps, x, weights = w)
#' (v1 - v0)/eps
#' weibull_root_jacobian(shape, x, weights = w)
#'
weibull_root_jacobian <- function(shape, x, weights = rep(1, length(x))) {
  # coefficients
  c1 <- sum(weights)
  c2 <- sum(weights*log(x))
  c3 <- sum(weights*x^shape)
  c4 <- sum(weights*x^shape*log(x))
  c5 <- sum(weights*x^shape*log(x)^2)
  # return derivative of root
  -c1/shape^2 - c1*c5/c3 + c1*c4^2/c3^2
}

