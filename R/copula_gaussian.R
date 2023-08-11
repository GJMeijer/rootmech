# Gaussian copula functions
# 07/08/2023 - G. J. Meijer


#' Fit gaussian copula using loglikelihood fitting
#'
#' @description
#' Find the correlation coefficient for the best gaussian copula
#' fit through given two-dimensional cumulative density data.
#'
#' @inheritParams copula_gaussian_multiplier
#' @param weights vector with weights for each observation
#' @param eps small value, so that the correlation coefficient is set on
#'   determined on the interval -(1 - eps) <= rho <= (1 - eps). rho >= 1
#'   and rho <= -1 will lead to numerical problems.
#' @return dataframe with weighted loglikelihood of the best fit (field
#'   `loglikelihood`) and the optimal value of the correlation coefficient
#'   (field `rho`)
#' @export
#' @examples
#' # generate data
#' xy <- bivariate_normal_random(1000, rho = -0.6)
#' px <- pnorm(xy[, 1])
#' py <- pnorm(xy[, 2])
#'
#' # fit correlation coefficient
#' ft <- copula_gaussian_fit(px, py)
#'
#' # confidence interval
#' dc <- bivariate_normal_predictioninterval(rho = ft$rho)
#'
#' # plot observations + prediction interval
#' plot(xy$x, xy$y)
#' lines(dc$x, dc$y)
#'
copula_gaussian_fit <- function(
    px,
    py,
    weights = rep(1, length(px)),
    eps = 1e-12
) {
  # do not use any observations with NA values in px and/or py
  i <- !(is.na(px) | is.na(py))
  px <- px[i]
  py <- py[i]
  # initial guess
  rho0 <- copula_gaussian_initialguess(px, py)
  # optimize
  sol <- stats::optim(
    rho0,
    copula_gaussian_likelihood,
    gr = copula_gaussian_likelihood_jacobian,
    px = px,
    py = py,
    weights = weights,
    method = "L-BFGS-B",
    lower = -1 + eps,
    upper = 1 - eps,
    control = list(fnscale = -1)
  )
  # get actual loglikelihood score
  x <- stats::qnorm(px, mean = 0, sd = 1)
  y <- stats::qnorm(py, mean = 0, sd = 1)
  dx <- stats::dnorm(x, mean = 0, sd = 1)
  dy <- stats::dnorm(y, mean = 0, sd = 1)
  cm <- copula_gaussian_multiplier(px, py, rho = sol$par)
  L <- sum(weights*(dx + dy + cm))
  # return best rho + loglikelihood value
  data.frame(
    loglikelihood = L,
    rho = sol$par
  )
}


#' Objective loglikelihood function for gaussian copula fitting
#'
#' @description
#' Returns scaled version of the loglikelihood of the gaussian copula fit.
#' There is a single fitting parameter: correlation coefficient `rho`.
#'
#' The loglikelihood value is given by:
#'
#'   L = sum(weights*(dx + dy + cm))
#'
#' where dx and dy are the probability densities of the marginal distributions,
#' and cm the copula multiplier. Since only cm is dependent on the value
#' of rho, we can instead maximise instead:
#'
#'   L = sum(weights*cm)
#'
#' @inheritParams copula_gaussian_fit
#' @param rho correlation coefficient to fit
#' @return weighted loglikeihood value
copula_gaussian_likelihood <- function(
    rho,
    px,
    py,
    weights = rep(1, length(px))
) {
  # log-transformed copula multipliers
  cm <- copula_gaussian_multiplier(px, py, rho = rho, log = TRUE)
  # loglikelihood (scaled)
  sum(weights*cm)
}


#' Jacobian of the function `copula_gaussian_likelihood()`
#'
#' @description
#' Returns the derivative of the (scalar) function
#' `copula_gaussian_likelihood()` with respect to input parameter `rho`
#'
#' @inheritParams copula_gaussian_likelihood
#' @return derivative of `copula_gaussian_likelihood()` with respect to `rho`
#' @examples
#' # Test: compare analytical with numerical jacobian
#' rho <- 0.6
#' px <- seq(0.4, 0.5, l = 10)
#' py <- seq(0.2, 0.4, l = 10)
#' weights <- seq(0.8, 1.3, l = 10)
#'
#' eps <- 1e-6
#' L1 <- copula_gaussian_likelihood(rho, px, py)
#' L2 <- copula_gaussian_likelihood(rho + eps, px, py)
#' J2 <- (L2 - L1)/eps
#' J <- copula_gaussian_likelihood_jacobian(rho, px, py)
#'
#' c(J, J2)
copula_gaussian_likelihood_jacobian <- function(
    rho,
    px,
    py,
    weights = rep(1, length(px))
) {
  # log-transformed copula multipliers - derivatives
  cm_jac <- copula_gaussian_multiplier_jacobian(px, py, rho = rho, log = TRUE)
  # sum weighted log-transformed probabilities
  sum(weights*cm_jac$rho)
}


#' Initial guess for gaussian copula fitting
#'
#' @description
#' Make an initial estimate for correlation coefficient rho, which is
#' used in the gaussian copula fitting function `copula_gaussian_fit()`
#'
#' @inheritParams copula_gaussian_fit
#' @return initial guess for correlation coefficient `rho`
copula_gaussian_initialguess <- function(px, py) {
  x <- stats::qnorm(px, mean = 0, sd = 1)
  y <- stats::qnorm(py, mean = 0, sd = 1)
  stats::cov(x, y)
}


#' Gaussian copula multiplier
#'
#' @description
#' Return gaussian copula multipliers based on known 2-d cumulative probability
#' density values and known correlation coefficient.
#'
#' The function is vectorised.
#'
#' @param px,py vectors with cumulative probability densities
#' @param rho correlation coefficient
#' @param log if `log = TRUE`, return the log-transformed copula multipliers
#'   instead
#' @return vector with gaussian copula multipliers for each pair of `px` and
#' `py`.
#' @examples
#' # input
#' rho <- 0.8
#' mu <- c(2, 3)
#' sd <- c(1, 1)
#'
#' # generate some data
#' x <- mu[1] + 3*sd[1]*seq(-1, 1, l = 101)
#' y <- mu[2] + 3*sd[2]*seq(-1, 1, l = 101)
#'
#' # grid
#' df <- expand.grid(x = x, y = y, KEEP.OUT.ATTRS = FALSE)
#' # probability density, and cumulative probability density
#' df$dx <- dnorm(xy$x, mean = mu[1], sd = sd[1])
#' df$dy <- dnorm(xy$y, mean = mu[2], sd = sd[2])
#' df$px <- pnorm(xy$x, mean = mu[1], sd = sd[1])
#' df$py <- pnorm(xy$y, mean = mu[2], sd = sd[2])
#'
#' # copula multiplier
#' df$cm <- copula_gaussian_multiplier(df$px, df$py, rho = rho)
#'
#' # matrix with probability
#' z <- matrix(with(df, cm*dx*dy), nrow = length(x), ncol = length(y))
#'
#' # plot
#' contour(x, y, z)
#'
copula_gaussian_multiplier <- function(px, py, rho = 0, log = FALSE) {
  a <- sqrt(2)*erfi(2*px - 1)
  b <- sqrt(2)*erfi(2*py - 1)
  if (log == TRUE) {
    -0.5*log(1 - rho^2) - 0.5*(rho^2*(a^2 + b^2) - 2*rho*a*b)/(1 - rho^2)
  } else {
    exp(-(rho^2*(a^2 + b^2) - 2*rho*a*b)/(2*(1 - rho^2)))/sqrt(1 - rho^2)
  }
}


#' Derivative of `copula_gaussian_multiplier()`
#'
#' @description
#' Calculate derivative of copula multiplier calculated using the function
#' `copula_gaussian_multiplier()` with respect to its input parameters `px`,
#' `py` and `rho`.
#'
#' The function is vectorised.
#'
#' @inheritParams copula_gaussian_multiplier
#' @return list with derivatives, with fieldnames `px`, `py` and `rho`.
#' @examples
#' px <- 0.1
#' py <- 0.2
#' rho <- -0.8
#' log <- FALSE
#'
#' cm0 <- copula_gaussian_multiplier(px, py, rho, log = log)
#' J <- copula_gaussian_multiplier_jacobian(px, py, rho, log = log)
#'
#' eps <- 1e-6
#' J2 <- (c(
#'   copula_gaussian_multiplier(px + eps, py, rho, log = log),
#'   copula_gaussian_multiplier(px, py + eps, rho, log = log),
#'   copula_gaussian_multiplier(px, py, rho + eps, log = log)
#' ) - cm0)/eps
#'
#' J
#' J2
#'
copula_gaussian_multiplier_jacobian <- function(px, py, rho = 0, log = FALSE) {
  # derivatives of a and b
  erfix <- erfi(2*px - 1)
  erfiy <- erfi(2*py - 1)
  a <- sqrt(2)*erfix
  b <- sqrt(2)*erfiy
  da_dpx <- sqrt(2*pi)*exp(erfix^2)
  db_dpy <- sqrt(2*pi)*exp(erfiy^2)
  # derivatives of multiplier (cm)
  if (log == TRUE) {
    dcm_da <- -(rho*(b - a*rho))/(rho^2 - 1)
    dcm_db <- -(rho*(a - b*rho))/(rho^2 - 1)
    dcm_drho <- (-rho^3 + a*b*rho^2 + (-a^2 - b^2 + 1)*rho + a*b)/(rho^2 - 1)^2
  } else {
    dcm_da <- -(exp((rho^2*(a^2 + b^2) - 2*a*b*rho)/(2*rho^2 - 2))*(-2*a*rho^2 + 2*b*rho))/((1 - rho^2)^(1/2)*(2*rho^2 - 2))
    dcm_db <- -(exp((rho^2*(a^2 + b^2) - 2*a*b*rho)/(2*rho^2 - 2))*(-2*b*rho^2 + 2*a*rho))/((1 - rho^2)^(1/2)*(2*rho^2 - 2))
    dcm_drho <- (exp((rho^2*(a^2 + b^2) - 2*a*b*rho)/(2*rho^2 - 2))*(-a^2*rho + a*b*rho^2 + a*b - b^2*rho - rho^3 + rho))/(1 - rho^2)^(5/2)
  }
  # return derivatives
  list(
    px = dcm_da*da_dpx,
    py = dcm_db*db_dpy,
    rho = dcm_drho
  )
}

