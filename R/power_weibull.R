# Functions for power-law fitting using Weibull probabilities
# 07/08/2023 - G. J. Meijer


#' Fit power law + weibull using loglikelihood fitting
#'
#' @description
#' Fit a power law distribution plus a Weibull shape parameter to a series
#' of test data (x, y), using the log-likelihood method.
#'
#' The variation in the data, defined as the ratio (y/y_fit) where y_fit is
#' predicted using a power law fit, is assumed to be Weibull distribution
#' with a mean of 1 and a to be fitted shape parameter.'
#'
#' To find the best fitting parameters for the power law, the probability of
#' all measurements is maximised.
#'
#' Measurements may be weighted according to their x-values in the form:
#' probability_weigthed = probability_{x,y}^weights, where `weights` is the
#' vector with weighting for each measured point. Given
#' that most biomechanical test relatively few thick roots, but thick roots
#' have a large effect on the calculated reinforcement, I suggest weighting
#' by the cross-sectional area, `weights = x^2`.
#'
#'
#' Function uses the function `optim()` from the package `stats` to conduct
#' the optimalisation.
#'
#' Initial parameter guesses for the log-likelihood fitting are estimated
#' using simple linear fits on transformed values of the input. These extra
#' initial fitting procedures ensure the `optim()` does not return inaccurate
#' values because the intial parameter guesses are too far removed from the
#' solution.
#'
#' @param x measured x-values (e.g. root diameters)
#' @param y measured y-values (e.g. root tensile strength)
#' @param weights weighting for each measurement. Default = 1, but a strong
#'   case can be made for weighting with `weights = x^2` because of the
#'   large effect of thick roots on root-reinforcement
#' @param guess initial guess for multiplier (log-transformed), power-law,
#'   and weibull shape parameter (log-transformed). If not defined, an initial
#'   guess is made using the function `power_weibull_initialguess()`
#' @return a list containing the fields
#'   * `loglikelihood`: the log-likelihood score of the fit
#'   * `multiplier`: fit power-law multiplier
#'   * `power`: fit power coefficient
#'   * `shape`: fit Weibull shape parameter
#' @examples
#' # input parameters
#' x <- seq(1, 10, l = 251)
#' alpha <- 5
#' beta <- -0.4
#' kappa <- 4
#' y <- alpha*x^beta*rweibull(length(x), kappa, 1/gamma(1 + 1/kappa))
#'
#' # fit
#' ft <- power_weibull_fit(x, y)
#' ft
#'
#' # plot
#' xp <- seq(min(x), max(x), l = 101)
#' yp <- ft$multiplier*xp^ft$power
#' plot(x, y)
#' lines(xp, yp, col = "red")
#' @export
#'
power_weibull_fit <- function(
    x, y,
    weights = rep(1, length(x)),
    guess = NULL
) {
  # do not use any observations with NA values in x and/or y
  i <- !(is.na(x) | is.na(y))
  x <- x[i]
  y <- y[i]
  # make guess if needed
  if (is.null(guess)) {
    guess <- power_weibull_initialguess(x, y)
  }
  # log-likelihood fit
  ft <- stats::optim(
    guess,
    power_weibull_likelihood,
    gr = power_weibull_likelihood_jacobian,
    method = "BFGS",
    lx = log(x), ly = log(y), weights = weights
  )
  # return
  data.frame(
    loglikelihood = -ft$value,
    multiplier = exp(ft$par[1])*gamma(1 + 1/exp(ft$par[3])),
    power = ft$par[2],
    shape = exp(ft$par[3])
  )
}


#' Likelihood function for power-law weibull loglikelihood fitting
#'
#' @description
#' Calculates the negative loglikelihood of the power-law + weibull shape
#' fit. Multiplier and weibull shape parameters are inputted as log-
#' transformed values to avoid issues with negative numbers.
#'
#' @param par fitting parameter vector: log(multiplier), power, log(shape).
#' @param lx,ly vectors with log-transformed x and y measurements
#' @param weights weighting for each measurement
#' @return negative log-likelihood score
#'
power_weibull_likelihood <- function(
    par,
    lx,
    ly,
    weights = rep(1, length(lx))
){
  # split param
  la <- par[1]
  b <- par[2]
  lk <- par[3]
  # likelihood
  (-lk*sum(weights) +
      la*exp(lk)*sum(weights) +
      b*exp(lk)*sum(weights*lx) -
      (exp(lk) - 1)*sum(weights*ly) +
      exp(-la*exp(lk))*sum(weights*exp(exp(lk)*(ly - b*lx))))
}


#' Derivative of function `power_weibull_likelihood()`
#'
#' @description
#' Generates the derivative of the results of the function
#' `power_weibull_likelihood()` with respect to its input argument `par`.
#'
#' Function is vectorised.
#'
#' @inheritParams power_weibull_likelihood
#' @return vector with derivatives with respect to `par`
#' @examples
#' # Compare analytical and numberical jacobians
#' x <- seq(1, 10, l = 51)
#' y <- 5*x^(0.5) + runif(length(x))
#' lx <- log(x)
#' ly <- log(y)
#' par <- c(log(3), -0.2, log(3.4))
#' weights <- x
#'
#' J <- power_weibull_likelihood_jacobian(
#'   par, lx, ly, weights = weights
#' )
#'
#' eps <- 1e-6
#' L <- power_weibull_likelihood(
#'   par, lx, ly, weights = weights
#' )
#' J2 <- rep(NA, 3)
#' for (i in 1:3) {
#'   dx <- rep(0, 3)
#'   dx[i] <- dx[i] + eps
#'   J2[i] <- (power_weibull_likelihood(
#'     par + dx, lx, ly, weights = weights) - L)/eps
#' }
#'
#' J
#' J2
#'
power_weibull_likelihood_jacobian <- function(
    par,
    lx,
    ly,
    weights = rep(1, length(lx))
) {
  # split param
  la <- par[1]
  b <- par[2]
  lk <- par[3]
  # derivatives
  dL_dla <- (
    exp(lk)*sum(weights) +
    -exp(lk - la*exp(lk))*sum(weights*exp(exp(lk)*(ly - b*lx)))
  )
  dL_db <- (
    exp(lk)*sum(weights*lx) +
    -exp(-la*exp(lk))*sum(weights*lx*exp(exp(lk)*(ly - b*lx) + lk))
  )
  dL_dlk <- (
    -sum(weights) +
    la*exp(lk)*sum(weights) +
    b*exp(lk)*sum(weights*lx) +
    -exp(lk)*sum(weights*ly) +
    -la*exp(lk - la*exp(lk))*sum(weights*exp(exp(lk)*(ly - b*lx))) +
    exp(-la*exp(lk))*sum(weights*(ly - b*lx)*exp(exp(lk)*(ly - b*lx) + lk))
  )
  # return
  c(dL_dla, dL_db, dL_dlk)
}


#' Create initial guess for weibull probability power-law fit
#'
#' @md
#' @description
#' Generate an initial guess for the power-law fit using weibull distributions
#' and loglikelihood fitting.
#'
#' The function makes a crude guess for the initial guess used in the function
#' `power_weibull()`, which consists of three numbers:
#'
#' * `log(alpha/gamma(1 + 1/kappa)) = log(alpha) - log(gamma(1 + 1/kappa))`
#' * `beta`
#' * `log(kappa)`
#'
#' where `alpha` and `beta` are the multiplier and power-law coefficient of
#' the best power-law fit, and `kappa` the shape parameter of the Weibull
#' distribution.
#' @param x,y arrays with measured values
#' @return vector with three elements, see Description
#' @examples
#' # input data
#' nr <- 51
#' dr <- seq(1, 10, l = nr)
#' kappa <- 10000
#' alpha <- 10
#' beta <- -0.2
#' tru <- alpha*dr^beta * rweibull(nr, kappa, 1/gamma(1 + 1/kappa))
#'
#' # initial guess
#' par <- power_weibull_initialguess(dr, tru)
#' # predict
#' xf <- seq(min(dr), max(dr), l = 251)
#' yf <- exp(par[1])*gamma(1 + 1/exp(par[3]))*xf^par[2]
#'
#' # plot fit
#' plot(dr, tru)
#' lines(xf, yf, col = "red")
#'
power_weibull_initialguess <- function(x, y) {
  # linear fit through log-transformed x, y data to obtain a and b for power-law fit
  ft1 <- stats::lm(log(y) ~ log(x))
  la <- ft1$coef[1]
  b <- ft1$coef[2]
  # get weibull shape parameter - crude guess using linear fitting on
  # log-transformed cumulative probability curve
  t <- y/exp(stats::predict(ft1))
  c <- seq(0.01, 0.99, l = length(t))
  ft2 <- stats::lm(log(-log(1 - c)) ~ log(sort(t)))
  k <- ft2$coef[2]
  # return
  as.vector(c(la - log(gamma(1 + 1/k)), b, log(k)))
}


