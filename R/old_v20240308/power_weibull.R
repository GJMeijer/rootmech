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
#' probability_weigthed = probability_{x,y}^weights, where `weights` is the
#' vector with weighting for each measured point. Given
#' that most biomechanical test relatively few thick roots, but thick roots
#' have a large effect on the calculated reinforcement, I suggest weighting
#' by the cross-sectional area, `weights = x^2`.
#'
#' The problem is solved by solving for the roots of the derivatives of the
#' loglikelihood function with respect to `power` and `shape`. The power-law
#' multiplier can be expressed in terms of those two parameters, thus
#' simplifying the solving method. Root solving is done using the
#' `rootSolve::multiroot()` function, with an analytically defined jacobian.
#'
#' Initial estimate for `power` is determined by linear fitting of the
#' log-log transformed data. Subsequently, a Weibull distribution is fitted
#' to the scaled data (scaled by `x^power`) using the function `weibull_fit()`.
#'
#' @param x measured x-values (e.g. root diameters). These values are assumed
#'   to be already normalised by a reference value to ensure a unitless
#'   parameter.
#' @param y measured y-values (e.g. root tensile strength)
#' @param weights weighting for each measurement. Default = 1, but a strong
#'   case can be made for weighting with `weights = x^2` because of the
#'   large effect of thick roots on root-reinforcement
#' @return a list containing the fields
#'   * `loglikelihood`: the log-likelihood score of the fit
#'   * `multiplier`: fit power-law multiplier
#'   * `power`: fit power coefficient
#'   * `shape`: fit Weibull shape parameter
#' @export
#' @examples
#' # input parameters
#' x <- seq(1, 10, l = 101)
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
#'
#' # compare likelihood scores
#' y2 <- ft$multiplier*x^ft$power
#' p <- dweibull(y/y2, ft$shape, 1/gamma(1 + 1/ft$shape))
#'
#' ft$loglikelihood
#' sum(log(p))
#'
power_weibull_fit <- function(
    x,
    y,
    weights = rep(1, length(x))
) {
  # initial guess
  ft0 <- stats::lm(log(y) ~ log(x), weights = weights)
  power0 <- as.numeric(ft0$coefficients[2])
  # remove power-law effect from y-values
  yp <- y/(x^power0)
  # fit Weibull distribution to adjusted y-data
  shape0 <- weibull_fit(yp)$shape
  # find beta and kappa using a root solve
  sol <- rootSolve::multiroot(
    power_weibull_root,
    start = c(power0, shape0),
    jacfunc = power_weibull_root_jacobian,
    x = x,
    y = y,
    weights = weights
  )
  # assign
  power <- sol$root[1]
  shape <- sol$root[2]
  # get Weibull scale multiplier
  alpha <- (sum(weights*y^shape*x^(-power*shape))/sum(weights))^(1/shape)
  # calculate loglikelihood
  lambda <- alpha*x^power
  logpi <- log(shape) - shape*log(lambda) + (shape - 1)*log(y) - (y/lambda)^shape
  # return dataframe
  data.frame(
    loglikelihood = sum(weights*logpi),
    multiplier = alpha*gamma(1 + 1/shape),
    power = power,
    shape = shape
  )
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
  # split parameters: beta = power coefficient, kappa = Weibull shape parameter
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
#' # check jacobian - compare to numerical derivative
#' x <- seq(1, 10, l = 101)
#' alpha <- 20
#' beta <- -2
#' kappa <- 6.5
#' y <- alpha*x^beta*rweibull(length(x), kappa, 1/gamma(1 + 1/kappa))
#'
#' par <- c(beta, kappa)
#' eps <- 1e-6
#' r0 <- power_weibull_root(par, x, y)
#' r1 <- power_weibull_root(par + c(eps, 0), x, y)
#' r2 <- power_weibull_root(par + c(0, eps), x, y)
#' (cbind(r1, r2) - r0)/eps
#' power_weibull_root_jacobian(par, x, y)
#'
power_weibull_root_jacobian <- function(par, x, y, weights = rep(1, length(x))) {
  # split parameters: beta = power coefficient, kappa = Weibull shape parameter
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


#' Generate prediction interval for log-transformed power-law fit
#'
#' @description
#' Generate a prediction interval with a certain amount of confidence it
#' contains a random new observation. Fit is generated by the function
#' `power_logtransform_fit()`.
#'
#' @inheritParams power_normal_predictioninterval
#' @param shape Weibull scale parameter
#' @param scale Weibull scale parameter
#' @return dataframe with values for x (`x`), the mean value for y (`y`), and
#'   the lower and upper bounds of the confidence interval (`lower`, `upper`)
#' @export
#' @examples
#' x <- seq(2, 10, l = 250)
#' y <- 10*x^-0.5 * rweibull(length(x), shape = 4, scale = 1/gamma(1 + 1/4))
#' ft <- power_weibull_fit(x, y)
#' df <- power_weibull_predictioninterval(
#'   min(x), max(x), ft$multiplier, ft$power, ft$shape)
#'
#' plot(x, y, ylim = c(0, max(y)))
#' lines(df$x, df$mean, col = "red")
#' lines(df$x, df$lower, col = "blue")
#' lines(df$x, df$upper, col = "blue")
#'
power_weibull_predictioninterval <- function(
    xmin,
    xmax,
    multiplier,
    power,
    shape,
    scale = 1/gamma(1 + 1/shape),
    confidence = 0.95,
    n = 101
) {
  # generate or x-values range
  df <- data.frame(x = seq(xmin, xmax, l = n))
  # generate fit
  df$mean <- multiplier*df$x^power
  # distance - log-transformed space
  m <- stats::qweibull(
    0.5 + 0.5*c(-1, 1)*confidence,
    shape = shape,
    scale = scale
  )
  # lower and upper prediction interval
  df$lower <- m[1]*df$mean
  df$upper <- m[2]*df$mean
  # return dataframe
  df
}
