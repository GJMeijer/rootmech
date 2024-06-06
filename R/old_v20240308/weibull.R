#' Fit weibull distribution to a series of observations
#'
#' @description
#' Fit a Weibull probability distribution to a series of observations. The
#' best fit is determined by using a root solve on the  derivative of the
#' loglikelihood function, using the default R function `uniroot()`
#'
#' @param x vector with observations
#' @param shape if not `NULL`, a fixed value for the shape parameter
#' @param scale if not `NULL`, a fixed value for the scale parameter
#' @param mu if not `NULL`, a fixed value for the mean of the Weibull
#'   distribution after scaling y values by x^power
#' @param weights a vector with weights for each observation
#' @param shape_range default `uniroot` interval. It is assumed that the actual
#'   shape parameter is somewhere in this range, otherwise the function will
#'   not work
#' @return dataframe with fields for `loglikelihood`, `shape` and `scale`
#' @export
#' @examples
#' shape <- 4
#' scale <- 8
#' x <- rweibull(1000, shape, scale)
#'
#' ft <- weibull_fit(x, shape = 3, mu = mean(x))
#'
#' hist(x, breaks = 25, freq = FALSE)
#' xp <- seq(0, max(x), l = 101)
#' yp <- dweibull(xp, ft$shape, ft$scale)
#' lines(xp, yp)
#'
weibull_fit <- function(
    x,
    shape = NULL,
    scale = NULL,
    mu = NULL,
    weights = rep(1, length(x)),
    shape_range = c(1, 10)
) {
  # remove NA values
  x <- x[!is.na(x)]
  # fit, depending on which parameters are fixed
  if (!is.numeric(shape)) {
    if (!is.numeric(scale)) {
      if (!is.numeric(mu)) {
        # nothing known - fit shape and scale
        shape <- stats::uniroot(
          function(k) {
            1/k -
              sum(weights*log(x)*x^k)/sum(weights*x^k) +
              sum(weights*log(x))/sum(weights)
          },
          shape_range,
          extendInt = "downX"
        )$root
        scale <- (sum(weights*x^shape)/sum(weights))^(1/shape)
      } else {
        # mean known - fit shape
        shape <- stats::uniroot(
          function(k) {
            digamma(1 + 1/k)/k*((gamma(1 + 1/k)/mu)^k*sum(weights*x^k) - sum(weights)) +
              sum(weights)/k - log(mu/gamma(1 + 1/k))*sum(weights) + sum(weights*log(x)) -
              sum(weights*log(x*gamma(1 + 1/k)/mu)*(x*gamma(1 + 1/k)/mu)^k)
          },
          shape_range,
          extendInt = "downX"
        )$root
        scale <- mu/gamma(1 + 1/shape)
      }
    } else {
      if (!is.numeric(mu)) {
        # scale known - fit shape
        shape <- stats::uniroot(
          function(k) {
            sum(weights)/k -
              log(scale)*sum(weights) +
              sum(weights*log(x)) -
              sum(weights*log(x/scale)*(x/scale)^k)
          },
          shape_range,
          extendInt = "downX"
        )$root
      } else {
        # scale and mean known - fit nothing but get shape
        shape <- stats::uniroot(
          function(k) scale*gamma(1 + 1/k) - mu,
          shape_range,
          extendInt = "upX"
        )$root
      }
    }
  } else {
    if (!is.numeric(scale)) {
      if (!is.numeric(mu)) {
        # shape known - fit scale
        scale <- (sum(weights*x^shape)/sum(weights))^(1/shape)
      } else {
        # shape and mean known - fit nothing
        scale <- mu/gamma(1 + 1/shape)
      }
    }
  }
  # calculate likelihood
  logp <- log(shape) - shape*log(scale) + (shape - 1)*log(x) - (x/scale)^shape
  loglikelihood <- sum(weights*logp)
  # return
  data.frame(
    loglikelihood = loglikelihood,
    shape = shape,
    scale = scale
  )
}


#' Fit a weibull survival function using non-linear least-squares regreggion
#'
#' @description
#' Fit a Weibull survival function to a series of `x` observations. The
#' Weibull parameters are found by minimising the sum of squared residuals,
#' using R's `stats::optim()` function.
#'
#' @param x vector with observations
#' @return dataframe with fitting parameters `shape` and `scale`
#' @export
#' @examples
#' # generate data and fit
#' x <- sort(stats::rweibull(25, 4, 6))
#' ft1 <- weibull_survival_nls_fit(x)
#'
#' # experimental data
#' n <- length(x)
#' y <- seq(1 - 1/(2*n), 1/(2*n), l = n)
#'
#' # likelihood fit
#' ft2 <- weibull_fit(x)
#'
#' # create prediction lines
#' xp <- seq(min(x), max(x), l = 101)
#' yp1 <- 1 - pweibull(xp, ft1$shape, ft1$scale)
#' yp2 <- 1 - pweibull(xp, ft2$shape, ft2$scale)
#'
#' # plot
#' plot(x, y)
#' lines(xp, yp1, col = "red")
#' lines(xp, yp2, col = "blue")
#'
weibull_survival_nls_fit <- function(x) {
  # remove NA and sort
  x <- sort(x[!is.na(x)])
  # y-positions
  y <- seq(1 - 1/(2*length(x)), 1/(2*length(x)), l = length(x))
  # optimise R^2
  sol <- stats::optim(
    weibull_survival_nls_initialguess(x),
    weibull_survival_nls_rss,
    gr = weibull_survival_nls_rss_jacobian,
    x = x, y = y,
    method = "L-BFGS-B",
    lower = c(0, 0),
    upper = c(Inf, Inf)
  )
  # calculate likelihood
  logp <- log(sol$par[1]) - sol$par[1]*log(sol$par[2]) + (sol$par[1] - 1)*log(x) - (x/sol$par[2])^sol$par[1]
  loglikelihood <- sum(logp)
  # return
  data.frame(
    loglikelihood = loglikelihood,
    shape = sol$par[1],
    scale = sol$par[2],
    rss = sol$value
  )
}


#' Initial guess for `weibull_survival_nls_fit()` function
#'
#' @description
#' Generate initial guess for `weibull_survival_nls_fit()` function. The
#' guess is made by linear fitting on the linearlised x versus y data
#'
#' @inheritParams weibull_survival_nls_fit
#' @param y cumulative values for each observation in `x`
#' @return vector with guess for shape and scale parameter
#'
weibull_survival_nls_initialguess <- function(
    x,
    y = seq(1 - 1/(2*length(x)), 1/(2*length(x)), l = length(x))
) {
  x2 <- sort(x)
  y2 <- log(-log(y))
  ft <- stats::lm(y2 ~ log(x2))
  as.vector(c(ft$coef[2], exp(-ft$coef[1]/ft$coef[2])))
}


#' Residual sum-of-squares for `weibull_survival_nls_fit()` optimalisation
#'
#' @description
#' Calculate the residual sum of squares for Weibull survival function
#' fitting.
#'
#' @inheritParams weibull_survival_nls_initialguess
#' @param par vector with fitting values (Weibull shape and scale parameters)
#' @return residual sum of squares (scalar)
#'
weibull_survival_nls_rss <- function(
    par,
    x,
    y = seq(1 - 1/(2*length(x)), 1/(2*length(x)), l = length(x))
) {
  yp <- exp(-(x/par[2])^par[1])
  sum((y - yp)^2)
}


#' Jacobian of `weibull_survival_nls_rss()`
#'
#' @description
#' Returns the derivative of the function `weibull_survival_nls_rss()` with
#' respect to fitting parameters `par`
#'
#' @inheritParams weibull_survival_nls_rss
#' @return vector with derivatives with respect to each element in `par`
#' @examples
#' # test derivative
#' par <- c(3, 5)
#' x <- rweibull(100, par[1], par[2])
#' eps <- 1e-6
#' z0 <- weibull_survival_nls_rss(par, x)
#' z1a <- weibull_survival_nls_rss(par + c(eps, 0), x)
#' z1b <- weibull_survival_nls_rss(par + c(0, eps), x)
#' (c(z1a, z1b) - z0)/eps
#' weibull_survival_nls_rss_jacobian(par, x)
#'
weibull_survival_nls_rss_jacobian <- function(
    par,
    x,
    y = seq(1 - 1/(2*length(x)), 1/(2*length(x)), l = length(x))
) {
  yp <- exp(-(x/par[2])^par[1])
  dyp_dpar1 <- -log(x/par[2]) * (x/par[2])^par[1] * yp
  dyp_dpar2 <- par[1]/par[2] * (x/par[2])^par[1] * yp
  c(sum(2*(yp - y)*dyp_dpar1), sum(2*(yp - y)*dyp_dpar2))
}

