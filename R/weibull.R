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
    shape_range = c(1, 100)
) {
  # remove NA values
  x <- x[!is.na(x)]
  # fit, depending on which parameters are fixed
  if (!is.numeric(shape)) {
    if (!is.numeric(scale)) {
      if (!is.numeric(mu)) {
        # nothing known - fit shape and scale
        shape <- uniroot(
          function(k) {
            1/k -
              sum(weights*log(x)*x^k)/sum(weights*x^k) +
              sum(weights*log(x))/sum(weights)
          },
          shape_range
        )$root
        scale <- (sum(weights*x^shape)/sum(weights))^(1/shape)
      } else {
        # mean known - fit shape
        shape <- uniroot(
          function(k) {
            digamma(1 + 1/k)/k*((gamma(1 + 1/k)/mu)^k*sum(weights*x^k) - sum(weights)) +
              sum(weights)/k - log(mu/gamma(1 + 1/k))*sum(weights) + sum(weights*log(x)) -
              sum(weights*log(x*gamma(1 + 1/k)/mu)*(x*gamma(1 + 1/k)/mu)^k)
          },
          shape_range
        )$root
        scale <- mu/gamma(1 + 1/shape)
      }
    } else {
      if (!is.numeric(mu)) {
        # scale known - fit shape
        shape <- uniroot(
          function(k) {
            sum(weights)/k -
              log(scale)*sum(weights) +
              sum(weights*log(x)) -
              sum(weights*log(x/scale)*(x/scale)^k)
          },
          shape_range
        )$root
      } else {
        # scale and mean known - fit nothing but get shape
        shape <- uniroot(
          function(k) scale*gamma(1 + 1/k) - mu,
          shape_range
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
