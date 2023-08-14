#' Fit power-law function using R^2
#'
#' @description
#' Fit a power law function between `x` and `y` in the form y = a*x^b.
#'
#' The fit is generated using the `stats::nlm()` non-linear least-squares
#' fitting function in R.
#'
#' @inheritParams power_weibull_fit
#' @return dataframe with R^2 value of best fit (field `r2`) as well as the
#'   multiplier (`multiplier`) and power coefficient (`power`)
#' @export
#' @examples
#' x <- seq(2, 10, l = 25)
#' y <- 10*x^-0.5 * rweibull(length(x), shape = 4, scale = 1/gamma(1 + 1/4))
#' ft <- power_r2_fit(x, y)
#'
#' xp <- seq(min(x), max(x), l = 101)
#' yp <- ft$multiplier*xp^ft$power
#'
#' plot(x, y)
#' lines(xp, yp)
#'
power_r2_fit <- function(
    x, y,
    weights = rep(1, length(x)),
    par0 = NULL
) {
  # generate initial guess if required
  if (is.null(par0)) {
    ft0 <- lm(log(y) ~ log(x), weights = weights)
    par0 <- as.vector(c(exp(coef(ft0)[1]), coef(ft0)[2]))
  }
  # non-linear least-squares fitting
  ft <- nls(
    y ~ a*x^b,
    start = list(a = par0[1], b = par0[2])
  )
  # get r2 value of fit
  yp <- predict(ft)
  r2 <- 1 - sum((yp - y)^2)/sum((y - mean(y))^2)
  # return dataframe with key data
  data.frame(
    r2 = r2,
    multiplier = coef(ft)[1],
    power = coef(ft)[2]
  )
}
