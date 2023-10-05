#' @examples
#' # generate some data
#' par <- c(0, 1)
#' n <- 100
#' x <- sort(rnorm(n, par[1], par[2]))
#'
#' # normal distribution assumption
#' anderson_darling(x, par)
#' cmstatr::anderson_darling_normal(x = x)
#' nortest::ad.test(x)
#'
#' # weibull assumptions
#' pr <- weibull_fit(x)
#' anderson_darling(x, c(pr$shape, pr$scale), distribution = "weibull")
#'
anderson_darling <- function(x, par, distribution = "normal") {
  # sort
  x <- sort(x)
  # Get cumulative values
  if (distribution == "normal") {
    cum <- stats::pnorm(x, mean = par[1], sd = par[2])
    cum2 <- stats::pnorm(-x, mean = par[1], sd = par[2])
  } else if (distribution == "weibull") {
    cum <- stats::pweibull(x, shape = par[1], scale = par[2])
  }
  # calculate unadjusted Anderson-Darling statistic
  n <- length(x)
  i <- seq(n)
  A <- -n - 1/n*sum((2*i - 1)*(log(cum) + log(rev(cum2))))
  # p-value, i.e. the likelihood that the data is from the assumed distribution
  p <- ifelse(
    A >= 0.60,
    p <- exp(1.2937 - 5.709*A + 0.0186*A^2),
    ifelse(
      A >= 0.34,
      p <- exp(0.9177 - 4.279*A - 1.38*A^2),
      ifelse(
        A >= 0.20,
        1 - exp(-8.318 + 42.796*A - 59.938*A^2),
        1 - exp(-13.436 + 101.14*A - 223.73*A^2)
      )
    )
  )
  A
}
