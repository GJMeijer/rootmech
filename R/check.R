#' Check if all points are perfectly on power-law line
#'
#' @md
#' @description
#' Checked by comparing co-linearity in log-log space
#'
#' @param x,y vectors of diameter and strength data
#' @param tol tolerance
#' @return list with fields:
#'   * `colinear`: boolean that is true when all points are colinear
#'   * `multiplier`: multiplier of power law fit
#'   * `exponent`: multiplier of power law fit
#' @export
#' @examples
#' x <- seq(3, 6, l = 5)
#' y <- 20*x^-0.5
#' check_perfectfit(x, y)
#'
#' y[2] <- y[2] + 0.1
#' check_perfectfit(x, y)
#'
check_perfectfit <- function(x, y, tol = .Machine$double.eps^0.5) {
  # unique pairs of x and y only
  df <- unique(data.frame(x = log(x), y = log(y)))
  # calculate derivatives
  dx <- diff(df$x)
  dy <- diff(df$y)
  # check
  if (nrow(df) == 1) {
    colinear <- TRUE
    exponent <- 0
    multiplier <- exp(df$y[1])
  } else if (nrow(df) == 2) {
    colinear <- TRUE
    exponent <- dy/dx
    multiplier <- exp(df$y[1] - exponent*df$x[1])
  } else {
    # calculate cross-products of vectors connecting subsequent data points
    cross_prod <- utils::head(dx, -1)*utils::tail(dy, -1) - utils::tail(dx, -1)*utils::head(dy, -1)
    if (length(dx) <= 1) {
      colinear <- TRUE
    } else {
      colinear <- all(abs(cross_prod) < tol)
    }
    exponent <- mean(dy/dx)
    multiplier <- mean(exp(utils::head(df$y, -1) - exponent*utils::head(df$x, -1)))
  }
  # return
  list(
    colinear = colinear,
    multiplier = multiplier,
    exponent = exponent
  )
}
