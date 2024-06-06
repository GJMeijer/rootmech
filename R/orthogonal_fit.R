#' Orthagonal linear fitting
#'
#' @description
#' Conduct an orthagonal linear fit, minimising the sum of squared distances
#' between datapoints and the linear fit
#'
#' @param x,y vectors with x and y data
#' @return vector with angle (in radians) and offset, describing the fitting
#'   line. The angle is defined as positive rotating from the x-axis to the
#'   y-axis, and the offset is defined as the distance between the origin
#'   and the fitting line
#' @export
#' @examples
#' x <- c(1, 2, 3, 4, 5)
#' y <- c(2, 3, 3, 5, 5)
#' ft <- orthagonal_fit(x, y)
#' xp <- seq(min(x), max(x), l = 11)
#' yp <- ft[2]/cos(ft[1]) + xp*tan(ft[1])
#' plot(x, y, xlim = c(0, 6), ylim = c(0, 6))
#' lines(xp, yp)
#'
orthagonal_fit <- function(x, y) {
  # number of points
  n <- length(x)
  # parameters
  a <- sum(x^2) - sum(x)*sum(x)/n
  b <- sum(y^2) - sum(y)*sum(y)/n
  c <- sum(x*y) - sum(x)*sum(y)/n
  # angle
  t <- 0.5*atan2(2*c, a - b)
  # offset
  offset <- mean(-x*sin(t) + y*cos(t))
  # return
  c(t, offset)
}
