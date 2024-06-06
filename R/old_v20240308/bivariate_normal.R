# Bivariate normal distribution functions
# 07/08/2023 - G. J. Meijer


#' Variance-covariance matrix for bivariate normal distribution
#'
#' @description
#' Generate variance-covariance matrix for a bivariate normal distribution,
#' based on known standard deviations and correlation coefficient.
#'
#' @param sd vector with standard deviations
#' @param rho correlation coefficient
#' @return 2*2 matrix
#' @export
#' @examples
#' bivariate_normal_variancematrix(sd = c(2, 3), rho = 0.4)
#'
bivariate_normal_variancematrix <- function(sd = c(1, 1), rho = 0.0) {
  matrix(
    c(sd[1]^2, rho*sd[1]*sd[2], rho*sd[1]*sd[2], sd[2]^2),
    nrow = 2,
    ncol = 2
  )
}


#' Generate random bivariate normal data
#'
#' @description
#' Generate random bivariate normally distributed data, with known means,
#' standard deviations and a correlation coefficient.
#'
#' @param n number of samples to draw
#' @param mu vector with means for x and y
#' @param sd vector with standard deviations for x and y
#' @param rho correlation coefficient for the bivariate normal distribution
#' @return dataframe with `x` and `y` fields for randomly generated data
#' @export
#' @examples
#' # example1: mean = 0, sd = 1
#' x <- bivariate_normal_random(100)
#' plot(x[, 1], x[, 2])
#'
#' # correlated data, with offset
#' df <- bivariate_normal_random(100, mu = c(3, 6), sd = c(1, 2), rho = -0.5)
#' plot(df$x, df$y)
#'
bivariate_normal_random <- function(
    n,
    mu = c(0, 0),
    sd = c(1, 1),
    rho = 0
) {
  # generate variance-covariance matrix
  Sigma <- bivariate_normal_variancematrix(sd = sd, rho = rho)
  # random data, without mean taken into account
  M <- matrix(stats::rnorm(n*2), ncol = 2) %*% chol(Sigma)
  # add means, and return
  data.frame(
    x = M[, 1] + mu[1],
    y = M[, 2] + mu[2]
  )
}


#' Generate prediction interval for bivariate normal distribution
#'
#' @description
#' Generate the coordinates of the prediction interval ellipsoid for a
#' bivariate normal distribution.
#'
#' @inheritParams bivariate_normal_random
#' @param confidence value of confidence required, centred around the mean (
#'   i.e. setting `confidence = 0.95` gives the ellipsoid that contains a random
#'   point with 95\% confidence, with tails of 2.5\%)
#' @param n number of points on ellipsoid
#' @return dataframe with coordinates `x` and `y` for each confidence level
#'   in `confidence`
#' @export
#' @examples
#' # parameters
#' mu <- c(2, 5)
#' sd <- c(1, 0.2)
#' rho <- 0.9
#' confidence <- 0.95
#'
#' # generate some random data
#' df <- bivariate_normal_random(100, mu = mu, sd = sd, rho = rho)
#' # generate confidence interval
#' dc <- bivariate_normal_predictioninterval(mu = mu, sd = sd, rho = rho)
#'
#' # plot
#' plot(df$x, df$y)
#' lines(dc$x, dc$y)
#'
bivariate_normal_predictioninterval <- function(
    mu = c(0, 0),
    sd = c(1, 1),
    rho = 0,
    confidence = 0.95,
    n = 360
) {
  # variance-covariance matrix
  Sigma <- bivariate_normal_variancematrix(sd = sd, rho = rho)
  # eigenvalue decomposition
  eig <- eigen(Sigma)
  # multiplier - chi-squared distribution
  chisq <- stats::qchisq(confidence, df = 2)
  # angles
  theta <- seq(0, 2*pi, l = n + 1)
  # combinations
  df <- expand_grid_df(
    data.frame(theta = theta),
    data.frame(confidence = confidence, chisq = chisq)
  )
  # bases of ellipsiod
  B1 <- sqrt(eig$values[1]*df$chisq)*cos(df$theta)
  B2 <- sqrt(eig$values[2]*df$chisq)*sin(df$theta)
  B <- matrix(c(B1, B2), ncol = 2, byrow = FALSE)
  # rotate ellipsoid using eigenvectors
  R <- B %*% t(eig$vectors)
  df$x <- R[, 1] + mu[1]
  df$y <- R[, 2] + mu[2]
  # return
  df[, c("confidence", "x", "y")]
}


#' Probability density for a bivariate normal distribution
#'
#' @inheritParams bivariate_normal_random
#' @param x n*2 matrix with observations
#' @param log if `log = TRUE`, return the log-transformed probability density
#'   rather than the probability density
#' @return vector with probability density for each observation
#' @export
#' @examples
#' # input
#' mu <- c(5, 2)
#' sd <- c(2, 1)
#' rho <- -0.6
#'
#' # generate grid
#' x <- mu[1] + 3*sd[1]*seq(-1, 1, l = 101)
#' y <- mu[2] + 3*sd[2]*seq(-1, 1, l = 101)
#' xy <- expand.grid(x = x, y = y, KEEP.OUT.ATTRS = FALSE)
#'
#' # calculate probability density
#' d <- matrix(
#'   bivariate_normal_density(xy, mu = mu, sd = sd, rho = rho),
#'   nrow = length(x),
#'   ncol = length(y)
#' )
#'
#' # plot contour plot
#' contour(x, y, d)
#'
bivariate_normal_density <- function(
    x,
    mu = c(0, 0),
    sd = c(1, 1),
    rho = 0,
    log = FALSE
) {
  # normalise data by mean and sds
  xn1 <- (x[, 1] - mu[1])/sd[1]
  xn2 <- (x[, 2] - mu[2])/sd[2]
  # temporary variable
  z <- xn1^2 + xn2^2 - 2*rho*xn1*xn2
  # return
  if (log == TRUE) {
    -z/(2*(1 - rho^2)) - log(2*pi) - log(sd[1]) - log(sd[2]) - 0.5*log(1 - rho^2)
  } else {
    exp(-z/(2*(1 - rho^2)))/(2*pi*sd[1]*sd[2]*sqrt(1 - rho^2))
  }
}


#' Major and minor axes for bivariate normal confidence ellipsoid
#'
#' @description
#' Generate coordinates to plot the minor and major axes of the
#' bivariate normal confidence ellipsoid
#'
#' @inheritParams bivariate_normal_random
#' @param confidence confidence value of the ellipsoid
#' @param n number of points to use for each axis
#' @return dataframe with axis type (fields `axis`, can be `major` or
#'   `minor`), and the coordiante fields `x` and `y`
#' @export
#' @examples
#' bivariate_normal_axes(rho = 0.5, n = 10)
#'
bivariate_normal_axes <- function(
    mu = c(0, 0),
    sd = c(1, 1),
    rho = 0,
    confidence = 0.95,
    n = 100
) {
  # multiplier - chi-squared distribution
  chisq <- stats::qchisq(confidence, df = 2)
  # major axis
  df_major <- data.frame(
    axis = "major",
    x = mu[1] + sqrt(chisq)*sd[1]*seq(-1, 1, l = n),
    y = mu[2] + sqrt(chisq)*rho*sd[2]*seq(-1, 1, l = n)
  )
  # minor axis
  df_minor <- data.frame(
    axis = "minor",
    x = mu[1] + sqrt(chisq)*rho*sd[1]*seq(-1, 1, l = n),
    y = mu[2] + sqrt(chisq)*sd[2]*seq(-1, 1, l = n)
  )
  # return
  rbind(df_major, df_minor)
}


#' Convert bivariate normal values to weibull
#'
#' @description
#' Convert bivariate normally distributed values to weibull-distributed
#' values. This is done by setting the cumulative densities of both
#' distributions the same
#'
#' @inheritParams bivariate_normal_random
#' @param x,y positions in bivariate normal coordinates
#' @param shape_x,shape_y Weibull shape parameters for x and y
#' @param scale_x,scale_y Weibull scale parameters for x and y
#' @return dataframe with Weibull positions `x` and `y`
#' @export
#' @examples
#' df_norm <- bivariate_normal_predictioninterval(rho = 0.5)
#' plot(df_norm$x, df_norm$y, "l")
#'
#' df_weibull <- bivariate_normal_toweibull(
#'   df_norm$x, df_norm$y, 4, 6, scale_y = 4)
#' plot(df_weibull$x, df_weibull$y, "l")
bivariate_normal_toweibull <- function(
    x,
    y,
    shape_x,
    shape_y,
    mu = c(0, 0),
    sd = c(1, 1),
    scale_x = 1/gamma(1 + 1/shape_x),
    scale_y = 1/gamma(1 + 1/shape_y)
) {
  # cumulative density function
  px <- stats::pnorm(x, mean = mu[1], sd = sd[1])
  py <- stats::pnorm(y, mean = mu[2], sd = sd[2])
  # convert to weibull, and return dataframe
  data.frame(
    x = stats::qweibull(px, shape = shape_x, scale = scale_x),
    y = stats::qweibull(py, shape = shape_y, scale = scale_y)
  )
}
