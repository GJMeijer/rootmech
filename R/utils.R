#' Determine rounded plot axis limits
#'
#' @description
#' Determine plot axis limits based on input values, preferred tick step sizes
#' and a maximum number of ticks
#'
#' @md
#' @param x values on axis that should be within the plot range
#' @param lim two-parameter vector with fixed axis limits. If `NA`, an axis
#'   limit can be freely chosen
#' @param step array with possible step sizes. These values are relative to
#'   the order of magnitude of the difference in input `x`, i.e. relative to
#'   `10^floor(log10(max(x) - min(x)))`
#' @param ticks the maximum number of ticks that is acceptable
#' @return a list with fields:
#'   * `lim`: two-parameter vector with min and max limit of the plot range
#'   * `breaks`: vector with plot breakpoints
#'
#' @export
#' @examples
#' # data
#' x <- 10*rnorm(100)
#'
#' # basic range
#' z <- round_range(x)
#' plot(x, seq(length(x)), xlim = z$lim, xaxs = "i", xaxt = "n")
#' axis(side = 1, at = z$breaks)
#'
#' # fix the lower bound
#' z <- round_range(x, lim = c(-100, NA))
#' plot(x, seq(length(x)), xlim = z$lim, xaxs = "i", xaxt = "n")
#' axis(side = 1, at = z$breaks)
#'
round_range <- function(
    x,
    lim = c(NA, NA),
    step = c(0.1, 0.2, 0.5, 1, 2, 5),
    ticks = 7
) {
  # get all limits
  lims <- c(
    min(min(x), lim[1], na.rm = TRUE),
    max(max(x), lim[2], na.rm = TRUE)
  )
  # order of magnitude of difference
  oom <- floor(log10(diff(lims)))
  # round scaled upper and lower values using potential step sizes (step)
  z0 <- floor((lims[1]/10^oom)/step)*step
  z1 <- ceiling((lims[2]/10^oom)/step)*step
  # calculate number of ticks with distance <step> required for each option
  ntick <- ceiling((z1 - z0)/step)
  ntick[ntick > ticks] <- 0
  # pick option that has the most ticks but number still smaller than <n>
  i <- which.max(ntick)
  # set range
  range <- c(z0[i], z1[i])*10^oom
  # substitute fixed values
  range[!is.na(lim)] <- lim[!is.na(lim)]
  # return range and ticks
  list(
    lim = range,
    breaks = seq(range[1], range[2], step[i]*10^oom)
  )
}


#' Round a double to a number of significant digits and convert to character
#'
#' @description
#' Round a number to a number of significant digits and convert to character
#' string. Always returns `digits` number of significant digits, also when
#' there are trailing zeros. if `round = TRUE`, the number is rounded to
#' `digits` number of decimals instead.
#'
#' @param x numeric value to convert to character string
#' @param digits number of significant digits
#' @param round if `TRUE`, round `x` to `digits` number of decimals. If
#'   `FALSE`, round to `digits` number of significant digits
#' @return character string
#' @export
#' @examples
#' numeric2character(pi, digits = 4)
#' numeric2character(3, digits = 2, round = TRUE)
#' numeric2character(sqrt(5)/100000, digits = 3)
#'
numeric2character <- function(x, digits = 3, round = FALSE) {
  if (round == TRUE) {
    xr <- round(x, digits = digits)
    digits2 <- digits + floor(log10(abs(x))) + 1
    purrr::map2_chr(
      xr, digits2,
      function(x, n) {
        if (n <= 0) {
          if (digits > 0) {
            paste0(c("0.", rep(0, digits)), collapse = "")
          } else {
            "0"
          }
        } else {
          gsub(
            "\\.$", "",
            formatC(
              x,
              digits = n,
              format = "fg",
              flag = "#"
            )
          )
        }
      }
    )
  } else {
    gsub(
      "\\.$", "",
      formatC(
        signif(x, digits = digits),
        digits = digits,
        format = "fg",
        flag = "#"
      )
    )
  }
}


#' Generate confidence ellipse for 2-dimensional data
#'
#' @md
#' @description
#' Generate coordinates for a x% confidence ellipse.
#'
#' @param x,y vectors with data
#' @param level confidence level
#' @param n number of points to generate on ellipse circumference
#' @return list with fields:
#'   * `ellipse`: dataframe with `x` and `y` coordinates of the
#'     confidence ellipse
#'   * `major` dataframe with `x` and `y` coordinates of the major
#'     axis within the ellipse
#'   * `minor` dataframe with `x` and `y` coordinates of the minor
#'     axis within the ellipse
#'   * `mean`: dataframe with `x` and `y` coordinates of the mean
#' @export
#' @examples
#' x <- rnorm(100, 2, 1)
#' y <- rnorm(100, 4, 0.5) + 0.5*x
#' df <- confidence_ellipse(x, y)
#' plot(x, y)
#' lines(df$ellipse$x, df$ellipse$y, col = "red")
#' lines(df$major$x, df$major$y, col = "red")
#' lines(df$minor$x, df$minor$y, col = "red")
#' lines(df$mean$x, df$mean$y, type = "p", col = "red")
#'
confidence_ellipse <- function(x, y, level = 0.95, n = 181) {
  # means
  mux <- mean(x)
  muy <- mean(y)
  # variance-covariance matrix
  sdx <- stats::sd(x)
  sdy <- stats::sd(y)
  corxy <- stats::cor(x, y)
  Sigma <- matrix(
    c(sdx^2, corxy*sdx*sdy, corxy*sdx*sdy, sdy^2),
    nrow = 2
  )
  # eigenvalues and vectors
  eig <- eigen(Sigma)
  # generate ellipsoid
  radii <- sqrt(eig$values)*sqrt(stats::qchisq(level, df = 2))
  # generate ellipsoid
  theta <- seq(-pi, pi, l = n)
  df1 <- data.frame(
    x = radii[1]*cos(theta),
    y = radii[2]*sin(theta)
  )
  # rotate ellipsoid, and then move
  angle <- atan2(eig$vectors[1, 2], eig$vectors[1, 1])
  df2 <- data.frame(
    x = mux + df1$x*cos(angle) + df1$y*sin(angle),
    y = muy - df1$x*sin(angle) + df1$y*cos(angle)
  )
  # major and minor axes
  df_major <- data.frame(
    x = mux + c(-radii[1], radii[1])*cos(angle),
    y = muy - c(-radii[1], radii[1])*sin(angle)
  )
  df_minor <- data.frame(
    x = mux + c(-radii[2], radii[2])*sin(angle),
    y = muy + c(-radii[2], radii[2])*cos(angle)
  )
  # return
  list(
    mean = data.frame(x = mux, y = muy),
    covariance = Sigma,
    ellipse = df2,
    major = df_major,
    minor = df_minor
  )
}
