#' Fit a power-law curve to a series of (x, y) data
#'
#' @description
#' Fit a power-law curve to a serries of weighted (x, y) data. Can use
#' different assumptions for the desciption of intra-diameter variation
#'
#' @md
#' @param x measured x-values (e.g. root diameters). These values are assumed
#'   to be already normalised by a reference value to ensure a unitless
#'   parameter.
#' @param y measured y-values (e.g. root tensile strength)
#' @param model model used to describe intra-diameter variation. Can be
#'   one of the following:
#'   * `normal_strength`: normal distribution - constant standard deviation
#'     with x (homoscedastic)
#'   * `normal_force`: normal distribution - constant standard deviation with
#'     x in terms of forces (y*(pi/4*x^2))
#'   * `normal_scaled`: normal distribution - standard deviation scales with
#'     mean
#'   * `normal_sdpower`: normal distribution - seperate power-law describin the
#'     mean and standard deviation
#'   * `lognormal`: lognormal distribution of ratio y/mean
#'   * `weibull`: weibull distribution; mean = power-law, constant shape
#'     parameter
#'   * `gamma`: gamma distribution; mean = power law, constant shape parameter
#'   * `logistic`: logistic distribution; mean = power-law, scale parameter
#'     scales with mean
#'   * `uniform`: uniform distribution; mean = power-law. width of uniform
#'     distribution scales with mean
#' @param weights vector with weighting for each observation (default: 1)
#' @param correction only for `model == "lognormal"`. If `TRUE`, means
#'   correction is applied. If `FALSE`, the result from linear regression
#'   is returned instead (not correctly implementing the mean of log etc)
#' @param method solving method. Can be `newton` for gradient-descent (using
#'   `rootSolve::multiroot()`), `bisection` for a bisection algorithm
#'   (using `stats::uniroot()`) or `chull` for a method involving convex
#'   hulls (using, `grDevices::chull()`, only available for
#'   `model == "uniform"`)). Not all methods are available for each model,
#'   see individual documentation for fitting functions for each model
#' @param range initial search range around best guess, when using a
#'   bisection algorithm (`method == "bisection`). `range` is added to the
#'   best guess to provide the bounds of the initial search window
#' @return list with fitting results. Will contain fields
#'   * `loglikelihood` fitted loglikelihood
#'   * `multiplier` fitted power-law multiplier
#'   * `exponent`: fitted power-law exponent
#'   * `scale`: scale parameter (for models `weibull`, `gamma` and
#'     `logistic`)
#'   * `sd`: standard deviation (for model `normal_strength`)
#'   * `sd_multiplier`: standard deviation power-law multiplier (for models
#'     `normal_sdpower`, `normal_scaled` and `normal_force`)
#'   * `sd_exponent`: standard deviation power-law exponent (for model
#'     `normal_sdpower` and `normal_force`)
#'   * `width`: width of uniform distribution, at x = 1 (for model `uniform`)
#' @export
#' @examples
#' # generate some data
#' y0 <- 20
#' beta <- -0.5
#' kappa <- 4
#' lambda <- 1/gamma(1 + 1/kappa)
#' x <- seq(1, 8, l = 51)
#' y <- y0*x^beta*rweibull(length(x), kappa, lambda)
#'
#' power_fit(x, y, model = "normal_strength")
#' power_fit(x, y, model = "normal_force")
#' power_fit(x, y, model = "normal_scaled")
#' power_fit(x, y, model = "normal_sdpower")
#' power_fit(x, y, model = "lognormal")
#' power_fit(x, y, model = "weibull")
#' power_fit(x, y, model = "gamma")
#' power_fit(x, y, model = "logistic")
#' power_fit(x, y, model = "uniform")
#'
power_fit <- function(
    x,
    y,
    model = "normal",
    weights = rep(1, length(x)),
    correction = TRUE,
    method = "bisection",
    range = c(-1, 1)
) {
  if (tolower(model) == "weibull") {
    # Weibull distribution - scale parameter scales with mean
    power_weibull_fit(
      x,
      y,
      weights = weights
    )
  } else if (tolower(model) %in% c("normal_strength", "normal")) {
    # Normal distribution - homoscedatic residuals in terms of strength
    power_normal_unscaled_fit(
      x,
      y,
      weights = weights,
      method = method,
      range = range
    )
  } else if (tolower(model) == "normal_force") {
    # Normal distribution - homoscedatic residuals in terms of force
    ft <- power_normal_unscaled_fit(
      x,
      y*(pi/4*x^2),
      weights = weights,
      method = method,
      range = range
    )
    ft$multiplier <- 4/pi*ft$multiplier
    ft$exponent <- ft$exponent - 2
    ft$sd_multiplier <- 4/pi*ft$sd
    ft$sd_exponent <- -2
    ft <- ft[-4]
  } else if (tolower(model) == "normal_scaled") {
    # Normal distribution - standard deviation scales with mean
    power_normal_scaled_fit(
      x,
      y,
      weights = weights,
      method = method,
      range = range
    )
  } else if (tolower(model) == "normal_sdpower") {
    # Normal distribution - seperate power-laws for mean and standard deviation
    power_normal_sdpower_fit(
      x,
      y,
      weights = weights
    )
  } else if (tolower(model) %in% c("lognormal", "log-normal")) {
    # Log-normal distribution - heteroscedatic log-strength residuals
    ft <- power_lognormal_fit(
      x,
      y,
      weights = weights
    )
    if (correction == FALSE) {
      ft$multiplier <- ft$multiplier/exp(0.5*ft$sigmaL^2)
    }
    ft
  } else if (tolower(model) == "gamma") {
    # Gamma distribution - scale parameter scales with mean
    power_gamma_fit(
      x, y,
      weights = weights,
      method = method,
      range = range
    )
  } else if (tolower(model) == "uniform") {
    # Uniform distribution - width scales with mean
    power_uniform_fit(
      x,
      y,
      weights = weights
    )
  } else if (tolower(model) == "logistic") {
    # Logistic distribution - scale parameter scales with mean
    power_logistic_fit(
      x,
      y,
      weights = weights
    )
  } else {
    stop("'model' not recognised")
  }
}

#' Prediction interval
#'
#' @description
#' A short description...
#'
#' @param x a
#' @param model a
#' @param multiplier a
#' @param exponent a
#' @param level a
#' @param ... a
#' @return bla
#' @export
#' @examples
#' # generate some data
#' y0 <- 20
#' beta <- -0.5
#' kappa <- 4
#' lambda <- 1/gamma(1 + 1/kappa)
#' x <- seq(1, 8, l = 101)
#' y <- y0*x^beta*rweibull(length(x), kappa, lambda)
#'
#' # weibull fit
#' ftw <- power_fit(x, y, "weibull")
#'
#' # generate prediction interval
#' dfw <- power_predictioninterval(
#'   seq(min(x), max(x), l = 251), "weibull",
#'   ftw$multiplier, ftw$exponent, shape = ftw$shape
#' )
#'
#' # plot
#' plot(x, y)
#' lines(dfw$x, dfw$y, col = "red")
#' lines(dfw$x, dfw$ymin, col = "blue")
#' lines(dfw$x, dfw$ymax, col = "blue")
#'
#' # normal fit
#' ftn <- power_fit(x, y, "normal_strength")
#'
#' # generate prediction interval
#' dfn <- power_predictioninterval(
#'   seq(min(x), max(x), l = 251), "normal_strength",
#'   ftn$multiplier, ftn$exponent, sd = ftn$sd
#' )
#'
#' # plot
#' plot(x, y)
#' lines(dfn$x, dfn$y, col = "red")
#' lines(dfn$x, dfn$ymin, col = "blue")
#' lines(dfn$x, dfn$ymax, col = "blue")
power_predictioninterval <- function(
  x,
  model,
  multiplier,
  exponent,
  level = 0.95,
  ...
) {
  # mean
  y <- multiplier*x^exponent
  # cumulative lower and upper
  cum <- 0.5 + 0.5*level*c(-1, 1)
  # get upper and lower prediction intervals
  if (tolower(model) == "weibull") {
    # Weibull distribution - scale parameter scales with mean
    f <- stats::qweibull(cum, shape, 1/gamma(1 + 1/shape))
    ymin <- y*f[1]
    ymax <- y*f[2]
  } else if (tolower(model) %in% c("normal_strength", "normal")) {
    # Normal distribution - homoscedatic residuals in terms of strength
    f <- stats::qnorm(cum, 0, sd)
    ymin <- y + f[1]
    ymax <- y + f[2]
  } else if (tolower(model) == "normal_force") {
    # Normal distribution - homoscedatic residuals in terms of force
    sd_exponent <- -2
    f <- stats::qnorm(cum, 0, 1)
    ymin <- y + f[1]*sd_multiplier*x^sd_exponent
    ymax <- y + f[2]*sd_multiplier*x^sd_exponent
  } else if (tolower(model) == "normal_scaled") {
    # Normal distribution - standard deviation scales with mean
    f <- stats::qnorm(cum, 1, sd_multiplier/multiplier)
    ymin <- y*f[1]
    ymax <- y*f[2]
  } else if (tolower(model) == "normal_sdpower") {
    # Normal distribution - seperate power-laws for mean and standard deviation
    f <- stats::qnorm(cum, 0, sd_multiplier)
    ymin <- y + f[1]*x^sd_exponent
    ymax <- y + f[2]*x^sd_exponent
  } else if (tolower(model) %in% c("lognormal", "log-normal")) {
    # Log-normal distribution - heteroscedatic log-strength residuals
    f <- stats::qlnorm(cum, -0.5*sigmaL^2, sigmaL)
    ymin <- y*f[1]
    ymax <- y*f[2]
  } else if (tolower(model) == "gamma") {
    # Gamma distribution - scale parameter scales with mean
    ymin <- stats::qgamma(cum[1], shape = shape, scale = 1/shape*multiplier*x^exponent)
    ymax <- stats::qgamma(cum[2], shape = shape, scale = 1/shape*multiplier*x^exponent)
  } else if (tolower(model) == "uniform") {
    # Uniform distribution - width scales with mean
    ymin <- (multiplier - 0.5*confidence*width)*x^exponent
    ymax <- (multiplier + 0.5*confidence*width)*x^exponent
  } else if (tolower(model) == "logistic") {
    # Logistic distribution - scale parameter scales with mean
    ymin <- stats::qlogis(cum[1], multiplier*x^exponent, scale*x^exponent)
    ymax <- stats::qlogis(cum[2], multiplier*x^exponent, scale*x^exponent)
  } else {
    stop("'model' not recognised")
  }
  # return dataframe
  data.frame(
    x = x,
    y = y,
    ymin = ymin,
    ymax = ymax
  )
}


#' Confidence interval
#'
#' @description
#' A short description...
#'
#' @param x,y a
#' @param model a
#' @param weights a
#' @param level a
#' @param xp a
#' @param multiplier a
#' @param exponent a
#' @param sd a
#' @param sd_multiplier a
#' @param sd_exponent a
#' @param shape a
#' @param scale a
#' @param sdlog a
#' @return bla
#' @export
#' @examples
#' # generate some data
#' y0 <- 20
#' beta <- -0.5
#' kappa <- 3
#' lambda <- 1/gamma(1 + 1/kappa)
#' x <- seq(1, 8, l = 51)
#' y <- y0*x^beta*rweibull(length(x), kappa, lambda)
#'
#' df <- power_confidenceinterval(x, y, "lognormal")
#'
#' plot(x, y)
#' lines(df$x, df$y, col = "red")
#' lines(df$x, df$ymin, col = "blue")
#' lines(df$x, df$ymax, col = "blue")
#'
#' library("ggplot2")
#' ggplot() +
#'   geom_point(aes(x = x, y = y)) +
#'   geom_smooth(aes(x = x, y = y), method = "lm") +
#'   scale_y_log10() +
#'   geom_ribbon(
#'     data = df,
#'     aes(x = x, ymin = ymin, ymax = ymax),
#'     alpha = 0.25, fill = "red"
#'   )
#'
power_confidenceinterval <- function(
    x,
    y,
    model,
    weights = rep(1, length(x)),
    level = 0.95,
    xp = seq(min(x, na.rm = TRUE), max(x, na.rm = TRUE), l = 101),
    multiplier = NULL,
    exponent = NULL,
    sd = NULL,
    sd_multiplier = NULL,
    sd_exponent = NULL,
    shape = NULL,
    scale = NULL,
    sdlog = NULL
) {
  # calculate 2nd derivative of loglikelihood function
  if (tolower(model) == "weibull") {
    # Weibull distribution - scale parameter scales with mean
    if (is.null(multiplier) | is.null(exponent) | is.null(shape)) {
      ft <- power_weibull_fit(x, y, weights = weights)
      multiplier <- ft$multiplier
      exponent <- ft$exponent
      shape <- ft$shape
    }
    d2logL_dpar2 <- power_weibull_loglikelihood(
      c(multiplier, exponent, shape),
      x, y, weights = weights,
      deriv = 2
    )
  } else if (tolower(model) %in% c("normal_strength", "normal")) {
    # Normal distribution - homoscedatic residuals in terms of strength
    f <- stats::qnorm(cum, 0, sd)
    ymin <- y + f[1]
    ymax <- y + f[2]
  } else if (tolower(model) == "normal_force") {
    # Normal distribution - homoscedatic residuals in terms of force
    sd_exponent <- -2
    f <- stats::qnorm(cum, 0, 1)
    ymin <- y + f[1]*sd_multiplier*x^sd_exponent
    ymax <- y + f[2]*sd_multiplier*x^sd_exponent
  } else if (tolower(model) == "normal_scaled") {
    # Normal distribution - standard deviation scales with mean
    f <- stats::qnorm(cum, 1, sd_multiplier/multiplier)
    ymin <- y*f[1]
    ymax <- y*f[2]
  } else if (tolower(model) == "normal_sdpower") {
    # Normal distribution - seperate power-laws for mean and standard deviation
    f <- stats::qnorm(cum, 0, sd_multiplier)
    ymin <- y + f[1]*x^sd_exponent
    ymax <- y + f[2]*x^sd_exponent
  } else if (tolower(model) %in% c("lognormal", "log-normal")) {
    # Log-normal distribution - heteroscedatic log-strength residuals
    if (is.null(multiplier) | is.null(exponent) | is.null(sdlog)) {
      ft <- power_lognormal_fit(x, y, weights = weights)
      multiplier <- ft$multiplier
      exponent <- ft$exponent
      sdlog <- ft$sdlog
    }
    d2logL_dpar2 <- power_lognormal_loglikelihood(
      c(multiplier, exponent, sdlog),
      x, y, weights = weights,
      deriv = 2
    )
  } else if (tolower(model) == "gamma") {
    # Gamma distribution - scale parameter scales with mean
    ymin <- stats::qgamma(cum[1], shape = shape, scale = 1/shape*multiplier*x^exponent)
    ymax <- stats::qgamma(cum[2], shape = shape, scale = 1/shape*multiplier*x^exponent)
  } else if (tolower(model) == "uniform") {
    # Uniform distribution - width scales with mean
    ymin <- (multiplier - 0.5*confidence*width)*x^exponent
    ymax <- (multiplier + 0.5*confidence*width)*x^exponent
  } else if (tolower(model) == "logistic") {
    # Logistic distribution - scale parameter scales with mean
    ymin <- stats::qlogis(cum[1], multiplier*x^exponent, scale*x^exponent)
    ymax <- stats::qlogis(cum[2], multiplier*x^exponent, scale*x^exponent)
  } else {
    stop("'model' not recognised")
  }
  # initialise dataframe
  df <- data.frame(x = xp)
  df$y <- multiplier*df$x^exponent
  # guess variance-covariance matrix from 2nd derivative of loglikelihood
  Sigma <- solve(-d2logL_dpar2)
  # derivatives of power-law function
  dyp_dmultiplier <- df$x^exponent
  dyp_dpower <- multiplier*log(df$x)*df$x^exponent
  # get confidence interval using delta method
  var <- (
    dyp_dmultiplier^2*Sigma[1, 1] +
      2*dyp_dmultiplier*dyp_dpower*Sigma[1, 2] +
      dyp_dpower^2*Sigma[2, 2]
  )
  # multiplier of standard deviation
  mult <- stats::qnorm(0.5 + 0.5*level)  # gaussian distributed
  # confidence interval
  df$ymin <- df$y - mult*sqrt(var)
  df$ymax <- df$y + mult*sqrt(var)
  # return dataframe
  df
}
