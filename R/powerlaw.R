#' Fit a power-law curve to a series of (x, y) data
#'
#' @description
#' Fit a power-law curve to a series of weighted (x, y) data. Can use
#' different assumptions for the description of intra-diameter variation
#'
#' @md
#' @param x measured x-values (e.g. root diameters). These values are assumed
#'   to be already normalised by a reference value to ensure a dimensionless
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
#'   * `normal_sdpower`: normal distribution - separate power-law describing the
#'     mean and standard deviation
#'   * `lognormal`: lognormal distribution of ratio y/mean
#'   * `lognormal_uncorrected`: lognormal distribution with no means correction
#'     (i.e. result of linear fitting on log-transformed data)
#'   * `weibull`: Weibull distribution; mean = power-law, constant shape
#'     parameter
#'   * `gamma`: gamma distribution; mean = power law, constant shape parameter
#'   * `logistic`: logistic distribution; mean = power-law, scale parameter
#'     scales with mean
#'   * `uniform`: uniform distribution; mean = power-law. width of uniform
#'     distribution scales with mean
#'   * `gumbel`: Gumbel distribution; mean = power law, scale parameter
#'     scales with mean
#' @param weights vector with weighting for each observation (default: 1)
#' @param method solving method. Can be
#'   * `newton` for gradient-descent (using
#'   `rootSolve::multiroot()`)
#'   * `bisection` for a bisection algorithm
#'   (using `stats::uniroot()`)
#'   * `chull` for a method involving convex
#'   hulls (using, `grDevices::chull()`. This method is only available for
#'   `model == "uniform"`.
#'
#'   Not all methods are available for each model,
#'   see individual documentation for fitting functions for each model (e.g.
#'   by typing `?powerlaw_fit_weibull` into the R console).
#' @param range initial search range around best guess, when using a
#'   bisection algorithm (`method == "bisection`). `range` is added to the
#'   best guess to provide the bounds of the initial search window
#' @param start optional starting values for fitting. If not defined, a guess
#'   is made
#' @return list with fitting results. Will contain fields
#'   * `loglikelihood` fitted loglikelihood
#'   * `multiplier` fitted power-law multiplier
#'   * `exponent`: fitted power-law exponent
#'   * `sd_multiplier`: standard deviation power-law multiplier (for models
#'     `normal_x`()
#'   * `sd_exponent`: standard deviation power-law exponent (for models
#'     `normal_x`)
#'   * `sdlog`: log-standard deviation (for model `lognormal`)
#'   * `scale`: scale parameter (for models `logistic`, and `gumbel`)
#'   * `shape`: shape parameter (for models `gamma` and `weibull`)
#'   * `width`: width of uniform distribution, at x = 1 (for model `uniform`)
#'   * `ks_distance`: Kolmogorov-Smirnov distance of best fit
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
#' powerlaw_fit(x, y, model = "gamma")
#' powerlaw_fit(x, y, model = "gumbel")
#' powerlaw_fit(x, y, model = "logistic")
#' powerlaw_fit(x, y, model = "lognormal")
#' powerlaw_fit(x, y, model = "normal_strength")
#' powerlaw_fit(x, y, model = "normal_force")
#' powerlaw_fit(x, y, model = "normal_scaled")
#' powerlaw_fit(x, y, model = "normal_sdpower")
#' powerlaw_fit(x, y, model = "uniform")
#' powerlaw_fit(x, y, model = "weibull")
#'
powerlaw_fit <- function(
    x,
    y,
    model = "normal",
    weights = rep(1, length(x)),
    method = "bisection",
    range = c(-1, 1),
    start = NULL
) {
  if (tolower(model) == "gamma") {
    # Gamma distribution - scale parameter scales with mean
    ft <- powerlaw_gamma_fit(
      x,
      y,
      weights = weights,
      method = method,
      range = range,
      start = start
    )
  } else if (tolower(model) == "gumbel") {
    ft <- powerlaw_gumbel_fit(
      x,
      y,
      weights = weights,
      method = method,
      range = range,
      start = start
    )
  } else if (tolower(model) == "logistic") {
    # Logistic distribution - scale parameter scales with mean
    ft <- powerlaw_logistic_fit(
      x,
      y,
      weights = weights,
      start = start
    )
  } else if (tolower(model) %in% c("lognormal", "lognormal_corrected")) {
    # Log-normal distribution (with correction) - heteroscedatic log-strength residuals
    ft <- powerlaw_lognormal_fit(
      x,
      y,
      weights = weights
    )
  } else if (tolower(model) %in% ("lognormal_uncorrected")) {
    # Log-normal distribution - no correction for mean
    ft <- powerlaw_lognormal_uncorrected_fit(
      x,
      y,
      weights = weights
    )
  } else if (tolower(model) %in% c("normal_strength", "normal")) {
    # Normal distribution - homoscedatic residuals in terms of strength
    ft <- powerlaw_normal_fit(
      x,
      y,
      sd_exponent = 0,
      weights = weights,
      method = method,
      range = range,
      start = start
    )
  } else if (tolower(model) == "normal_force") {
    # Normal distribution - homoscedatic residuals in terms of force
    ft <- powerlaw_normal_fit(
      x,
      y,
      sd_exponent = -2,
      weights = weights,
      method = method,
      range = range,
      start = start
    )
  } else if (tolower(model) == "normal_scaled") {
    # Normal distribution - standard deviation scales with mean
    ft <- powerlaw_normal_fit(
      x,
      y,
      sd_exponent = "scaled",
      weights = weights,
      method = method,
      range = range,
      start = start
    )
  } else if (tolower(model) == "normal_sdpower") {
    # Normal distribution - separate power-laws for mean and standard deviation
    ft <- powerlaw_normal_fit(
      x,
      y,
      sd_exponent = NULL,
      weights = weights,
      method = method,
      range = range,
      start = start
    )
  } else if (tolower(model) == "uniform") {
    # Uniform distribution - width scales with mean
    ft <- powerlaw_uniform_fit(
      x,
      y,
      weights = weights
    )
  } else if (tolower(model) == "weibull") {
    # Weibull distribution - scale parameter scales with mean
    ft <- powerlaw_weibull_fit(
      x,
      y,
      weights = weights,
      start = start
    )
  } else {
    stop("'model' not recognised")
  }
  # add KS distance
  ft$ks_distance <- powerlaw_ks(
    x,
    y,
    model,
    ft$multiplier,
    ft$exponent,
    sd_multiplier = if ("sd_multiplier" %in% names(ft)) ft$sd_multiplier else NULL,
    sd_exponent = if ("sd_exponent" %in% names(ft)) ft$sd_exponent else NULL,
    sdlog = if ("sdlog" %in% names(ft)) ft$sdlog else NULL,
    scale = if ("scale" %in% names(ft)) ft$scale else NULL,
    shape = if ("shape" %in% names(ft)) ft$shape else NULL,
    width = if ("width" %in% names(ft)) ft$width else NULL,
    weights = weights
  )$ks_distance
  # return
  ft
}


#' Generate prediction intervals for power law fits
#'
#' @md
#' @description
#' Calculate the upper and lower limits of the prediction interval for a
#' power law fit.
#'
#' @param x a-values at which to generate the interval
#' @param model fitting model, see function `powerlaw_fit()`
#' @param multiplier fitted power law multiplier
#' @param exponent fitted power law exponent
#' @param sd_multiplier fitted standard deviation multiplier
#' @param sd_exponent fitted standard deviation exponent
#' @param sdlog fitted lognormal standard deviation
#' @param shape fitted shape parameter
#' @param scale fitted scale parameter
#' @param width fitted width parameter
#' @param level significance level
#' @param gamma Euler-Mascheroni constant
#' @return dataframe with fields:
#'   * `x`: x-coordinates
#'   * `y`: average (power law)
#'   * `ymin`: lower limit of the prediction interval
#'   * `ymax`: upper limit of the prediction interval
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
#' ftw <- powerlaw_fit(x, y, "weibull")
#' # generate prediction interval
#' dfw <- powerlaw_predictioninterval(
#'   seq(min(x), max(x), l = 251), "weibull",
#'   ftw$multiplier, ftw$exponent, shape = ftw$shape
#' )
#' # plot
#' plot(x, y)
#' lines(dfw$x, dfw$y, col = "red")
#' lines(dfw$x, dfw$ymin, col = "blue")
#' lines(dfw$x, dfw$ymax, col = "blue")
#'
#' # normal fit
#' ftn <- powerlaw_fit(x, y, "normal_strength")
#' # generate prediction interval
#' dfn <- powerlaw_predictioninterval(
#'   seq(min(x), max(x), l = 251), "normal_strength",
#'   ftn$multiplier, ftn$exponent, sd_multiplier = ftn$sd_multiplier
#' )
#' # plot
#' plot(x, y)
#' lines(dfn$x, dfn$y, col = "red")
#' lines(dfn$x, dfn$ymin, col = "blue")
#' lines(dfn$x, dfn$ymax, col = "blue")
#'
#' # gumbel fit
#' ftg <- powerlaw_fit(x, y, "gumbel")
#' # generate prediction interval
#' dfg <- powerlaw_predictioninterval(
#'   seq(min(x), max(x), l = 251), "gumbel",
#'   ftg$multiplier, ftg$exponent, scale = ftg$scale
#' )
#' # plot
#' plot(x, y)
#' lines(dfg$x, dfg$y, col = "red")
#' lines(dfg$x, dfg$ymin, col = "blue")
#' lines(dfg$x, dfg$ymax, col = "blue")
#'
powerlaw_predictioninterval <- function(
  x,
  model,
  multiplier,
  exponent,
  sd_multiplier = NA,
  sd_exponent = 0,
  sdlog = NA,
  shape = NA,
  scale = NA,
  width = NA,
  level = 0.95,
  gamma = 0.57721566490153286060651209008240243104215933593992
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
    f <- stats::qnorm(cum, 0, sd_multiplier)
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
    # Normal distribution - separate power-laws for mean and standard deviation
    f <- stats::qnorm(cum, 0, sd_multiplier)
    ymin <- y + f[1]*x^sd_exponent
    ymax <- y + f[2]*x^sd_exponent
  } else if (tolower(model) %in% c("lognormal", "log-normal")) {
    # Log-normal distribution - heteroscedatic log-strength residuals
    f <- stats::qlnorm(cum, -0.5*sdlog^2, sdlog)
    ymin <- y*f[1]
    ymax <- y*f[2]
  } else if (tolower(model) %in% c("lognormal_uncorrected")) {
    f <- stats::qlnorm(cum, 0, sdlog)
    ymin <- y*f[1]
    ymax <- y*f[2]
  } else if (tolower(model) == "gamma") {
    # Gamma distribution - scale parameter scales with mean
    ymin <- stats::qgamma(cum[1], shape = shape, scale = 1/shape*multiplier*x^exponent)
    ymax <- stats::qgamma(cum[2], shape = shape, scale = 1/shape*multiplier*x^exponent)
  } else if (tolower(model) == "uniform") {
    # Uniform distribution - width scales with mean
    ymin <- y - 0.5*level*width*x^exponent
    ymax <- y + 0.5*level*width*x^exponent
  } else if (tolower(model) == "logistic") {
    # Logistic distribution - scale parameter scales with mean
    ymin <- stats::qlogis(cum[1], multiplier*x^exponent, scale*x^exponent)
    ymax <- stats::qlogis(cum[2], multiplier*x^exponent, scale*x^exponent)
  } else if (tolower(model) == "gumbel") {
    # Gumbel distribution
    f <- multiplier - gamma*scale - scale*log(-log(cum))
    ymin <- f[1]*x^exponent
    ymax <- f[2]*x^exponent
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


#' Estimate variance-covariance matrix for power law fit
#'
#' @description
#' Estimate the variance-covariance matrix for a power law fit,
#' either using bootstrapping or the loglikelihood estimation
#'
#' @inheritParams powerlaw_predictioninterval
#' @inheritParams powerlaw_fit
#' @param method if `fisher`, the covariance matrix is estimated as the
#'   inverse of the negative second partial derivative of the loglikelihood
#'   function with respect to fitting parameters. If `bootstrap`, the covariance
#'   matrix is estimated using bootstrapping
#' @param n number of bootstrap samples
#' @return the variance-covariance matrix, with names rows and columns
#'   according to the fitting parameter
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
#' ft <- powerlaw_fit(x, y, "weibull")
#' powerlaw_covariancematrix(
#'   x, y,
#'   "weibull",
#'   method = "fisher",
#'   multiplier = ft$multiplier,
#'   exponent = ft$exponent,
#'   shape = ft$shape
#' )
#'
#' ft <- powerlaw_fit(x, y, "lognormal")
#' powerlaw_covariancematrix(
#'   x, y,
#'   "lognormal",
#'   method = "fisher",
#'   multiplier = ft$multiplier,
#'   exponent = ft$exponent,
#'   sdlog = ft$sdlog
#' )
#'
#' ft <- powerlaw_fit(x, y, "lognormal_uncorrected")
#' powerlaw_covariancematrix(
#'   x, y,
#'   "lognormal_uncorrected",
#'   method = "fisher",
#'   multiplier = ft$multiplier,
#'   exponent = ft$exponent,
#'   sdlog = ft$sdlog
#' )
#'
powerlaw_covariancematrix <- function(
    x,
    y,
    model,
    method = "fisher",
    multiplier = NULL,
    exponent = NULL,
    sd_multiplier = NULL,
    sd_exponent = NULL,
    sdlog = NULL,
    shape = NULL,
    scale = NULL,
    width = NULL,
    weights = rep(1, length(x)),
    n = 500
) {
  # change method to bootstrap, in case of uniform
  if ((method == "fisher") & (model == "uniform")) {
    method <- "bootstrap"
  }
  # names for results fields
  if (tolower(model) %in% c("gamma", "weibull")) {
    fields <- c("multiplier", "exponent", "shape")
  } else if (tolower(model) %in% c("gumbel", "logistic")) {
    fields <- c("multiplier", "exponent", "scale")
  } else if (tolower(model) %in% c("lognormal", "lognormal_corrected", "lognormal_uncorrected")) {
    fields <- c("multiplier", "exponent", "sdlog")
  } else if (tolower(model) %in% c("normal_strength", "normal", "normal_force", "normal_scaled")) {
    fields <- c("multiplier", "exponent", "sd_multiplier")
  } else if (tolower(model) == "normal_sdpower") {
    fields <- c("multiplier", "exponent", "sd_multiplier", "sd_exponent")
  } else if (tolower(model) == "uniform") {
    fields <- c("multiplier", "exponent", "width")
  } else {
    stop("'model' not recognised")
  }
  # calculate variance-covariance matrix
  if (tolower(method) == "fisher") {
    # fisher information - inverse of negative 2nd derivative of loglikelihood ####
    if (tolower(model) == "gamma") {
      # Gamma distribution - scale parameter scales with mean
      J2 <- powerlaw_gamma_loglikelihood(
        c(multiplier, exponent, shape),
        x, y, weights = weights,
        deriv = 2
      )
    } else if (tolower(model) == "gumbel") {
      # Gumbel distribution
      J2 <- powerlaw_gumbel_loglikelihood(
        c(multiplier, exponent, scale),
        x, y, weights = weights,
        deriv = 2
      )
    } else if (tolower(model) == "logistic") {
      # Logistic distribution - scale parameter scales with mean
      J2 <- powerlaw_logistic_loglikelihood(
        c(multiplier, exponent, scale),
        x, y, weights = weights,
        deriv = 2
      )
    } else if (tolower(model) %in% c("lognormal")) {
      # Log-normal distribution - heteroscedatic log-strength residuals
      J2 <- powerlaw_lognormal_loglikelihood(
        c(multiplier, exponent, sdlog),
        x, y, weights = weights,
        deriv = 2
      )
    } else if (tolower(model) %in% c("lognormal_uncorrected")) {
      # log-normal distribution (without means-correction)
      J2 <- powerlaw_lognormal_uncorrected_loglikelihood(
        c(multiplier, exponent, sdlog),
        x, y, weights = weights,
        deriv = 2
      )
    } else if (tolower(model) %in% c("normal_strength", "normal")) {
      # Normal distribution - homoscedatic residuals in terms of strength
      J2 <- powerlaw_normal_freebeta_loglikelihood(
        c(multiplier, exponent, sd_multiplier),
        x, y, weights = weights,
        delta = 0,
        deriv = 2
      )
    } else if (tolower(model) == "normal_force") {
      # Normal distribution - homoscedatic residuals in terms of force
      J2 <- powerlaw_normal_freebeta_loglikelihood(
        c(multiplier, exponent, sd_multiplier),
        x, y, weights = weights,
        delta = -2,
        deriv = 2
      )
    } else if (tolower(model) == "normal_scaled") {
      # Normal distribution - standard deviation scales with mean
      J2 <- powerlaw_normal_linkedbetadelta_loglikelihood(
        c(multiplier, exponent, sd_multiplier),
        x, y, weights = weights,
        deriv = 2
      )
    } else if (tolower(model) == "normal_sdpower") {
      # Normal distribution - separate power-laws for mean and standard deviation
      J2 <- powerlaw_normal_freebetadelta_loglikelihood(
        c(multiplier, exponent, sd_multiplier, sd_exponent),
        x, y, weights = weights,
        deriv = 2
      )
    } else if (tolower(model) == "uniform") {
      # Uniform distribution - width scales with mean
      # undefined - bootstrap used - change method in first lines of function
    } else if (tolower(model) == "weibull") {
      # Weibull distribution - scale parameter scales with mean
      J2 <- powerlaw_weibull_loglikelihood(
        c(multiplier, exponent, shape),
        x, y, weights = weights,
        deriv = 2
      )
    } else {
      stop("'model' not recognised")
    }
    # return variance-covariance matrix
    Sigma <- solve(-J2)
  } else if (tolower(method == "bootstrap")) {
    # bootstrap ####
    # generate list of simulations
    S <- matrix(NA, n, length(fields))
    for (i in 1:n) {
      j <- sample(seq(length(x)), length(x), replace = TRUE)
      S[i, ] <- unlist(powerlaw_fit(
        x[j], y[j],
        model = model,
        weights = weights[j]
      )[fields])
    }
    # calculate variance-covariance matrix
    Sigma <- stats::cov(S)
  }
  # assign names to matrix
  rownames(Sigma) <- fields
  colnames(Sigma) <- fields
  # return
  Sigma
}


#' Generate confidence interval for power-law fit
#'
#' @md
#' @description
#' Generate confidence intervals for power law fits, based on the
#' first delta method, and a calculate variance-covariance matrix
#'
#' @inheritParams powerlaw_predictioninterval
#' @param Sigma Variance-covariance matrix, with row and column names
#'   indicating the variable names
#' @return dataframe with fields:
#'   * `x`: x-coordinates
#'   * `y`: average (power law)
#'   * `ymin`: lower limit of confidence interval
#'   * `ymax`: upper limit of confidence interval
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
#' ftw <- powerlaw_fit(x, y, "weibull")
#' # generate covariance matrix
#' Sw <- powerlaw_covariancematrix(
#'   x, y, "weibull",
#'   multiplier = ftw$multiplier,
#'   exponent = ftw$exponent,
#'   shape = ftw$shape
#' )
#' Cw <- powerlaw_confidenceinterval(
#'   seq(min(x), max(x), l = 101),
#'   ftw$multiplier,
#'   ftw$exponent,
#'   Sw
#' )
#' # plot
#' plot(x, y)
#' lines(Cw$x, Cw$y, col = "red")
#' lines(Cw$x, Cw$ymin, col = "blue")
#' lines(Cw$x, Cw$ymax, col = "blue")
#'
#' # gamma fit
#' ftg <- powerlaw_fit(x, y, "gamma")
#' # generate covariance matrix
#' Sg <- powerlaw_covariancematrix(
#'   x, y, "gamma",
#'   multiplier = ftg$multiplier,
#'   exponent = ftg$exponent,
#'   shape = ftg$shape
#' )
#' Cg <- powerlaw_confidenceinterval(
#'   seq(min(x), max(x), l = 101),
#'   ftw$multiplier,
#'   ftw$exponent,
#'   Sg
#' )
#' # plot
#' plot(x, y)
#' lines(Cg$x, Cg$y, col = "red")
#' lines(Cg$x, Cg$ymin, col = "blue")
#' lines(Cg$x, Cg$ymax, col = "blue")
#'
#'
powerlaw_confidenceinterval <- function(
    x,
    multiplier,
    exponent,
    Sigma,
    level = 0.95
) {
  # derivatives of power-law function
  dy_dmultiplier <- x^exponent
  dy_dpower <- multiplier*log(x)*x^exponent
  # get confidence interval using delta method
  var <- (
    dy_dmultiplier^2*Sigma["multiplier", "multiplier"] +
      2*dy_dmultiplier*dy_dpower*Sigma["multiplier", "exponent"] +
      dy_dpower^2*Sigma["exponent", "exponent"]
  )
  # multiplier of standard deviation
  mult <- stats::qnorm(0.5 + 0.5*level)  # gaussian distributed
  # confidence interval
  y <- multiplier*x^exponent
  ymin <- y - mult*sqrt(var)
  ymax <- y + mult*sqrt(var)
  # return dataframe
  data.frame(x = x, y = y, ymin = ymin, ymax = ymax)
}


#' Kolmogorov-Smirnov of power-law fit
#'
#' @description
#' Calculate Kolmogorov-Smirnov (KS) distance for power-law fit.
#'
#' Also generates dataframes to plot measured and fitted cumulative density
#' results.
#'
#' @md
#' @inheritParams powerlaw_predictioninterval
#' @inheritParams powerlaw_fit
#' @param n number of points on fitted cumulative trace
#' @param gamma Euler-Mascheroni constant
#' @return list with fields:
#'   * `ks_distance`: KS-distance
#'   * `data`: dataframe (with `x` and `y` fields) for cumulative density
#'     distribution of data
#'   * `fit`: dataframe (with `x` and `y` fields) for cumulative density
#'     distribution of fit
#'   * `difference`: dataframe (with `x` and `y` fields) of locations at
#'     maximum distance (KS-distance) between the data and fitted curves#'
#' @export
#' @examples
#' # generate some data
#' y0 <- 20
#' beta <- -0.5
#' kappa <- 4
#' lambda <- 1/gamma(1 + 1/kappa)
#' x <- seq(1, 8, l = 51)
#' y <- y0*x^beta*rweibull(length(x), kappa, lambda)
#' w <- x
#' w <- rep(1, length(x))
#'
#' # gamma
#' ft <- powerlaw_fit(x, y, "gamma", weights = w)
#' ks <- powerlaw_ks(
#'   x, y, "gamma", ft$multiplier, ft$exponent, shape = ft$shape, weights = w
#' )
#' # gumbel
#' ft <- powerlaw_fit(x, y, "gumbel")
#' ks <- powerlaw_ks(
#'   x, y, "gumbel", ft$multiplier, ft$exponent, scale = ft$scale
#' )
#' # logistic
#' ft <- powerlaw_fit(x, y, "logistic")
#' ks <- powerlaw_ks(
#'   x, y, "logistic", ft$multiplier, ft$exponent, scale = ft$scale
#' )
#' # lognormal
#' ft <- powerlaw_fit(x, y, "lognormal")
#' ks <- powerlaw_ks(
#'   x, y, "lognormal", ft$multiplier, ft$exponent, sdlog = ft$sdlog
#' )
#' # normal
#' ft <- powerlaw_fit(x, y, "normal_strength")
#' ks <- powerlaw_ks(
#'   x, y, "normal_strength", ft$multiplier, ft$exponent,
#'   sd_multiplier = ft$sd_multiplier
#' )
#' # weibull
#' ft <- powerlaw_fit(x, y, "weibull")
#' ks <- powerlaw_ks(
#'   x, y, "weibull", ft$multiplier, ft$exponent, shape = ft$shape
#' )
#' # uniform
#' ft <- powerlaw_fit(x, y, "uniform")
#' ks <- powerlaw_ks(
#'   x, y, "uniform", ft$multiplier, ft$exponent, width = ft$width
#' )
#'
#' # plot
#' plot(ks$data$x, ks$data$y, "l")
#' lines(ks$fit$x, ks$fit$y, col = "blue")
#' lines(ks$difference$x, ks$difference$y, col = "red")
#'
powerlaw_ks <- function(
    x,
    y,
    model,
    multiplier,
    exponent,
    sd_multiplier = NA,
    sd_exponent = 0,
    sdlog = NA,
    shape = NA,
    scale = NA,
    width = NA,
    n = 101,
    gamma = 0.57721566490153286060651209008240243104215933593992,
    weights = rep(1, length(x))
) {
  # change model name to lowercase
  model <- tolower(model)
  # fits
  if (model == "gamma") {
    # gamma ####
    yn <- y/(multiplier*x^exponent)
    or <- order(yn)
    yn <- yn[or]
    C2 <- rep(stats::pgamma(yn, shape = shape, scale = 1/shape), each = 2)
    xp <- seq(min(yn, na.rm = TRUE), max(yn, na.rm = TRUE), l = n)
    yp <- stats::pgamma(xp, shape, scale = 1/shape)
  } else if (model == "gumbel") {
    # gumbel ####
    yn <- (y - multiplier*x^exponent)/(scale*x^exponent)
    or <- order(yn)
    yn <- yn[or]
    C2 <- rep(exp(-exp(-(yn + gamma))), each = 2)
    xp <- seq(min(yn, na.rm = TRUE), max(yn, na.rm = TRUE), l = n)
    yp <- exp(-exp(-(xp + gamma)))
  } else if (model %in% c("logistic")) {
    # logistic ####
    yn <- (y - multiplier*x^exponent)/(scale*x^exponent)
    or <- order(yn)
    yn <- yn[or]
    C2 <- rep(stats::plogis(yn, location = 0, scale = 1), each = 2)
    xp <- seq(min(yn, na.rm = TRUE), max(yn, na.rm = TRUE), l = n)
    yp <- stats::plogis(xp, location = 0, scale = 1)
  } else if (model %in% c("lognormal", "lognormal_corrected")) {
    # lognormal ####
    yn <- y/(multiplier*x^exponent)
    or <- order(yn)
    yn <- yn[or]
    C2 <- rep(stats::plnorm(yn, -0.5*sdlog^2, sdlog), each = 2)
    xp <- seq(min(yn, na.rm = TRUE), max(yn, na.rm = TRUE), l = n)
    yp <- stats::plnorm(xp, -0.5*sdlog^2, sdlog)
  } else if (model %in% c("lognormal_uncorrected")) {
    # lognormal (uncorrected) ####
    yn <- y/(multiplier*x^exponent)
    or <- order(yn)
    yn <- yn[or]
    C2 <- rep(stats::plnorm(yn, 0, sdlog), each = 2)
    xp <- seq(min(yn, na.rm = TRUE), max(yn, na.rm = TRUE), l = n)
    yp <- stats::plnorm(xp, 0, sdlog)
  } else if (model %in% c("normal", "normal_strength", "normal_force", "normal_scaled", "normal_sdpower")) {
    # normal ####
    if (model %in% c("normal", "normal_strength")) {
      sd_exponent <- 0
    } else if (model == "normal_force") {
      sd_exponent <- -2
    } else if (model == "normal_scaled") {
      sd_exponent <- exponent
    }
    yn <- (y - multiplier*x^exponent)/(sd_multiplier*x^sd_exponent)
    or <- order(yn)
    yn <- yn[or]
    C2 <- rep(stats::pnorm(yn, 0, 1), each = 2)
    xp <- seq(min(yn, na.rm = TRUE), max(yn, na.rm = TRUE), l = n)
    yp <- stats::pnorm(xp, 0, 1)
  } else if (model == "uniform") {
    # uniform ####
    yn <- y/(multiplier*x^exponent)
    or <- order(yn)
    yn <- yn[or]
    C2 <- rep(stats::punif(yn, 1 - 0.5*width/multiplier, 1 + 0.5*width/multiplier), each = 2)
    xp <- seq(1 - 0.5*width/multiplier, 1 + 0.5*width/multiplier, l = n)
    yp <- stats::punif(xp, 1 - 0.5*width/multiplier, 1 + 0.5*width/multiplier)
  } else if (model == "weibull") {
    # weibull ####
    yn <- y/(multiplier*x^exponent)
    or <- order(yn)
    yn <- yn[or]
    C2 <- rep(stats::pweibull(yn, shape, 1/gamma(1 + 1/shape)), each = 2)
    xp <- seq(min(yn, na.rm = TRUE), max(yn, na.rm = TRUE), l = n)
    yp <- stats::pweibull(xp, shape, 1/gamma(1 + 1/shape))
  }
  # real cumulative
  w <- (weights[or])/sum(weights)
  nx <- length(x)
  C1 <- rep(c(0, cumsum(w)), each = 2)[2:(2*nx + 1)]
  # distances between experimental and fitted cumulative curves
  distances <- abs(C2 - C1)
  i_max <- which.max(distances)
  # return list
  list(
    ks_distance = distances[i_max],
    data = data.frame(
      x = rep(yn, each = 2),
      y = C1
    ),
    fit = data.frame(x = xp, y = yp),
    difference = data.frame(
      x = rep(yn[floor((i_max + 1)/2)], 2),
      y = c(C1[i_max], C2[i_max])
    )
  )
}


#' Fit a power-law curve to a series of (x, y) data - iterative weighting
#'
#' @md
#' @description
#' Fit a power-law curve to a series of weighted (x, y) data. Can use
#' different assumptions for the description of intra-diameter variation.
#'
#' This function is a wrapper around the function `powerlaw_fit()`, designed
#' to deal with cases where there is very unequal weighting between various
#' observations. This may cause the initial guess, as used in a
#' gradient-descent root solving method, to be too far of the target, which
#' will inhibit convergence and may cause the algorithm to break.
#'
#' This function gradually finds the solution by scaling the weighting factors,
#' and using the results from a fit to inform the 'next' fit with altered
#' weighting factors.
#'
#' There are two loops:
#' * in loop 1, the weighting is gradually reduced:
#'
#'     w = w^(frac^n)
#'
#'   where `frac` is a user-defined scaling fraction, and 'n' a counter which
#'   increases with every iteration by a factor 'dn'
#' * loop 2 starts once a stable solution is found in loop 1. The factor 'n'
#'   is now reduced by 'dn' every step. The results from the previous analysis
#'   are used as an initial guess for the next. If no stable solution is
#'   found, dn is halved and a new guess is made. This process is continued
#'   until n has reduced to 0 again and a stable solution is obtained
#'
#' @inheritParams powerlaw_fit
#' @param n_max maximum number of iterations per loop. If exceeded, the
#'   function does not return anything
#' @param frac the weight scaling factor, see description
#' @param dn increment in n, see description
#' @return same output as `powerlaw_fit()`. Returns `NULL` if too many
#'   iterations are required
#' @export
#' @examples
#' y0 <- 20
#' beta <- -0.5
#' shape <- 4
#' x <- seq(1, 10, l = 51)
#' y <- y0*x^beta*stats::rweibull(length(x), shape, 1/gamma(1 + 1/shape))
#' y[x <= 3] <- 10*y[x <= 3]
#' weights <- x^1.1
#'
#' model <- "logistic"
#' ft1 <- powerlaw_fit_iterativeweighting(x, y, model, weights = weights)
#'
#' xp <- seq(min(x), max(x), l = 101)
#' yp <- ft1$multiplier*xp^ft1$exponent
#' plot(x, y)
#' lines(xp, yp, col = "red")
#'
powerlaw_fit_iterativeweighting <- function(
    x,
    y,
    model = "normal",
    weights = rep(1, length(x)),
    method = "bisection",
    range = c(-1, 1),
    frac = 0.5,
    n_max = 100,
    dn = 1
) {
  # initial value of n
  n <- 0
  # Loop 1: decrease weight scaling until stable solution found
  for (i in 1:n_max) {
    # current weighting
    w <- weights^(frac^n)
    # find solution - return NA with function throws an error
    sol <- tryCatch(
      powerlaw_fit(
        x, y, model = model, weights = w,
        method = method, range = range
      ),
      warning = function(w) list(loglikelihood = NA)
    )
    if (is.na(sol$loglikelihood)) {
      # decrease weighting
      n <- n + dn
    } else {
      # solution found -> break loop
      break
    }
  }
  # Loop 2: increase again, and use initial guess from last step as starting point
  if (n > 0) {
    n <- n - dn
    for (i in 1:n_max) {
      # find initial guess, based on previous solution
      if (tolower(model) %in% c("gamma", "weibull")) {
        start <- c(sol$exponent, sol$shape)
      } else if (tolower(model) == "gumbel") {
        start <- c(sol$exponent, sol$scale)
      } else if (tolower(model) == "logistic") {
        start <- c(sol$multiplier, sol$exponent, sol$scale)
      } else if (tolower(model) %in% c("normal", "normal_strength", "normal_force", "normal_scaled")) {
        start <- sol$exponent
      } else if (tolower(model) == "normal_sdpower") {
        start <- c(sol$exponent, sol$sd_exponent)
      }
      # current weighting
      w <- weights^(frac^n)
      # find solution - return NA with function throws an error
      sol_new <- tryCatch(
        powerlaw_fit(
          x, y, model = model, weights = w,
          method = method, range = range,
          start = start
        ),
        warning = function(w) list(loglikelihood = NA)
      )
      if (is.na(sol_new$loglikelihood)) {
        # no solution found - decrease weighting again
        dn <- 0.5*dn
        n <- n + dn
      } else {
        # solution found, increase weighting and update current solution
        n <- n - dn
        sol <- sol_new
      }
      if (n < 0) {
        # final solution reached -> break loop
        n <- n + dn
        break
      }
    }
  }
  # return solution
  if (n <= 0) {
    sol
  } else {
    NULL
  }
}
