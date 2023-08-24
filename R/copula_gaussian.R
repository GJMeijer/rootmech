# Gaussian copula functions
# 07/08/2023 - G. J. Meijer


#' Fit gaussian copula using loglikelihood fitting
#'
#' @description
#' Find the correlation coefficient for the best gaussian copula
#' fit through given two-dimensional cumulative density data.
#'
#' @inheritParams copula_gaussian_multiplier
#' @param weights vector with weights for each observation
#' @param eps small value, so that the correlation coefficient is set on
#'   determined on the interval -(1 - eps) <= rho <= (1 - eps). rho >= 1
#'   and rho <= -1 will lead to numerical problems.
#' @return dataframe with weighted loglikelihood of the best fit (field
#'   `loglikelihood`) and the optimal value of the correlation coefficient
#'   (field `rho`)
#' @export
#' @examples
#' # generate data
#' xy <- bivariate_normal_random(1000, rho = -0.6)
#' px <- pnorm(xy[, 1])
#' py <- pnorm(xy[, 2])
#'
#' # fit correlation coefficient
#' ft <- copula_gaussian_fit(px, py)
#'
#' # confidence interval
#' dc <- bivariate_normal_predictioninterval(rho = ft$rho)
#'
#' # plot observations + prediction interval
#' plot(xy$x, xy$y)
#' lines(dc$x, dc$y)
#'
copula_gaussian_fit <- function(
    px,
    py,
    weights = rep(1, length(px)),
    eps = 1e-12
) {
  # do not use any observations with NA values in px and/or py
  i <- !(is.na(px) | is.na(py))
  px <- px[i]
  py <- py[i]
  # initial guess
  rho0 <- copula_gaussian_initialguess(px, py)
  # optimize
  sol <- stats::optim(
    rho0,
    copula_gaussian_likelihood,
    gr = copula_gaussian_likelihood_jacobian,
    px = px,
    py = py,
    weights = weights,
    method = "L-BFGS-B",
    lower = -1 + eps,
    upper = 1 - eps,
    control = list(fnscale = -1)
  )
  # get actual loglikelihood score
  x <- stats::qnorm(px, mean = 0, sd = 1)
  y <- stats::qnorm(py, mean = 0, sd = 1)
  dx <- stats::dnorm(x, mean = 0, sd = 1)
  dy <- stats::dnorm(y, mean = 0, sd = 1)
  cm <- copula_gaussian_multiplier(px, py, rho = sol$par)
  L <- sum(weights*(dx + dy + cm))
  # return best rho + loglikelihood value
  data.frame(
    loglikelihood = L,
    rho = sol$par
  )
}


#' GGplot for copula fitting
#'
#' @description
#' Generate a plot, using ggplot, of the copula fit between root strain to
#' failure and root tensile strength, both normalised by the power-law fitted
#' predictions.
#'
#' A Gaussian copula is used, and the marginal distributions are both
#' Weibull-distributed with known shape parameter and a mean of 1.
#'
#' Prediction intervals can be generated.
#'
#' @importFrom rlang .data
#' @inheritParams power_weibull_plot
#' @param x,y arrays with x and y-observations.
#' @param shape_x,shape_y Weibull shape parameters for parameters `x` and `y`
#' @param rho correlation coefficient for the Gaussian copula
#' @param ellipsoid_axes if `TRUE`, plot the major and minor axis of the
#'   Gaussian prediction ellipsoid
#' @param ellipsoid_axes_linetype array with linetypes for major and minor
#'   axis of the prediction ellipsoid
#' @return ggplot object
#' @export
#' @examples
#' # input
#' shape_x <- 4
#' shape_y <- 6
#' rho <- -0.4
#'
#' # create random data
#' df <- biomech_random(rep(1, 100), 1, 0, shape_x, 1, 0, shape_y, rho = rho)
#'
#' # annotations
#' ann <- create_annotations(rho, prefix = "rho==", parse = TRUE)
#'
#' # plot
#' copula_gaussian_plot(df$epsru, df$tru, shape_x, shape_y, rho = rho,
#'   annotations = ann)
#'
copula_gaussian_plot <- function(
    x,
    y,
    shape_x,
    shape_y,
    rho = 0,
    confidence = c(0.5, 0.95),
    annotations = NULL,
    xlab = expression(epsilon[r*","*u]/epsilon[r*","*u*","*fit]),
    ylab = expression(t[r*","*u]/t[r*","*u*","*fit]),
    xlim = c(NA, NA),
    ylim = c(NA, NA),
    ticks = 7,
    ellipsoid_axes = TRUE,
    ellipsoid_axes_linetype = c(1, 2),
    settings = plot_settings()
) {
  # initiate plot
  plt <- ggplot2::ggplot() +
    ggplot2::theme_bw() +
    ggplot2::xlab(xlab) +
    ggplot2::ylab(ylab)
  # generate prediction intervals
  if (!is.null(rho) & length(confidence) >= 1) {
    df_ell <- bivariate_normal_predictioninterval(rho = rho, confidence = confidence)
    df_conf <- bivariate_normal_toweibull(df_ell$x, df_ell$y, shape_x, shape_y)
    df_conf$confidence <- df_ell$confidence
    # add prediction intervals to plot
    plt <- plt + ggplot2::geom_polygon(
      data = df_conf,
      ggplot2::aes(x = .data$x, y = .data$y, group = as.factor(.data$confidence)),
      alpha = settings$alpha_fit/length(confidence),
      fill = settings$fill_fit
    )
    # calculate axes
    if (ellipsoid_axes == TRUE) {
      df_ax <- bivariate_normal_axes(rho = rho, confidence = max(confidence))
      df_axes <- bivariate_normal_toweibull(df_ax$x, df_ax$y, shape_x, shape_y)
      df_axes$axis <- df_ax$axis
      # add axes to plot
      plt <- plt + ggplot2::geom_path(
        data = df_axes,
        ggplot2::aes(x = .data$x, y = .data$y, linetype = as.factor(.data$axis)),
        color = settings$color_fit,
        linewidth = settings$linewidth_fit,
        show.legend = FALSE
      ) +
        ggplot2::scale_linetype_manual(values = ellipsoid_axes_linetype)
    }
    # determine axes limits
    xlims <- round_range(c(x, df_conf$x), lim = xlim, ticks = ticks)
    ylims <- round_range(c(y, df_conf$y), lim = ylim, ticks = ticks)
  } else {
    # determine axes limits
    xlims <- round_range(x, lim = xlim, ticks = ticks)
    ylims <- round_range(y, lim = ylim, ticks = ticks)
  }
  # add points to plot
  plt <- plt + ggplot2::geom_point(
    data = data.frame(x = x, y = y),
    ggplot2::aes(x = .data$x, y = .data$y),
    color = settings$color_meas,
    shape = settings$shape_meas,
    size = settings$size_meas
  )
  # add labels to plot
  if (!is.null(annotations)) {
    annotations$xplot <- xlims$lim[1] + diff(xlims$lim)*annotations$x
    annotations$yplot <- ylims$lim[1] + diff(ylims$lim)*annotations$y
    plt <- plt + ggplot2::geom_text(
      data = annotations,
      ggplot2::aes(x = .data$xplot, y = .data$yplot, label = .data$label, hjust = .data$hjust, vjust = .data$vjust),
      parse = annotations$parse[1],
      color = settings$color_ann,
      size = settings$size_ann
    )
  }
  # change axis limits
  plt <- plt +
    ggplot2::scale_x_continuous(breaks = xlims$breaks) +
    ggplot2::scale_y_continuous(breaks = ylims$breaks) +
    ggplot2::coord_cartesian(xlim = xlims$lim, ylim = ylims$lim, expand = FALSE)
  # return plot
  plt
}


#' Objective loglikelihood function for gaussian copula fitting
#'
#' @description
#' Returns scaled version of the loglikelihood of the gaussian copula fit.
#' There is a single fitting parameter: correlation coefficient `rho`.
#'
#' The loglikelihood value is given by:
#'
#'   L = sum(weights*(dx + dy + cm))
#'
#' where dx and dy are the probability densities of the marginal distributions,
#' and cm the copula multiplier. Since only cm is dependent on the value
#' of rho, we can instead maximise instead:
#'
#'   L = sum(weights*cm)
#'
#' @inheritParams copula_gaussian_fit
#' @param rho correlation coefficient to fit
#' @return weighted loglikeihood value
#' @export
copula_gaussian_likelihood <- function(
    rho,
    px,
    py,
    weights = rep(1, length(px))
) {
  # log-transformed copula multipliers
  cm <- copula_gaussian_multiplier(px, py, rho = rho, log = TRUE)
  # loglikelihood (scaled)
  sum(weights*cm)
}


#' Jacobian of the function `copula_gaussian_likelihood()`
#'
#' @description
#' Returns the derivative of the (scalar) function
#' `copula_gaussian_likelihood()` with respect to input parameter `rho`
#'
#' @inheritParams copula_gaussian_likelihood
#' @return derivative of `copula_gaussian_likelihood()` with respect to `rho`
#' @export
#' @examples
#' # Test: compare analytical with numerical jacobian
#' rho <- 0.6
#' px <- seq(0.4, 0.5, l = 10)
#' py <- seq(0.2, 0.4, l = 10)
#' weights <- seq(0.8, 1.3, l = 10)
#'
#' eps <- 1e-6
#' L1 <- copula_gaussian_likelihood(rho, px, py)
#' L2 <- copula_gaussian_likelihood(rho + eps, px, py)
#' J2 <- (L2 - L1)/eps
#' J <- copula_gaussian_likelihood_jacobian(rho, px, py)
#'
#' c(J, J2)
#'
copula_gaussian_likelihood_jacobian <- function(
    rho,
    px,
    py,
    weights = rep(1, length(px))
) {
  # log-transformed copula multipliers - derivatives
  cm_jac <- copula_gaussian_multiplier_jacobian(px, py, rho = rho, log = TRUE)
  # sum weighted log-transformed probabilities
  sum(weights*cm_jac$rho)
}


#' Initial guess for gaussian copula fitting
#'
#' @description
#' Make an initial estimate for correlation coefficient rho, which is
#' used in the gaussian copula fitting function `copula_gaussian_fit()`
#'
#' @inheritParams copula_gaussian_fit
#' @return initial guess for correlation coefficient `rho`
copula_gaussian_initialguess <- function(px, py) {
  x <- stats::qnorm(px, mean = 0, sd = 1)
  y <- stats::qnorm(py, mean = 0, sd = 1)
  stats::cov(x, y)
}


#' Gaussian copula multiplier
#'
#' @description
#' Return gaussian copula multipliers based on known 2-d cumulative probability
#' density values and known correlation coefficient.
#'
#' The function is vectorised.
#'
#' @param px,py vectors with cumulative probability densities
#' @param rho correlation coefficient
#' @param log if `log = TRUE`, return the log-transformed copula multipliers
#'   instead
#' @return vector with gaussian copula multipliers for each pair of `px` and
#' `py`.
#' @export
#' @examples
#' # input
#' rho <- 0.8
#' mu <- c(2, 3)
#' sd <- c(1, 1)
#'
#' # generate some data
#' x <- mu[1] + 3*sd[1]*seq(-1, 1, l = 101)
#' y <- mu[2] + 3*sd[2]*seq(-1, 1, l = 101)
#'
#' # grid
#' df <- expand.grid(x = x, y = y, KEEP.OUT.ATTRS = FALSE)
#' # probability density, and cumulative probability density
#' df$dx <- dnorm(df$x, mean = mu[1], sd = sd[1])
#' df$dy <- dnorm(df$y, mean = mu[2], sd = sd[2])
#' df$px <- pnorm(df$x, mean = mu[1], sd = sd[1])
#' df$py <- pnorm(df$y, mean = mu[2], sd = sd[2])
#'
#' # copula multiplier
#' df$cm <- copula_gaussian_multiplier(df$px, df$py, rho = rho)
#'
#' # matrix with probability
#' z <- matrix(with(df, cm*dx*dy), nrow = length(x), ncol = length(y))
#'
#' # plot
#' contour(x, y, z)
#'
copula_gaussian_multiplier <- function(px, py, rho = 0, log = FALSE) {
  a <- sqrt(2)*erfi(2*px - 1)
  b <- sqrt(2)*erfi(2*py - 1)
  if (log == TRUE) {
    -0.5*log(1 - rho^2) - 0.5*(rho^2*(a^2 + b^2) - 2*rho*a*b)/(1 - rho^2)
  } else {
    exp(-(rho^2*(a^2 + b^2) - 2*rho*a*b)/(2*(1 - rho^2)))/sqrt(1 - rho^2)
  }
}


#' Derivative of `copula_gaussian_multiplier()`
#'
#' @description
#' Calculate derivative of copula multiplier calculated using the function
#' `copula_gaussian_multiplier()` with respect to its input parameters `px`,
#' `py` and `rho`.
#'
#' The function is vectorised.
#'
#' @inheritParams copula_gaussian_multiplier
#' @return list with derivatives, with fieldnames `px`, `py` and `rho`.
#' @export
#' @examples
#' px <- 0.1
#' py <- 0.2
#' rho <- -0.8
#' log <- FALSE
#'
#' cm0 <- copula_gaussian_multiplier(px, py, rho, log = log)
#' J <- copula_gaussian_multiplier_jacobian(px, py, rho, log = log)
#'
#' eps <- 1e-6
#' J2 <- (c(
#'   copula_gaussian_multiplier(px + eps, py, rho, log = log),
#'   copula_gaussian_multiplier(px, py + eps, rho, log = log),
#'   copula_gaussian_multiplier(px, py, rho + eps, log = log)
#' ) - cm0)/eps
#'
#' J
#' J2
#'
copula_gaussian_multiplier_jacobian <- function(px, py, rho = 0, log = FALSE) {
  # derivatives of a and b
  erfix <- erfi(2*px - 1)
  erfiy <- erfi(2*py - 1)
  a <- sqrt(2)*erfix
  b <- sqrt(2)*erfiy
  da_dpx <- sqrt(2*pi)*exp(erfix^2)
  db_dpy <- sqrt(2*pi)*exp(erfiy^2)
  # derivatives of multiplier (cm)
  if (log == TRUE) {
    dcm_da <- -(rho*(b - a*rho))/(rho^2 - 1)
    dcm_db <- -(rho*(a - b*rho))/(rho^2 - 1)
    dcm_drho <- (-rho^3 + a*b*rho^2 + (-a^2 - b^2 + 1)*rho + a*b)/(rho^2 - 1)^2
  } else {
    dcm_da <- -(exp((rho^2*(a^2 + b^2) - 2*a*b*rho)/(2*rho^2 - 2))*(-2*a*rho^2 + 2*b*rho))/((1 - rho^2)^(1/2)*(2*rho^2 - 2))
    dcm_db <- -(exp((rho^2*(a^2 + b^2) - 2*a*b*rho)/(2*rho^2 - 2))*(-2*b*rho^2 + 2*a*rho))/((1 - rho^2)^(1/2)*(2*rho^2 - 2))
    dcm_drho <- (exp((rho^2*(a^2 + b^2) - 2*a*b*rho)/(2*rho^2 - 2))*(-a^2*rho + a*b*rho^2 + a*b - b^2*rho - rho^3 + rho))/(1 - rho^2)^(5/2)
  }
  # return derivatives
  list(
    px = dcm_da*da_dpx,
    py = dcm_db*db_dpy,
    rho = dcm_drho
  )
}

