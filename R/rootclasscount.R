#' Multisegment power-law fitting of root class count
#'
#' @description
#' Fits a multisegment, C^0-continuous power-law probability distribution to
#' root diameter class data. Each class has an lower (`xmin`) and upper (`xmax`)
#' diameter limit, and `y` number of roots within (alternatively, `y` could
#' be root length, for example acquired with WinRhizo measurements).
#'
#' Fitting is performed by maximising the likelihood.
#'
#' @param xmin,xmax vectors with lower and upper limits for each root diameter
#'   class. `xmax > xmin >= TRUE` for each element.
#' @param y vector with number of roots (or total root length) per class
#' @param ns number of power-law segments to use in fit
#' @param n0 number of extra optional segments breakpoint positions to include
#'   in generating an initial guess
#' @param fixed vector with `2*ns + 1` elements, with values of the
#'   fitting parameter (`ns + 1` breakpoints, `ns` power-law coefficients)
#'   that should be fixed and therefore not fitted. Values that are `NA` will
#'   be fitted while any non-NA will be held constant
#' @param guess initial guess for the fitting parameters (a vector with
#'   `ns + 1` x-breakpoints, and `ns` power-law coefficients). If not
#'   provided, an initial guess is made according to the function
#'   `rootclass_initialguess()`,
#' @param weights_multiplier vector with fitting multiplication weights for
#'   each class
#' @param weights_power vector with fitting power-law coefficients for each
#'   class
#' @return a list containing the loglikelihood score of the best fit
#'   (`loglikelihood`) and a dataframe (`par`) with power-law coefficients
#'   for each segment
#' @export
#' @examples
#' no <- 9
#' xc <- seq(2, 6.5, l = no + 1)
#' xmin <- xc[1:no]
#' xmax <- xc[2:(no + 1)]
#' y <- (8 - (0.5*(xmin + xmax) - 5)^2)
#'
#' ns <- 2
#'
#' ft1 <- rootclass_fit(xmin, xmax, y, ns)
#' ft2 <- rootclass_fit(xmin, xmax, y, ns, fixed = c(1, NA, 7, NA, -3))
#'
#' rootclass_cumulative_plot(ft1$par, xmin, xmax, y)
#' rootclass_cumulative_plot(ft2$par, xmin, xmax, y)
rootclass_fit <- function(
    xmin,
    xmax,
    y,
    ns,
    n0 = 0,
    fixed = rep(NA, 2*ns + 1),
    guess = NULL,
    weights_multiplier = rep(1, length(xmin)),
    weights_power = rep(0, length(xmin))
) {
  # expand length of weights vectors to same length as limits
  weights_multiplier <- rep(weights_multiplier, length(xmin))[1:length(xmin)]
  weights_power <- rep(weights_power, length(xmin))[1:length(xmin)]
  # initial guess
  if (is.null(guess)) {
    guess <- rootclass_initialguess(xmin, xmax, y, ns, fixed = fixed, n0 = n0)
  }
  # constraints
  constr <- rootclass_constraints(ns, xmin, xmax, y, fixed = fixed)
  # fit
  sol <- stats::constrOptim(
    guess,
    rootclass_loglikelihood,
    grad = rootclass_loglikelihood_jacobian,
    ui = constr$ui,
    ci = constr$ci,
    method = "BFGS",
    control = list(fnscale = -1),
    xmin = xmin,
    xmax = xmax,
    y = y,
    fixed = fixed,
    weights_multiplier = weights_multiplier,
    weights_power = weights_power
  )
  # generate list with outputs
  parall <- fixed
  parall[is.na(fixed)] <- sol$par
  xb <- parall[1:(ns + 1)]
  b <- parall[(ns + 2):(2*ns + 1)]
  a <- rootclass_multipliers(xb, b)
  pti <- rootclass_probfit(xb, a, b)
  list(
    loglikelihood = sol$value,
    par = data.frame(
      xmin = xb[1:ns],
      xmax = xb[2:(ns + 1)],
      multiplier = a/sum(pti),
      power = b,
      total = pti/sum(pti)
    )
  )
}


#' Plot observed and fitted cumulative density of root class distributions
#'
#' @description
#' Plot the observed and multisegment power-law fit for the cumulative density
#' distribution of variable `x`, for root diameter class count data
#'
#' Plots are generated using the ggplot2 package.
#'
#' @importFrom rlang .data
#' @inheritParams rootcount_cumulative_plot
#' @param xmin,xmax vector with lower and upper limits of diameter classes
#' @param y vector with number of roots (or root length) in each class
#' @param par dataframe with fitting values per segment. For more information,
#'   see documentation in output argument in function `rootclass_fit()`
#' @param n number of points to use for fitting line in each segment
#' @return ggplot object
#' @export
#' @examples
#' no <- 9
#' xc <- seq(2, 6.8, l = no + 1)
#' xmin <- xc[1:no]
#' xmax <- xc[2:(no + 1)]
#' y <- (8 - (0.5*(xmin + xmax) - 5)^2)
#'
#' ns <- 2
#' par <- rootclass_fit(xmin, xmax, y, ns)$par
#'
#' rootclass_cumulative_plot(par, xmin, xmax, y)
#'
rootclass_cumulative_plot <- function(
    par,
    xmin,
    xmax,
    y,
    n = 101,
    xlab = expression("Root diameter"~d[r]~"[mm]"),
    ylab = "Cumulative probability density [-]",
    xlim = c(0, NA),
    ylim = c(0, 1),
    ticks = 7,
    settings = plot_settings()
) {
  # cumulative trace - observations
  df1 <- rootclass_cumulative_observations(xmin, xmax, y)
  # cumulative trace - fit
  df2 <- rootcount_cumulative_fitted(par, n = n)
  # axes limits
  xlims <- round_range(c(xmin, xmax, par$xmin), lim = xlim, ticks = ticks)
  ylims <- round_range(c(0, 1), lim = ylim, ticks = ticks)
  # plot
  ggplot2::ggplot() +
    ggplot2::theme_bw() +
    ggplot2::geom_path(
      data = df1,
      ggplot2::aes(x = .data$x, y = .data$y),
      color = settings$color_meas,
      linetype = settings$linetype_meas,
      linewidth = settings$linewidth_meas
    ) +
    ggplot2::geom_point(
      data = df1,
      ggplot2::aes(x = .data$x, y = .data$y),
      color = settings$color_meas,
      size = settings$size_meas
    ) +
    ggplot2::geom_path(
      data = df2,
      ggplot2::aes(x = .data$x, y = .data$y),
      color = settings$color_fit,
      linetype = settings$linetype_fit,
      linewidth = settings$linewidth_fit
    ) +
    ggplot2::xlab(xlab) +
    ggplot2::ylab(ylab) +
    ggplot2::scale_x_continuous(breaks = xlims$breaks) +
    ggplot2::scale_y_continuous(breaks = ylims$breaks) +
    ggplot2::coord_cartesian(
      xlim = xlims$lim,
      ylim = ylims$lim,
      expand = FALSE
    )
}


#' Plot observed and fitted root probability distributions
#'
#' @description
#' Plot the observed and multisegment power-law fit for the probability density
#' distribution of binned root count data.
#'
#' Plots are generated using the ggplot2 package.
#'
#' @importFrom rlang .data
#' @inheritParams rootclass_cumulative_plot
#' @return ggplot object
#' @export
#' @examples
#' no <- 9
#' xc <- seq(2, 6.8, l = no + 1)
#' xmin <- xc[1:no]
#' xmax <- xc[2:(no + 1)]
#' y <- (8 - (0.5*(xmin + xmax) - 5)^2)
#'
#' ns <- 2
#' par <- rootclass_fit(xmin, xmax, y, ns)$par
#'
#' rootclass_density_plot(par, xmin, xmax, y)
rootclass_density_plot <- function(
    par,
    xmin,
    xmax,
    y,
    n = 101,
    xlab = expression("Root diameter"~d[r]~"[mm]"),
    ylab = "Probability density [-]",
    xlim = c(0, NA),
    ylim = c(0, NA),
    ticks = 7,
    settings = plot_settings()
) {
  # observed data
  df1 <- data.frame(xmin = xmin, xmax = xmax, y = y)
  df1$p <- with(df1, y/sum(y)/(xmax - xmin))
  # fitted density curve
  df2 <- rootcount_density_fitted(par, n = n)
  df2 <- data.frame(
    bundle_id = c(1, df2$bundle_id, max(df2$bundle_id)),
    x = c(df2$x[1], df2$x, utils::tail(df2$x, 1)),
    y = c(0, df2$y, 0)
  )
  # axes limits
  xlims <- round_range(c(df1$xmin, df1$xmax, df2$x), lim = xlim, ticks = ticks)
  ylims <- round_range(c(df1$p, df2$y), lim = ylim, ticks = ticks)
  # plot
  ggplot2::ggplot() +
    ggplot2::theme_bw() +
    ggplot2::geom_rect(
      data = df1,
      ggplot2::aes(xmin = .data$xmin, xmax = .data$xmax, ymin = 0, ymax = .data$p),
      color = settings$color_meas,
      fill = settings$fill_meas,
      alpha = settings$alpha_meas
    ) +
    ggplot2::geom_path(
      data = df2,
      ggplot2::aes(x = .data$x, y = .data$y),
      color = settings$color_fit,
      linetype = settings$linetype_fit,
      linewidth = settings$linewidth_fit
    ) +
    ggplot2::xlab(xlab) +
    ggplot2::ylab(ylab) +
    ggplot2::scale_x_continuous(breaks = xlims$breaks) +
    ggplot2::scale_y_continuous(breaks = ylims$breaks) +
    ggplot2::coord_cartesian(
      xlim = xlims$lim,
      ylim = ylims$lim,
      expand = FALSE
    )
}


#' Initial guess for `rootclass_fit()` function
#'
#' @description
#' Generate a simple initial guess for the `rootclass_fit()` function.
#'
#' Breakpoints are distributed evenly across the domain determined by the
#' midpoints of the smallest and largest root class. Power-law power
#' coefficients are simply taken as 0
#'
#' @inheritParams rootclass_fit
#' @param frac parameter for estimating outer breakpoints, if not fixed. The
#'   breakpoint is estimated as the smallest value of x-values that are at
#'   `frac` normalised position (0 to 1) within their class. For the last
#'   breakpoint position, this is the position `1 - frac`.
#' @return vector with initial guess
#' @export
#' @examples
#' no <- 9
#' xc <- seq(2, 7, l = no + 1)
#' xmin <- xc[1:no]
#' xmax <- xc[2:(no + 1)]
#' y <- (8 - (0.5*(xmin + xmax) - 5)^2)
#' ns <- 2
#' rootclass_initialguess(xmin, xmax, y, ns)
#'
rootclass_initialguess <- function(
    xmin,
    xmax,
    y,
    ns,
    n0 = 0,
    fixed = rep(NA, 2*ns + 1),
    frac = 0.25,
    weights_multiplier = rep(1, length(xmin)),
    weights_power = rep(0, length(xmin))
) {
  # estimate fitting domain
  fixed_xb <- fixed[1:(ns + 1)]
  if (is.na(fixed_xb[1])) {
    x1 <- min(xmin + frac*(xmax - xmin))
  } else {
    x1 <- fixed_xb[1]
  }
  if (is.na(fixed_xb[ns + 1])) {
    x2 <- max(xmax - frac*(xmax - xmin))
  } else {
    x2 <- fixed_xb[ns + 1]
  }
  # breakpoint selection
  xb <- cbind(
    x1,
    breakpoint_options(ns - 1, x1, x2, fixed = fixed_xb[2:ns], n0 = n0),
    x2,
    deparse.level = 0
  )
  # number of unfixed power-coefficients
  nb_unfixed <- sum(is.na(fixed[(ns+2):(2*ns + 1)]))
  # fit for each of the breakpoint options
  Lopts <- apply(
    xb,
    1,
    function(xbi) {
      rootclass_fit(
        xmin, xmax, y, ns,
        weights_multiplier = weights_multiplier,
        weights_power = weights_power,
        guess = rep(0, nb_unfixed),
        fixed = c(xbi, fixed[(ns + 2):(2*ns + 1)])
      )
    }
  )
  # find option with minimum loglikelihood score
  i <- which.max(sapply(Lopts, function(x) x$loglikelihood))
  # vector with all initial parameters
  parall <- c(
    Lopts[[i]]$par$xmin,
    Lopts[[i]]$par$xmax[ns],
    Lopts[[i]]$par$power
  )
  # return non-fixed parameters
  parall[is.na(fixed)]
}


#' Generate contraints for `rootclass_fit()` optimalisation
#'
#' @description
#' Generate constraint matrix `ui` and constraint vector `ci` to use in
#' the fitting function `rootclass_fit()`, which uses the `stats::constrOptim()`
#' optimalisation method (see `stats::constrOptim()` for more detail).
#'
#' The following constraints are used:
#' - Smallest breakpoint larger than 0
#' - Smallest breakpoint smaller than the smallest upper class limit in `xmax`
#'   that contains roots (`y > 0`)
#' - Breakpoints vector in increasing order
#' - Largest breakpoint larger than the largest lower class limit in `xmin`
#'   that contains roots (`y > 0`)
#'
#' @inheritParams rootclass_fit
#' @return a list with the constraint matrix (`ui`) and vector (`ci`)
#' @export
#' @examples
#' ns <- 3
#' xmin <- seq(2, 5)
#' xmax <- seq(3, 6)
#' y <- c(1, 4)
#' rootclass_constraints(ns, xmin, xmax, y)
#'
rootclass_constraints <- function(ns, xmin, xmax, y, fixed = rep(NA, 2*ns + 1)) {
  # generate matrix
  ui <- matrix(
    c(1, -1, rep(c(-1, rep(0, ns + 2), 1), ns), 1, rep(0, ns*(ns + 3))),
    nrow = ns + 3,
    ncol = 2*ns + 1,
    byrow = FALSE
  )
  # generate vector
  ci <- c(
    0,
    -min(xmax[y > 0], na.rm = TRUE),
    rep(0, ns),
    max(xmin[y > 0], na.rm = TRUE)
  )
  # account for fixed vectors
  # add known values to coefficients in ci - and remove from matrix
  free <- is.na(fixed)
  ci <- ci - (as.matrix(ui[, !free], ncol = sum(!free)) %*% fixed[!free])
  ui <- as.matrix(ui[, free], ncol = sum(free))
  # return list
  list(ui = ui, ci = ci)
}


#' Loglikelihood score of multi-segment power-law fit for root count data
#'
#' @inheritParams rootclass_fit
#' @param par vector with fitting parameters. For `ns` segments, the vector will
#'   consist of `ns + 1` consecutive x-breakpoints, followed by `ns` power-law
#'   cower coefficients
#' @return loglikelihood value (scalar)
#' @export
#' @examples
#' xb <- c(1.5, 3.5, 6.5, 7.5, 8.5)
#' b <- seq(-2, 2, l = length(xb) - 1)
#' par <- c(xb, b)
#' xmin <- c(2,3,4)
#' xmax <- c(2.5,4,6)
#' y <- c(3, 4, 5)
#'
#' rootclass_loglikelihood(par, xmin, xmax, y)
#'
rootclass_loglikelihood <- function(
    par,
    xmin, xmax, y,
    fixed = rep(NA, length(par)),
    weights_multiplier = rep(1, length(xmin)),
    weights_power = rep(0, length(xmin))
) {
  # generate vector with all parameters, including fixed
  parall <- fixed
  parall[is.na(fixed)] <- par
  # split parameters
  ns <- (length(parall) - 1)/2
  no <- length(xmin)
  xb <- parall[1:(ns + 1)]
  b <- parall[(ns + 2):(2*ns + 1)]
  # multipliers
  a <- rootclass_multipliers(xb, b)
  # unscaled probability in observations
  pi <- rootclass_probobs(xb, a, b, xmin, xmax)
  # unscaled total probability
  pt <- sum(rootclass_probfit(xb, a, b))
  # generate weights
  pw <- rootclass_probobs(
    xb, a, b, xmin, xmax,
    weights_multiplier = weights_multiplier, weights_power = weights_power
  )
  weights <- y*pw/pi
  # loglikelihood
  sum(weights*(log(pi) - log(pt)))
}


#' Derivative of function `rootclass_loglikelihood()`
#'
#' @description
#' Generates the derivative of the results of the function
#' `rootclass_loglikelihood()` with respect to fitting value input
#' `par`.
#'
#' Function is vectorised.
#'
#' @inheritParams rootclass_loglikelihood
#' @return a scalar containing the derivatives for each value in fitting vector
#'   `par`
#' @export
#' @examples
#' xb <- c(1.5, 3.5, 5.5)
#' b <- seq(-0.5, 0.5, l = length(xb) - 1)
#' par <- c(xb, b)
#' xmin <- c(2, 3, 4)
#' xmax <- c(2.5, 4, 6)
#' y <- c(3, 4, 5)
#' weights_multiplier <- runif(length(xmin), 2, 5)
#' weights_power <- runif(length(xmin), -2, 2)
#'
#' L <- rootclass_loglikelihood(
#'   par, xmin, xmax, y, weights_multiplier = weights_multiplier,
#'   weights_power = weights_power
#' )
#' J <- rootclass_loglikelihood_jacobian(
#'   par, xmin, xmax, y, weights_multiplier = weights_multiplier,
#'   weights_power = weights_power
#' )
#'
#' eps <- 1e-6
#' J2 <- rep(0, length(par))
#' for (i in 1:length(par)) {
#'   dx <- rep(0, length(par))
#'   dx[i] <- eps
#'   J2[i] <- (rootclass_loglikelihood(par + dx, xmin, xmax, y,
#'     weights_multiplier = weights_multiplier,
#'     weights_power = weights_power
#'    ) - L)/eps
#' }
#'
#' J
#' J2
#'
rootclass_loglikelihood_jacobian <- function(
    par,
    xmin,
    xmax,
    y,
    fixed = rep(NA, length(par)),
    weights_multiplier = rep(1, length(xmin)),
    weights_power = rep(0, length(xmin))
) {
  # generate vector with all parameters, including fixed
  parall <- fixed
  parall[is.na(fixed)] <- par
  # split parameters
  ns <- (length(parall) - 1)/2
  no <- length(xmin)
  xb <- parall[1:(ns + 1)]
  b <- parall[(ns + 2):(2*ns + 1)]
  # multipliers
  a <- rootclass_multipliers(xb, b)
  Ja <- rootclass_multipliers_jacobian(xb, b)
  # unscaled probability in observations
  pi <- rootclass_probobs(xb, a, b, xmin, xmax)
  Jpi <- rootclass_probobs_jacobian(xb, a, b, xmin, xmax)
  pi_par <- cbind(
    Jpi$xb + Jpi$a %*% Ja$xb,
    Jpi$b + Jpi$a %*% Ja$b
  )
  # unscaled total probability
  pt <- sum(rootclass_probfit(xb, a, b))
  Jpt <- rootclass_probfit_jacobian(xb, a, b)
  pt_par <- c(
    colSums(Jpt$xb + Jpt$a %*% Ja$xb),
    colSums(Jpt$b + Jpt$a %*% Ja$b)
  )
  # weights
  pw <- rootclass_probobs(xb, a, b, xmin, xmax, weights_multiplier = weights_multiplier, weights_power = weights_power)
  Jpw <- rootclass_probobs_jacobian(xb, a, b, xmin, xmax, weights_multiplier = weights_multiplier, weights_power = weights_power)
  pw_par <- cbind(
    Jpw$xb + Jpw$a %*% Ja$xb,
    Jpw$b + Jpw$a %*% Ja$b
  )
  weights <- y*pw/pi
  weights_par <- y*(pw_par/pi - pw*pi_par/pi^2)
  # loglikelihood derivative
  J <- colSums(weights_par*(log(pi) - log(pt)) + weights*(pi_par/pi)) - sum(weights)*pt_par/pt
  # only return derivative with respect to non-fixed values
  J[is.na(fixed)]
}


#' Scaling multipliers for power-laws for each segment
#'
#' @description
#' Calculates the power-law multiplication factors for each segment. These
#' are yet unscaled (so the total probability != 1), and instead assume
#' that the multiplier for the first segment is equal to 1.
#'
#' Correct probability scaling is done within the function
#' `rootclass_loglikelihood()` instead.
#'
#' @param xb vector with x-breakpoints (length = number of segments + 1)
#' @param b vector with power law coefficients for each segment
#' @return a vector with power-law multipliers for each fitted segment
#' @export
#' @examples
#' # check if correct, by checking if probability curve is C^0 continuous
#' xb <- c(1.5, 3.5, 6.5, 7.5, 8.5)
#' b <- seq(-2, 2, l = length(xb) - 1)
#' a <- rootclass_multipliers(xb, b)
#' df <- data.frame(a = a, b = b, x1 = xb[1:(length(xb) - 1)],
#'   x2 = xb[2:length(xb)])
#' df2 <- expand_grid_df(data.frame(s = seq(0, 1, l = 25)), df)
#' df2$x <- with(df2, x1 + s*(x2 - x1))
#' df2$y <- with(df2, a*x^b)
#' plot(df2$x, df2$y, "l")
#'
rootclass_multipliers <- function(xb, b) {
  c(1, cumprod(xb[2:(length(xb) - 1)]^(-diff(b))))
}


#' Derivative of function `rootclass_multipliers()`
#'
#' @description
#' Generates the derivative of the results of the function
#' `rootclass_multipliers()` with respect to its input parameters `xb` and `b`.
#'
#' Function is vectorised.
#'
#' @inheritParams rootclass_multipliers
#' @return a list containing the derivatives for each input parameter. Has
#'   fields `xb` and `b`.
#' @export
#' @examples
#' xb <- seq(1.5, 7.5, l = 4)
#' b <- seq(-2.5, 2, l = length(xb) - 1)
#' a <- rootclass_multipliers(xb, b)
#' J <- rootclass_multipliers_jacobian(xb, b)
#'
#' Jxb <- matrix(0, length(a), length(xb))
#' eps <- 1e-4
#' for (i in 1:length(xb)) {
#'   dx <- rep(0, length(xb))
#'   dx[i] <- eps
#'   a2 <- rootclass_multipliers(xb + dx, b)
#'   Jxb[, i] <- (a2 - a)/eps
#' }
#'
#' Jb <- matrix(0, length(a), length(b))
#' eps <- 1e-4
#' for (i in 1:length(b)) {
#'   dx <- rep(0, length(b))
#'   dx[i] <- eps
#'   a2 <- rootclass_multipliers(xb, b + dx)
#'   Jb[, i] <- (a2 - a)/eps
#' }
#'
#' J
#' Jxb
#' Jb
#'
rootclass_multipliers_jacobian <- function(xb, b) {
  # basic
  ns <- length(b)
  a <- rootclass_multipliers(xb, b)
  # derivative with respect to xb
  if (ns == 1) {
    list(
      xb = matrix(0, nrow = 1, ncol = 2),
      b = matrix(0, nrow = 1, ncol = 1)
    )
  } else {
    # derivatives with respect to xb
    Jxb <- outer(c(0, a[2:ns]), c(0, -diff(b)/xb[2:(length(xb) - 1)], 0), FUN = "*")
    Jxb[upper.tri(Jxb, diag = FALSE)] <- 0
    # derivative with respect to b
    M1 <- outer(c(0, a[2:ns]), log(xb[2:ns]))
    M1[upper.tri(M1, diag = TRUE)] <- 0
    # return
    list(
      xb = Jxb,
      b = cbind(M1, 0) - cbind(0, M1)
    )
  }
}


#' Total probability in each measured root class
#'
#' @description
#' Calculate the total probability within each measured root class, based on
#' the fitted power-law segment. These values are not yet scaled for the total
#' probability.
#'
#' @inheritParams rootclass_fit
#' @inheritParams rootclass_multipliers
#' @param a power-law multipliers for each segment, as calculated with the
#'   function `rootclass_multipliers()`
#' @return vector with (unscaled) total probabilities
#' @export
#' @examples
#' xb <- c(1.5, 3.5, 6.5, 7.5, 8.5)
#' a <- seq(2, 4, l = length(xb) - 1)
#' b <- seq(-2, 2, l = length(xb) - 1)
#' xmin <- c(2, 3, 4)
#' xmax <- c(2.5, 4, 20)
#' rootclass_probobs(xb, a, b, xmin, xmax)
#'
rootclass_probobs <- function(
    xb,
    a,
    b,
    xmin,
    xmax,
    weights_multiplier = rep(1, length(xmin)),
    weights_power = rep(0, length(xmin))
) {
  # number of segments + observed intervals
  ns <- length(b)
  no <- length(xmin)
  # integration limits for each segment/observation combination - matrix
  L <- pmin(outer(xmin, xb, FUN = "pmax"), xmax)
  L1 <- L[, 1:ns]
  L2 <- L[, 2:(ns + 1)]
  # matrix forms for a and b
  Ma <- outer(weights_multiplier, a)
  Mb <- outer(weights_power, b, FUN = "+")
  # integrate each combination
  tc <- power_integrate(Mb, L1, L2, multiplier = Ma)
  # return unweighted probabilty per observation
  rowSums(tc)
}


#' Derivative of function `rootclass_probobs()`
#'
#' @description
#' Generates the derivative of the results of the function
#' `rootclass_probobs()` with respect to input parameters `xb`, `a` and `b`.
#'
#' Function is vectorised.
#'
#' @inheritParams rootclass_probobs
#' @return a list containing the derivatives for each input parameter. Has
#'   fields `xb`, `a` and `b`.
#' @export
#' @examples
#' xb <- c(2.2, 3.5, 6.5, 7.5, 8.5)
#' a <- seq(2, 4, l = length(xb) - 1)
#' b <- seq(-2, 2, l = length(xb) - 1)
#' xmin <- c(1, 3, 4)
#' xmax <- c(2.5, 7, 20)
#' weights_multiplier <- runif(length(xmin), 3, 8)
#' weights_power <- runif(length(xmin), -2, 2)
#' p <- rootclass_probobs(xb, a, b, xmin, xmax,
#'   weights_multiplier = weights_multiplier, weights_power = weights_power)
#' J <- rootclass_probobs_jacobian(xb, a, b, xmin, xmax,
#'   weights_multiplier = weights_multiplier, weights_power = weights_power)
#'
#' eps <- 1e-6
#' Jxb <- matrix(0, nrow = length(p), ncol = length(xb))
#' for (i in 1:length(xb)) {
#'   dx <- rep(0, length(xb))
#'   dx[i] <- eps
#'   Jxb[, i] <- (rootclass_probobs(xb + dx, a, b, xmin, xmax,
#'     weights_multiplier = weights_multiplier,
#'     weights_power = weights_power) - p)/eps
#' }
#' Ja <- matrix(0, nrow = length(p), ncol = length(a))
#' for (i in 1:length(a)) {
#'   dx <- rep(0, length(a))
#'   dx[i] <- eps
#'   Ja[, i] <- (rootclass_probobs(
#'     xb, a + dx, b, xmin, xmax, weights_multiplier = weights_multiplier,
#'     weights_power = weights_power
#'   ) - p)/eps
#' }
#' Jb <- matrix(0, nrow = length(p), ncol = length(b))
#' for (i in 1:length(b)) {
#'   dx <- rep(0, length(b))
#'   dx[i] <- eps
#'   Jb[, i] <- (rootclass_probobs(
#'     xb, a, b + dx, xmin, xmax, weights_multiplier = weights_multiplier,
#'     weights_power = weights_power
#'   ) - p)/eps
#' }
#'
#' J
#' Jxb
#' Ja
#' Jb
#'
rootclass_probobs_jacobian <- function(
    xb,
    a,
    b,
    xmin,
    xmax,
    weights_multiplier = rep(1, length(xmin)),
    weights_power = rep(0, length(xmin))
) {
  # number of segments + observed intervals
  ns <- length(b)
  no <- length(xmin)
  # integration limits for each segment/observation combination - matrix
  L <- pmin(outer(xmin, xb, FUN = "pmax"), xmax)
  L1 <- L[, 1:ns]
  L2 <- L[, 2:(ns + 1)]
  # matrix forms for a and b
  Ma <- outer(weights_multiplier, a)
  Mb <- outer(weights_power, b, FUN = "+")
  # integrate each combination
  tc <- power_integrate(Mb, L1, L2, multiplier = Ma)
  Jtc <- power_integrate_jacobian(Mb, L1, L2, multiplier = Ma)
  # mask - breakpoint inside segment?
  Mi <- outer(xmin, xb, FUN = "<=") & outer(xmax, xb, FUN = ">=")
  tmp1 <- Jtc$lower
  tmp1[!Mi[, 1:ns]] <- 0
  tmp2 <- Jtc$upper
  tmp2[!Mi[, 2:(ns + 1)]] <- 0
  dtc_dxb <- cbind(tmp1, 0) + cbind(0, tmp2)
  # return
  list(
    xb = dtc_dxb,
    a = weights_multiplier*Jtc$multiplier,
    b = Jtc$power
  )
}


#' Total probability in fitting segment
#'
#' @description
#' Calculate the total probability within each fitted power-law segment.
#' These values are not yet scaled for the total probability.
#'
#' @inheritParams rootclass_probobs
#' @return vector with (unscaled) total probabilities
#' @export
#' @examples
#' xb <- c(1.5, 3.5, 6.5, 7.5, 8.5)
#' a <- seq(2, 4, l = length(xb) - 1)
#' b <- seq(-2, 2, l = length(xb) - 1)
#' rootclass_probfit(xb, a, b)
#'
rootclass_probfit <- function(xb, a, b) {
  ns <- length(b)
  power_integrate(b, xb[1:ns], xb[2:(ns + 1)], multiplier = a)
}


#' Derivative of function `rootclass_probfit()`
#'
#' @description
#' Generates the derivative of the results of the function
#' `rootclass_probfit()` with respect to input parameters `xb`, `a` and `b`.
#'
#' Function is vectorised.
#'
#' @inheritParams rootclass_probfit
#' @return a list containing the derivatives for each input parameter. Has
#'   fields `xb`, `a` and `b`.
#' @export
#' @examples
#' xb <- c(1.5, 3.5, 6.5, 7.5, 8.5)
#' a <- seq(2, 4, l = length(xb) - 1)
#' b <- seq(-2, 2, l = length(xb) - 1)
#' p <- rootclass_probfit(xb, a, b)
#' J <- rootclass_probfit_jacobian(xb, a, b)
#'
#' eps <- 1e-6
#' Jxb <- matrix(0, length(p), length(xb))
#' for (i in 1:length(xb)) {
#'   dx <- rep(0, length(xb))
#'   dx[i] <- eps
#'   Jxb[, i] <- (rootclass_probfit(xb + dx, a, b) - p)/eps
#' }
#' Ja <- matrix(0, length(p), length(a))
#' for (i in 1:length(a)) {
#'   dx <- rep(0, length(a))
#'   dx[i] <- eps
#'   Ja[, i] <- (rootclass_probfit(xb, a + dx, b) - p)/eps
#' }
#' Jb <- matrix(0, length(p), length(b))
#' for (i in 1:length(b)) {
#'   dx <- rep(0, length(b))
#'   dx[i] <- eps
#'   Jb[, i] <- (rootclass_probfit(xb, a, b + dx) - p)/eps
#' }
#'
#' J
#' Jxb
#' Ja
#' Jb
#'
rootclass_probfit_jacobian <- function(xb, a, b) {
  ns <- length(b)
  Jtf <- power_integrate_jacobian(b, xb[1:ns], xb[2:(ns + 1)], multiplier = a)
  if (ns == 1) {
    list(
      xb = c(Jtf$lower, Jtf$upper),
      a = Jtf$multiplier,
      b = Jtf$power
    )
  } else {
    list(
      xb = cbind(diag(Jtf$lower), 0) + cbind(0, diag(Jtf$upper)),
      a = diag(Jtf$multiplier),
      b = diag(Jtf$power)
    )
  }
}


#' Cumulative probability density curve for root class count observations
#'
#' @description
#' Generate a observation-value versus cumulative probability density
#' curve for a series of root counts and root diameter classes
#'
#' @inheritParams rootclass_fit
#' @return a dataframe with x-values (`x`) and the cumulative probability
#'   density (`y`)
#' @export
#' @examples
#' no <- 9
#' xc <- seq(2, 7, l = no + 1)
#' xmin <- xc[1:no]
#' xmax <- xc[2:(no + 1)]
#' y <- (8 - (0.5*(xmin + xmax) - 5)^2)
#'
#' rootclass_cumulative_observations(xmin, xmax, y)
#'
rootclass_cumulative_observations <- function(xmin, xmax, y) {
  # dataframe with all classes
  tmp <- data.frame(xmin = xmin, xmax = xmax, total = y, gradient = y/(xmax - xmin))
  # vector with all unique breakpoints
  xb <- sort(unique(c(xmin, xmax)))
  nb <- length(xb)
  # for each domain in `xb` (columns), determine if observation class is present in range
  t <- outer(tmp$xmin, xb[1:(nb - 1)], ">=") & outer(tmp$xmax, xb[2:nb], "<=")
  # gradient within each domain
  gr <- colSums(tmp$gradient*t)
  # return dataframe
  data.frame(
    x = xb,
    y = cumsum(c(0, gr*diff(xb)))/sum(y)
  )
}

