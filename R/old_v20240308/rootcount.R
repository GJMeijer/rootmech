#' Fit a multisegment power-law probability distribution
#'
#' @description
#' Fit a C^0 continuous, multisegment power-law probability distribution to a
#' set of observations `x`. Observations can be weighted if required.
#'
#' The fit is sensitive to the initial guess. An initial guess is generated
#' by the function `rootcount_initialguess()`, see documentation for more
#' details. The function generates `ns - 1 + n0` optinal breakpoints and fits
#' with fixed breakpoint positions. The combination of breakpoints that yields
#' the best likelihood is then used as an initial guess for the 'main' fitting
#' process.
#'
#' @param x vector with observations
#' @param ns number of segments to use
#' @param xmin,xmax minimum and maximum values of `x` to use in the fit.
#'   By default, `xmin = min(x)` and `xmax = max(x)` since this will maximise
#'   the probability. These values may however be overwritten if required.
#' @param n0 number of extra optional segments breakpoint positions to include
#'   in generating an initial guess
#' @param fixed vector with `2*ns - 1` elements, with values of the
#'   fitting parameter (`ns - 1` breakpoints, `ns` power-law coefficients)
#'   that should be fixed and therefore not fitted. Values that are `NA` will
#'   be fitted while any non-NA will be held constant
#' @param weights vector with weights for each observation. The individual
#'   probabilities are raised to the power `weights`. By default, all
#'   observations are weighted equally (`weights = 1`)
#' @param guess vector with user-defined initial guess for th vector of fitting
#'   parameters (`ns - 1` breakpoints, followed by `ns` power coefficients).
#'   If not defined, a simple guess is made
#' @return a list with the loglikelihood (field `loglikelihood`) and a
#'   dataframe (field `par`) with the fitting results for each segment.
#'   Each segment is defined using the lower and upper x-limit (`xmin`, `xmax`),
#'   power-law multiplier (`multiplier`) and power-law coefficient (`power`)
#'   and the total probability in the segment (`total`).
#' @export
#' @examples
#' # test fit
#' x <- rweibull(100, 4, 6)
#' rootcount_fit(x, 2)
#'
#' # multiple segments
#' x <- c(
#'   rweibull(100, shape = 4, scale = 1),
#'   rweibull(50, shape = 6, scale = 3),
#'   rweibull(10, shape = 12, scale = 6)
#' )
#' ft1 <- rootcount_fit(x, 3)
#' rootcount_cumulative_plot(x, ft1$par)
#' # multiple segments, with multiple initial guesses
#' ft2 <- rootcount_fit(x, 3, n0 = 5)
#' rootcount_cumulative_plot(x, ft2$par)
#'
#' # test fixed argument
#' ft3 <- rootcount_fit(x, 3, fixed = c(3, 4, 0.1, NA, 2.5))
#' rootcount_cumulative_plot(x, ft3$par)
#'
rootcount_fit <- function(
    x,
    ns,
    xmin = min(x, na.rm = TRUE),
    xmax = max(x, na.rm = TRUE),
    n0 = 0,
    fixed = rep(NA, 2*ns - 1),
    weights = rep(1, length(x)),
    guess = NULL
) {
  # remove NA values from array
  x <- x[!is.na(x)]
  # return errors
  if (any(x > xmax)) stop("Observations beyond domain, some x>xmax")
  if (any(x < xmin)) stop("Observations beyond domain, some x<xmin")
  # generate
  if (any(is.na(fixed))) {
    # fit - only if not all parameters already fixed
    # initial guess
    if (is.null(guess)) {
      guess <- rootcount_initialguess(
        x,
        ns,
        xmin = xmin,
        xmax = xmax,
        n0 = n0,
        fixed = fixed,
        weights = weights
      )
    } else {
      if (length(guess) != sum(is.na(fixed))) {
        stop("Length of initial guess vector `guess` not compatible with number of
             `NA` values in `fixed`argument")
      }
    }
    # cases
    if (ns == 1) {
      # single segment
      sol <- stats::optim(
        guess,
        rootcount_loglikelihood_single,
        gr = rootcount_loglikelihood_single_jacobian,
        method = "BFGS",
        control = list(fnscale = -1),
        x = x,
        weights = weights,
        xmin = xmin,
        xmax = xmax
      )
      # parameter vectors
      L <- sol$value
      xb <- NULL
      b <- sol$par
      pt <- 1
    } else {
      # get constraints
      constr <- rootcount_constraints(ns, xmin, xmax, fixed = fixed)
      # optimize
      sol <- stats::constrOptim(
        guess,
        rootcount_loglikelihood,
        grad = rootcount_loglikelihood_jacobian,
        constr$ui,
        constr$ci,
        method = "BFGS",
        control = list(fnscale = -1),
        x = x,
        fixed = fixed,
        weights = weights,
        xmin = xmin,
        xmax = xmax
      )
      # parameter vectors
      L <- sol$value
      parall <- rootcount_constructpar(sol$par, fixed = fixed)
      xb <- parall[1:(ns - 1)]
      b <- parall[ns:(2*ns - 1)]
      pt <- rootcount_multiplier(xb, b, xmin, xmax)
    }
  } else {
    # all parameters known/fixed
    if (ns == 1) {
      L <- rootcount_loglikelihood_single(
        fixed, x, xmin = xmin, xmax = xmax, weights = weights)
      xb <- NULL
      b <- fixed
      pt <- 1
    } else {
      L <- rootcount_loglikelihood(
        fixed, x, xmin = xmin, xmax = xmax, weights = weights,
        fixed = rep(NA, 2*ns - 1))
      xb <- fixed[1:(ns - 1)]
      b <- fixed[ns:(2*ns - 1)]
      pt <- rootcount_multiplier(xb, b, xmin, xmax)
    }
  }
  # generate output
  out <- list(
    loglikelihood = L,
    par = data.frame(
      xmin = c(xmin, xb),
      xmax = c(xb, xmax),
      power = b,
      total = pt/sum(pt)
    )
  )
  # calculate multipliers
  out$par$multiplier <- with(out$par, ifelse(
    is_near(power, -1),
    total/log(xmax/xmin),
    total*(1 + power)/(xmax^(power + 1) - xmin^(power + 1))
  ))
  # return
  out
}


#' Plot observed and fitted cumulative density of root distributions
#'
#' @description
#' Plot the observed and multisegment power-law fit for the cumulative density
#' distribution of variable `x`.
#'
#' Plots are generated using the ggplot2 package.
#'
#' @importFrom rlang .data
#' @inheritParams power_weibull_plot
#' @param x vector with observations
#' @param par dataframe with fitting values per segment. For more information,
#'   see documentation in output argument in function `rootcount_fit()`
#' @param n number of points to use for fitting line in each segment
#' @return ggplot object
#' @export
#' @examples
#' x <- c(
#'   rweibull(25, shape = 4, scale = 1),
#'   rweibull(25, shape = 6, scale = 3),
#'   rweibull(25, shape = 12, scale = 6)
#' )
#'
#' ft <- rootcount_fit(x, 6)
#'
#' rootcount_cumulative_plot(x, ft$par)
rootcount_cumulative_plot <- function(
    x,
    par,
    n = 101,
    xlab = expression("Root diameter"~d[r]~"[mm]"),
    ylab = "Cumulative probability density [-]",
    xlim = c(0, NA),
    ylim = c(0, 1),
    ticks = 7,
    settings = plot_settings()
) {
  # cumulative trace - observations
  df1 <- rootcount_cumulative_observations(x)
  # cumulative trace - fit
  df2 <- rootcount_cumulative_fitted(par, n = n)
  # axes limits
  xlims <- round_range(c(x, par$xmin), lim = xlim, ticks = ticks)
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
#' distribution of variable `x`. The observed data is binned using `bins`
#' number of equal-width bins, and per bin the average probabilty density is
#' calculated.
#'
#' Plots are generated using the ggplot2 package.
#'
#' @importFrom rlang .data
#' @inheritParams rootcount_cumulative_plot
#' @param bins number of equally spaced bins to use to plot observed data
#' @return ggplot object
#' @export
#' @examples
#' x <- c(
#'   rweibull(50, shape = 4, scale = 1),
#'   rweibull(50, shape = 6, scale = 3),
#'   rweibull(50, shape = 12, scale = 6)
#' )
#'
#' ft <- rootcount_fit(x, 6)
#'
#' rootcount_density_plot(x, ft$par, bins = 30)
#' rootcount_cumulative_plot(x, ft$par)
rootcount_density_plot <- function(
    x,
    par,
    bins = 20,
    n = 101,
    xlab = expression("Root diameter"~d[r]~"[mm]"),
    ylab = "Probability density [-]",
    xlim = c(0, NA),
    ylim = c(0, NA),
    ticks = 7,
    settings = plot_settings()
) {
  # cut - to get estimation of max density in bar
  ct <- floor((x - min(x))/(max(x) - min(x))*bins) + 1
  ct[ct > bins] <- bins
  tb <- table(ct)
  df1 <- data.frame(x = min(x) + (max(x) - min(x))*seq(0.5/bins, 1 - 0.5/bins, l = bins), n = 0)
  df1$n[as.integer(row.names(tb))] <- as.vector(tb)
  df1$y <- df1$n/length(x)/(max(x) - min(x))*bins
  # fitted density curve
  df2 <- rootcount_density_fitted(par, n = n)
  df2 <- data.frame(
    bundle_id = c(1, df2$bundle_id, max(df2$bundle_id)),
    x = c(df2$x[1], df2$x, utils::tail(df2$x, 1)),
    y = c(0, df2$y, 0)
  )
  # axes limits
  xlims <- round_range(c(df1$x, df2$x), lim = xlim, ticks = ticks)
  ylims <- round_range(c(df1$y, df2$y), lim = ylim, ticks = ticks)
  # plot
  ggplot2::ggplot() +
    ggplot2::theme_bw() +
    ggplot2::geom_col(
      data = df1,
      ggplot2::aes(x = .data$x, y = .data$y),
      width = (max(x) - min(x))/bins,
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


#' Initial guess for multisegment power-law fitting
#'
#' @description
#' Generates a single initial guess for multisegment power-law fitting.
#' The number of optional breakpoints is `ns - 1 + n0`, which are equally
#' spaced along the domain `xmin` to `xmax`. A fit is made for each of the
#' combinations of the breakpoints (breakpoints are fixed in this analysis, so
#' only power-laws are fitted), assuming all power-law coefficients are equal
#' to zero.
#'
#' The combinations of breakpoints with the largest likelihood score is
#' returned. This will serve as a decent initial guess for the root count
#' fitting function, in which the exact locations of the breakpoints may also
#' vary.
#'
#' Larger values of `n0` will probably result in a better guess. This however
#' also results in longer runtimes, as more options need to be evaluated.
#'
#' @inheritParams rootcount_fit
#' @return a vector with an initial guess
#' @export
#' @examples
#' rootcount_initialguess(seq(1.2, 5.6, l = 50), 3)
#' rootcount_initialguess(seq(1.2, 5.6, l = 50), 3, n0 = 3)
#'
rootcount_initialguess <- function(
    x,
    ns,
    xmin = min(x, na.rm = TRUE),
    xmax = max(x, na.rm = TRUE),
    n0 = 0,
    fixed = rep(NA, 2*ns - 1),
    weights = rep(1, length(x))
) {
  # single segment - only fitting coefficient is power-coefficient
  if (ns == 1) {
    if (is.na(fixed[1])) {
      par0 <- 0
    } else {
      par0 <- fixed[1]
    }
    # multiple segments
  } else {
    # matrix with breakpoints options
    fixed_xb <- fixed[1:(ns - 1)]
    if (any(is.na(fixed_xb))) {
      # some breakpoints unknown - generate options
      xb <- breakpoint_options(
        ns - 1,
        xmin,
        xmax,
        fixed = fixed_xb,
        n0 = n0
      )
      xb <- matrix(xb[, is.na(fixed_xb)], nrow = nrow(xb))
      # number of unfixed power-coefficients
      nb_unfixed <- sum(is.na(fixed[ns:(2*ns - 1)]))
      # fit for each of the breakpoint options
      Lopts <- apply(
        xb,
        1,
        function(xbi) {
          rootcount_fit(
            x, ns,
            xmin = xmin, xmax = xmax,
            weights = weights,
            guess = rep(0, nb_unfixed),
            fixed = c(xbi, fixed[ns:(2*ns - 1)])
          )
        }
      )
      # find option with minimum loglikelihood score
      i <- which.max(sapply(Lopts, function(x) x$loglikelihood))
      # vector with all initial parameters
      parall <- c(Lopts[[i]]$par$xmin[2:ns], Lopts[[i]]$par$power)
    } else {
      # all breakpoints fixed - return simple guess - all non-fixed powers = 0
      parall <- rep(0, 2*ns + 1)
    }
    # return vector with non-fixed values
    par0 <- parall[is.na(fixed)]
  }
  # return guess
  par0
}


#' Combine fitting parameters and known parameters into single vector
#'
#' @description
#' Generate a single vector with fitting parameters and known ('fixed')
#' parameters for root count fitting.
#'
#' @param par vector with fitting parameters
#' @param fixed vector with known parameters. `NA` values in this vector
#'   indicate unknwon values (i.e., those in `par`).
#' @return vector with all coefficients
#' @export
#' @examples
#' rootcount_constructpar(c(1,2,3), fixed = c(NA, 5, NA, NA, 6))
#'
rootcount_constructpar <- function(par, fixed = rep(NA, length(par))) {
  fixed[is.na(fixed)] <- par
  fixed
}


#' log-transformed probability for multi-segment power-law
#'
#' @description
#' Returns log-transformed probabilities for observations x. Probabilities are
#' scaled by an unknown multiplier, which is seperately determined (function
#' `rootcount_multiplier()`) based on the fact that the total probability
#' must equal one.
#'
#' @param logxb vector with log-transformed values of the x-breakpoints between
#'   segments
#' @param b vector with power-law power coefficients for each segment
#' @param logx vector with log-transformed observations x
#' @return vector with scaled log-transformed probabilities for each
#'   observation in `x`
#' @export
#' @examples
#' x <- seq(1, 4, l = 14)
#' b <- c(1, 2, 3)
#' xb <- c(2.2, 3.6)
#' rootcount_logprob(log(xb), b, log(x))
#'
rootcount_logprob <- function(logxb, b, logx) {
  # logx - logxb, matrix form
  M <- outer(logx, logxb, FUN = "-")
  # Heavyside function
  H <- (M >= 0)
  # difference in b
  db <- diff(b)
  # weighted log-transformed probability
  b[1]*logx + colSums(db*t(M*H))
}


#' Derivative of function `rootcount_logprob()`
#'
#' @description
#' Generates the derivative of the results of the function `rootcount_logprob()`
#' with respect to its input parameters `logxb` and `b`.
#'
#' Function is vectorised.
#'
#' @inheritParams rootcount_logprob
#' @return a list containing the derivatives for each input parameter. Has
#'   fields `logxb` and `b`.
#' @export
#' @examples
#' # Test jacobian by comparing to numberical solution
#' logxb <- c(1.25, 1.35)
#' b <- c(2.5, 3.1, -2.5)
#' logx <- seq(1.2, 1.5, l = 9)
#'
#' p0 <- rootcount_logprob(logxb, b, logx)
#' J <- rootcount_logprob_jacobian(logxb, b, logx)
#'
#' eps <- 1e-6
#' Jb <- matrix(0, nrow = length(logx), ncol = length(b))
#' for (i in 1:length(b)) {
#'   dx <- rep(0, length(b))
#'   dx[i] <- eps
#'   p1 <- rootcount_logprob(logxb, b + dx, logx)
#'   Jb[, i] <- (p1 - p0)/eps
#' }
#'
#' Jlogxb <- matrix(0, nrow = length(logx), ncol = length(logxb))
#' for (i in 1:length(logxb)) {
#'   dx <- rep(0, length(logxb))
#'   dx[i] <- eps
#'   p1 <- rootcount_logprob(logxb + dx, b, logx)
#'   Jlogxb[, i] <- (p1 - p0)/eps
#' }
#'
#' J
#' Jlogxb
#' Jb
#'
rootcount_logprob_jacobian <- function(logxb, b, logx) {
  # number of segments and observations
  ns <- length(b)
  nx <- length(logx)
  # logx - logxb, matrix form
  M <- outer(logx, logxb, FUN = "-")
  # Heavyside function
  H <- (M >= 0)
  # derivatives with respect to b
  Jb <- matrix(0, nrow = nx, ncol = ns)
  Jb[, 1] <- Jb[, 1] + logx
  Jb[, 2:ns] <- Jb[, 2:ns] + M*H
  Jb[, 1:(ns - 1)] <- Jb[, 1:(ns - 1)] - M*H
  # derivatives with respect to logxb
  Jlogxb <- t(-diff(b)*t(H))
  # return list with derivatives
  list(
    logxb = Jlogxb,
    b = Jb
  )
}


#' Scaling multipliers for individual probabilities
#'
#' @description
#' Calculates the scaling factor that should be applied to individual
#' probabilities to ensure that the total probability equals one.
#'
#' Returns a vector, which, when summed, gives the value of the multiplier.
#'
#' Probabilities  calculated using `rootcount_logprob()` should be divided
#' by this multiplier to get the actual probability.
#'
#' @inheritParams rootcount_fit
#' @param xb vector with x-breakpoints between segments
#' @param b vector with power-law power coefficients for each segment
#' @return vector with multiplier contributions of each segment, which, when
#'   summed give the probability multiplier
#' @export
#' @examples
#' # TEST - does probability sum to one?
#' # parameters
#' x <- seq(1.5, 6.5, l = 1000)
#' xb <- c(2, 3.2, 5)
#' b <- c(1.1, 3, -2, -0.4)
#' xmin <- min(x)
#' xmax <- max(x)
#'
#' # individual probabilities
#' pip <- exp(rootcount_logprob(log(xb), b, log(x)))
#' # integral prod
#' pt <- sum(rootcount_multiplier(xb, b, xmin, xmax))
#' # quick estmation of integral
#' sum(pip/pt)*diff(x)[1]
#'
rootcount_multiplier <- function(xb, b, xmin, xmax) {
  # number of segments
  ns <- length(b)
  # boundaries
  g <- c(xmin, xb, xmax)
  # integrals
  I <- power_integrate(b, g[1:ns], g[2:(ns + 1)])
  # multiplication constants
  m <- cumprod(c(1, xb^(-diff(b))))
  # return total per segment
  m*I
}


#' Derivative of function `rootcount_multiplier()`
#'
#' @description
#' Generates the derivative of the results of the function
#' `rootcount_multiplier()` with respect to its input parameters `xb` and `b`.
#'
#' Function is vectorised.
#'
#' @inheritParams rootcount_multiplier
#' @return a list containing the derivatives for each input parameter. Has
#'   fields `xb` and `b`.
#' @export
#' @examples
#' # Test jacobian by comparing to numberical solution
#' xb <- c(3, 4)
#' b <- c(-1.2, -0.2, 1.3)
#' xmin <- 2.5
#' xmax <- 5.5
#'
#' t <- rootcount_multiplier(xb, b, xmin, xmax)
#' J <- rootcount_multiplier_jacobian(xb, b, xmin, xmax)
#'
#' eps <- 1e-6
#' Jxb <- matrix(0, nrow = length(b), ncol = length(xb))
#' for (i in 1:length(xb)) {
#'   dx <- rep(0, length(xb))
#'   dx[i] <- eps
#'   Jxb[, i] <- (rootcount_multiplier(xb + dx, b, xmin, xmax) - t)/eps
#' }
#'
#' Jb <- matrix(0, nrow = length(b), ncol = length(b))
#' for (i in 1:length(b)) {
#'   dx <- rep(0, length(b))
#'   dx[i] <- eps
#'   Jb[, i] <- (rootcount_multiplier(xb, b + dx, xmin, xmax) - t)/eps
#' }
#'
#' J
#' Jxb
#' Jb
rootcount_multiplier_jacobian <- function(xb, b, xmin, xmax) {
  # number of segments
  ns <- length(b)
  # boundaries
  g <- c(xmin, xb, xmax)
  # integrals
  I <- power_integrate(b, g[1:ns], g[2:(ns + 1)])
  JI <- power_integrate_jacobian(b, g[1:ns], g[2:(ns + 1)])
  if (ns == 2) {
    I_xb <- rbind(JI$upper[1:(ns - 1)], JI$lower[2:ns])
  } else {
    I_xb <- rbind(0, diag(JI$lower[2:ns])) + rbind(diag(JI$upper[1:(ns - 1)]), 0)
  }
  I_b <- diag(JI$power)
  # multiplication constants
  k <- xb^(-diff(b))
  m <- c(1, cumprod(k))
  if (ns == 2) {
    m_xb <- rbind(0, -diff(b)*xb^(-diff(b) - 1))
    m_b <- rbind(0, c(1, -1)*log(xb)*xb^(-diff(b)))
  } else {
    k_xb <- diag(-diff(b)*xb^(-diff(b) - 1))
    k_b <- cbind(diag(log(xb)*k), 0) + cbind(0, -diag(log(xb)*k))
    m_k_temp <- outer(cumprod(k), k, FUN = "/")
    m_k <- m_k_temp*lower.tri(m_k_temp, diag = TRUE)
    m_xb <- rbind(0, m_k %*% k_xb)
    m_b <- rbind(0, m_k %*% k_b)
  }
  # return total per segment
  list(
    xb = diag(I) %*% m_xb + diag(m) %*% I_xb,
    b = diag(I) %*% m_b + diag(m) %*% I_b
  )
}


#' Loglikelihood score of multi-segment power-law fit
#'
#' @description
#' Return the (weighted) loglikelihood score for a multi-segment
#' power-law fit.
#'
#' @inheritParams rootcount_fit
#' @param par vector with fitting parameters. For n segments, the vector will
#'   consist of n-1 consecutive x-breakpoints, followed by n power-law
#'   cower coefficients
#' @return scalar with weighted loglikelihood value
#' @export
#' @examples
#' # parameters
#' x <- seq(1.5, 6.5, l = 25)
#' xb <- c(2, 3.2, 5)
#' b <- c(1.1, 3, -2, -0.4)
#'
#' # loglikelihood
#' rootcount_loglikelihood(c(xb, b), x)
#'
rootcount_loglikelihood <- function(
    par,
    x,
    xmin = min(x, na.rm = TRUE),
    xmax = max(x, na.rm = TRUE),
    weights = rep(1, length(x)),
    fixed = rep(NA, length(par))
) {
  # generate single vector
  parall <- rootcount_constructpar(par, fixed = fixed)
  # split input
  ns <- (length(parall) + 1)/2
  xb <- parall[1:(ns - 1)]
  b <- parall[ns:(2*ns - 1)]
  # individual probabilities
  lpi <- rootcount_logprob(log(xb), b, log(x))
  # log of sums of integrals per segment
  lpt <- log(sum(rootcount_multiplier(xb, b, xmin, xmax)))
  # return loglikeihood
  sum(weights*lpi) - sum(weights)*lpt
}


#' Derivative of function `rootcount_loglikelihood()`
#'
#' @description
#' Generates the derivative of the results of the function
#' `rootcount_loglikelihood()` with respect to fitting value input
#' `par`.
#'
#' Function is vectorised.
#'
#' @inheritParams rootcount_loglikelihood
#' @return a scalar containing the derivatives for each value in fitting vector
#'   `par`
#' @export
#' @examples
#' # Test jacobian by comparing to numberical solution
#' x <- seq(1.4, 6.9, l = 25)
#' xb <- c(2, 3.2, 5.1)
#' b <- c(1.1, 2.6, -2.1, -0.4)
#' par <- c(xb, b)
#' w <- seq(0.8, 1.2, l = length(x))
#'
#' L <- rootcount_loglikelihood(par, x, weights = w)
#' J <- rootcount_loglikelihood_jacobian(par, x, weights = w)
#'
#' eps <- 1e-6
#' J2 <- rep(0, length(par))
#' for (i in 1:length(par)) {
#'   dx <- rep(0, length(par))
#'   dx[i] <- eps
#'   J2[i] <- (rootcount_loglikelihood(par + dx, x, weights = w) - L)/eps
#' }
#'
#' J
#' J2
#'
rootcount_loglikelihood_jacobian <- function(
    par,
    x,
    xmin = min(x, na.rm = TRUE),
    xmax = max(x, na.rm = TRUE),
    weights = rep(1, length(x)),
    fixed = rep(NA, length(par))
) {
  # generate single vector
  parall <- rootcount_constructpar(par, fixed = fixed)
  # split input
  ns <- (length(parall) + 1)/2
  xb <- parall[1:(ns - 1)]
  b <- parall[ns:(2*ns - 1)]
  # individual probabilities
  Jlpi <- rootcount_logprob_jacobian(log(xb), b, log(x))
  # log of sums of integrals per segment
  pti <- rootcount_multiplier(xb, b, xmin, xmax)
  Jpti <- rootcount_multiplier_jacobian(xb, b, xmin, xmax)
  # return loglikeihood
  c(
    colSums(weights*t(1/xb*t(Jlpi$logxb))) -
      sum(weights)*colSums(Jpti$xb)/sum(pti),
    colSums(weights*Jlpi$b) -
      sum(weights)*colSums(Jpti$b)/sum(pti)
  )[is.na(fixed)]
}


#' Generate contraints for `rootcount_fit()` optimalisation
#'
#' @description
#' Generate constraint matrix `ui` and constraint vector `ci` to use in
#' the fitting function `rootcount_fit()`. This is used when there are more
#' than two segments, in which case the `stats::constrOptim()` function is used
#' which required these constraints.
#'
#' See `stats::constrOptim()` for more detail.
#'
#' @inheritParams rootcount_fit
#' @return a list with the constraint matrix (`ui`) and vector (`ci`)
#' @export
#' @examples
#' rootcount_constraints(4, 2, 5)
#'
rootcount_constraints <- function(ns, xmin, xmax, fixed = rep(NA, 2*ns - 1)) {
  # generate matrix
  ui <- matrix(
    c(1, -1, rep(c(rep(0, ns - 1), 1, -1), ns - 2), rep(0, ns^2)),
    nrow = ns,
    ncol = 2*ns - 1,
    byrow = FALSE
  )
  # generate vector
  ci <- c(xmin, rep(0, ns - 2), -xmax)
  # account for fixed vectors
  # add known values to coefficients in ci - and remove from matrix
  free <- is.na(fixed)
  ci <- ci - (as.matrix(ui[, !free], ncol = sum(!free)) %*% fixed[!free])
  ui <- as.matrix(ui[, free], ncol = sum(free))
  # return
  list(ui = ui, ci = ci)
}


#' Loglikelihood for fitting a single segment
#'
#' @description
#' Calculate the loglikelihood score for a single-segment power-law
#' fit. In this case, the only fitting parameter is a single power-law power
#' coefficient.
#'
#' @inheritParams rootcount_loglikelihood
#' @return scalar with weighted loglikelihood value
#' @export
#' @examples
#' rootcount_loglikelihood_single(2, seq(2, 4, l = 10))
#'
rootcount_loglikelihood_single <- function(
    par,
    x,
    xmin = min(x, na.rm = TRUE),
    xmax = max(x, na.rm = TRUE),
    weights = rep(1, length(x))
) {
  # individual probabilities
  lpi <- par*log(x)
  # total of integral
  pt <- power_integrate(par, xmin, xmax)
  # likelihood
  sum(weights*lpi) - sum(weights)*log(pt)
}


#' Derivative of function `rootcount_loglikelihood_single()`
#'
#' @description
#' Generates the derivative of the results of the function
#' `rootcount_loglikelihood_single()` with respect to fitting value input
#' `par`.
#'
#' Function is vectorised.
#'
#' @inheritParams rootcount_loglikelihood_single
#' @return a scalar containing the derivatives for each value in fitting vector
#'   `par`
#' @export
#' @examples
#' # Test jacobian by comparing to numberical solution
#' x <- seq(2.2, 6.5, l = 51)
#' par <- 2.2
#' w <- seq(0.8, 1.4, l = length(x))
#'
#' J <- rootcount_loglikelihood_single_jacobian(par, x, weights = w)
#'
#' eps <- 1e-6
#' J2 <- c(rootcount_loglikelihood_single(par + eps, x, weights = w) -
#'   rootcount_loglikelihood_single(par, x, weights = w))/eps
#'
#' J
#' J2
#'
rootcount_loglikelihood_single_jacobian <- function(
    par,
    x,
    xmin = min(x, na.rm = TRUE),
    xmax = max(x, na.rm = TRUE),
    weights = rep(1, length(x))
) {
  # derivative: individual probabilities
  lpi_par <- log(x)
  # derivative: total of integral
  pt <- power_integrate(par, xmin, xmax)
  pt_par <- power_integrate_jacobian(par, xmin, xmax)$power
  # likelihood
  sum(weights*lpi_par) - sum(weights)*pt_par/pt
}


#' Generate a cumulative probability density curve for observations
#'
#' @description
#' Generate a observation-value versus cumulative probability density
#' curve for a series of observations
#'
#' @inheritParams rootcount_fit
#' @return a dataframe with x-values (`x`) and the cumulative probability
#'   density (`y`)
#' @export
#' @examples
#' rootcount_cumulative_observations(seq(2, 5, l = 10))
#'
rootcount_cumulative_observations <- function(x) {
  data.frame(
    x = rep(sort(x), each = 2),
    y = rep(seq(0, 1, l = length(x) + 1), each = 2)[2:(2*length(x) + 1)]
  )
}


#' Generate a cumulative probability density curve for multisegment fit
#'
#' @description
#' Generate a observation versus cumulative probability density
#' curve for a multisegment power-law fit, as generated by the
#' fitting function `rootcount_fit()`.
#'
#' @param par dataframe with fitting parameters per power-law segment. This
#'   is the dataframe returned by the function `rootcount_fit` as the field
#'   `par` in the outputted list
#' @param n number of points to use on curve for each segment
#' @return a dataframe with x-values (`x`), the cumulative probability
#' density (`y`), and an identifier for each bundle (`bundle_id`, integers
#' 1, 2, 3 etc)
#' @export
#' @examples
#' ft <- rootcount_fit(seq(2, 6, l = 51), 2)
#' rootcount_cumulative_fitted(ft$par, n = 10)
#'
rootcount_cumulative_fitted <- function(par, n = 101) {
  # start and end values
  par$bundle_id <- seq(1, nrow(par))
  par$ymax <- cumsum(par$total)
  par$ymin <- c(0, par$ymax)[1:nrow(par)]
  # get all x-coordinates
  df <- expand_grid_df(data.frame(s = seq(0, 1, l = n)), par)
  df$x <- df$xmin + (df$xmax - df$xmin)*df$s
  # get y-coordinates
  df$y <- with(df, ymin + (ymax - ymin)*ifelse(
    is_near(power, -1),
    log(x/xmin)/log(xmax/xmin),
    (x^(power + 1) - xmin^(power + 1))/(xmax^(power + 1) - xmin^(power + 1))
  ))
  # return
  df[, c("bundle_id", "x", "y")]
}


#' Generate a probability density curve for multisegment fit
#'
#' @description
#' Generate a observation versus probability density
#' curve for a multisegment power-law fit, as generated by the
#' fitting function `rootcount_fit()`.
#'
#' @param par dataframe with fitting parameters per power-law segment. This
#'   is the dataframe returned by the function `rootcount_fit` as the field
#'   `par` in the outputted list
#' @param n number of points to use on curve for each segment
#' @return a dataframe with x-values (`x`), the cumulative probability
#' density (`y`), and an identifier for each bundle (`bundle_id`, integers
#' 1, 2, 3 etc)
#' @export
#' @examples
#' ft <- rootcount_fit(seq(2, 6, l = 51), 2)
#' df <- rootcount_density_fitted(ft$par, n = 101)
#' plot(df$x, df$y, "l")
#'
rootcount_density_fitted <- function(par, n = 101) {
  # assign bundle id
  par$bundle_id <- seq(1, nrow(par))
  # get multiplication coefficients
  par$alpha <- with(par, ifelse(
    is_near(power, -1),
    total/log(xmax/xmin),
    total*(1 + power)/(xmax^(power + 1) - xmin^(power + 1))
  ))
  # get all x-coordinates
  df <- expand_grid_df(data.frame(s = seq(0, 1, l = n)), par)
  df$x <- df$xmin + (df$xmax - df$xmin)*df$s
  # get y-coordinates
  df$y <- with(df, alpha*x^power)
  # return
  df[, c("bundle_id", "x", "y")]
}


#' Predict fitted cumulative probability density for observations
#'
#' @description
#' Calculate the fitted cumulative probability density for observations `x`,
#' given the fit results `par` generated by the function `rootcount_fit()`
#' or `rootclass_fit()`.
#'
#' @inheritParams rootcount_density_fitted
#' @param x array with observations
#' @return array with cumulative probability densities
#' @export
#' @examples
#' x <- rweibull(100, 4, 6)
#' par <- rootcount_fit(x, 2)$par
#'
#' xp <- seq(min(x), max(x), l = 251)
#' yp <- rootcount_cumulative_prediction(xp, par)
#'
#' plot(xp, yp, "l")
#'
rootcount_cumulative_prediction <- function(x, par) {
  # all breakpoints
  xb <- c(-Inf, par$xmin[2:nrow(par)], Inf)
  # find domain for each observation
  i <- cut(x, xb, labels = FALSE)
  # cumulative total
  totalcum <- c(0, cumsum(par$total))
  # cumulative probability densities
  with(par, totalcum[i] + total[i]*ifelse(
    is_near(power[i], -1),
    log(x/xmin[i])/log(xmax[i]/xmin[i]),
    (x^(power[i] + 1) - xmin[i]^(power[i] + 1))/(xmax[i]^(power[i] + 1) - xmin[i]^(power[i] + 1))
  ))
}


#' Predict fitted probability density for observations
#'
#' @description
#' Calculate the fitted probability density for observations `x`,
#' given the fit results `par` generated by the function `rootcount_fit()`
#'
#' @inheritParams rootcount_density_fitted
#' @param x array with observations
#' @return array with probability densities
#' @export
#' @examples
#' x <- rweibull(100, 4, 6)
#' par <- rootcount_fit(x, 2)$par
#'
#' xp <- seq(min(x), max(x), l = 251)
#' yp <- rootcount_density_prediction(xp, par)
#'
#' plot(xp, yp, "l")
#'
rootcount_density_prediction <- function(x, par) {
  # all breakpoints
  xb <- c(-Inf, par$xmin[2:nrow(par)], Inf)
  # find domain for each observation
  i <- cut(x, xb, labels = FALSE)
  # probability densities
  par$multiplier[i]*x^par$power[i]
}


#' Generate random root diameters based on multisegment power-law probabilities
#'
#' @description
#' Generate random root diameter observations based on a multisegment power-law
#' fit result.
#'
#' @inheritParams rootcount_density_fitted
#' @param n number of random observations required
#' @return a vector with randomly drawn observations
#' @export
#' @examples
#' # generate some data and fit
#' x1 <- rweibull(100, shape = 4, scale = 6)
#' par <- rootcount_fit(x1, 2)$par
#'
#' # draw random values from fit
#' x2 <- rootcount_random(1000, par)
#'
#' # plot cumulative data
#' c1 <- seq(0, 1, l = length(x1))
#' c2 <- seq(0, 1, l = length(x2))
#' plot(sort(x1), c1)
#' lines(sort(x2), c2, "l", col = "red")
#'
rootcount_random <- function(n, par) {
  # draw uniform distribution -> cumulative
  p <- stats::runif(n, min = 0, max = 1)
  # find segment
  i <- cut(p, c(-Inf, cumsum(par$total)), labels = FALSE)
  # cumulative probability at the start of eachs segment
  pc <- c(0, cumsum(par$total))
  # convert power-law cumulative to observations
  with(par, ifelse(
    is_near(power[i], -1),
    xmin[i]*(xmax[i]/xmin[i])^((p - pc[i])/total[i]),
    ((p - pc[i])/total[i]*(xmax[i]^(power[i] + 1) - xmin[i]^(power[i] + 1)) + xmin[i]^(power[i] + 1))^(1/(power[i] + 1))
  ))
}
