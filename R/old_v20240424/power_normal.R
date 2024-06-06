# Fit power-law with normal distributions of residuals
# 26/03/2024 - G. J. Meijer


#' Power law fitting with normally distributed residuals
#'
#' @description
#' Fit a power-law curve to a series of `x`, `y` data. The mean of the
#' data follows a power-law in the form:
#'
#'   mu = multiplier*x^exponent
#'
#' The residuals are assumed normally distributed, with a mean of zero and
#' a standard deviation sd in the form:
#'
#'   sd = sd_multiplier*x^(sd_exponent)
#'
#' When doing non-linear least squares fitting, `sd_exponent = 0` (default
#' setting).
#'
#' The problem is solved by maximising the (weight) likelihood function.
#' Individual probabilities for each observation are given by
#'
#'   p = 1/(sd*sqrt(2*pi))*exp(-0.5*((y - mu)/sd)^2)
#'
#' and the weighted logliklihood:
#'
#'   log L = sum w*log(p)
#'
#' where `w` is the weighting for each observation (default, w = 1).
#'
#' The problem can be written in terms of unknown `exponent` and `sd_exponent`
#' only, through assuming that at the maximum likelihood the partial
#' derivatives of the likelihood function are equal to zero.
#'
#' Unknowns are found using a root solving algorithm using gradient descent,
#' as implemented in the function `rootSolve::multiroot()`.
#'
#' @param x,y vectors with x and y observations. Must be the same length, and
#'   contain numerical values only
#' @param weights weighting parameter for each observation
#' @param sd_exponent if a numeric value, the power coefficients for standard
#'   deviations is fixed at that value. If `sd_exponent = NA`, it is assumed to
#'   be an independent parameter that must be fitted. In all other cases, it is
#'   assumed that `sd_exponent = exponent`.
#' @param method choose `newton` for gradient descent solving, using the
#'   `rootSolve::multiroot()` function, or `bisection` for bisection root
#'   solving algorithm using the `stats::uniroot()` function.
#' @param range two-value array to add to best guess, to define initial
#'   interval for bisection algorithm
#' @return a list with the loglikelihood of the fit (`loglikelihood`), the
#'   multiplier and power coefficient for the mean (`multiplier`, `exponent`)
#'   and the multiplier and power coefficient for the standard deviation
#'   (`sd_multiplier`, `sd_exponent`). The returned likelihood is defined as
#'   the loglikelihood of (sd_multiplier/multiplier)*(y - mu)/sd, i.e. the
#'   likelihood of (y/multiplier - x^exponent)/(x^sd_exponent)
#' @examples
#' # parameters
#' multiplier <- 50
#' exponent <- -0.3
#' sd_multiplier <- 15
#' sd_exponent <- -0.6
#'
#' # generate data
#' n <- 51
#' x <- seq(2, 10, l = n)
#' y <- abs(multiplier*x^exponent +
#'   sd_multiplier*stats::rnorm(length(x))*x^sd_exponent)
#'
#' # fit
#' ft1 <- power_normal_fit(x, y, sd_exponent = 0)
#' ft2 <- power_normal_fit(x, y, sd_exponent = NA)
#' ft3 <- power_normal_fit(x, y, sd_exponent = "scaled")
#'
#' # prediction
#' xp <- seq(min(x), max(x), l = 251)
#' yp1 <- ft1$multiplier*xp^ft1$exponent
#' yp2 <- ft2$multiplier*xp^ft2$exponent
#' yp3 <- ft3$multiplier*xp^ft3$exponent
#'
#' # plot
#' plot(x, y)
#' lines(xp, yp1, col = "red")
#' lines(xp, yp2, col = "blue")
#' lines(xp, yp3, col = "darkgreen")
#'
power_normal_fit <- function(
    x,
    y,
    sd_exponent = 0,
    weights = rep(1, length(x)),
    method = "newton",
    range = c(-1, 1)
) {
  # fit powers
  if (is.numeric(sd_exponent)) {
    # FIT WITH KNOWN SD POWER COEFFICIENT
    # initial guess
    exponent0 <- as.numeric(
      stats::lm(log(y) ~ log(x), weights = weights)$coef[2]
    )
    # fit
    if (method == "newton") {
      exponent <- rootSolve::multiroot(
        power_normal_root_fixd,
        start = exponent0,
        jacfunc = power_normal_root_fixd_jacobian,
        x = x,
        y = y,
        sd_exponent = sd_exponent,
        weights = weights
      )$root
    } else if (method == "bisection") {
      exponent <- stats::uniroot(
        power_normal_root_fixd,
        exponent0 + range,
        extendInt = "downX",
        x = x,
        y = y,
        sd_exponent = sd_exponent,
        weights = weights
      )$root
    } else {
      stop("`method` not recognised")
    }
  } else if (is.na(sd_exponent)) {
    # FIT BOTH POWER COEFFICIENTS
    # initial guess - solve with assumed power coefficients
    exponent0 <- as.numeric(
      stats::lm(log(y) ~ log(x), weights = weights)$coef[2]
    )
    if (method == "newton") {
      sd_exponent0 <- rootSolve::multiroot(
        power_normal_root_fixb,
        start = exponent0,
        jacfunc = power_normal_root_fixb_jacobian,
        x = x,
        y = y,
        exponent = exponent0,
        weights = weights
      )$root
    } else if (method == "bisection") {
      sd_exponent0 <- stats::uniroot(
        power_normal_root_fixd,
        exponent0 + range,
        extendInt = "downX",
        x = x,
        y = y,
        exponent = exponent0,
        weights = weights
      )$root
    } else {
      stop("`method` not recognised")
    }
    # fit powers
    sol <- rootSolve::multiroot(
      power_normal_root_freebd,
      start = c(exponent0, sd_exponent0),
      jacfunc = power_normal_root_freebd_jacobian,
      x = x,
      y = y,
      weights = weights
    )
    exponent <- sol$root[1]
    sd_exponent <- sol$root[2]
  } else {
    ## FIT ASSUMING BOTH POWERS ARE EQUAL
    # initial guess
    exponent0 <- as.numeric(
      stats::lm(log(y) ~ log(x), weights = weights)$coef[2]
    )
    # fit
    if (method == "newton") {
      exponent <- rootSolve::multiroot(
        power_normal_root_linkedbd,
        start = exponent0,
        jacfunc = power_normal_root_linkedbd_jacobian,
        x = x,
        y = y,
        weights = weights
      )$root
    } else if (method == "bisection") {
      exponent <- stats::uniroot(
        power_normal_root_linkedbd,
        exponent0 + range,
        extendInt = "downX",
        x = x,
        y = y,
        weights = weights
      )$root
    } else {
      stop("`method` not recognised")
    }
    sd_exponent <- exponent
  }
  # calculate other parameters
  c1 <- sum(weights)
  c3 <- sum(weights*x^(2*exponent - 2*sd_exponent))
  c7 <- sum(weights*y*x^(exponent - 2*sd_exponent))
  c11 <- sum(weights*y^2*x^(-2*sd_exponent))
  multiplier <- c7/c3
  sd_multiplier <- sqrt((c11 - multiplier*c7)/c1)
  # likelihood calculations
  yp <- (y/multiplier - x^exponent)/(x^sd_exponent)
  logpi <- stats::dnorm(yp, mean = 0, sd = sd_multiplier/multiplier, log = TRUE)
  # return
  list(
    loglikelihood = sum(weights*logpi),
    multiplier = multiplier,
    exponent = exponent,
    sd_multiplier = sd_multiplier,
    sd_exponent = sd_exponent
  )
}


#' Roots for power law + normal residuals - fixed sd_exponent
#'
#' @description
#' Function that returns roots to solved for power-law fitting with normally
#' distributed residuals, assuming a fixed and known value for the exponent for
#' the standard deviations
#'
#' @inheritParamts power_normal_fit
#' @param exponent unknown power coefficient
#' @return root
#'
power_normal_root_fixd <- function(
    exponent,
    x,
    y,
    sd_exponent = 0,
    weights = rep(1, length(x))
) {
  # coefficients
  c1 <- sum(weights)
  c3 <- sum(weights*x^(2*exponent - 2*sd_exponent))
  c4 <- sum(weights*x^(2*exponent - 2*sd_exponent)*log(x))
  c7 <- sum(weights*y*x^(exponent - 2*sd_exponent))
  c8 <- sum(weights*y*x^(exponent - 2*sd_exponent)*log(x))
  c11 <- sum(weights*y^2*x^(-2*sd_exponent))
  # multipliers and their derivatives
  a <- c7/c3
  g <- (c11 - a*c7)/c1
  da_dd <- 2*(a*c4 - c8)/c3
  dg_db <- a*c3/c1*da_dd
  # return root
  -0.5*c1/g*dg_db
}


#' Jacobian of roots for power law + normal residuals - fixed sd_exponent
#'
#' @description
#' Function that returns the jacobian of the roots to solve for power-law
#' fitting with normally distributed residuals, assuming a fixed and known
#' value for the exponent for standard deviations
#'
#' @inheritParamts power_normal_root_fixd
#' @return derivative of root as defined by function `power_normal_root_fixd()`
#' @examples
#' # parameters
#' multiplier <- 48
#' exponent <- -0.32
#' sd_multiplier <- 5.2
#' sd_exponent <- 0.2
#'
#' # generate data
#' n <- 251
#' x <- seq(2, 10, l = n)
#' y <- abs(multiplier*x^exponent + sd_multiplier*stats::rnorm(length(x))*x^sd_exponent)
#' w <- runif(n)
#'
#' # check derivative - compare analytical to numerical
#' eps <- 1e-6
#' v0 <- power_normal_root_fixd(exponent, x, y, sd_exponent = sd_exponent, weights = w)
#' v1 <- power_normal_root_fixd(exponent + eps, x, y, sd_exponent = sd_exponent, weights = w)
#' (v1 - v0)/eps
#' power_normal_root_fixd_jacobian(exponent, x, y, sd_exponent = sd_exponent, weights = w)
#'
power_normal_root_fixd_jacobian <- function(
    exponent,
    x,
    y,
    sd_exponent = 0,
    weights = rep(1, length(x))
) {
  # coefficients
  c1 <- sum(weights)
  c3 <- sum(weights*x^(2*exponent - 2*sd_exponent))
  c4 <- sum(weights*x^(2*exponent - 2*sd_exponent)*log(x))
  c5 <- sum(weights*x^(2*exponent - 2*sd_exponent)*log(x)^2)
  c7 <- sum(weights*y*x^(exponent - 2*sd_exponent))
  c8 <- sum(weights*y*x^(exponent - 2*sd_exponent)*log(x))
  c9 <- sum(weights*y*x^(exponent - 2*sd_exponent)*log(x)^2)
  c11 <- sum(weights*y^2*x^(-2*sd_exponent))
  # multipliers and their derivatives
  a <- c7/c3
  g <- (c11 - a*c7)/c1
  da_db <- (c8 - 2*a*c4)/c3
  da_dd <- 2*(a*c4 - c8)/c3
  dg_db <- a*c3/c1*da_dd
  d2g_db2 <- -2*c3/c1*da_db^2 + 2/c1*(2*a^2*c5 - a*c9)
  # return derivative of root with respect to `power`
  0.5*c1*(dg_db/g)^2 - 0.5*c1/g*d2g_db2
}


#' Roots for power law + normal residuals - free powers
#'
#' @description
#' Function that returns roots to solved for power-law fitting with normally
#' distributed residuals, assuming both the exponents for the mean
#' and standard deviation are independent and unknown.
#'
#' @inheritParamts power_normal_fit
#' @param par vector with unknown power coefficients `exponent` and `sd_exponent`
#' @return roots of both root equations
#'
power_normal_root_freebd <- function(
    par,
    x,
    y,
    weights = rep(1, length(x))
) {
  # split parameters
  exponent <- par[1]
  sd_exponent <- par[2]
  # coefficients
  c1 <- sum(weights)
  c2 <- sum(weights*log(x))
  c3 <- sum(weights*x^(2*exponent - 2*sd_exponent))
  c4 <- sum(weights*x^(2*exponent - 2*sd_exponent)*log(x))
  c7 <- sum(weights*y*x^(exponent - 2*sd_exponent))
  c8 <- sum(weights*y*x^(exponent - 2*sd_exponent)*log(x))
  c11 <- sum(weights*y^2*x^(-2*sd_exponent))
  c12 <- sum(weights*y^2*x^(-2*sd_exponent)*log(x))
  # multipliers and their derivatives
  a <- c7/c3
  g <- (c11 - a*c7)/c1
  da_dd <- 2*(a*c4 - c8)/c3
  dg_db <- a*c3/c1*da_dd
  dg_dd <- 2/c1*(-a^2*c4 + 2*a*c8 - c12)
  # roots
  dL_db <- -0.5*c1/g*dg_db
  dL_dd <- -0.5*c1/g*dg_dd - c2
  # return
  c(dL_db, dL_dd)
}


#' Jacobian of roots for power law + normal residuals - free powers
#'
#' @description
#' Function that returns the jacobian of the roots to solve for power-law
#' fitting with normally distributed residuals, assuming both the exponents
#' for the mean and standard deviation are independent
#' and unknown.
#'
#' @inheritParamts power_normal_root_freebd
#' @return matrix with partial derivatives of roots as defined by
#'   function `power_normal_root_freebd()`
#' @examples
#' # parameters
#' multiplier <- 48
#' exponent <- -0.32
#' sd_multiplier <- 5.2
#' sd_exponent <- 0.2
#'
#' # generate data
#' n <- 251
#' x <- seq(2, 10, l = n)
#' y <- abs(multiplier*x^exponent + sd_multiplier*stats::rnorm(length(x))*x^sd_exponent)
#' w <- runif(n)
#'
#' # check derivative by comparing analytical and numerical solutions
#' eps <- 1e-6
#' par <- c(exponent, sd_exponent)
#' v0 <- power_normal_root_freebd(par, x, y, weights = w)
#' v1 <- power_normal_root_freebd(par + c(eps, 0), x, y, weights = w)
#' v2 <- power_normal_root_freebd(par + c(0, eps), x, y, weights = w)
#' (cbind(v1, v2) - v0)/eps
#' power_normal_root_freebd_jacobian(par, x, y, weights = w)
#'
power_normal_root_freebd_jacobian <- function(
    par,
    x,
    y,
    weights = rep(1, length(x))
) {
  # split parameters
  exponent <- par[1]
  sd_exponent <- par[2]
  # coefficients
  c1 <- sum(weights)
  c3 <- sum(weights*x^(2*exponent - 2*sd_exponent))
  c4 <- sum(weights*x^(2*exponent - 2*sd_exponent)*log(x))
  c5 <- sum(weights*x^(2*exponent - 2*sd_exponent)*log(x)^2)
  c7 <- sum(weights*y*x^(exponent - 2*sd_exponent))
  c8 <- sum(weights*y*x^(exponent - 2*sd_exponent)*log(x))
  c9 <- sum(weights*y*x^(exponent - 2*sd_exponent)*log(x)^2)
  c11 <- sum(weights*y^2*x^(-2*sd_exponent))
  c12 <- sum(weights*y^2*x^(-2*sd_exponent)*log(x))
  c13 <- sum(weights*y^2*x^(-2*sd_exponent)*log(x)^2)
  # multipliers and their derivatives
  a <- c7/c3
  g <- (c11 - a*c7)/c1
  da_db <- (c8 - 2*a*c4)/c3
  da_dd <- 2*(a*c4 - c8)/c3
  dg_db <- a*c3/c1*da_dd
  dg_dd <- 2/c1*(-a^2*c4 + 2*a*c8 - c12)
  d2g_db2 <- -2*c3/c1*da_db^2 + 2/c1*(2*a^2*c5 - a*c9)
  d2g_dbdd <- -2*c3/c1*da_db*da_dd - 4/c1*(a^2*c5 - a*c9)
  d2g_dd2 <- 2*c3/c1*da_dd^2 + 4/c1*(a^2*c5 - 2*a*c9 + c13)
  # derivatives of roots
  d2L_db2 <- 0.5*c1*(dg_db/g)^2 - 0.5*c1/g*d2g_db2
  d2L_dbdd <- 0.5*c1/g^2*dg_db*dg_dd - 0.5*c1/g*d2g_dbdd
  d2L_dd2 <- 0.5*c1*(dg_dd/g)^2 - 0.5*c1/g*d2g_dd2
  # return
  matrix(c(d2L_db2, d2L_dbdd, d2L_dbdd, d2L_dd2), nrow = 2, ncol = 2)
}


#' Roots for power law + normal residuals - linked powers
#'
#' @description
#' Function that returns roots to solved for power-law fitting with normally
#' distributed residuals, assuming the exponent for the standard deviation is
#' equal to the exponent for the mean.
#'
#' @inheritParamts power_normal_fit
#' @param exponent unknown exponent
#' @return root
#'
power_normal_root_linkedbd <- function(
    exponent,
    x,
    y,
    weights = rep(1, length(x))
) {
  # coefficients
  c1 <- sum(weights)
  c2 <- sum(weights*log(x))
  c7 <- sum(weights*y*x^(-exponent))
  c8 <- sum(weights*y*x^(-exponent)*log(x))
  c11 <- sum(weights*y^2*x^(-2*exponent))
  c12 <- sum(weights*y^2*x^(-2*exponent)*log(x))
  # multipliers and their derivatives
  a <- c7/c1
  g <- (c11 - a*c7)/c1
  dg_db <- 2*(a*c8 - c12)/c1
  # return root
  -0.5*c1/g*dg_db - c2
}


#' Jacobian of roots for power law + normal residuals - linked powers
#'
#' @description
#' Function that returns the jacobian of the roots to solve for power-law
#' fitting with normally distributed residuals, assuming the exponent for
#' the standard deviation is equal to the exponent for the mean.
#'
#' @inheritParamts power_normal_root_linkedbd
#' @return derivative of root as defined by function
#'   `power_normal_root_linkedbd()`
#' @examples
#' # parameters
#' multiplier <- 48
#' exponent <- -0.32
#' sd_multiplier <- 5.2
#' sd_exponent <- exponent
#'
#' # generate data
#' n <- 251
#' x <- seq(2, 10, l = n)
#' y <- abs(multiplier*x^exponent + sd_multiplier*stats::rnorm(length(x))*x^sd_exponent)
#' w <- runif(n)
#'
#' # check derivative - compare analytical and numerical solutions
#' eps <- 1e-6
#' v0 <- power_normal_root_linkedbd(exponent, x, y, weights = w)
#' v1 <- power_normal_root_linkedbd(exponent + eps, x, y, weights = w)
#' (v1 - v0)/eps
#' power_normal_root_linkedbd_jacobian(exponent, x, y, weights = w)
#'
power_normal_root_linkedbd_jacobian <- function(
    exponent,
    x,
    y,
    weights = rep(1, length(x))
) {
  # coefficients
  c1 <- sum(weights)
  c2 <- sum(weights*log(x))
  c7 <- sum(weights*y*x^(-exponent))
  c8 <- sum(weights*y*x^(-exponent)*log(x))
  c9 <- sum(weights*y*x^(-exponent)*log(x)^2)
  c11 <- sum(weights*y^2*x^(-2*exponent))
  c12 <- sum(weights*y^2*x^(-2*exponent)*log(x))
  c13 <- sum(weights*y^2*x^(-2*exponent)*log(x)^2)
  # multipliers and their derivatives
  a <- c7/c1
  g <- (c11 - a*c7)/c1
  dg_db <- 2*(a*c8 - c12)/c1
  d2g_db2 <- 4*c13/c1 - 2*(c8^2 + c7*c9)/c1^2
  # return derivative of root
  0.5*c1*(dg_db/g)^2 - 0.5*c1/g*d2g_db2
}


#' Roots for power law + normal residuals - fixed sd_exponent
#'
#' @description
#' Function that returns roots to solved for power-law fitting with normally
#' distributed residuals, assuming a fixed and known value for the exponent
#' for means
#'
#' @inheritParamts power_normal_fit
#' @param exponent unknown power coefficient
#' @return root
#'
power_normal_root_fixb <- function(
    sd_exponent,
    x,
    y,
    exponent = 0,
    weights = rep(1, length(x))
) {
  # coefficients
  c1 <- sum(weights)
  c2 <- sum(weights*log(x))
  c3 <- sum(weights*x^(2*exponent - 2*sd_exponent))
  c4 <- sum(weights*x^(2*exponent - 2*sd_exponent)*log(x))
  c7 <- sum(weights*y*x^(exponent - 2*sd_exponent))
  c8 <- sum(weights*y*x^(exponent - 2*sd_exponent)*log(x))
  c11 <- sum(weights*y^2*x^(-2*sd_exponent))
  c12 <- sum(weights*y^2*x^(-2*sd_exponent)*log(x))
  # multipliers and their derivatives
  a <- c7/c3
  g <- (c11 - a*c7)/c1
  dg_dd <- 2/c1*(-a^2*c4 + 2*a*c8 - c12)
  # roots
  -0.5*c1/g*dg_dd - c2
}


#' Jacobian of roots for power law + normal residuals - fixed sd_exponent
#'
#' @description
#' Function that returns the jacobian of the roots to solve for power-law
#' fitting with normally distributed residuals, assuming a fixed and known
#' value for the exponent for standard deviations
#'
#' @inheritParamts power_normal_root_fixd
#' @return derivative of root as defined by function `power_normal_root_fixd()`
#' @examples
#' # parameters
#' multiplier <- 48
#' exponent <- -0.32
#' sd_multiplier <- 5.2
#' sd_exponent <- 0.2
#'
#' # generate data
#' n <- 251
#' x <- seq(2, 10, l = n)
#' y <- abs(multiplier*x^exponent + sd_multiplier*stats::rnorm(length(x))*x^sd_exponent)
#' w <- runif(n)
#'
#' # check derivative - compare analytical to numerical
#' eps <- 1e-6
#' v0 <- power_normal_root_fixb(sd_exponent, x, y, exponent = exponent, weights = w)
#' v1 <- power_normal_root_fixb(sd_exponent + eps, x, y, exponent = exponent, weights = w)
#' (v1 - v0)/eps
#' power_normal_root_fixb_jacobian(sd_exponent, x, y, exponent = exponent, weights = w)
#'
power_normal_root_fixb_jacobian <- function(
    sd_exponent,
    x,
    y,
    exponent = 0,
    weights = rep(1, length(x))
) {
  # coefficients
  c1 <- sum(weights)
  c3 <- sum(weights*x^(2*exponent - 2*sd_exponent))
  c4 <- sum(weights*x^(2*exponent - 2*sd_exponent)*log(x))
  c5 <- sum(weights*x^(2*exponent - 2*sd_exponent)*log(x)^2)
  c7 <- sum(weights*y*x^(exponent - 2*sd_exponent))
  c8 <- sum(weights*y*x^(exponent - 2*sd_exponent)*log(x))
  c9 <- sum(weights*y*x^(exponent - 2*sd_exponent)*log(x)^2)
  c11 <- sum(weights*y^2*x^(-2*sd_exponent))
  c12 <- sum(weights*y^2*x^(-2*sd_exponent)*log(x))
  c13 <- sum(weights*y^2*x^(-2*sd_exponent)*log(x)^2)
  # multipliers and their derivatives
  a <- c7/c3
  g <- (c11 - a*c7)/c1
  da_dd <- 2*(a*c4 - c8)/c3
  dg_dd <- 2/c1*(-a^2*c4 + 2*a*c8 - c12)
  d2g_dd2 <- 2*c3/c1*da_dd^2 + 4/c1*(a^2*c5 - 2*a*c9 + c13)
  # return derivatives of root
  0.5*c1*(dg_dd/g)^2 - 0.5*c1/g*d2g_dd2
}



#' test <- function(power, sd_exponent, x, y, weights) {
#'   c1 <- sum(weights)
#'   c3 <- sum(weights*x^(2*power - 2*sd_exponent))
#'   c7 <- sum(weights*y*x^(power - 2*sd_exponent))
#'   c11 <- sum(weights*y^2*x^(-2*sd_exponent))
#'   a <- c7/c3
#'   g <- (c11 - a*c7)/c1
#'   c(a, g)
#' }
#' #' @examples
#' #' # parameters
#' #' multiplier <- 50
#' #' exponent <- -0.3
#' #' sd_multiplier <- 15
#' #' sd_exponent <- -0.6
#' #'
#' #' # generate data
#' #' n <- 51
#' #' x <- seq(2, 10, l = n)
#' #' y <- abs(multiplier*x^exponent + sd_multiplier*stats::rnorm(length(x))*x^sd_exponent)
#' #' w <- stats::runif(n)
#' #'
#' #' # test
#' #' eps <- 1e-6
#' #' v0 <- test(exponent, sd_exponent, x, y, w)
#' #' v1 <- test(exponent + eps, sd_exponent, x, y, w)
#' #' v2 <- test(exponent, sd_exponent + eps, x, y, w)
#' #' t((cbind(v1, v2) - v0)/eps)
#' #' test_jac(exponent, sd_exponent, x, y, w)
#' #'
#' test_jac <- function(exponent, sd_exponent, x, y, weights) {
#'   c1 <- sum(weights)
#'   c3 <- sum(weights*x^(2*exponent - 2*sd_exponent))
#'   c4 <- sum(weights*x^(2*exponent - 2*sd_exponent)*log(x))
#'   c5 <- sum(weights*x^(2*exponent - 2*sd_exponent)*log(x)^2)
#'   c7 <- sum(weights*y*x^(exponent - 2*sd_exponent))
#'   c8 <- sum(weights*y*x^(exponent - 2*sd_exponent)*log(x))
#'   c9 <- sum(weights*y*x^(exponent - 2*sd_exponent)*log(x)^2)
#'   c11 <- sum(weights*y^2*x^(-2*sd_exponent))
#'   c12 <- sum(weights*y^2*x^(-2*sd_exponent)*log(x))
#'   a <- c7/c3
#'   g <- (c11 - a*c7)/c1
#'   da_db <- (c8 - 2*a*c4)/c3
#'   da_dd <- 2*(a*c4 - c8)/c3
#'   dg_db <- a*c3/c1*da_dd
#'   dg_dd <- 2/c1*(-a^2*c4 + 2*a*c8 - c12)
#'   matrix(c(da_db, da_dd, dg_db, dg_dd), nrow = 2, byrow = TRUE)
#' }
#' #' @examples
#' #' # parameters
#' #' multiplier <- 50
#' #' exponent <- -0.9
#' #' sd_multiplier <- 15
#' #' sd_exponent <- -0.6
#' #'
#' #' # generate data
#' #' n <- 51
#' #' x <- seq(2, 10, l = n)
#' #' y <- abs(multiplier*x^exponent + sd_multiplier*stats::rnorm(length(x))*x^sd_exponent)
#' #' w <- stats::runif(n)
#' #'
#' #' # test
#' #' eps <- 1e-6
#' #' v0 <- test_jac(exponent, sd_exponent, x, y, w)[2,]
#' #' v1 <- test_jac(exponent + eps, sd_exponent, x, y, w)[2,]
#' #' v2 <- test_jac(exponent, sd_exponent + eps, x, y, w)[2,]
#' #' t(cbind(v1, v2) - v0)/eps
#' #' test_jac2(exponent, sd_exponent, x, y, w)
#' #'
#' test_jac2 <- function(exponent, sd_exponent, x, y, weights) {
#'   c1 <- sum(weights)
#'   c3 <- sum(weights*x^(2*exponent - 2*sd_exponent))
#'   c4 <- sum(weights*x^(2*exponent - 2*sd_exponent)*log(x))
#'   c5 <- sum(weights*x^(2*exponent - 2*sd_exponent)*log(x)^2)
#'   c7 <- sum(weights*y*x^(exponent - 2*sd_exponent))
#'   c8 <- sum(weights*y*x^(exponent - 2*sd_exponent)*log(x))
#'   c9 <- sum(weights*y*x^(exponent - 2*sd_exponent)*log(x)^2)
#'   c11 <- sum(weights*y^2*x^(-2*sd_exponent))
#'   c12 <- sum(weights*y^2*x^(-2*sd_exponent)*log(x))
#'   c13 <- sum(weights*y^2*x^(-2*sd_exponent)*log(x)^2)
#'   a <- c7/c3
#'   g <- (c11 - a*c7)/c1
#'   da_db <- (c8 - 2*a*c4)/c3
#'   da_dd <- 2*(a*c4 - c8)/c3
#'   d2g_db2 <- -2*c3/c1*(da_db)^2 + 2/c1*(2*a^2*c5 - a*c9)
#'   d2g_dbdd <- -2*c3/c1*da_db*da_dd - 4/c1*(a^2*c5 - a*c9)
#'   d2g_dd2 <- 2*c3/c1*da_dd^2 + 4/c1*(a^2*c5 - 2*a*c9 + c13)
#'   matrix(c(d2g_db2, d2g_dbdd, d2g_dbdd, d2g_dd2), nrow = 2, byrow = TRUE)
#' }


#' Covariance matrix of fitting parameters for `power_normal_fit()`
#'
#' @description
#' Estimate the variance-covariance matrix for power law fitting parameters, based
#' on Fisher information matrix, for normally distributed residuals
#'
#' @inheritParams power_normal_fit
#' @param multiplier,exponent power-law multiplier and power coefficient for the
#'   mean
#' @param sd_multiplier,sd_exponent power-law multiplier and power coefficient for
#'   the standard deviation
#' @return 2*2 matrix (for parameters `multiplier` and `exponent`)
#' @examples
#' #' # parameters
#' multiplier <- 50
#' exponent <- -0.3
#' sd_multiplier <- 15
#' sd_exponent <- -0.6
#'
#' # generate data
#' n <- 51
#' x <- seq(2, 10, l = n)
#' y <- abs(multiplier*x^exponent +
#'   sd_multiplier*stats::rnorm(length(x))*x^sd_exponent)
#'
#' # covariance matrix
#' power_normal_covariancematrix(
#'   x, y,
#'   multiplier, exponent,
#'   sd_multiplier, sd_exponent
#' )
#'
power_normal_covariancematrix <- function(
    x,
    y,
    multiplier,
    exponent,
    sd_multiplier,
    sd_exponent,
    weights = rep(1, length(x))
) {
  # coefficients
  c3 <- sum(weights*x^(2*exponent - 2*sd_exponent))
  c4 <- sum(weights*x^(2*exponent - 2*sd_exponent)*log(x))
  c5 <- sum(weights*x^(2*exponent - 2*sd_exponent)*log(x)^2)
  c8 <- sum(weights*y*x^(exponent - 2*sd_exponent)*log(x))
  c9 <- sum(weights*y*x^(exponent - 2*sd_exponent)*log(x)^2)
  # second derivatives of log-transformed probabilities
  dL2_da2 <- -c3/sd_multiplier^2
  dL2_dadb <- (c8 - 2*multiplier*c4)/sd_multiplier^2
  dL2_db2 <- (multiplier*c9 - 2*multiplier^2*c5)/sd_multiplier^2
  # Fisher information matrix - multiplier & power
  fisher <- -matrix(c(dL2_da2, dL2_dadb, dL2_dadb, dL2_db2), nrow = 2)
  # return matrix
  solve(fisher)
}


#' Generate prediction interval for normal fit
#'
#' @description
#' Generate prediction interval with specified level of confidence, based on
#' power law fit with gaussian residuals
#'
#' @inheritParams power_normal_covariancematrix
#' @param x root diameters at which to predict interval
#' @param level confidence level (fraction)
#' @return dataframe with field for diameter (`x`), average power law strength
#'   (`y`) and the lower and upper bound of the prediction interval (`ymin`,
#'   `ymax`)
#' @examples
#' x <- seq(1, 10, l = 10)
#' power_normal_predictioninterval(x, 20, -0.5, 2, -0.4)
#'
power_normal_predictioninterval <- function(
    x,
    multiplier,
    exponent,
    sd_multiplier,
    sd_exponent,
    level = 0.95
) {
  y <- multiplier*x^exponent
  sd <- sd_multiplier*x^sd_exponent
  offset <- stats::qnorm(0.5 - 0.5*level, 0, sd)
  data.frame(x = x, y = y, ymin = y - offset, ymax = y + offset)
}
