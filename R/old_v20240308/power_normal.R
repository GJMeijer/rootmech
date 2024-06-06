#' Power law fitting with normally distributed residuals
#'
#' @description
#' Fit a power-law curve to a series of `x`, `y` data. The mean of the
#' data follows a power-law in the form:
#'
#'   mu = multiplier*x^power
#'
#' The residuals are assumed normally distributed, with a mean of zero and
#' a standard deviation sd in the form:
#'
#'   sd = sd_multiplier*x^(sd_power)
#'
#' When doing non-linear least squares fitting, `sd_power = 0` (default
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
#' The problem can be written in terms of unknown `power` and `sd_power`
#' only, through assuming that at the maximum likelihood the partial
#' derivatives of the likelihood function are equal to zero.
#'
#' Unknowns are found using a root solving algorithm, either using the
#' bisection method (`method = "uniroot"`) or gradient descent
#' (`method = "newton"`). If (`method = "optim")` unknowns are found by
#' maximising the loglikelihood function using the `stats::optim()`
#' algorithm.
#'
#' @param x,y vectors with x and y observations. Must be the same length, and
#'   contain numerical values only
#' @param weights weighting parameter for each observation
#' @param multiplier fixed value of multiplier for the mean. Will be fitted if
#'   not explicitly defined (`NA`).
#' @param power fixed value of power coefficient for the mean. Will be fitted if
#'   not explicitly defined (`NA`).
#' @param sd_multiplier fixed value of multiplier for the standard deviation.
#'   Will be fitted if not explicitly defined (`NA`).
#' @param sd_power fixed value of multiplier for the standard deviation.
#'   Will be fitted if not explicitly defined (`NA`). Default setting is
#'   `sd_power = 0`, i.e. the standard deviation is set as constant with `x`.
#'   Can also be set to `fixed` in which case the coefficient if assumed the
#'   same value as `power`.
#' @param method `newton` for gradient descent root solving. `uniroot` for
#'   bisection root solving. `optim` for loglikelihood optimisations.
#'   `newton` can only be used when one of the power coefficients (`power`
#'   or `sd_power`) is known.
#' @param interval initial search interval for root bisection method, when
#'   method `uniroot` is used
#' @param extendInt interval extension setting for `uniroot` method
#' @param guess initial guess for power law coefficient(s) when method `newton`
#'   is used. If not defined, a decent guess is made based on simple fitting
#'   (see function `power_normal_initialguess()`)
#' @return a dataframe with the loglikelihood of the fit (`loglikelihood`), the
#'   multiplier and power coefficient for the mean (`multiplier`, `power`) and
#'   the multiplier and power coefficient for the standard deviation
#'   (`sd_multiplier`, `sd_power`)
#' @export
#' @examples
#' beta <- -0.3
#' delta <- -1
#' x <- seq(2, 10, l = 100)
#' y <- abs(50*x^beta + 40*rnorm(length(x))*x^delta)
#'
#' ft <- power_normal_fit(x, y, power = -0.3, sd_power = -1)
#' ft <- power_normal_fit(x, y, sd_power = -1)
#' ft <- power_normal_fit(x, y, sd_power = NULL)
#' ft <- power_normal_fit(x, y, sd_power = "fixed")
#' ft <- power_normal_fit(x, y)
#' ft <- power_normal_fit(x, y, sd_power = NA)
#'
#' xp <- seq(min(x), max(x), l = 101)
#' yp <- ft$multiplier*xp^ft$power
#' range <- 1*ft$sd_multiplier*xp^ft$sd_power
#' plot(x, y, ylim = c(0, max(c(yp + range, y))))
#' lines(xp, yp, col = "red")
#' lines(xp, yp - range, col = "red")
#' lines(xp, yp + range, col = "red")
#'
power_normal_fit <- function(
    x,
    y,
    weights = rep(1, length(x)),
    multiplier = NA,
    power = NA,
    sd_multiplier = NA,
    sd_power = 0,
    method = "newton",
    interval = c(-2, 2),
    extendInt = "yes",
    guess = NULL
) {
  # get fitting case
  case <- power_normal_case(beta = power, delta = sd_power)
  # errors
  if ((method == "uniroot") & (case == 6)) {
    warning("uniroot method cannot be used when both power coefficients are unknown. Method changed to `newton`")
    method <- "newton"
  }
  # no fitting required - all power law coefficients known
  if (case == 1) {
    val <- c(power, sd_power)
  } else if (case == 2) {
    val <- c(power, power)
  } else {
    # fitting required - make initial guess
    if (!is.numeric(guess) & (method %in% c("newton", "optim"))) {
      guess <- power_normal_initialguess(
        x, y, weights = weights,
        case = case,
        multiplier = multiplier, power = power,
        sd_multiplier = sd_multiplier, sd_power = sd_power
      )
    }
    # find solutions for power coefficient(s)
    if (method == "newton") {
      val <- rootSolve::multiroot(
        power_normal_root,
        guess,
        jacfunc = power_normal_root_jacobian,
        jactype = "fullusr",
        x = x,
        y = y,
        case = case,
        alpha = multiplier,
        beta = power,
        gamma = sd_multiplier,
        delta = sd_power,
        weights = weights,
      )$root
    } else if (method == "uniroot") {
      val <- stats::uniroot(
        power_normal_root,
        interval = interval,
        extendInt = extendInt,
        x = x,
        y = y,
        case = case,
        weights = weights,
        alpha = multiplier,
        beta = power,
        gamma = sd_multiplier,
        delta = sd_power
      )$root
    } else if (method == "optim") {
      val <- stats::optim(
        guess,
        power_normal_loglikelihood,
        gr = power_normal_loglikelihood_jacobian,
        method = "BFGS",
        x = x,
        y = y,
        case = case,
        weights = weights,
        alpha = multiplier,
        beta = power,
        gamma = sd_multiplier,
        delta = sd_power,
        control = list(fnscale = -1)
      )$par
    }
  }
  # assign parameters
  t1 <- power_normal_betadelta(val, case = case, beta = power, delta = sd_power)
  C <- power_normal_coef(t1[1], t1[2], x, y, weights = weights)
  t2 <- power_normal_alphagammasq(C, alpha = multiplier, gamma = sd_multiplier)
  # calculate likelihood
  mu <- t2[1]*x^t1[1]
  sigma <- sqrt(t2[2])*x^t1[2]
  logpi <- -log(sigma) - 0.5*log(2*pi) - 0.5*((y - mu)/sigma)^2
  # return dataframe
  data.frame(
    loglikelihood = sum(weights*logpi),
    multiplier = t2[1],
    power = t1[1],
    sd_multiplier = sqrt(t2[2]),
    sd_power = t1[2]
  )
}


#' Detect power law (normally distributed residuals) fitting case
#'
#' @md
#' @description
#' Detect type of power-law fitting (using normally distributed residuals).
#' Can have values
#'
#' * 1: `power` known, `sd_power` known
#' * 2: `power` known, `sd_power` = `power`
#' * 3: `power` known, `sd_power` unknown
#' * 4: `power` unknown, `sd_power` known
#' * 5: `power` unknown, `sd_power` = `power`
#' * 6: `power` unknown, `sd_power` unknown
#'
#' @param beta mean power coefficient (`power`)
#' @param delta standard deviation power coefficient (`sd_power`)
#' @return integer, from 1 to 6, see description
#'
power_normal_case <- function(
    beta = NA,
    delta = NA
) {
  if (is.numeric(beta)) {
    if (is.numeric(delta)) {
      1 # beta and gamma both known
    } else if (is.character(delta)) {
      2 # beta known, delta = beta
    } else {
      3 # beta known, delta unknown
    }
  } else {
    if (is.numeric(delta)) {
      4 # beta unknown, delta known
    } else if (is.character(delta)) {
      5 # beta = delta, unknown
    } else {
      6 # beta unknown, delta unknown
    }
  }
}


power_normal_initialguess <- function(
    x,
    y,
    case = 4,
    multiplier = NA,
    power = NA,
    sd_multiplier = NA,
    sd_power = 0,
    weights = rep(1, length(x)),
    range = c(-2, 2)
) {
  # get mean - power and multiplier
  if (case %in% c(4, 5, 6)) {
    if (!is.numeric(multiplier)) {
      ft1 <- stats::lm(log(y) ~ log(x), weights = weights)
      sd <- sqrt(sum(ft1$residuals^2)/length(x))
      power <- as.numeric(ft1$coefficients[2])
      multiplier <- as.numeric(exp(ft1$coefficients[1] + 0.5*sd^2))
    } else {
      ft1 <- stats::lm(log(y/multiplier) ~ 0 + log(x), weights = weights)
      power <- as.numeric(ft1$coefficients[1])
    }
  } else if (!is.numeric(multiplier)) {
    multiplier <- mean(y/x^power)
  }
  # guess for delta required
  if (case %in% c(3, 6)) {
    res <- y - (multiplier*x^power)
    sd_power <- power_normal_fit(
      x, res, weights = weights,
      multiplier = 0, power = 0,
      sd_multiplier = sd_multiplier, sd_power = sd_power,
      method = "uniroot",
      interval = power + c(-1, 1)
    )$sd_power
  } else {
    sd_power <- NULL
  }
  # return guess
  if (case %in% c(4, 5, 6)) {
    c(power, sd_power)
  } else {
    c(sd_power)
  }
}


#' Intermediate coefficients used in function `power_normal_fit()`
#'
#' @description
#' Defines eight intermediate variables that are used in the root solving
#' procedure used in the function `power_normal_fit()` in order
#' to obtain fitting coefficients.
#'
#' @inheritParams power_normal_fit
#' @param beta power coefficient of the fit
#' @param delta power coefficient of the standard deviation of residuals
#' @return vector with coefficients
#'
power_normal_coef <- function(
    beta,
    delta,
    x,
    y,
    weights = rep(1, length(x))
) {
  c(
    sum(weights),
    sum(weights*log(x)),
    sum(weights*y^2*x^(-2*delta)),
    sum(weights*y^2*x^(-2*delta)*log(x)),
    sum(weights*y*x^(beta - 2*delta)),
    sum(weights*y*x^(beta - 2*delta)*log(x)),
    sum(weights*x^(2*beta - 2*delta)),
    sum(weights*x^(2*beta - 2*delta)*log(x))
  )
}


#' Jacobian of `power_normal_coef()`
#'
#' @description
#' Return the derivative of intermediate coefficients (determined using the
#' function `power_normal_coef()`) with respect to input arguments
#' `beta` and `delta`.
#'
#' @inheritParams power_normal_coef
#' @return a 8*2 matrix (rows = coefficients, columns = derivatives, i.e. beta
#'   and delta)
#' @examples
#' # Test derivative - compare to simple numerical solution
#' beta <- -0.6
#' delta <- 0.2
#' x <- seq(2, 10, l = 25)
#' y <- 10*x^-beta * rweibull(length(x), shape = 4, scale = 1/gamma(1 + 1/4))
#' w <- runif(length(x))
#'
#' eps <- 1e-6
#' C <- power_normal_coef(beta, delta, x, y, weights = w)
#' C1 <- power_normal_coef(beta + eps, delta, x, y, weights = w)
#' C2 <- power_normal_coef(beta, delta + eps, x, y, weights = w)
#' (cbind(C1, C2) - C)/eps
#' power_normal_coef_jacobian(beta, delta, x, y, weights = w)
#'
power_normal_coef_jacobian <- function(
    beta,
    delta,
    x,
    y,
    weights = rep(1, length(x))
) {
  cbind(
    c(
      0,
      0,
      0,
      0,
      sum(weights*y*x^(beta - 2*delta)*log(x)),
      sum(weights*y*x^(beta - 2*delta)*log(x)^2),
      2*sum(weights*x^(2*beta - 2*delta)*log(x)),
      2*sum(weights*x^(2*beta - 2*delta)*log(x)^2)
    ),
    c(
      0,
      0,
      -2*sum(weights*y^2*x^(-2*delta)*log(x)),
      -2*sum(weights*y^2*x^(-2*delta)*log(x)^2),
      -2*sum(weights*y*x^(beta - 2*delta)*log(x)),
      -2*sum(weights*y*x^(beta - 2*delta)*log(x)^2),
      -2*sum(weights*x^(2*beta - 2*delta)*log(x)),
      -2*sum(weights*x^(2*beta - 2*delta)*log(x)^2)
    )
  )
}


#' Combine fitted and fixed power coefficients
#'
#' @description
#' Generate a vector of the two power coefficients (`beta` or `power`, and
#' `delta` or `sd_power`) from values to fit and fixed values
#'
#' @inheritParams power_normal_coef
#' @param val vector with powers to fit fit, in root solving
#' @param case fitting case integer, see function `power_normal_case()`
#' @return 2-value vector with both power coefficients
#'
power_normal_betadelta <- function(val, beta = NA, delta = 0, case = 4) {
  if (case == 2) {
    delta <- beta
  } else if (case == 3) {
    delta <- val
  } else if (case == 4) {
    beta <- val
  } else if (case == 5) {
    beta <- val
    delta <- val
  } else if (case == 6) {
    beta <- val[1]
    delta <- val[2]
  }
  c(beta, delta)
}


#' Calculate power-law multipliers from known power coefficients
#'
#' @description
#' Generate a vector of the two power-law multipliers (`alpha` or `multiplier`,
#' and `gamma` or `sd_multiplier`) from known power coefficients.
#'
#' @param C intermediate coefficients, defined by function `power_normal_coef()`
#' @param alpha fixed value of mean multiplier, if set
#' @param gamma fixed value of standard deviation multiplier, if set
#' @return 2-value vector with both multipliers coefficients. The squared value
#'   of the standard deviation multiplier is returned
#'
power_normal_alphagammasq <- function(C, alpha = NA, gamma = NA) {
  if (!is.numeric(alpha)) {
    alpha <- C[5]/C[7]
  }
  if (!is.numeric(gamma)) {
    gammasq <- (C[3] - 2*alpha*C[5] + alpha^2*C[7])/C[1]
  } else {
    gammasq <- gamma^2
  }
  c(alpha, gammasq)
}


#' Jacobian of function `power_normal_alphagammasq()`
#'
#' @description
#' Returns the derivatives of `alpha` and squared `gamma` with respect to
#' coefficients `C`, as defined by function `power_normal_alphagammasq()`
#'
#' @inheritParams power_normal_alphagammasq
#' @return 2 by 8 jacobian matrix with derivatives of `alpha` and `gamma^2`
#'
power_normal_alphagammasq_jacobian <- function(C, alpha = NA, gamma = NA) {
  if (!is.numeric(alpha)) {
    alpha <- C[5]/C[7]
    dalpha_dC <- c(0, 0, 0, 0, 1/C[7], 0, -C[5]/C[7]^2, 0)
  } else {
    dalpha_dC <- rep(0, 8)
  }
  if (!is.numeric(gamma)) {
    gammasq <- (C[3] - 2*alpha*C[5] + alpha^2*C[7])/C[1]
    dgammasq_dalpha <- c(-2*C[5] + 2*alpha*C[7])/C[1]
    dgammasq_dC <- c(
      -(C[3] - 2*alpha*C[5] + alpha^2*C[7])/C[1]^2,
      0,
      1/C[1],
      0,
      -2*alpha/C[1],
      0,
      1*alpha^2/C[1],
      0
    )
  } else {
    dgammasq_dalpha <- 0
    dgammasq_dC <- rep(0, 8)
  }
  rbind(
    dalpha_dC,
    dgammasq_dC + dgammasq_dalpha %*% dalpha_dC
  )
}


#' Roots to solve in function `power_normal_fit()`
#'
#' @description
#' Returns the root values that need to be solved in order to solve the
#' power-law fitting problem with normally distributed standard deviations
#' (function `power_normal_fit`)
#'
#' @inheritParams power_normal_fit
#' @param val vector with power coefficients that need to be solved for
#' @param case fitting case integer, see function `power_normal_case()`
#' @return a one or two parameter vector (depending on the length of `val`,
#'   i.e. the fitting case) that need to root-solved
#'
power_normal_root <- function(
    val,
    x,
    y,
    case = 4,
    alpha = NA,
    beta = NA,
    gamma = NA,
    delta = 0,
    weights = rep(1, length(x))
) {
  # get beta and delta
  t1 <- power_normal_betadelta(val, case = case, beta = beta, delta = delta)
  # calculate coefficients
  C <- power_normal_coef(t1[1], t1[2], x, y, weights = weights)
  # get alpha and gamma_squared
  t2 <- power_normal_alphagammasq(C, alpha = alpha, gamma = gamma)
  # define roots
  rbeta <- t2[1]*C[6] - t2[1]^2*C[8]
  rdelta <- -C[2]*t2[2] + C[4] - 2*t2[1]*C[6] + t2[1]^2*C[8]
  # return correct root
  if (case == 3) {
    rdelta
  } else if (case == 4) {
    rbeta
  } else if (case == 5) {
    rbeta + rdelta
  } else if (case == 6) {
    c(rbeta, rdelta)
  }
}


#' Jacobian of function `power_normal_root()`
#'
#' @description
#' Returns the derivatives of the function `power_normal_root()` (the roots
#' that need to be solved for) with respect to fitting parameters `val`
#'
#' @inheritParams power_normal_root
#' @return jacobian matrix (square, with size `len(val)`#'
#' @examples
#' alpha <- NULL
#' beta <- -0.3
#' gamma <- NULL
#' delta <- -1
#' x <- seq(2, 10, l = 100)
#' y <- abs(50*x^beta + 20*rnorm(length(x))*x^delta)
#' y <- abs(y)
#' eps <- 1e-6
#'
#' # case 3 - beta known, delta unknown
#' val <- delta
#' r0 <- power_normal_root(val, x, y, case = 3, beta = beta)
#' r1 <- power_normal_root(val + eps, x, y, case = 3, beta = beta)
#' (r1 - r0)/eps
#' power_normal_root_jacobian(val, x, y, case = 3, beta = beta)
#'
#' # case 4 - beta unknown, delta known
#' val <- beta
#' r0 <- power_normal_root(val, x, y, case = 4, delta = delta)
#' r1 <- power_normal_root(val + eps, x, y, case = 4, delta = delta)
#' (r1 - r0)/eps
#' power_normal_root_jacobian(val, x, y, case = 4, delta = delta)
#'
#' # case 5 - beta unknown, delta = beta
#' val <- beta
#' r0 <- power_normal_root(val, x, y, case = 5)
#' r1 <- power_normal_root(val + eps, x, y, case = 5)
#' (r1 - r0)/eps
#' power_normal_root_jacobian(val, x, y, case = 5)
#'
#' # case 6 - both beta and delta unknown
#' val <- c(beta, delta)
#' r0 <- power_normal_root(val, x, y, case = 6)
#' r1 <- power_normal_root(val + c(eps, 0), x, y, case = 6)
#' r2 <- power_normal_root(val + c(0, eps), x, y, case = 6)
#' (cbind(r1, r2) - r0)/eps
#' power_normal_root_jacobian(val, x, y, case = 6)
#'
power_normal_root_jacobian <- function(
    val,
    x,
    y,
    case = 4,
    alpha = NA,
    beta = NA,
    gamma = NA,
    delta = 0,
    weights = rep(1, length(x))
) {
  # get beta and delta
  t1 <- power_normal_betadelta(val, case = case, beta = beta, delta = delta)
  # calculate coefficients
  C <- power_normal_coef(t1[1], t1[2], x, y, weights = weights)
  dC_dt1 <- power_normal_coef_jacobian(t1[1], t1[2], x, y, weights = weights)
  # get alpha and gamma_squared
  t2 <- power_normal_alphagammasq(C, alpha = alpha, gamma = gamma)
  dt2_dC <- power_normal_alphagammasq_jacobian(C, alpha = alpha, gamma = gamma)
  dt2_dt1 <- dt2_dC %*% dC_dt1
  # define roots
  rbeta <- t2[1]*C[6] - t2[1]^2*C[8]
  drbeta_dt2 <- c(C[6] - 2*t2[1]*C[8], 0)
  drbeta_dC <- c(0, 0, 0, 0, 0, t2[1], 0, -t2[1]^2)
  drbeta_dt1 <- (drbeta_dC + drbeta_dt2 %*% dt2_dC) %*% dC_dt1
  rdelta <- -C[2]*t2[2] + C[4] - 2*t2[1]*C[6] + t2[1]^2*C[8]
  drdelta_dt2 <- c(-2*C[6] + 2*t2[1]*C[8], -C[2])
  drdelta_dC <- c(0, -t2[2], 0, 1, 0, -2*t2[1], 0, t2[1]^2)
  drdelta_dt1 <- (drdelta_dC + drdelta_dt2 %*% dt2_dC) %*% dC_dt1
  # return correct root
  if (case == 3) {
    as.matrix(drdelta_dt1[2])
  } else if (case == 4) {
    as.matrix(drbeta_dt1[1])
  } else if (case == 5) {
    as.matrix(sum(drbeta_dt1 + drdelta_dt1))
  } else if (case == 6) {
    rbind(drbeta_dt1, drdelta_dt1)
  }
}


power_normal_loglikelihood <- function(
    val,
    x,
    y,
    case = 4,
    alpha = NA,
    beta = NA,
    gamma = NA,
    delta = 0,
    weights = rep(1, length(x))
) {
  # get beta and delta
  t1 <- power_normal_betadelta(val, case = case, beta = beta, delta = delta)
  # calculate coefficients
  C <- power_normal_coef(t1[1], t1[2], x, y, weights = weights)
  # get alpha and gamma_squared
  t2 <- power_normal_alphagammasq(C, alpha = alpha, gamma = gamma)
  # calculate loglikelihood
  L <- -0.5*log(t2[2])*C[1] - t1[2]*C[2] - 0.5*log(2*pi)*C[1] -
    0.5*C[3]/t2[2] + t2[1]*C[5]/t2[2] - 0.5*t2[1]^2*C[7]/t2[2]
  # return
  L
}


#' @examples
#' alpha <- NULL
#' beta <- -0.3
#' gamma <- NULL
#' delta <- -1
#' x <- seq(2, 10, l = 25)
#' y <- abs(50*x^beta + 20*rnorm(length(x))*x^delta)
#' y <- abs(y)
#' eps <- 1e-6
#'
#' # case 3 - beta known, delta unknown
#' val <- delta
#' r0 <- power_normal_loglikelihood(val, x, y, case = 3, beta = beta)
#' r1 <- power_normal_loglikelihood(val + eps, x, y, case = 3, beta = beta)
#' (r1 - r0)/eps
#' power_normal_loglikelihood_jacobian(val, x, y, case = 3, beta = beta)
#'
#' # case 4 - beta unknown, delta known
#' val <- beta
#' r0 <- power_normal_loglikelihood(val, x, y, case = 4, delta = delta)
#' r1 <- power_normal_loglikelihood(val + eps, x, y, case = 4, delta = delta)
#' (r1 - r0)/eps
#' power_normal_loglikelihood_jacobian(val, x, y, case = 4, delta = delta)
#'
#' # case 5 - beta unknown, delta = beta
#' val <- beta
#' r0 <- power_normal_loglikelihood(val, x, y, case = 5)
#' r1 <- power_normal_loglikelihood(val + eps, x, y, case = 5)
#' (r1 - r0)/eps
#' power_normal_loglikelihood_jacobian(val, x, y, case = 5)
#'
#' # case 6 - both beta and delta unknown
#' val <- c(beta, delta)
#' r0 <- power_normal_loglikelihood(val, x, y, case = 6)
#' r1 <- power_normal_loglikelihood(val + c(eps, 0), x, y, case = 6)
#' r2 <- power_normal_loglikelihood(val + c(0, eps), x, y, case = 6)
#' (c(r1, r2) - r0)/eps
#' power_normal_loglikelihood_jacobian(val, x, y, case = 6)
#'
power_normal_loglikelihood_jacobian <- function(
    val,
    x,
    y,
    case = 4,
    alpha = NA,
    beta = NA,
    gamma = NA,
    delta = 0,
    weights = rep(1, length(x))
) {
  # get beta and delta
  t1 <- power_normal_betadelta(val, case = case, beta = beta, delta = delta)
  # calculate coefficients
  C <- power_normal_coef(t1[1], t1[2], x, y, weights = weights)
  dC_dt1 <- power_normal_coef_jacobian(t1[1], t1[2], x, y, weights = weights)
  # get alpha and gamma_squared
  t2 <- power_normal_alphagammasq(C, alpha = alpha, gamma = gamma)
  dt2_dC <- power_normal_alphagammasq_jacobian(C, alpha = alpha, gamma = gamma)
  # derivative of loglikelihood function
  L <- -0.5*log(t2[2])*C[1] - t1[2]*C[2] - 0.5*log(2*pi)*C[1] -
    0.5*C[3]/t2[2] + t2[1]*C[5]/t2[2] - 0.5*t2[1]^2*C[7]/t2[2]
  dL_dt1 <- c(
    0,
    -C[2]
  )
  dL_dt2 <- c(
    (C[5] - t2[1]*C[7])/t2[2],
    -0.5*C[1]/t2[2] + (0.5*C[3] - t2[1]*C[5] + 0.5*t2[1]^2*C[7])/(t2[2]^2)
  )
  dL_dC <- c(
    -0.5*log(t2[2]) - 0.5*log(2*pi),
    -t1[2],
    -0.5/t2[2],
    0,
    t2[1]/t2[2],
    0,
    -0.5*t2[1]^2/t2[2],
    0
  )
  # derivative with respect to t1
  jac <- dL_dt1 + (dL_dt2 %*% dt2_dC + dL_dC) %*% dC_dt1
  # return
  if (case == 3) {
    jac[2]
  } else if (case == 4) {
    jac[1]
  } else if (case == 5) {
    jac[1] + jac[2]
  } else if (case == 6) {
    as.vector(jac)
  }
}

