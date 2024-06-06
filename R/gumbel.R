#' Gumbel likelihood fitting
#'
#' @description
#' Determine best fitting Gumbel distribution using weighted loglikelihood
#' fitting.
#'
#' @md
#' @inheritParams weibull_fit
#' @param scale_min minimum value of scale parameter. Set to a small value,
#'   since by definition scale > 0
#' @return a list containing the fields
#'   * `loglikelihood`: the log-likelihood score of the fit
#'   * `location`: fitted gumbel location parameter
#'   * `scale`: fitted gumbel scale parameter
#' @export
#' @examples
#' mu <- 20
#' theta <- 3
#' x <- rgumbel(20, mu, theta)
#' gumbel_fit(x)
#'
gumbel_fit <- function(
    x,
    weights = rep(1, length(x)),
    method = "bisection",
    range = c(-1, 1),
    start = NULL,
    scale_min = 1e-1
) {
  # initial guess
  if (is.null(start)) {
    start <- gumbel_initialguess(x, weights = weights)
  }
  # fit scale parameter
  if (method == "newton") {
    theta <- rootSolve::multiroot(
      gumbel_root,
      start,
      jacfunc = gumbel_root_jacobian,
      x = x,
      weights = weights
    )$root
  } else if (method == "bisection") {
    theta <- stats::uniroot(
      gumbel_root,
      pmax(scale_min, start + range),
      extendInt = "downX",
      x = x,
      weights = weights
    )$root
  } else {
    stop("`method` not recognised")
  }
  # calculate location parameter
  c1 <- sum(weights)
  c3 <- sum(weights*exp(-x/theta))
  mu <- theta*log(c1/c3)
  # calculate loglikelihood
  logL <- gumbel_loglikelihood(c(mu, theta), x, weights = weights)
  # return list
  list(
    loglikelihood = logL,
    location = mu,
    scale = theta
  )
}

#' Initial guess for Gumbel distribution fitting
#'
#' @description
#' Return an initial guess for the Gumbel scale parameter, to be used as a
#' starting point in function `gumbel_fit()`.
#'
#' The estimate is made using linear regression on the linearised cumulative
#' density function, i.e by fitting
#'
#'   log(-log(1 - P)) = mu/scale - y/scale
#'
#' where P is the cumulative density for each value x, and mu the
#' location parameter of the Gumbel distribution
#'
#' @export
#' @return initial guess for Gumbel shape parameter
#' @keywords internal
#'
gumbel_initialguess <- function(x, weights = rep(1, length(x))) {
  if (FALSE) {
    # guess based on linearising cumulative density
    n <- length(x)
    C <- seq(0.5/n, 1 - 0.5/n, l = n)
    yf <- log(-log(C))
    xf <- sort(x)
    ft1 <- stats::lm(yf ~ xf)
    as.numeric(-1/ft1$coef[2])
  } else {
    # guess from 1st method of moments - lognormal fit
    muL <- sum(weights*log(x))/sum(weights)
    sdL <- sqrt(sum(weights*(log(x) - muL)^2)/sum(weights))
    mean <- exp(muL + 0.5*sdL^2)
    var <- (exp(sdL^2) - 1)*exp(2*muL + sdL^2)
    sqrt(6*var)/pi
  }
}


#' Calculate Gumbel loglikelihood
#'
#' @description
#' Calculate the weighted loglikelihood for a gumbel fit.
#' Can also be used to calculate the first or second partial
#' derivative with respect to fitting parameters.
#'
#' @inheritParams gumbel_fit
#' @param par vector with fitting parameters (location and scale parameter)
#' @param deriv order of partial derivative requested
#' @return loglikelihood, or its partial derivatives to order `deriv`
#' @keywords internal
#'
gumbel_loglikelihood <- function(
    par,
    x,
    weights = rep(1, length(x)),
    deriv = 0
) {
  # unpack input
  mu <- par[1]
  theta <- par[2]
  # coefficients
  c1 <- sum(weights)
  c2 <- sum(weights*x)
  c3 <- sum(weights*exp(-x/theta))
  c4 <- sum(weights*x*exp(-x/theta))
  c5 <- sum(weights*x^2*exp(-x/theta))
  # loglikelihood
  if (deriv == 0) {
    c1*(mu/theta - log(theta)) - c2/theta - c3*exp(mu/theta)
  } else if (deriv == 1) {
    dlogL_dmu <- c1/theta - c3/theta*(exp(mu/theta))
    dlogL_dtheta <- -c1*(mu/theta^2 + 1/theta) +
      c2/theta^2 +
      (c3*mu - c4)/theta^2*exp(mu/theta)
    c(dlogL_dmu, dlogL_dtheta)
  } else if (deriv == 2) {
    d2logL_dmu2 <- -c3/theta^2*exp(mu/theta)
    d2logL_dmudtheta <- 1/theta^2*(c3*(1 + mu/theta) - c4/theta)*exp(mu/theta) -
      c1/theta^2
    d2logL_dtheta2 <- c1/theta^2 +
      2/theta^3*(c1*mu - c2) -
      1/theta^3*(c3*mu - c4)*exp(mu/theta)*(2 + mu/theta) +
      1/theta^4*(c4*mu - c5)*exp(mu/theta)
    matrix(
      c(d2logL_dmu2, d2logL_dmudtheta, d2logL_dmudtheta, d2logL_dtheta2),
      nrow = 2
    )
  }
}


#' Root to solve for Gumbel fitting
#'
#' @description
#' Root equation to solve for Gumbel fitting
#'
#' @inheritParams gumbel_fit
#' @param theta shape parameter
#' @return root
#' @keywords internal
#'
gumbel_root <- function(theta, x, weights = rep(1, length(x))) {
  # coefficients
  c1 <- sum(weights)
  c2 <- sum(weights*x)
  c3 <- sum(weights*exp(-x/theta))
  c4 <- sum(weights*x*exp(-x/theta))
  # root
  1/theta^2*(c2 - c1*c4/c3) - c1/theta
}


#' Jacobian of root to solve for Gumbel fitting
#'
#' @description
#' Jacobian of root equation to solve for Gumbel fitting.
#'
#' @inheritParams gumbel_root
#' @return derivative of function `gumbel_root()` with respect to input
#'   argument `theta`
#' @keywords internal
#'
gumbel_root_jacobian <- function(theta, x, weights = rep(1, length(x))) {
  # coefficients
  c1 <- sum(weights)
  c2 <- sum(weights*x)
  c3 <- sum(weights*exp(-x/theta))
  c4 <- sum(weights*x*exp(-x/theta))
  c5 <- sum(weights*x^2*exp(-x/theta))
  # root
  c1/(c3*theta^4)*(c4^2/c3 - c5) - 2/theta^3*(c2 - c1*c4/c3) + c1/theta^2
}


#' Gumbel probability distribution density
#'
#' @description
#' Calculate Gumbel probability density
#'
#' @param x vector with x-values
#' @param location Gumbel location parameter
#' @param scale Gumbel scale parameter
#' @param weights vector with weights for each x-value
#' @param log if `TRUE`, return log-transformed probability densities
#' @return vector with probability density for each value in x
#' @export
#' @examples
#' dgumbel(seq(1, 3, l = 5), location = 2, scale = 2)
#'
dgumbel <- function(
  x,
  location = 0,
  scale = 1,
  weights = rep(1, length(x)),
  log = FALSE
) {
  z <- (x - location)/scale
  if (log == FALSE) {
    p <- 1/scale*exp(-z - exp(-z))
    p^weights
  } else {
    logp <- -log(scale) - z - exp(-z)
    weights*logp
  }
}


#' Generate random Gumbel values
#'
#' @description
#' Draw random values from a Gumbel probability distribution
#'
#' @param n number of samples to draw
#' @param location Gumbel location parameter
#' @param scale Gumbel scale parameter
#' @return vector with samples
#' @export
#' @examples
#' hist(rgumbel(1000, location = 20, scale = 2))
#'
rgumbel <- function(
  n,
  location = 0,
  scale = 1
) {
  y <- stats::runif(n, 0, 1)
  location - scale*log(-log(y))
}

