#' Fit power-law assuming Weibull distributed residuals
#'
#' @description
#' Fit a power-law regression to a series of (x, y) data using weighted
#' loglikelihood optimalisation.
#'
#' The power law curve describes the mean of the weibull distribution. The shape
#' parameter of the Weibull distribution is constant with x
#'
#' @inheritParams powerlaw_fit
#' @param start (optimal) initial guess for power-law exponent and weibull shape
#'   parameter. If not defined, an educated guess is made using the function
#'   `powerlaw_weibull_initialguess()`.
#' @return a list containing the fields
#'   * `loglikelihood`: the log-likelihood score of the fit
#'   * `multiplier`: fitted power-law multiplier
#'   * `exponent`: fitted power-law exponent
#'   * `shape`: fitted weibull shape parameter
#' @export
#' @examples
#' #' generate some data
#' y0 <- 20
#' beta <- -0.5
#' kappa <- 10
#' lambda <- 1/gamma(1 + 1/kappa)
#' x <- seq(1, 8, l = 101)
#' y <- y0*x^beta*rweibull(length(x), kappa, lambda)
#' weights <- runif(length(x), 0.8, 1.2)
#'
#' # fit
#' ft <- powerlaw_weibull_fit(x, y, weights = weights)
#'
#' # plot
#' xp <- seq(min(x), max(x), l = 101)
#' yp <- ft$multiplier*xp^ft$exponent
#' plot(x, y)
#' lines(xp, yp, col = "red")
#'
powerlaw_weibull_fit <- function(
    x,
    y,
    weights = rep(1, length(x)),
    start = NULL
) {
  # initial guess
  if (is.null(start)) {
    start <- powerlaw_weibull_initialguess(x, y, weights = weights)
  }
  # fit using Newton-Raphson
  sol <- rootSolve::multiroot(
    powerlaw_weibull_root,
    start,
    jacfunc = powerlaw_weibull_root_jacobian,
    x = x,
    y = y,
    weights = weights
  )
  beta <- sol$root[1]
  kappa <- sol$root[2]
  # calculate multiplier
  c1 <- sum(weights)
  c4 <- sum(weights*x^(-beta*kappa)*y^kappa)
  y0 <- gamma(1 + 1/kappa)*(c4/c1)^(1/kappa)
  # loglikelihood
  logL <- powerlaw_weibull_loglikelihood(c(y0, beta, kappa), x, y, weights = weights)
  # return list
  list(
    loglikelihood = logL,
    multiplier = y0,
    exponent = beta,
    shape = kappa
  )
}


#' Calculate power-law weibull loglikelihood
#'
#' @description
#' Calculate the weighted loglikelihood for a power-law fit with weibull
#' residuals. Can also be used to calculate the first or second partial
#' derivative with respect to fitting parameters.
#'
#' @inheritParams powerlaw_weibull_fit
#' @param par vector with fitting parameters (power-law multiplier, power-law
#'   exponent, and weibull shape parameter)
#' @param deriv order of partial derivative requested
#' @return loglikelihood, or its partial derivatives to order `deriv`
#' @keywords internal
#'
powerlaw_weibull_loglikelihood <- function(
    par,
    x,
    y,
    weights = rep(1, length(x)),
    deriv = 0
) {
  # unpack input parameters
  y0 <- par[1]
  beta <- par[2]
  kappa <- par[3]
  # coefficients
  c1 <- sum(weights)
  c2 <- sum(weights*log(x))
  c3 <- sum(weights*log(y))
  c4 <- sum(weights*x^(-beta*kappa)*y^kappa)
  c5 <- sum(weights*x^(-beta*kappa)*y^kappa*log(x))
  c6 <- sum(weights*x^(-beta*kappa)*y^kappa*log(y))
  c7 <- sum(weights*x^(-beta*kappa)*y^kappa*log(x)^2)
  c8 <- sum(weights*x^(-beta*kappa)*y^kappa*log(x)*log(y))
  c9 <- sum(weights*x^(-beta*kappa)*y^kappa*log(y)^2)
  # gamma functions
  g <- gamma(1 + 1/kappa)
  p <- digamma(1 + 1/kappa)
  q <- psigamma(1 + 1/kappa, deriv = 1)
  # loglikelihood
  if (deriv == 0) {
    (log(kappa) - kappa*log(y0) + kappa*log(g))*c1 -
      beta*kappa*c2 +
      (kappa - 1)*c3 -
      (g/y0)^kappa*c4
  } else if (deriv == 1) {
    c(
      -c1*kappa/y0 + c4*kappa*g^kappa*y0^(-kappa - 1),
      -kappa*c2 + kappa*(g/y0)^kappa*c5,
      c1*(1/kappa + log(g/y0) - p/kappa) - beta*c2 + c3 -
        (g/y0)^kappa*(c4*(log(g/y0) - p/kappa) - beta*c5 + c6)
    )
  } else if (deriv == 2) {
    d2logL_dy02 <- c1*kappa/y0^2 - c4*kappa*(kappa + 1)*g^kappa*y0^(-kappa - 2)
    d2logL_dy0dbeta <- -c5*kappa^2*g^kappa*y0^(-kappa - 1)
    d2logL_dy0dkappa <- -c1/y0 + g^kappa*y0^(-kappa - 1)*
      (c4*(1 + kappa*log(g/y0) - p) + kappa*(c6 - beta*c5))
    d2logL_dbeta2 <- -kappa^2*(g/y0)^kappa*c7
    d2logL_dbetadkappa <- -c2 + (g/y0)^kappa*
      (c5*(1 + kappa*log(g/y0) - p) + kappa*(c8 - beta*c7))
    d2logL_dkappa2 <- c1/kappa^2*(q/kappa - 1) - (g/y0)^kappa*(
      2*(log(g/y0) - p/kappa)*(c6 - beta*c5) +
        (log(g/y0) - p/kappa)^2*c4 +
        (c4*q/kappa^3 + beta^2*c7 - 2*beta*c8 + c9)
    )
    matrix(
      c(
        d2logL_dy02, d2logL_dy0dbeta, d2logL_dy0dkappa,
        d2logL_dy0dbeta, d2logL_dbeta2, d2logL_dbetadkappa,
        d2logL_dy0dkappa, d2logL_dbetadkappa, d2logL_dkappa2
      ),
      nrow = 3
    )
  }
}


#' Generate initial guess for `powerlaw_weibull_fit()`
#'
#' @description
#' Generate an intiial guess for power-law exponent and weibull shape parameter,
#' to be used in function `powerlaw_weibull_fit()`.
#'
#' The power-law exponent is estimated using linear regression on
#' log-transformed x and y data. The weibull shape parameter is subsequently
#' estimated from weibull fitting of the scaled data (y/x^beta)
#'
#' @inheritParams powerlaw_weibull_fit
#' @return estimate for power-law exponent and weibull shape parameter
#' @keywords internal
#'
powerlaw_weibull_initialguess <- function(x, y, weights = rep(1, length(x))) {
  # guess exponent from linear regression on log-data
  ft0 <- stats::lm(log(y) ~ log(x), weights = weights)
  beta <- as.numeric(ft0$coef[2])
  # linear regression on linearised cumulative weibull density
  kappa <- weibull_fit(y/x^beta, weights = weights)$shape
  # return beta and kappa
  c(beta, kappa)
}


#' Root to solve for power-law Weibull fitting
#'
#' @description
#' Root equation to solve for power-law Weibull fitting
#'
#' @inheritParams powerlaw_weibull_fit
#' @param par fitting parameter (power-law exponent, and Weibull shape
#'   parameter)
#' @return two-parameter vector with roots
#' @keywords internal
#'
powerlaw_weibull_root <- function(
    par,
    x,
    y,
    weights = rep(1, length(x))
) {
  # unpack input parameters
  beta <- par[1]
  kappa <- par[2]
  # coefficients
  c1 <- sum(weights)
  c2 <- sum(weights*log(x))
  c3 <- sum(weights*log(y))
  c4 <- sum(weights*x^(-beta*kappa)*y^kappa)
  c5 <- sum(weights*x^(-beta*kappa)*y^kappa*log(x))
  c6 <- sum(weights*x^(-beta*kappa)*y^kappa*log(y))
  # roots
  dlogL_dbeta <- kappa*c1*c5/c4 - kappa*c2
  dlogL_dkappa <- (1/kappa + (beta*c5 - c6)/c4)*c1 - beta*c2 + c3
  # return
  c(dlogL_dbeta, dlogL_dkappa)
}


#' Jacobian of root to solve for power-law Weibull fitting
#'
#' @description
#' Jacobian of root equation to solve for power-law Weibull fitting.
#'
#' @inheritParams powerlaw_weibull_root
#' @return derivative of function `powerlaw_weibull_root()` with respect to input
#'   argument `par`
#' @keywords internal
#'
powerlaw_weibull_root_jacobian <- function(
    par,
    x,
    y,
    weights = rep(1, length(x))
) {
  # unpack input parameters
  beta <- par[1]
  kappa <- par[2]
  # coefficients
  c1 <- sum(weights)
  c2 <- sum(weights*log(x))
  c3 <- sum(weights*log(y))
  c4 <- sum(weights*x^(-beta*kappa)*y^kappa)
  c5 <- sum(weights*x^(-beta*kappa)*y^kappa*log(x))
  c6 <- sum(weights*x^(-beta*kappa)*y^kappa*log(y))
  c7 <- sum(weights*x^(-beta*kappa)*y^kappa*log(x)^2)
  c8 <- sum(weights*x^(-beta*kappa)*y^kappa*log(x)*log(y))
  c9 <- sum(weights*x^(-beta*kappa)*y^kappa*log(y)^2)
  # derivatives
  d2logL_dbeta2 <- c1*kappa^2/c4*(c5^2/c4 - c7)
  d2logL_dbetadkappa <- c1/c4*(c5 - beta*kappa*c7 + kappa*c8 + kappa*c5/c4*(beta*c5 - c6)) - c2
  d2logL_dkappa2 <- (-1/kappa^2 - (beta^2*c7 - 2*beta*c8 + c9)/c4 + (beta*c5 - c6)^2/c4^2)*c1
  # return matrix
  matrix(
    c(d2logL_dbeta2, d2logL_dbetadkappa, d2logL_dbetadkappa, d2logL_dkappa2),
    nrow = 2
  )
}

