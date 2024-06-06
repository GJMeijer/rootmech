#' @examples
#' #' generate some data
#' y0 <- 20
#' beta <- -0.5
#' kappa <- 10
#' lambda <- 1/gamma(1 + 1/kappa)
#' x <- seq(1, 8, l = 100001)
#' y <- y0*x^beta*rweibull(length(x), kappa, lambda)
#' weights <- runif(length(x), 0.8, 1.2)
#'
#' # fit
#' ft <- power_weibull_fit(x, y, weights = weights)
#'
#' # plot
#' xp <- seq(min(x), max(x), l = 101)
#' yp <- ft$multiplier*xp^ft$exponent
#' plot(x, y)
#' lines(xp, yp, col = "red")
#'
power_weibull_fit <- function(
    x,
    y,
    start = NULL,
    weights = rep(1, length(x))
) {
  # initial guess
  if (is.null(start)) {
    ft0 <- stats::lm(log(y) ~ log(x), weights = weights)
    beta0 <- as.numeric(ft0$coefficients[2])
    kappa0 <- weibull_fit(y/x^beta0)$shape
    start <- c(beta0, kappa0)
  }
  # fit using Newton-Raphson
  sol <- rootSolve::multiroot(
    power_weibull_root,
    start,
    jacfunc = power_weibull_root_jacobian,
    x = x,
    y = y,
    weights = weights
  )
  beta <- sol$root[1]
  kappa <- sol$root[2]
  # calculate multiplier
  c1 <- sum(weights)
  c4 <- sum(weights*x^(-beta*kappa)*y^kappa)
  y0 <- gamma(1 + 1/kappa)*(kappa*c4/((kappa - 1)*c1))^(1/kappa)
  # loglikelihood
  logL <- power_weibull_loglikelihood(sol$root, x, y, weights = weights)
  power_weibull_loglikelihood(c(y0, beta, kappa), x, y, weights)
  # return list
  list(
    loglikelihood = logL,
    multiplier = y0,
    exponent = beta,
    shape = kappa
  )
}


power_weibull_loglikelihood <- function(
    par,
    x,
    y,
    weights = rep(1, length(x))
) {
  # two parameter -> calculate multiplier from equations
  if (length(par) == 2) {
    # unpack input parameters
    beta <- par[1]
    kappa <- par[2]
    # coefficients
    c1 <- sum(weights)
    c2 <- sum(weights*log(x))
    c3 <- sum(weights*log(y))
    c4 <- sum(weights*x^(-beta*kappa)*y^kappa)
    # loglikelihood
    (log(kappa) + log(gamma(1 + 1/kappa)))*c1 +
      (kappa - 1)/kappa*(log(kappa - 1) + log(c1) - log(kappa) - log(c4) - 1)*c1 +
      (kappa - 1)*c3 - beta*(kappa - 1)*c2
  # three parameters
  } else if (length(par) == 3) {
    # unpack input parameters
    y0 <- par[1]
    beta <- par[2]
    kappa <- par[3]
    # coefficients
    c1 <- sum(weights)
    c2 <- sum(weights*log(x))
    c3 <- sum(weights*log(y))
    c4 <- sum(weights*x^(-beta*kappa)*y^kappa)
    # loglikelihood
    c1*log(kappa) + c1*kappa*log(gamma(1 + 1/kappa)) + c3*(kappa - 1) -
      c1*(kappa - 1)*log(y0) - c2*beta*(kappa - 1) -
      c4*(gamma(1 + 1/kappa)/y0)^kappa
  }
}


#' @examples
#' # input parameters
#' y0 <- 20
#' beta <- -0.5
#' kappa <- 10
#' lambda <- 1/gamma(1 + 1/kappa)
#' x <- seq(1, 8, l = 51)
#' y <- y0*x^beta*rweibull(length(x), kappa, lambda)
#' weights <- runif(length(x), 0.8, 1.2)
#'
#' # first derivative
#' par <- c(y0, beta, kappa)
#' eps <- 1e-6
#' L0 <- power_weibull_loglikelihood(par, x, y, weights = weights)
#' L1 <- power_weibull_loglikelihood(par + c(eps, 0, 0), x, y, weights = weights)
#' L2 <- power_weibull_loglikelihood(par + c(0, eps, 0), x, y, weights = weights)
#' L3 <- power_weibull_loglikelihood(par + c(0, 0, eps), x, y, weights = weights)
#' (c(L1, L2, L3) - L0)/eps
#' power_weibull_loglikelihood_jacobian(par, x, y, weights = weights, deriv = 1)
#'
#' # 2nd derivative
#' J0 <- power_weibull_loglikelihood_jacobian(par, x, y, weights = weights, deriv = 1)
#' J1 <- power_weibull_loglikelihood_jacobian(par + c(eps, 0, 0), x, y, weights = weights, deriv = 1)
#' J2 <- power_weibull_loglikelihood_jacobian(par + c(0, eps, 0), x, y, weights = weights, deriv = 1)
#' J3 <- power_weibull_loglikelihood_jacobian(par + c(0, 0, eps), x, y, weights = weights, deriv = 1)
#' (cbind(J1, J2, J3) - J0)/eps
#' power_weibull_loglikelihood_jacobian(par, x, y, weights = weights, deriv = 2)
#'
power_weibull_loglikelihood_jacobian <- function(
    par,
    x,
    y,
    weights = rep(1, length(x)),
    deriv = 1
) {
  if (length(par) == 3) {
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
    # first derivative
    if (deriv == 1) {
      dlogL_dy0 <- -c1*(kappa - 1)/y0 + c4*kappa*g^kappa*y0^(-kappa - 1)
      dlogL_dbeta <- -c2*(kappa - 1) + kappa*c5*(g/y0)^kappa
      dlogL_dkappa <- c1/kappa*(1 - p) + c1*log(g/y0) + c3 - beta*c2 -
        (g/y0)^kappa*(c4*log(g/y0) - c4*p/kappa - beta*c5 + c6)
      c(dlogL_dy0, dlogL_dbeta, dlogL_dkappa)
    # second derivative
    } else if (deriv == 2) {
      d2logL_dy02 <- c1*(kappa - 1)/y0^2 - kappa*(kappa + 1)*c4*g^kappa*y0^(-kappa - 2)
      d2logL_dy0dbeta <- -c5*kappa^2*g^kappa*y0^(-kappa - 1)
      d2logL_dy0dkappa <- -c1/y0 + g^kappa*y0^(-kappa - 1)*
        (c4*(1 - p) + c4*kappa*log(g/y0) - kappa*(beta*c5 - c6))
      d2logL_dbeta2 <- -kappa^2*(g/y0)^kappa*c7
      d2logL_dbetadkappa <- -c2 + (g/y0)^kappa*
        (c5*(1 - p) + kappa*c5*log(g/y0) - kappa*(beta*c7 - c8))
      d2logL_dkappa2 <- -c1/kappa^2 + c1/kappa^3*q -
        (g/y0)^kappa*(log(g/y0) - p/kappa)*(c4*log(g/y0) - c4*p/kappa - beta*c5 + c6) -
        (g/y0)^kappa*((beta*c5 - c6)*(p/kappa - log(g/y0)) + c4*q/kappa^3 +
                        beta^2*c7 - 2*beta*c8 + c9)
      matrix(c(
        d2logL_dy02, d2logL_dy0dbeta, d2logL_dy0dkappa,
        d2logL_dy0dbeta, d2logL_dbeta2, d2logL_dbetadkappa,
        d2logL_dy0dkappa, d2logL_dbetadkappa, d2logL_dkappa2),
        nrow = 3
      )
    }
  }
}


#' @examples
#' # test root - compare analytical function to numerical approximation
#' y0 <- 20
#' beta <- -0.5
#' kappa <- 10
#' x <- seq(1, 8, l = 51)
#' y <- y0*x^beta*rweibull(length(x), kappa, lambda)
#' weights <- runif(length(x), 0.8, 1.2)
#'
#' eps <- 1e-6
#' f0 <- power_weibull_loglikelihood(c(beta, kappa), x, y, weights = weights)
#' f1 <- power_weibull_loglikelihood(c(beta + eps, kappa), x, y, weights = weights)
#' f2 <- power_weibull_loglikelihood(c(beta, kappa + eps), x, y, weights = weights)
#' (c(f1, f2) - f0)/eps
#' power_weibull_root(c(beta, kappa), x, y, weights = weights)
#'
power_weibull_root <- function(
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
  dlogL_dbeta <- (kappa - 1)*(c1*c5/c4 - c2)
  dlogL_dkappa <- c1/kappa + c1*(kappa - 1)/(c4*kappa)*(beta*c5 - c6) +
    c1/(kappa^2)*(log(kappa - 1) + log(c1) - log(kappa) - log(c4) - digamma(1 + 1/kappa)) +
    c3 - beta*c2
  # return
  c(dlogL_dbeta, dlogL_dkappa)
}


#' @examples
#' # test jacobian - compare analytical function to numerical approximation
#' y0 <- 20
#' beta <- -0.5
#' kappa <- 10
#' x <- seq(1, 8, l = 51)
#' y <- y0*x^beta*rweibull(length(x), kappa, lambda)
#' weights <- runif(length(x), 0.8, 1.2)
#'
#' j0 <- power_weibull_root(par, x, y, weights = weights)
#' j1 <- power_weibull_root(par + c(eps, 0), x, y, weights = weights)
#' j2 <- power_weibull_root(par + c(0, eps), x, y, weights = weights)
#' (cbind(j1, j2) - j0)/eps
#' power_weibull_root_jacobian(par, x, y, weights = weights)
#'
power_weibull_root_jacobian <- function(
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
  d2logL_dbeta2 <- kappa*(kappa - 1)*c1/c4*(c5^2/c4 - c7)
  d2logL_dbetadkappa <- c1*c5/c4 - c2 + beta*(kappa - 1)*c1/c4*(c5^2/c4 - c7) +
    (kappa - 1)*c1/c4*(c8 - c5*c6/c4)
  d2logL_dkappa2 <- c1*(kappa - 1)/(c4*kappa)*((beta*c5 - c6)^2/c4 - beta^2*c7 + 2*beta*c8 - c9) +
    2*c1/kappa^2*((2 - kappa)/(2*kappa - 2) + (beta*c5 - c6)/c4) +
    -2*c1/kappa^3*(log(kappa - 1) + log(c1) - log(kappa) - log(c4) - digamma(1 + 1/kappa) + 0.5) +
    c1/kappa^4*psigamma(1 + 1/kappa, deriv = 1)
  # return matrix
  matrix(
    c(d2logL_dbeta2, d2logL_dbetadkappa, d2logL_dbetadkappa, d2logL_dkappa2),
    nrow = 2
  )
}


#' @examples
#' # generate some data
#' y0 <- 20
#' beta <- -0.5
#' kappa <- 4
#' lambda <- 1/gamma(1 + 1/kappa)
#' x <- seq(1, 8, l = 51)
#' y <- y0*x^beta*rweibull(length(x), kappa, lambda)
#' weights <- runif(length(x), 0.8, 1.2)
#'
#' # fit
#' ft <- power_weibull_fit(x, y, weights = weights)
#'
#' # Covariance
#' power_weibull_covariance(
#'   ft$multiplier, ft$exponent, ft$shape,
#'   x, y, weights = weights
#' )
#'
power_weibull_covariance <- function(
    multiplier,
    exponent,
    shape,
    x,
    y,
    weights = weights
) {
  # 2nd derivative of loglikelihood
  J <- power_weibull_loglikelihood_jacobian(
    c(multiplier, exponent, shape),
    x, y, weights = weights,
    deriv = 2
  )
  # inverse of fisher information (negative second derivative of log(L))
  solve(-J)
}
