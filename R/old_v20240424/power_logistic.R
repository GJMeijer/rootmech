#' @examples
#' y0 <- 20
#' beta <- -0.5
#' s <- 0.1
#' x <- seq(1, 6, l = 51)
#' y <- abs(y0*x^beta*rlogis(length(x), 1, s))
#' w <- runif(length(x), 0.8, 1.2)
#'
#' ft <- power_logistic_fit(x, y, weights = w)
#' xp <- seq(min(x), max(x), l = 101)
#' yp <- ft$multiplier*xp^ft$exponent
#'
#' plot(x, y)
#' lines(xp, yp, col = "red")
#'
power_logistic_fit <- function(x, y, start = NULL, weights = rep(1, length(x))) {
  # initial guess
  if (is.null(start)) {
    start <- power_logistic_initialguess(x, y, weights = weights)
  }
  # solve
  sol <- rootSolve::multiroot(
    function(par) {
      power_logistic_loglikelihood(par, x, y, weights = weights, deriv = 1)
    },
    start,
    jacfunc = function(par) {
      power_logistic_loglikelihood(par, x, y, weights = weights, deriv = 2)
    }
  )
  # loglikelihood
  logL <- power_logistic_loglikelihood(sol$root, x, y, weights = weights)
  # return
  list(
    loglikelihood = logL,
    multiplier = sol$root[1],
    exponent = sol$root[2],
    shape = sol$root[3]
  )
}


power_logistic_initialguess <- function(x, y, weights = rep(1, length(x))) {
  # guess exponent from linear regression on log-data
  ft0 <- stats::lm(log(y) ~ log(x), weights = weights)
  beta <- ft0$coef[2]
  # guess for multiplier and shape, by fitting linearised CDF of (y/x^beta)
  yr <- sort(y/x^beta)
  yp <- seq(0.5/length(x), 1 - 0.5/length(x), l = length(x))
  ypl <- atanh(2*yp - 1)
  ft1 <- stats::lm(ypl ~ yr)
  s <- -0.5/ft1$coef[1]
  y0 <- 0.5/(s*ft1$coef[2])
  # return
  as.numeric(c(y0, beta, s))
}


#' @examples
#' # parameters and data
#' y0 <- 20
#' beta <- -0.5
#' s <- 0.1
#' x <- seq(2, 8, l = 51)
#' y <- y0*x^beta*rlogis(length(x), 1, s)
#' w <- runif(length(x), 0.8, 1.2)
#'
#' # check loglikelihood function
#' par <- c(y0, beta, s)
#' power_logistic_loglikelihood(par, x, y, weights = w, deriv = 0)
#' yp <- y0*x^beta
#' sum(w*stats::dlogis(y/yp, 1, s, log = TRUE))
#'
#' # test first derivative
#' eps <- 1e-6
#' f0 <- power_logistic_loglikelihood(par, x, y, weights = w, deriv = 0)
#' f1 <- power_logistic_loglikelihood(par + c(eps, 0, 0), x, y, weights = w, deriv = 0)
#' f2 <- power_logistic_loglikelihood(par + c(0, eps, 0), x, y, weights = w, deriv = 0)
#' f3 <- power_logistic_loglikelihood(par + c(0, 0, eps), x, y, weights = w, deriv = 0)
#' (c(f1, f2, f3) - f0)/eps
#' power_logistic_loglikelihood(par, x, y, weights = w, deriv = 1)
#'
#' # test second derivative
#' f0 <- power_logistic_loglikelihood(par, x, y, weights = w, deriv = 1)
#' f1 <- power_logistic_loglikelihood(par + c(eps, 0, 0), x, y, weights = w, deriv = 1)
#' f2 <- power_logistic_loglikelihood(par + c(0, eps, 0), x, y, weights = w, deriv = 1)
#' f3 <- power_logistic_loglikelihood(par + c(0, 0, eps), x, y, weights = w, deriv = 1)
#' (cbind(f1, f2, f3) - f0)/eps
#' power_logistic_loglikelihood(par, x, y, weights = w, deriv = 2)
#'
power_logistic_loglikelihood <- function(
    par,
    x,
    y,
    weights = rep(1, length(x)),
    deriv = 0
    ) {
  # split input
  y0 <- par[1]
  beta <- par[2]
  s <- par[3]
  # weighted loglikelihood
  if (deriv == 0) {
    -log(4*s)*sum(weights) -
      2*sum(weights*log(cosh((y/(y0*x^beta) - 1)/(2*s))))
  } else if (deriv == 1) {
    # intermediate parameters
    zeta <- tanh((y/(y0*x^beta) - 1)/(2*s))
    # derivatives of loglikelihood
    dlogL_dy0 <- 1/(s*y0^2)*sum(weights*zeta*y*x^(-beta))
    dlogL_dbeta <- 1/(s*y0)*sum(weights*zeta*y*log(x)*x^(-beta))
    dlogL_ds <- -sum(weights)/s + 1/(s^2)*sum(weights*zeta*(y/(y0*x^beta) - 1))
    # return
    c(dlogL_dy0, dlogL_dbeta, dlogL_ds)
  } else if (deriv == 2) {
    # intermediate parameters
    zeta <- tanh((y/(y0*x^beta) - 1)/(2*s))
    eta <- cosh((y/(y0*x^beta) - 1)/(2*s))^(-2)
    # derivatives of zeta
    dzeta_dy0 <- -eta*y/(2*s*y0^2*x^beta)
    dzeta_dbeta <- -eta*y*log(x)/(2*s*y0*x^beta)
    dzeta_ds <- -eta/(2*s^2)*(y/(y0*x^beta) - 1)
    # derivatives of loglikelihood
    d2logL_dy02 <- -2/(s*y0^3)*sum(weights*zeta*y/x^beta) +
      1/(s*y0^2)*sum(weights*y/x^beta*dzeta_dy0)
    d2logL_dy0dbeta <- -1/(s*y0^2)*sum(weights*zeta*y*log(x)/x^beta) +
      1/(s*y0^2)*sum(weights*y/x^beta*dzeta_dbeta)
    d2logL_dy0ds <- -1/(s^2*y0^2)*sum(weights*zeta*y/x^beta) +
      1/(s*y0^2)*sum(weights*y/x^beta*dzeta_ds)
    d2logL_dbeta2 <- -1/(s*y0)*sum(weights*zeta*y*log(x)^2/x^beta) +
      1/(s*y0)*sum(weights*y*log(x)/x^beta*dzeta_dbeta)
    d2logL_dbetads <- -1/(s^2*y0)*sum(weights*zeta*y*log(x)/x^beta) +
      1/(s*y0)*sum(weights*y*log(x)/x^beta*dzeta_ds)
    #d2logL_dbetads <- -1/(s*y0)*sum(weights*zeta*y*log(x)^2/x^beta) +
    #  1/s^2*sum(weights*(y/(y0*x^beta) - 1)*dzeta_ds)
    d2logL_ds2 <- sum(weights)/s^2 -
      2/s^3*sum(weights*zeta*(y/(y0*x^beta) - 1)) +
      1/s^2*sum(weights*(y/(y0*x^beta) - 1)*dzeta_ds)
    # return
    matrix(
      c(
        d2logL_dy02, d2logL_dy0dbeta, d2logL_dy0ds,
        d2logL_dy0dbeta, d2logL_dbeta2, d2logL_dbetads,
        d2logL_dy0ds, d2logL_dbetads, d2logL_ds2
      ),
      nrow = 3
    )
  }
}


#' @examples
#' # generate some data
#' y0 <- 20
#' beta <- -0.5
#' shape <- 0.1
#' x <- seq(1, 8, l = 51)
#' y <- y0*x^beta*rlogis(length(x), 1, shape)
#' weights <- runif(length(x), 0.8, 1.2)
#'
#' # fit
#' ft <- power_logistic_fit(x, y, weights = weights)
#'
#' # Covariance
#' power_logistic_covariancematrix(
#'   ft$multiplier, ft$exponent, ft$shape,
#'   x, y, weights = weights
#' )
#'
power_logistic_covariancematrix <- function(
    y0,
    beta,
    s,
    x,
    y,
    weights = weights
) {
  # 2nd derivative of loglikelihood
  J <- power_logistic_loglikelihood(
    c(y0, beta, s),
    x, y, weights = weights,
    deriv = 2
  )
  # inverse of fisher information (negative second derivative of log(L))
  solve(-J)
}
