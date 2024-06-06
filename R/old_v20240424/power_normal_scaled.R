#' @examples
#' # parameters
#' y0 <- 20
#' beta <- -0.5
#' sigma <- 0.2
#' x <- seq(1, 8, l = 50001)
#' y <- y0*x^beta*stats::rnorm(length(x), 1, sigma)
#' weights <- stats::runif(length(x), 0.9, 1.1)
#'
#' # fit
#' ft <- power_normal_scaled_fit(x, y, weights = weights)
#' xp <- seq(min(x), max(x), l = 251)
#' yp <- ft$multiplier*xp^ft$exponent
#'
#' # plot data and best fit
#' plot(x, y)
#' lines(xp, yp, col = "red")
#'
power_normal_scaled_fit <- function(
    x,
    y,
    weights = rep(1, length(x)),
    method = "newton",
    range = c(-1, 1),
    start = NULL
) {
  # initial guess
  if (is.null(start)) {
    ft0 <- stats::lm(log(y) ~ log(x), weights = weights)
    start <- as.numeric(ft0$coef[2])
  }
  # solve
  if (method == "newton") {
    beta <- rootSolve::multiroot(
      power_normal_scaled_root,
      start,
      jacfunc = power_normal_scaled_root_jacobian,
      x = x,
      y = y,
      weights = weights
    )$root
  } else if (method == "bisection") {
    beta <- stats::uniroot(
      power_normal_scaled_root,
      start + range,
      extendInt = "downX",
      x = x,
      y = y,
      weights = weights
    )$root
  } else {
    stop("`method` not recognised")
  }
  # multiplier <y0> and standard deviation <sigma>
  c1 <- sum(weights)
  c3 <- sum(weights*x^(-beta)*y)
  c6 <- sum(weights*x^(-2*beta)*y^2)
  y0 <- c6/c3
  sigma <- sqrt(1 - c3^2/(c1*c6))
  # loglikelihood
  logL <- power_normal_scaled_loglikelihood(
    c(y0, beta, sigma),
    x, y, weights
  )
  # return list
  list(
    loglikelihood = logL,
    multiplier = y0,
    exponent = beta,
    sd = sigma
  )
}


#' @examples
#' # define parameters
#' y0 <- 20
#' beta <- -0.5
#' sigma <- 0.2
#' x <- seq(1, 8, l = 51)
#' y <- y0*x^beta*stats::rnorm(length(x), 1, sigma)
#' w <- stats::runif(length(x), 0.8, 1.2)
#' par <- c(y0, beta, sigma)
#'
#' # check loglikelihood
#' yp <- y0*x^beta
#' sum(w*stats::dnorm(y/yp, 1, sigma, log = TRUE))
#' power_normal_scaled_loglikelihood(par, x, y, weights = w)
#'
#' # test first derivative
#' eps <- 1e-6
#' f0 <- power_normal_scaled_loglikelihood(par, x, y, weights = w, deriv = 0)
#' f1 <- power_normal_scaled_loglikelihood(par + c(eps, 0, 0), x, y, weights = w, deriv = 0)
#' f2 <- power_normal_scaled_loglikelihood(par + c(0, eps, 0), x, y, weights = w, deriv = 0)
#' f3 <- power_normal_scaled_loglikelihood(par + c(0, 0, eps), x, y, weights = w, deriv = 0)
#' (c(f1, f2, f3) - f0)/eps
#' power_normal_scaled_loglikelihood(par, x, y, weights = w, deriv = 1)
#'
#' # test second derivative
#' f0 <- power_normal_scaled_loglikelihood(par, x, y, weights = w, deriv = 1)
#' f1 <- power_normal_scaled_loglikelihood(par + c(eps, 0, 0), x, y, weights = w, deriv = 1)
#' f2 <- power_normal_scaled_loglikelihood(par + c(0, eps, 0), x, y, weights = w, deriv = 1)
#' f3 <- power_normal_scaled_loglikelihood(par + c(0, 0, eps), x, y, weights = w, deriv = 1)
#' (cbind(f1, f2, f3) - f0)/eps
#' power_normal_scaled_loglikelihood(par, x, y, weights = w, deriv = 2)
#'
power_normal_scaled_loglikelihood <- function(
    par,
    x,
    y,
    weights = rep(1, length(x)),
    deriv = 0
) {
  # split input
  y0 <- par[1]
  beta <- par[2]
  sigma <- par[3]
  # coefficients
  c1 <- sum(weights)
  c2 <- sum(weights*log(x))
  c3 <- sum(weights*x^(-beta)*y)
  c4 <- sum(weights*x^(-beta)*y*log(x))
  c5 <- sum(weights*x^(-beta)*y*log(x)^2)
  c6 <- sum(weights*x^(-2*beta)*y^2)
  c7 <- sum(weights*x^(-2*beta)*y^2*log(x))
  c8 <- sum(weights*x^(-2*beta)*y^2*log(x)^2)
  # loglikelihood
  if (deriv == 0) {
    -c1*(log(sigma) + 0.5*log(2*pi) + 0.5/sigma^2) +
      0.5/sigma^2*(2*c3/y0 - c6/y0^2)
  } else if (deriv == 1) {
    c(
      1/(sigma*y0)^2*(c6/y0 - c3),
      1/(sigma^2)*(c7/y0^2 - c4/y0),
      (c6/y0^2 - 2*c3/y0 + c1)/sigma^3 - c1/sigma
    )
  } else if (deriv == 2) {
    d2logL_dy02 <- 1/(sigma^2*y0^3)*(2*c3 - 3*c6/y0)
    d2logL_dy0dbeta <- 1/(sigma^2*y0^2)*(c4 - 2*c7/y0)
    d2logL_dy0dsigma <- 2/(sigma^3*y0^2)*(c3 - c6/y0)
    d2logL_dbeta2 <- 1/(sigma^2)*(c5/y0 - 2*c8/y0^2)
    d2logL_dbetadsigma <- 2/(sigma^3)*(c4/y0 - c7/y0^2)
    d2logL_dsigma2 <- c1/sigma^2 - 3/sigma^4*(c6/y0^2 - 2*c3/y0 + c1)
    matrix(
      c(
        d2logL_dy02, d2logL_dy0dbeta, d2logL_dy0dsigma,
        d2logL_dy0dbeta, d2logL_dbeta2, d2logL_dbetadsigma,
        d2logL_dy0dsigma, d2logL_dbetadsigma, d2logL_dsigma2
      ),
      nrow = 3
    )
  }
}


power_normal_scaled_root <- function(
    beta,
    x,
    y,
    weights = rep(1, length(x))
) {
  # calculate coefficients
  c3 <- sum(weights*x^(-beta)*y)
  c4 <- sum(weights*x^(-beta)*y*log(x))
  c6 <- sum(weights*x^(-2*beta)*y^2)
  c7 <- sum(weights*x^(-2*beta)*y^2*log(x))
  # root
  c3^2*c7/c6^2 - c3*c4/c6
}


#' @examples
#' # define parameters
#' y0 <- 20
#' beta <- -0.5
#' sigma <- 0.2
#' x <- seq(1, 8, l = 51)
#' y <- y0*x^beta*stats::rnorm(length(x), 1, sigma)
#' w <- stats::runif(length(x), 0.8, 1.2)
#'
#' eps <- 1e-6
#' f0 <- power_normal_scaled_root(beta, x, y, weights = w)
#' f1 <- power_normal_scaled_root(beta + eps, x, y, weights = w)
#' (f1 - f0)/eps
#' power_normal_scaled_root_jacobian(beta + eps, x, y, weights = w)
#'
power_normal_scaled_root_jacobian <- function(
    beta,
    x,
    y,
    weights = rep(1, length(x))
) {
  # calculate coefficients
  c1 <- sum(weights)
  c2 <- sum(weights*log(x))
  c3 <- sum(weights*x^(-beta)*y)
  c4 <- sum(weights*x^(-beta)*y*log(x))
  c5 <- sum(weights*x^(-beta)*y*log(x)^2)
  c6 <- sum(weights*x^(-2*beta)*y^2)
  c7 <- sum(weights*x^(-2*beta)*y^2*log(x))
  c8 <- sum(weights*x^(-2*beta)*y^2*log(x)^2)
  # jacobian of root
  (c4^2 + c3*c5)/c6 +
    -2*c3*(2*c4*c7 + c3*c8)/c6^2 +
    4*c3^2*c7^2/c6^3
}
