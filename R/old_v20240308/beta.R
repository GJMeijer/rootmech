#' Fit beta distribution to series of data
#'
#' @description
#' Fit a beta distribution to a series of data, using loglikelihood
#' optimalisation
#'
#' @param x vector with observations
#' @param weights weighting factor for each observations
#' @param guess (optional) user-defined guess for parameters
#' @return dataframe with fields `loglikelihood`, `alpha` and `beta`
#' @export
#' @examples
#' x <- rbeta(500, 2, 0.5)
#' ft <- beta_fit(x)
#' xp <- seq(0, 1, l = 101)
#' yp <- dbeta(xp, ft[1], ft[2])
#'
#' hist(x, freq = FALSE)
#' lines(xp, yp)
#'
beta_fit <- function(x, weights = rep(1, length(x)), guess = NULL) {
  # generate guess
  if (!(is.numeric(guess) & (length(guess) == 2))) {
    guess <- beta_loglikelihood_initialguess(x, weights = weights)
  }
  # find best fit
  sol <- stats::optim(
    guess,
    beta_loglikelihood,
    gr = beta_loglikelihood_jacobian,
    x = x,
    weights = weights,
    method = "L-BFGS-B",
    lower = c(1e-6, 1e-6),
    upper = c(Inf, Inf),
    control = list(fnscale = -1)
  )
  # return results
  data.frame(
    loglikelihood = sol$value,
    alpha = sol$par[1],
    beta = sol$par[2]
  )
}


#' Initial guess for beta distribution fit
#'
#' @description
#' Generate a (crude) initial guess for beta distribution fit
#'
#' @inheritParams beta_fit
#' @return two-parameter vector with initial guess
#'
beta_loglikelihood_initialguess <- function(x, weights = rep(1, length(x))) {
  c(0.5, 0.5)
}


#' Loglikelihood function for beta distribution
#'
#' @description
#' Returns the (weighted) loglikelihood score for a beta distribution fit
#'
#' @inheritParams beta_fit
#' @param val vector with guess for the two beta distribution parameters
#' @return loglikelihood score
#'
beta_loglikelihood <- function(val, x, weights = rep(1, length(x))) {
  c1 <- sum(weights)
  c2 <- sum(weights*log(x))
  c3 <- sum(weights*log(1 - x))
  -lbeta(val[1], val[2])*c1 + (val[1] - 1)*c2 + (val[2] - 1)*c3
}


#' Jacobian of `beta_loglikelihood()` function
#'
#' @description
#' Returns the derivative of the outout from `beta_loglikelihood()` function
#' with respect to input argument `val`
#'
#' @inheritParams beta_loglikelihood
#' @return vector with derivatives with respect to `val`#'
#' @examples
#' x <- runif(25)
#' var <- c(0.4, 0.7)
#'
#' eps <- 1e-6
#' z0 <- beta_loglikelihood(var, x)
#' z1 <- beta_loglikelihood(var + c(eps, 0), x)
#' z2 <- beta_loglikelihood(var + c(0, eps), x)
#' (c(z1, z2) - z0)/eps
#' beta_loglikelihood_jacobian(var, x)
#'
beta_loglikelihood_jacobian <- function(val, x, weights = rep(1, length(x))) {
  c1 <- sum(weights)
  c2 <- sum(weights*log(x))
  c3 <- sum(weights*log(1 - x))
  c(
    (digamma(val[1] + val[2]) - digamma(val[1]))*c1 + c2,
    (digamma(val[1] + val[2]) - digamma(val[2]))*c1 + c3
  )
}

