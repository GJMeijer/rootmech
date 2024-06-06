#' Anderson-Darling test for normally distributed data
#'
#' @description
#' Anderson-Darling test on (weighted) normally distributed data.
#'
#' Following the procedure outlines on page 8-41 of the book
#' "Composite Materials Handbook, Volume 1 - Polymer Matrix Composites -
#' Guidelines for Characterization of Structural Materials (CMH-17)"
#'
#' This follows the method used by the `cmstatr::anderson_darling_normal)`
#' method, with the exception that data may be weigted when obtaining the mean
#' and standard deviation of the data
#'
#' @param x vector with observations
#' @param weights weighting for individual observations in `x`
#' @return list with values for the Anderson-Darling test statistic (`A`), the
#'   adjusted statistic (adjusted for sample size, `A_star`) and the
#'   likelihood we are incorrect in assuming the sample is NOT normally
#'   distributed
#' @export
#' @examples
#' x <- stats::rnorm(1000, 5, 2)
#' anderson_darling_normal(x)
#'
anderson_darling_normal <- function(x, weights = rep(1, length(x))) {
  # order and complete cases
  weights <- weights[complete.cases(x)]
  x <- x[complete.cases(x)]
  ord <- order(x)
  weights <- weights[ord]
  x <- x[ord]
  n <- length(x)
  # weighted mean and standard deviation
  mu <- sum(weights*x)/sum(weights)
  sd <- sqrt(sum(weights*(x - mu)^2)*n/((n - 1)*sum(weights)))
  # calculate Anderson-Darling statistic
  logp1 <- stats::pnorm(x, mean = mu, sd = sd, log.p = TRUE)
  logp2 <- stats::pnorm(x, mean = mu, sd = sd, log.p = TRUE, lower.tail = FALSE)
  A <- mean((1 - 2*seq(1, n))*(logp1 + rev(logp2))) - n
  # adjusted statistic
  A_star <- (1 + 4/n - 25/n^2)*A
  # calculate p-value
  p <- 1/(1 + exp(-0.48 + 0.78*log(A_star) + 4.58*A_star))
  # return
  list(A = A, A_star = A_star, p = p)
}


#' Anderson-Darling test for Weibull distributed data
#'
#' @description
#' Anderson-Darling test on (weighted) Weibull distributed data.
#'
#' Following the procedure outlines on page 8-41 to 8-42 of the book
#' "Composite Materials Handbook, Volume 1 - Polymer Matrix Composites -
#' Guidelines for Characterization of Structural Materials (CMH-17)"
#'
#' This follows the method used by the `cmstatr::anderson_darling_weibull)`
#' method, with the exception that data may be weigted when obtaining the shape
#' and scale parameters for the data
#'
#' @param x vector with observations
#' @param weights weighting for individual observations in `x`
#' @return list with values for the Anderson-Darling test statistic (`A`), the
#'   adjusted statistic (adjusted for sample size, `A_star`) and the
#'   likelihood we are incorrect in assuming the sample is NOT normally
#'   distributed
#' @export
#' @examples
#' x <- stats::rweibull(1000, 5, 2)
#' anderson_darling_weibull(x)
#'
anderson_darling_weibull <- function(x, weights = rep(1, length(x))) {
  # order observations, and only use complete cases
  weights <- weights[complete.cases(x)]
  x <- x[complete.cases(x)]
  ord <- order(x)
  weights <- weights[ord]
  x <- x[ord]
  n <- length(x)
  # fit weibull distribution
  ft <- weibull_fit(x, weights = weights)
  # estimate z
  z <- (x/ft$scale)^ft$shape
  # calculate Anderson-darling statistic
  A <- mean((1 - 2*seq(1, n))*(log(1 - exp(-z)) - rev(z))) - n
  # adjusted statistic
  A_star <- (1 + 0.2/sqrt(n))*A
  # calculate p-value
  p <- 1/(1 + exp(-0.10 + 1.24*log(A_star) + 4.48*A_star))
  # return
  list(A = A, A_star = A_star, p = p)
}

