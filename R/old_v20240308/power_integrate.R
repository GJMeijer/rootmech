#' Integrate a power-law curve
#'
#' @description
#' Integrate a power law curve in the form:
#'
#'   y(x) = a*x^b
#'
#' between x=x1 and x=x2.
#'
#' Function is vectorised.
#'
#' @param power vector with powers
#' @param lower,upper vectors with lower and upper integration limits
#' @param multiplier (optional) vector with multipliers
#' @return vector with integrands
#' @export
#' @examples
#' power_integrate(2.1, 2, 5, multiplier = 3)
#' stats::integrate(function(x) 3*x^2.1, 2, 5)
#'
power_integrate <- function(power, lower, upper, multiplier = 1) {
  ifelse(
    is_near(power, -1),
    multiplier*log(upper/lower),
    multiplier*(upper^(power + 1) - lower^(power + 1))/(power + 1)
  )
}


#' Derivative of function `power_integrate()`
#'
#' @description
#' Generates the derivative of the results of the function `power_integrate()`
#' with respect to its input parameters `power`, `lower`, `upper` and
#' `multiplier`.
#'
#' Function is vectorised.
#'
#' @inheritParams power_integrate
#' @return a list containing the derivatives for each input parameter. Has
#'   fields `power`, `lower` and `upper` and `multiplier`.
#' @export
#' @examples
#' power <- -0.4
#' lower <- 1.5
#' upper <- 3.6
#' multiplier <- 2.1
#'
#' power_integrate_jacobian(power, lower, upper, multiplier = multiplier)
#'
#' I0 <- power_integrate(power, lower, upper, multiplier = multiplier)
#' eps <- 1e-6
#' (power_integrate(power + eps, lower, upper, multiplier = multiplier) - I0)/eps
#' (power_integrate(power, lower + eps, upper, multiplier = multiplier) - I0)/eps
#' (power_integrate(power, lower, upper + eps, multiplier = multiplier) - I0)/eps
#' (power_integrate(power, lower, upper, multiplier = multiplier + eps) - I0)/eps
#'
power_integrate_jacobian <- function(power, lower, upper, multiplier = 1) {
  list(
    power = multiplier*ifelse(
      is_near(power, -1),
      log(upper)^2/2 - log(lower)^2/2,
      (lower^(power + 1) - upper^(power + 1))/(power + 1)^2 -
        (lower^(power + 1)*log(lower) - upper^(power + 1)*log(upper))/(power + 1)
    ),
    lower = -multiplier*lower^power,
    upper = multiplier*upper^power,
    multiplier = power_integrate(power, lower, upper, multiplier = 1)
  )
}
