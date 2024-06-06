# Error functions
# 07/08/2023 - G. J. Meijer


#' Error function
#'
#' @description
#' Calculate the value of the error function
#'
#' @param x input values for the error function
#' @return erf(x)
#' @export
#' @examples
#' erf(0.5)
#'
erf <- function(x) {
  stats::pchisq(2*x^2, 1)*sign(x)
}


#' Inverse error function
#'
#' @description
#' Calculate the value of the inverse error function
#'
#' @param x input values for the inverse error function
#' @return erf^{-1}(x)
#' @export
#' @examples
#' erfi(0.5)
#'
erfi  <- function(x) {
  sqrt(stats::qchisq(abs(x), 1)/2)*sign(x)
}
