#' Compare two numeric vectors
#'
#' @description
#' This is a safe way of comparing if two vectors of floating
#' point numbers are (pairwise) equal. This is safer than
#' using `==``, because it has a built in tolerance. This function
#' is based on `dplyr::near()`.
#'
#' @param x,y Numeric vectors to compare
#' @param tol tolerance of comparison (optional)
#' @return a logical array (same size as input arrays) with
#'   element-wise comparison of `x` and `y``
#' @examples
#' sqrt(2)^2 == 2
#' is_near(sqrt(2)^2, 2)
#'
is_near <- function (x, y, tol = .Machine$double.eps^0.5) {
  abs(x - y) < tol
}


#' expand.grid function for dataframes
#'
#' @description
#' Generate all combinations of rows in a series of data.frames. This is the
#' dataframe equivalent of the standard R function `expand.grid()`, which only
#' works on vectors.
#'
#' The code is taken from Stack Overflow:
#' https://stackoverflow.com/questions/11693599/alternative-to-expand-grid-for-data-frames
#'
#' Iterations loop fastest through dataframes specified first.
#'
#' @param ... input dataframes
#' @return dataframe with all combinations of rows
#' @examples
#' df1 <- data.frame(a = c(1, 2), b = c(3, 4))
#' df2 <- data.frame(d = c(1, 2, 3), e = c(4, 5, 6), f = c(7, 8, 9))
#' expand.grid.df(df1, df2)
#'
expand.grid.df <- function(...) {
  Reduce(function(...) merge(..., by=NULL), list(...))
}
