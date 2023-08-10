is_near <- function (x, y, tol = .Machine$double.eps^0.5) {
  abs(x - y) < tol
}


expand.grid.df <- function(...) {
  Reduce(function(...) merge(..., by=NULL), list(...))
}
