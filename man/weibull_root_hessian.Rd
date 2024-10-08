% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/weibull.R
\name{weibull_root_hessian}
\alias{weibull_root_hessian}
\title{Hessian of root to solve for Weibull fitting}
\usage{
weibull_root_hessian(shape, x, weights = rep(1, length(x)))
}
\arguments{
\item{shape}{shape parameter}

\item{x}{array with observations}

\item{weights}{weighting factor for each observation in \code{x}}
}
\value{
second derivative of function `weibull_root()` with respect to input
  argument `shape`
}
\description{
Jacobian of root equation to solve for Weibull fitting.
}
\keyword{internal}
