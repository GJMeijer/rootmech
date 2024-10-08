% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/powerlaw_gamma.R
\name{powerlaw_gamma_root_beta_jacobian}
\alias{powerlaw_gamma_root_beta_jacobian}
\title{Jacobian of function `powerlaw_gamma_root_beta()`}
\usage{
powerlaw_gamma_root_beta_jacobian(beta, x, y, weights = rep(1, length(x)))
}
\arguments{
\item{beta}{power-law exponent to solve for}

\item{x}{measured x-values (e.g. root diameters). These values are assumed
to be already normalised by a reference value to ensure a dimensionless
parameter.}

\item{y}{measured y-values (e.g. root tensile strength)}

\item{weights}{vector with weighting for each observation (default: 1)}
}
\value{
derivative of `powerlaw_gamma_root_beta()` with respect to `beta`
}
\description{
Returns the derivative of the function `powerlaw_gamma_root_beta()`
with respect to input argument `beta`
}
\keyword{internal}
