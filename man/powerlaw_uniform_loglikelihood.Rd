% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/powerlaw_uniform.R
\name{powerlaw_uniform_loglikelihood}
\alias{powerlaw_uniform_loglikelihood}
\title{Calculate power-law normal (uniform) loglikelihood}
\usage{
powerlaw_uniform_loglikelihood(
  beta,
  x,
  y,
  weights = rep(1, length(x)),
  deriv = 0
)
}
\arguments{
\item{beta}{power-law exponent}

\item{x}{measured x-values (e.g. root diameters). These values are assumed
to be already normalised by a reference value to ensure a dimensionless
parameter.}

\item{y}{measured y-values (e.g. root tensile strength)}

\item{weights}{vector with weighting for each observation (default: 1)}

\item{deriv}{order of partial derivative requested}
}
\value{
loglikelihood, or its partial derivatives to order `deriv`
}
\description{
Calculate the weighted loglikelihood for a power-law fit with uniform
residuals. Can also be used to calculate the first partial
derivative with respect to fitting parameters.
}
\keyword{internal}
