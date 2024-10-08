% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/powerlaw_gamma.R
\name{powerlaw_gamma_root_beta}
\alias{powerlaw_gamma_root_beta}
\title{Root solve exponent in gamma power-law}
\usage{
powerlaw_gamma_root_beta(beta, x, y, weights = rep(1, length(x)))
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
value of root function to solve
}
\description{
Root solve equation to obtain the exponent in power-law fitting
with gamma-distribution
}
\keyword{internal}
