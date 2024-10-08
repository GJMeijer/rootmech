% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/powerlaw.R
\name{powerlaw_covariancematrix}
\alias{powerlaw_covariancematrix}
\title{Estimate variance-covariance matrix for power law fit}
\usage{
powerlaw_covariancematrix(
  x,
  y,
  model,
  method = "fisher",
  multiplier = NULL,
  exponent = NULL,
  sd_multiplier = NULL,
  sd_exponent = NULL,
  sdlog = NULL,
  shape = NULL,
  scale = NULL,
  width = NULL,
  weights = rep(1, length(x)),
  n = 500
)
}
\arguments{
\item{x}{a-values at which to generate the interval}

\item{y}{measured y-values (e.g. root tensile strength)}

\item{model}{fitting model, see function \code{powerlaw_fit()}}

\item{method}{if `fisher`, the covariance matrix is estimated as the
inverse of the negative second partial derivative of the loglikelihood
function with respect to fitting parameters. If `bootstrap`, the covariance
matrix is estimated using bootstrapping}

\item{multiplier}{fitted power law multiplier}

\item{exponent}{fitted power law exponent}

\item{sd_multiplier}{fitted standard deviation multiplier}

\item{sd_exponent}{fitted standard deviation exponent}

\item{sdlog}{fitted lognormal standard deviation}

\item{shape}{fitted shape parameter}

\item{scale}{fitted scale parameter}

\item{width}{fitted width parameter}

\item{weights}{vector with weighting for each observation (default: 1)}

\item{n}{number of bootstrap samples}
}
\value{
the variance-covariance matrix, with names rows and columns
  according to the fitting parameter
}
\description{
Estimate the variance-covariance matrix for a power law fit,
either using bootstrapping or the loglikelihood estimation
}
\examples{
# generate some data
y0 <- 20
beta <- -0.5
kappa <- 4
lambda <- 1/gamma(1 + 1/kappa)
x <- seq(1, 8, l = 51)
y <- y0*x^beta*rweibull(length(x), kappa, lambda)

ft <- powerlaw_fit(x, y, "weibull")
powerlaw_covariancematrix(
  x, y,
  "weibull",
  method = "fisher",
  multiplier = ft$multiplier,
  exponent = ft$exponent,
  shape = ft$shape
)

ft <- powerlaw_fit(x, y, "lognormal")
powerlaw_covariancematrix(
  x, y,
  "lognormal",
  method = "fisher",
  multiplier = ft$multiplier,
  exponent = ft$exponent,
  sdlog = ft$sdlog
)

ft <- powerlaw_fit(x, y, "lognormal_uncorrected")
powerlaw_covariancematrix(
  x, y,
  "lognormal_uncorrected",
  method = "fisher",
  multiplier = ft$multiplier,
  exponent = ft$exponent,
  sdlog = ft$sdlog
)

}
