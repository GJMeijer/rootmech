
#' @examples
#' # generate data
#' y0 <- 20
#' beta <- -0.5
#' sigmaL <- 0.2
#' muL <- log(y0) - 0.5*sigmaL^2
#' x <- seq(1, 8, l = 10001)
#' y <- x^beta*rlnorm(length(x), muL, sigmaL)
#'
#' # fit
#' ft <- power_lognormal_fit(x, y)
#' ft
#' xp <- seq(min(x), max(x), l = 101)
#' yp <- ft$multiplier*xp^ft$exponent
#'
#' # plot
#' plot(x, y)
#' lines(xp, yp, col = "red")
#'
power_lognormal_fit <- function(x, y, weights = rep(1, length(x))) {
  # coefficients
  c1 <- sum(weights)
  c2 <- sum(weights*log(x))
  c3 <- sum(weights*log(y))
  c4 <- sum(weights*log(x)^2)
  c5 <- sum(weights*log(x)*log(y))
  c6 <- sum(weights*log(y)^2)
  # exponennt
  beta <- (c1*c5 - c2*c3)/(c1*c4 - c2^2)
  #sigmaL <- ((c1^2*c4 - c1*c2^2 + (c1*(- c2^2 + c1*c4)*(c4*c1^2 - c1*c2^2 + 10*c1*c5^2 - 10*c4*c6*c1 + 10*c6*c2^2 - 20*c2*c3*c5 + 10*c4*c3^2))^(1/2))/(5*c1*(- c2^2 + c1*c4)))^(1/2)

  #sigmaL <- (-(c1*c2^2 - c1^2*c4 + (c1*(- c2^2 + c1*c4)*(c4*c1^2 - c1*c2^2 + 10*c1*c5^2 - 10*c4*c6*c1 + 10*c6*c2^2 - 20*c2*c3*c5 + 10*c4*c3^2))^(1/2))/(5*c1*(- c2^2 + c1*c4)))^(1/2)
  #((c1^2*c4 - c1*c2^2 + (c1*(- c2^2 + c1*c4)*(c4*c1^2 - c1*c2^2 + 10*c1*c5^2 - 10*c4*c6*c1 + 10*c6*c2^2 - 20*c2*c3*c5 + 10*c4*c3^2))^(1/2))/(5*c1*(- c2^2 + c1*c4)))^(1/2)

  sigmaL <- -(c6*c2^2 - 2*c2*c3*c5 + c4*c3^2 + c1*c5^2 - c1*c4*c6)/((- c2^2 + c1*c4)*(-(c1*(c6*c2^2 - 2*c2*c3*c5 + c4*c3^2 + c1*c5^2 - c1*c4*c6))/(- c2^2 + c1*c4))^(1/2))



  y0 <- exp((c3 - beta*c2)/c1 + 0.5*sigmaL^2)
  # return list
  list(
    loglikelihood = 0,
    multiplier = y0,
    exponent = beta,
    sigmaL = sigmaL
  )
}


x <- seq(1, 10, l = 5001)
y0 <- 20
beta <- -0.5
sigmaL <- 0.2
muL <- log(y0) - 0.5*sigmaL^2
weights <- rep(1, length(x))

y <- x^beta*rlnorm(length(x), muL, sigmaL)
plot(x, y)

c1 <- sum(weights)
c2 <- sum(weights*log(x))
c3 <- sum(weights*log(y))
c4 <- sum(weights*log(x)^2)
c5 <- sum(weights*log(x)*log(y))
c6 <- sum(weights*log(y)^2)

-c1*(log(sigmaL) + 0.5*log(2*pi)) +
  beta*c2 - c3 -
  (beta^2*c4 - 2*beta*c5 + c6 + (c3 - beta*c2)*(sigmaL^2 - 2*log(y0)) + c1*(log(y0)^2 - sigmaL^2*log(y0) + 0.25*sigmaL^4))/(2*sigmaL^2)

sum(weights*dlnorm(y/x^beta, muL, sigmaL, log = TRUE))

(c3 - beta*c2)/c1 + sigmaL^2/2


logp <- log(y0) + beta*log(x) - log(y) - log(sigmaL) - 0.5*log(2*pi) -
  (log(y) - log(y0) - beta*log(x) + 0.5*sigmaL^2)^2/(2*sigmaL^2)
sum(weights*logp)

exp(3/2*sigmaL^2 + (c3 - beta*c2)/c1)

#exp(sigmaL^2*((c3 - beta*c2)/c1 + 1.5))


eps <- 1e-6
f1 <- c1*(log(y0) - log(sigmaL) - 0.5*log(2*pi)) +
  beta*c2 - c3 -
  (beta^2*c4 - 2*beta*c5 + c6 + (c3 - beta*c2)*(sigmaL^2 - 2*log(y0)) + c1*(log(y0)^2 - sigmaL^2*log(y0) + 0.25*sigmaL^4))/(2*sigmaL^2)
f2 <- c1*(log(y0 + eps) - log(sigmaL) - 0.5*log(2*pi)) +
  beta*c2 - c3 -
  (beta^2*c4 - 2*beta*c5 + c6 + (c3 - beta*c2)*(sigmaL^2 - 2*log(y0 + eps)) + c1*(log(y0 + eps)^2 - sigmaL^2*log(y0 + eps) + 0.25*sigmaL^4))/(2*sigmaL^2)
(f2 - f1)/eps

3*c1/(2*y0) + (c3 - beta*c2 - c1*log(y0))/(y0*sigmaL^2)

(c1*c5 - c2*c3)/(c1*c4 - c2^2)

((c1^2*c4 - c1*c2^2 + (c1*(- c2^2 + c1*c4)*(c4*c1^2 - c1*c2^2 + 10*c1*c5^2 - 10*c4*c6*c1 + 10*c6*c2^2 - 20*c2*c3*c5 + 10*c4*c3^2))^(1/2))/(5*c1*(- c2^2 + c1*c4)))^(1/2)
(-(c1*c2^2 - c1^2*c4 + (c1*(- c2^2 + c1*c4)*(c4*c1^2 - c1*c2^2 + 10*c1*c5^2 - 10*c4*c6*c1 + 10*c6*c2^2 - 20*c2*c3*c5 + 10*c4*c3^2))^(1/2))/(5*c1*(- c2^2 + c1*c4)))^(1/2)
