y0 <- 20
beta <- -0.5
kappa <- 10
lambda <- 1/gamma(1 + 1/kappa)

x <- seq(1, 8, l = 1000)
y <- y0*x^beta*rweibull(length(x), kappa, lambda)
weights <- rep(1, length(x))

c1 <- sum(weights)
c2 <- sum(weights*log(x))
c3 <- sum(weights*log(y))
c4 <- sum(weights*x^(-beta*kappa)*y^kappa)
c5 <- sum(weights*x^(-beta*kappa)*y^kappa*log(x))
c6 <- sum(weights*x^(-beta*kappa)*y^kappa*log(y))
c7 <- sum(weights*x^(-beta*kappa)*y^kappa*log(x)^2)
c8 <- sum(weights*x^(-beta*kappa)*y^kappa*log(x)*log(y))
c9 <- sum(weights*x^(-beta*kappa)*y^kappa*log(y)^2)

zeta <- gamma(1 + 1/kappa)/y0
logL <- (log(kappa) + log(gamma(1 + 1/kappa)) + (kappa - 1)*log(zeta))*c1 + (kappa - 1)*c3 - beta*(kappa - 1)*c2 - zeta^kappa*c4
ybar <- y0*x^beta
lambda <- 1/gamma(1 + 1/kappa)
logL2 <- sum(weights*dweibull(y/ybar, kappa, lambda, log = TRUE))
c(logL, logL2)

zeta <- ((kappa - 1)*c1/(kappa*c4))^(1/kappa)
logL <- (log(kappa) + log(gamma(1 + 1/kappa)))*c1 + (kappa - 1)/kappa*(log(kappa - 1) + log(c1) - log(kappa) - log(c4) - 1)*c1 + (kappa - 1)*c3 - beta*(kappa - 1)*c2
y0 <- gamma(1 + 1/kappa)/zeta
ybar <- y0*x^beta
logL2 <- sum(weights*dweibull(y/ybar, kappa, lambda, log = TRUE))
c(logL, logL2)

dlogLdbeta <- (kappa - 1)*(c1*c5/c4 - c2)
logL <- (log(kappa) + log(gamma(1 + 1/kappa)))*c1 + (kappa - 1)/kappa*(log(kappa - 1) + log(c1) - log(kappa) - log(c4) - 1)*c1 + (kappa - 1)*c3 - beta*(kappa - 1)*c2
eps <- 1e-6
betaeps <- beta + eps
c4eps <- sum(weights*x^(-betaeps*kappa)*y^kappa)
logLeps <- (log(kappa) + log(gamma(1 + 1/kappa)))*c1 + (kappa - 1)/kappa*(log(kappa - 1) + log(c1) - log(kappa) - log(c4eps) - 1)*c1 + (kappa - 1)*c3 - betaeps*(kappa - 1)*c2
dlogLdbetaeps <- (logLeps - logL)/eps
c(dlogLdbeta, dlogLdbetaeps)

dlogLdkappa <- c1/kappa + c1*(kappa - 1)/kappa*((beta*c5 - c6)/c4) + c1/kappa^2*(log(kappa-1) + log(c1) - log(kappa) - log(c4) - digamma(1+1/kappa)) + c3 - beta*c2
eps <- 1e-6
kappaeps <- kappa + eps
c4kappa <- sum(weights*x^(-beta*kappaeps)*y^kappaeps)
logLeps <- (log(kappaeps) + log(gamma(1 + 1/kappaeps)))*c1 + (kappaeps - 1)/kappaeps*(log(kappaeps - 1) + log(c1) - log(kappaeps) - log(c4kappa) - 1)*c1 + (kappaeps - 1)*c3 - beta*(kappaeps - 1)*c2
dlogLdkappaeps <- (logLeps - logL)/eps
c(dlogLdkappa, dlogLdkappaeps)

y0test <- gamma(1 + 1/kappa)*(kappa*c4/((kappa - 1)*c1))^(1/kappa)
y0test

power_weibull_loglikelihood <- function(par, x, y, weights = rep(1, length(x))) {
  # unpack input parameters
  beta <- par[1]
  kappa <- par[2]
  # coefficients
  c1 <- sum(weights)
  c2 <- sum(weights*log(x))
  c3 <- sum(weights*log(y))
  c4 <- sum(weights*x^(-beta*kappa)*y^kappa)
  # loglikelihood
  (log(kappa) + log(gamma(1 + 1/kappa)))*c1 + 
    (kappa - 1)/kappa*(log(kappa - 1) + log(c1) - log(kappa) - log(c4) - 1)*c1 + 
    (kappa - 1)*c3 - beta*(kappa - 1)*c2
  # test loglikelihood 
  #zeta <- ((kappa - 1)*c1/(kappa*c4))^(1/kappa)
  #y0 <- gamma(1 + 1/kappa)/zeta
  #ypred <- y0*x^beta
  #sum(weights*dweibull(y/ypred, kappa, 1/gamma(1 + 1/kappa), log = TRUE))
}

power_weibull_root <- function(par, x, y, weights = rep(1, length(x))) {
  # unpack input parameters
  beta <- par[1]
  kappa <- par[2]
  # coefficients
  c1 <- sum(weights)
  c2 <- sum(weights*log(x))
  c3 <- sum(weights*log(y))
  c4 <- sum(weights*x^(-beta*kappa)*y^kappa)
  c5 <- sum(weights*x^(-beta*kappa)*y^kappa*log(x))
  c6 <- sum(weights*x^(-beta*kappa)*y^kappa*log(y))
  # roots
  dlogL_dbeta <- (kappa - 1)*(c1*c5/c4 - c2)
  dlogL_dkappa <- c1/kappa + c1*(kappa - 1)/(c4*kappa)*(beta*c5 - c6) + 
    c1/(kappa^2)*(log(kappa - 1) + log(c1) - log(kappa) - log(c4) - digamma(1 + 1/kappa)) + 
    c3 - beta*c2
  # return
  c(dlogL_dbeta, dlogL_dkappa)
}

power_weibull_root_jacobian <- function(par, x, y, weights = rep(1, length(x))) {
  # unpack input parameters
  beta <- par[1]
  kappa <- par[2]
  # coefficients
  c1 <- sum(weights)
  c2 <- sum(weights*log(x))
  c3 <- sum(weights*log(y))
  c4 <- sum(weights*x^(-beta*kappa)*y^kappa)
  c5 <- sum(weights*x^(-beta*kappa)*y^kappa*log(x))
  c6 <- sum(weights*x^(-beta*kappa)*y^kappa*log(y))
  c7 <- sum(weights*x^(-beta*kappa)*y^kappa*log(x)^2)
  c8 <- sum(weights*x^(-beta*kappa)*y^kappa*log(x)*log(y))
  c9 <- sum(weights*x^(-beta*kappa)*y^kappa*log(y)^2)
  # derivatives
  d2logL_dbeta2 <- kappa*(kappa - 1)*c1/c4*(c5^2/c4 - c7)
  d2logL_dbetadkappa <- c1*c5/c4 - c2 + beta*(kappa - 1)*c1/c4*(c5^2/c4 - c7) + 
    (kappa - 1)*c1/c4*(c8 - c5*c6/c4)
  d2logL_dkappa2 <- c1*(kappa - 1)/(c4*kappa)*((beta*c5 - c6)^2/c4 - beta^2*c7 + 2*beta*c8 - c9) + 
    2*c1/kappa^2*((2 - kappa)/(2*kappa - 2) + (beta*c5 - c6)/c4) + 
    -2*c1/kappa^3*(log(kappa - 1) + log(c1) - log(kappa) - log(c4) - digamma(1 + 1/kappa) + 0.5) + 
    c1/kappa^4*psigamma(1 + 1/kappa, deriv = 1)
  # return matrix
  matrix(
    c(d2logL_dbeta2, d2logL_dbetadkappa, d2logL_dbetadkappa, d2logL_dkappa2),
    nrow = 2
  )
}

power_weibull <- function(x, y, weights = rep(1, length(x))) {
  # initial guess
  ft0 <- stats::lm(log(y) ~ log(x), weights = weights)
  beta0 <- as.numeric(ft0$coefficients[2])
  kappa0 <- weibull_fit(y/x^beta0)$shape
  # fit using Newton-Raphson
  sol <- rootSolve::multiroot(
    power_weibull_root,
    c(beta0, kappa0),
    jacfunc = power_weibull_root_jacobian,
    x = x,
    y = y,
    weights = weights
  )
  beta <- sol$root[1]
  kappa <- sol$root[2]
  # calculate multiplier
  c1 <- sum(weights)
  c4 <- sum(weights*x^(-beta*kappa)*y^kappa)
  y0 <- gamma(1 + 1/kappa)*(kappa*c4/((kappa - 1)*c1))^(1/kappa)
  # loglikelihood
  logL <- power_weibull_loglikelihood(sol$root, x, y, weights = weights)
  # return list
  list(
    loglikelihood = logL,
    multiplier = y0,
    exponent = beta,
    shape = kappa
  )
}

sol <- rootSolve::multiroot(
  power_weibull_root,
  jacfunc = power_weibull_root_jacobian,
  c(beta, kappa),
  x = x, 
  y = y, 
  weights = weights
)
beta_fit <- sol$root[1]
kappa_fit <- sol$root[2]
c1 <- sum(weights)
c4 <- sum(weights*x^(-beta_fit*kappa_fit)*y^kappa_fit)
y0_fit <- gamma(1 + 1/kappa_fit)*((kappa_fit*c4)/((kappa_fit - 1)*c1))^(1/kappa_fit)

xp <- seq(min(x), max(x), l = 250)
yp <- y0_fit*xp^beta_fit
plot(x, y)
lines(xp, yp, col = "red")
c(y0_fit, beta_fit, kappa_fit)
ypred <- y0_fit*x^beta_fit

library("tidyverse")
df <- expand_grid(
  beta = beta_fit + 0.1*seq(-1, 1, l = 51),
  kappa = kappa_fit + 1*seq(-1, 1, l = 51)
) %>%
  mutate(logL = purrr::map2_dbl(
    beta, kappa,
    function(b, k) {power_weibull_loglikelihood(c(b, k), x, y, weights = weights)}
  ))
ggplot(df, aes(x = beta, y = kappa, z = logL, fill = logL)) + 
  geom_tile() +
  geom_contour()


par <- c(beta, kappa) #sol$root
eps <- 1e-6
f0 <- power_weibull_loglikelihood(par, x, y, weights = weights)
f1 <- power_weibull_loglikelihood(par + c(eps, 0), x, y, weights = weights)
f2 <- power_weibull_loglikelihood(par + c(0, eps), x, y, weights = weights)
(c(f1, f2) - f0)/eps
power_weibull_root(par, x, y, weights = weights)

j0 <- power_weibull_root(par, x, y, weights = weights)
j1 <- power_weibull_root(par + c(eps, 0), x, y, weights = weights)
j2 <- power_weibull_root(par + c(0, eps), x, y, weights = weights)
(cbind(j1, j2) - j0)/eps
power_weibull_root_jacobian(par, x, y, weights = weights)
