# data ####
y0 <- 20
beta <- -0.5
sigma0 <- 2
delta <- -0.75
x <- seq(1, 8, l = 101)
mu <- y0*x^beta
sigma <- sigma0*x^delta
y <- stats::rnorm(length(x), mu, sigma)
w <- stats::runif(length(x), 0.8, 1.2)

# test powerlaw_normal_freebetadelta_loglikelihood - loglikelihood ####
par <- c(y0, beta, sigma0, delta)
mu <- y0*x^beta
sigma <- sigma0*x^delta
logL1 <- sum(w*stats::dnorm(y, mu, sigma, log = TRUE))
logL2 <- powerlaw_normal_freebetadelta_loglikelihood(par, x, y, weights = w)
testthat::expect_equal(logL1, logL2)

# test first derivative - powerlaw_normal_freebetadelta_loglikelihood ####
eps <- 1e-9
f0 <- powerlaw_normal_freebetadelta_loglikelihood(par, x, y, weights = w)
f1 <- powerlaw_normal_freebetadelta_loglikelihood(par + c(eps, 0, 0, 0), x, y, weights = w)
f2 <- powerlaw_normal_freebetadelta_loglikelihood(par + c(0, eps, 0, 0), x, y, weights = w)
f3 <- powerlaw_normal_freebetadelta_loglikelihood(par + c(0, 0, eps, 0), x, y, weights = w)
f4 <- powerlaw_normal_freebetadelta_loglikelihood(par + c(0, 0, 0, eps), x, y, weights = w)
J_num <- (c(f1, f2, f3, f4) - f0)/eps
J_ana <- powerlaw_normal_freebetadelta_loglikelihood(par, x, y, weights = w, deriv = 1)
testthat::expect_equal(J_num, J_ana, tolerance = 0.001)

# test second derivative - powerlaw_normal_freebetadelta_loglikelihood ####
eps <- 1e-9
f0 <- powerlaw_normal_freebetadelta_loglikelihood(par, x, y, weights = w, deriv = 1)
f1 <- powerlaw_normal_freebetadelta_loglikelihood(par + c(eps, 0, 0, 0), x, y, weights = w, deriv = 1)
f2 <- powerlaw_normal_freebetadelta_loglikelihood(par + c(0, eps, 0, 0), x, y, weights = w, deriv = 1)
f3 <- powerlaw_normal_freebetadelta_loglikelihood(par + c(0, 0, eps, 0), x, y, weights = w, deriv = 1)
f4 <- powerlaw_normal_freebetadelta_loglikelihood(par + c(0, 0, 0, eps), x, y, weights = w, deriv = 1)
J_num <- (cbind(f1, f2, f3, f4, deparse.level = 0) - f0)/eps
J_ana <- powerlaw_normal_freebetadelta_loglikelihood(par, x, y, weights = w, deriv = 2)
testthat::expect_equal(J_num, J_ana, tolerance = 0.001)

# test derivative - powerlaw_normal_freebetadelta_root ####
eps <- 1e-9
par <- c(beta, delta)
f0 <- powerlaw_normal_freebetadelta_root(par, x, y, weights = w)
f1 <- powerlaw_normal_freebetadelta_root(par + c(eps, 0), x, y, weights = w)
f2 <- powerlaw_normal_freebetadelta_root(par + c(0, eps), x, y, weights = w)
J_num <- (cbind(f1, f2, deparse.level = 0) - f0)/eps
J_ana <- powerlaw_normal_freebetadelta_root_jacobian(par, x, y, weights = w)
testthat::expect_equal(J_num, J_ana, tolerance = 0.001)

# test powerlaw_normal_freebeta_loglikelihood - loglikelihood ####
par <- c(y0, beta, sigma0)
mu <- y0*x^beta
sigma <- sigma0*x^delta
logL1 <- sum(w*stats::dnorm(y, mu, sigma, log = TRUE))
logL2 <- powerlaw_normal_freebeta_loglikelihood(par, x, y, delta = delta, weights = w)
testthat::expect_equal(logL1, logL2)

# test first derivative - powerlaw_normal_freebeta_loglikelihood ####
eps <- 1e-9
f0 <- powerlaw_normal_freebeta_loglikelihood(par, x, y, delta = delta, weights = w)
f1 <- powerlaw_normal_freebeta_loglikelihood(par + c(eps, 0, 0), x, y, delta = delta, weights = w)
f2 <- powerlaw_normal_freebeta_loglikelihood(par + c(0, eps, 0), x, y, delta = delta, weights = w)
f3 <- powerlaw_normal_freebeta_loglikelihood(par + c(0, 0, eps), x, y, delta = delta, weights = w)
J_num <- (c(f1, f2, f3) - f0)/eps
J_ana <- powerlaw_normal_freebeta_loglikelihood(par, x, y, weights = w, delta = delta, deriv = 1)
testthat::expect_equal(J_num, J_ana, tolerance = 0.001)

# test second derivative - powerlaw_normal_freebeta_loglikelihood ####
eps <- 1e-9
f0 <- powerlaw_normal_freebeta_loglikelihood(par, x, y, delta = delta, weights = w, deriv = 1)
f1 <- powerlaw_normal_freebeta_loglikelihood(par + c(eps, 0, 0), x, y, delta = delta, weights = w, deriv = 1)
f2 <- powerlaw_normal_freebeta_loglikelihood(par + c(0, eps, 0), x, y, delta = delta, weights = w, deriv = 1)
f3 <- powerlaw_normal_freebeta_loglikelihood(par + c(0, 0, eps), x, y, delta = delta, weights = w, deriv = 1)
J_num <- (cbind(f1, f2, f3, deparse.level = 0) - f0)/eps
J_ana <- powerlaw_normal_freebeta_loglikelihood(par, x, y, weights = w, delta = delta, deriv = 2)
testthat::expect_equal(J_num, J_ana, tolerance = 0.001)

# test derivative - powerlaw_normal_freebeta_root ####
eps <- 1e-9
f0 <- powerlaw_normal_freebeta_root(beta, x, y, delta = delta, weights = w)
f1 <- powerlaw_normal_freebeta_root(beta + eps, x, y, delta = delta, weights = w)
J_num <- (f1 - f0)/eps
J_ana <- powerlaw_normal_freebeta_root_jacobian(beta, x, y, delta = delta, weights = w)
testthat::expect_equal(J_num, J_ana, tolerance = 0.001)

# test powerlaw_normal_freedelta_loglikelihood - loglikelihood ####
par <- c(y0, sigma0, delta)
mu <- y0*x^beta
sigma <- sigma0*x^delta
logL1 <- sum(w*stats::dnorm(y, mu, sigma, log = TRUE))
logL2 <- powerlaw_normal_freedelta_loglikelihood(par, x, y, beta = beta, weights = w)
testthat::expect_equal(logL1, logL2)

# test first derivative - powerlaw_normal_freedelta_loglikelihood ####
eps <- 1e-9
f0 <- powerlaw_normal_freedelta_loglikelihood(par, x, y, beta = beta, weights = w)
f1 <- powerlaw_normal_freedelta_loglikelihood(par + c(eps, 0, 0), x, y, beta = beta, weights = w)
f2 <- powerlaw_normal_freedelta_loglikelihood(par + c(0, eps, 0), x, y, beta = beta, weights = w)
f3 <- powerlaw_normal_freedelta_loglikelihood(par + c(0, 0, eps), x, y, beta = beta, weights = w)
J_num <- (c(f1, f2, f3) - f0)/eps
J_ana <- powerlaw_normal_freedelta_loglikelihood(par, x, y, weights = w, beta = beta, deriv = 1)
testthat::expect_equal(J_num, J_ana, tolerance = 0.001)

# test second derivative - powerlaw_normal_freedelta_loglikelihood ####
eps <- 1e-9
f0 <- powerlaw_normal_freedelta_loglikelihood(par, x, y, beta = beta, weights = w, deriv = 1)
f1 <- powerlaw_normal_freedelta_loglikelihood(par + c(eps, 0, 0), x, y, beta = beta, weights = w, deriv = 1)
f2 <- powerlaw_normal_freedelta_loglikelihood(par + c(0, eps, 0), x, y, beta = beta, weights = w, deriv = 1)
f3 <- powerlaw_normal_freedelta_loglikelihood(par + c(0, 0, eps), x, y, beta = beta, weights = w, deriv = 1)
J_num <- (cbind(f1, f2, f3, deparse.level = 0) - f0)/eps
J_ana <- powerlaw_normal_freedelta_loglikelihood(par, x, y, weights = w, beta = beta, deriv = 2)
testthat::expect_equal(J_num, J_ana, tolerance = 0.001)

# test derivative - powerlaw_normal_freebeta_root ####
eps <- 1e-6
f0 <- powerlaw_normal_freedelta_root(delta, x, y, beta = beta, weights = w)
f1 <- powerlaw_normal_freedelta_root(delta + eps, x, y, beta = beta, weights = w)
J_num <- (f1 - f0)/eps
J_ana <- powerlaw_normal_freedelta_root_jacobian(delta, x, y, beta = beta, weights = w)
testthat::expect_equal(J_num, J_ana, tolerance = 0.001)

# test powerlaw_normal_linkedbetadelta_loglikelihood - loglikelihood ####
par <- c(y0, beta, sigma0)
mu <- y0*x^beta
sigma <- sigma0*x^beta
logL1 <- sum(w*stats::dnorm(y, mu, sigma, log = TRUE))
logL2 <- powerlaw_normal_linkedbetadelta_loglikelihood(par, x, y, weights = w)
testthat::expect_equal(logL1, logL2)

# test first derivative - powerlaw_normal_linkedbetadelta_loglikelihood ####
eps <- 1e-9
f0 <- powerlaw_normal_linkedbetadelta_loglikelihood(par, x, y, weights = w)
f1 <- powerlaw_normal_linkedbetadelta_loglikelihood(par + c(eps, 0, 0), x, y, weights = w)
f2 <- powerlaw_normal_linkedbetadelta_loglikelihood(par + c(0, eps, 0), x, y, weights = w)
f3 <- powerlaw_normal_linkedbetadelta_loglikelihood(par + c(0, 0, eps), x, y, weights = w)
J_num <- (c(f1, f2, f3) - f0)/eps
J_ana <- powerlaw_normal_linkedbetadelta_loglikelihood(par, x, y, weights = w, deriv = 1)
testthat::expect_equal(J_num, J_ana, tolerance = 0.001)

# test second derivative - powerlaw_normal_linkedbetadelta_loglikelihood ####
eps <- 1e-9
f0 <- powerlaw_normal_linkedbetadelta_loglikelihood(par, x, y, weights = w, deriv = 1)
f1 <- powerlaw_normal_linkedbetadelta_loglikelihood(par + c(eps, 0, 0), x, y, weights = w, deriv = 1)
f2 <- powerlaw_normal_linkedbetadelta_loglikelihood(par + c(0, eps, 0), x, y, weights = w, deriv = 1)
f3 <- powerlaw_normal_linkedbetadelta_loglikelihood(par + c(0, 0, eps), x, y, weights = w, deriv = 1)
J_num <- (cbind(f1, f2, f3, deparse.level = 0) - f0)/eps
J_ana <- powerlaw_normal_linkedbetadelta_loglikelihood(par, x, y, weights = w, deriv = 2)
testthat::expect_equal(J_num, J_ana, tolerance = 0.001)

# test derivative - powerlaw_normal_linkedbetadelta_root_jacobian ####
eps <- 1e-6
f0 <- powerlaw_normal_linkedbetadelta_root(beta, x, y, weights = w)
f1 <- powerlaw_normal_linkedbetadelta_root(beta + eps, x, y, weights = w)
J_num <- (f1 - f0)/eps
J_ana <- powerlaw_normal_linkedbetadelta_root_jacobian(beta, x, y, weights = w)
testthat::expect_equal(J_num, J_ana, tolerance = 0.001)
