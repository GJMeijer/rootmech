# data ####
x <- seq(1, 7, l = 101)
y0 <- 25
beta <- -0.5
k <- 4
y <- y0*x^beta*stats::rgamma(length(x), shape = k, scale = 1/k)
w <- stats::runif(length(x), 0.8, 1.2)

# test loglikelihood ####
logL1 <- sum(w*stats::dgamma(y, k, scale = y0*x^beta/k, log = TRUE))
logL2 <- powerlaw_gamma_loglikelihood(c(y0, beta, k), x, y, weights = w)
testthat::expect_equal(logL1, logL2)

# test first derivative of loglikelihood ####
eps <- 1e-9
par <- c(y0, beta, k)
f0 <- powerlaw_gamma_loglikelihood(par, x, y, weights = w)
f1 <- powerlaw_gamma_loglikelihood(par + c(eps, 0, 0), x, y, weights = w)
f2 <- powerlaw_gamma_loglikelihood(par + c(0, eps, 0), x, y, weights = w)
f3 <- powerlaw_gamma_loglikelihood(par + c(0, 0, eps), x, y, weights = w)
J_num <- (c(f1, f2, f3) - f0)/eps
J_ana <- powerlaw_gamma_loglikelihood(par, x, y, weights = w, deriv = 1)
testthat::expect_equal(J_num, J_ana, tolerance = 0.001)

# test second derivative of loglikelihood ####
f0 <- powerlaw_gamma_loglikelihood(par, x, y, weights = w, deriv = 1)
f1 <- powerlaw_gamma_loglikelihood(par + c(eps, 0, 0), x, y, weights = w, deriv = 1)
f2 <- powerlaw_gamma_loglikelihood(par + c(0, eps, 0), x, y, weights = w, deriv = 1)
f3 <- powerlaw_gamma_loglikelihood(par + c(0, 0, eps), x, y, weights = w, deriv = 1)
J_num <- (cbind(f1, f2, f3, deparse.level = 0) - f0)/eps
J_ana <- powerlaw_gamma_loglikelihood(par, x, y, weights = w, deriv = 2)
testthat::expect_equal(J_num, J_ana, tolerance = 0.001)

# derivative of powerlaw_gamma_root_beta ####
eps <- 1e-9
f0 <- powerlaw_gamma_root_beta(beta, x, y, weights = w)
f1 <- powerlaw_gamma_root_beta(beta + eps, x, y, weights = w)
J_num <- (f1 - f0)/eps
J_ana <- powerlaw_gamma_root_beta_jacobian(beta, x, y, weights = w)
testthat::expect_equal(J_num, J_ana, tolerance = 0.001)

# derivative of powerlaw_gamma_root_k ####
eps <- 1e-9
f0 <- powerlaw_gamma_root_k(k, beta, x, y, weights = w)
f1 <- powerlaw_gamma_root_k(k + eps, beta, x, y, weights = w)
J_num <- (f1 - f0)/eps
J_ana <- powerlaw_gamma_root_k_jacobian(k, beta, x, y, weights = w)
testthat::expect_equal(J_num, J_ana, tolerance = 0.001)
