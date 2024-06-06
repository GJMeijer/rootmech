# data ####
y0 <- 20
beta <- -0.5
theta0 <- 4
x <- seq(1, 6, l = 100)
gamma <- 0.5772156649
theta <- theta0*x^beta
mu <- (y0 - gamma*theta0)*x^beta
w <- stats::runif(length(x), 0.9, 1.1)
y <- rgumbel(length(x), mu, theta)
par <- c(y0, beta, theta0)

# test loglikelihood ####
logL1 <- sum(w*dgumbel(y, mu, theta, log = TRUE))
logL2 <- powerlaw_gumbel_loglikelihood(par, x, y, weights = w)
testthat::expect_equal(logL1, logL2)

# test first derivative of loglikelihood ####
eps <- 1e-9
f0 <- powerlaw_gumbel_loglikelihood(par, x, y, weights = w)
f1 <- powerlaw_gumbel_loglikelihood(par + c(eps, 0, 0), x, y, weights = w)
f2 <- powerlaw_gumbel_loglikelihood(par + c(0, eps, 0), x, y, weights = w)
f3 <- powerlaw_gumbel_loglikelihood(par + c(0, 0, eps), x, y, weights = w)
J_num <- (c(f1, f2, f3) - f0)/eps
J_ana <- powerlaw_gumbel_loglikelihood(par, x, y, weights = w, deriv = 1)
testthat::expect_equal(J_num, J_ana, tolerance = 0.001)

# test second derivative of loglikelihood ####
f0 <- powerlaw_gumbel_loglikelihood(par, x, y, weights = w, deriv = 1)
f1 <- powerlaw_gumbel_loglikelihood(par + c(eps, 0, 0), x, y, weights = w, deriv = 1)
f2 <- powerlaw_gumbel_loglikelihood(par + c(0, eps, 0), x, y, weights = w, deriv = 1)
f3 <- powerlaw_gumbel_loglikelihood(par + c(0, 0, eps), x, y, weights = w, deriv = 1)
J_num <- (cbind(f1, f2, f3, deparse.level = 0) - f0)/eps
J_ana <- powerlaw_gumbel_loglikelihood(par, x, y, weights = w, deriv = 2)
testthat::expect_equal(J_num, J_ana, tolerance = 0.001)

# derivative of powerlaw_gumbel_root ####
eps <- 1e-9
f0 <- powerlaw_gumbel_root(c(beta, theta0), x, y, weights = w)
f1 <- powerlaw_gumbel_root(c(beta + eps, theta0), x, y, weights = w)
f2 <- powerlaw_gumbel_root(c(beta, theta0 + eps), x, y, weights = w)
J_num <- (cbind(f1, f2, deparse.level = 0) - f0)/eps
J_ana <- powerlaw_gumbel_root_jacobian(c(beta, theta0), x, y, weights = w)
testthat::expect_equal(J_num, J_ana, tolerance = 0.001)
