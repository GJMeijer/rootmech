# generate some data ####
x <- stats::rweibull(50, shape = 4, scale = 2)
mu <- 2
theta <- 0.5
par <- c(mu, theta)
w <- stats::runif(length(x), 0.9, 1.1)

# check if loglikelihood gives correct answer ####
L1 <- gumbel_loglikelihood(par, x, weights = w)
L2 <- sum(w*dgumbel(x, mu, theta, log = TRUE))
testthat::expect_equal(L1, L2)

# test first derivative of loglikelihoood function
eps <- 1e-6
f0 <- gumbel_loglikelihood(par, x, weights = w, deriv = 0)
f1 <- gumbel_loglikelihood(par + c(eps, 0), x, weights = w, deriv = 0)
f2 <- gumbel_loglikelihood(par + c(0, eps), x, weights = w, deriv = 0)
J_num <- (c(f1, f2) - f0)/eps
J_ana <- gumbel_loglikelihood(par, x, weights = w, deriv = 1)
testthat::expect_equal(J_num, J_ana, tolerance = 0.001)

# test second derivative of loglikelihood function
f0 <- gumbel_loglikelihood(par, x, weights = w, deriv = 1)
f1 <- gumbel_loglikelihood(par + c(eps, 0), x, weights = w, deriv = 1)
f2 <- gumbel_loglikelihood(par + c(0, eps), x, weights = w, deriv = 1)
J_num <- ((cbind(f1, f2, deparse.level = 0) - f0)/eps)
J_ana <- gumbel_loglikelihood(par, x, weights = w, deriv = 2)
testthat::expect_equal(J_num, J_ana, tolerance = 0.001)

# check if root is first partial derivative of loglikelihood ####
# numerical derivative
eps <- 1e-6
f0 <- gumbel_loglikelihood(par, x, weights = w)
f1 <- gumbel_loglikelihood(par + c(eps, 0), x, weights = w)
f2 <- gumbel_loglikelihood(par + c(0, eps), x, weights = w)
J_num <- (c(f1, f2) - f0)/eps
# analytical derivative
J_ana <- gumbel_loglikelihood(par, x, deriv = 1, weights = w)
# test
testthat::expect_equal(J_num, J_ana, tolerance = 0.001)

# test derivative of gumbel_root ####
# numerical derivative
eps <- 1e-6
v0 <- gumbel_root(theta, x, weights = w)
v1 <- gumbel_root(theta + eps, x, weights = w)
(v1 - v0)/eps
J_num <- gumbel_root_jacobian(theta, x, weights = w)
# analytical derivative
J_ana <- gumbel_root_jacobian(theta, x, weights = w)
# test
testthat::expect_equal(J_num, J_ana, tolerance = 0.001)
