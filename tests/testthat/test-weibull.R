# data ####
shape <- 4
scale <- 2
x <- stats::rweibull(50, shape, scale)
w <- stats::runif(length(x))

# check derivative of weibull root solving ####
eps <- 1e-6
v0 <- weibull_root(shape, x, weights = w)
v1 <- weibull_root(shape + eps, x, weights = w)
J_num <- (v1 - v0)/eps
J_ana <- weibull_root_jacobian(shape, x, weights = w)
testthat::expect_equal(J_num, J_ana, tolerance = 0.001)

# test loglikelihood
L0 <- weibull_loglikelihood(c(shape, scale), x, weights = w)
L1 <- sum(w*stats::dweibull(x, shape, scale, log = TRUE))
testthat::expect_equal(L0, L1)

# test first derivative of loglikelihood
f0 <- weibull_loglikelihood(c(shape, scale), x, weights = w, deriv = 0)
f1 <- weibull_loglikelihood(c(shape + eps, scale), x, weights = w, deriv = 0)
f2 <- weibull_loglikelihood(c(shape, scale + eps), x, weights = w, deriv = 0)
J_num <- (c(f1, f2) - f0)/eps
J_ana <- weibull_loglikelihood(c(shape, scale), x, weights = w, deriv = 1)
testthat::expect_equal(J_num, J_ana, tolerance = 0.001)

# test second derivative of loglikelihood
f0 <- weibull_loglikelihood(c(shape, scale), x, weights = w, deriv = 1)
f1 <- weibull_loglikelihood(c(shape + eps, scale), x, weights = w, deriv = 1)
f2 <- weibull_loglikelihood(c(shape, scale + eps), x, weights = w, deriv = 1)
J_num <- (cbind(f1, f2, deparse.level = 0) - f0)/eps
J_ana <- weibull_loglikelihood(c(shape, scale), x, weights = w, deriv = 2)
testthat::expect_equal(J_num, J_ana, tolerance = 0.001)
