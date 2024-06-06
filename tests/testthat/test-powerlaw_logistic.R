# data ####
y0 <- 20
beta <- -0.5
s0 <- 4
x <- seq(1, 6, l = 101)
y <- abs(stats::rlogis(length(x), y0*x^beta, s0*x^beta))
w <- runif(length(x), 0.8, 1.2)
par <- c(y0, beta, s0)

# test loglikelihood ####
logL1 <- sum(w*stats::dlogis(y, y0*x^beta, s0*x^beta, log = TRUE))
logL2 <- powerlaw_logistic_loglikelihood(par, x, y, weights = w, deriv = 0)
testthat::expect_equal(logL1, logL2)

# test first derivative of loglikelihood ####
eps <- 1e-9
f0 <- powerlaw_logistic_loglikelihood(par, x, y, weights = w)
f1 <- powerlaw_logistic_loglikelihood(par + c(eps, 0, 0), x, y, weights = w)
f2 <- powerlaw_logistic_loglikelihood(par + c(0, eps, 0), x, y, weights = w)
f3 <- powerlaw_logistic_loglikelihood(par + c(0, 0, eps), x, y, weights = w)
J_num <- (c(f1, f2, f3) - f0)/eps
J_ana <- powerlaw_logistic_loglikelihood(par, x, y, weights = w, deriv = 1)
testthat::expect_equal(J_num, J_ana, tolerance = 0.001)

# test second derivative of loglikelihood ####
f0 <- powerlaw_logistic_loglikelihood(par, x, y, weights = w, deriv = 1)
f1 <- powerlaw_logistic_loglikelihood(par + c(eps, 0, 0), x, y, weights = w, deriv = 1)
f2 <- powerlaw_logistic_loglikelihood(par + c(0, eps, 0), x, y, weights = w, deriv = 1)
f3 <- powerlaw_logistic_loglikelihood(par + c(0, 0, eps), x, y, weights = w, deriv = 1)
J_num <- (cbind(f1, f2, f3, deparse.level = 0) - f0)/eps
J_ana <- powerlaw_logistic_loglikelihood(par, x, y, weights = w, deriv = 2)
testthat::expect_equal(J_num, J_ana, tolerance = 0.001)
