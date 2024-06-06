# data ####
y0 <- 20
beta <- -0.5
sdlog <- 0.2
x <- seq(1, 7, l = 51)
y <- y0*x^beta*rlnorm(length(x), -0.5*sdlog^2, sdlog)
w <- runif(length(x), 0.8, 1.2)
par <- c(y0, beta, sdlog)

# test loglikelihood ####
logL1 <- sum(w*stats::dlnorm(y, log(y0) + beta*log(x) - 0.5*sdlog^2, sdlog, log = TRUE))
logL2 <- powerlaw_lognormal_loglikelihood(par, x, y, weights = w)
testthat::expect_equal(logL1, logL2)

# test first derivative of loglikelihood ####
eps <- 1e-9
f0 <- powerlaw_lognormal_loglikelihood(par, x, y, weights = w)
f1 <- powerlaw_lognormal_loglikelihood(par + c(eps, 0, 0), x, y, weights = w)
f2 <- powerlaw_lognormal_loglikelihood(par + c(0, eps, 0), x, y, weights = w)
f3 <- powerlaw_lognormal_loglikelihood(par + c(0, 0, eps), x, y, weights = w)
J_num <- (c(f1, f2, f3) - f0)/eps
J_ana <- powerlaw_lognormal_loglikelihood(par, x, y, weights = w, deriv = 1)
testthat::expect_equal(J_num, J_ana, tolerance = 0.001)

# test second derivative of loglikelihood ####
f0 <- powerlaw_lognormal_loglikelihood(par, x, y, weights = w, deriv = 1)
f1 <- powerlaw_lognormal_loglikelihood(par + c(eps, 0, 0), x, y, weights = w, deriv = 1)
f2 <- powerlaw_lognormal_loglikelihood(par + c(0, eps, 0), x, y, weights = w, deriv = 1)
f3 <- powerlaw_lognormal_loglikelihood(par + c(0, 0, eps), x, y, weights = w, deriv = 1)
J_num <- (cbind(f1, f2, f3, deparse.level = 0) - f0)/eps
J_ana <- powerlaw_lognormal_loglikelihood(par, x, y, weights = w, deriv = 2)
testthat::expect_equal(J_num, J_ana, tolerance = 0.001)

# compare fit results to linear regression ####
ft <- powerlaw_lognormal_fit(x, y, weights = w)
par1 <- c(ft$multiplier, ft$exponent, ft$sdlog)
ftl <- stats::lm(log(y) ~ log(x), weights = w)
betaL <- as.numeric(ftl$coef[2])
alphaL <- as.numeric(ftl$coef[1])
sdL <- sqrt(1/length(x)*sum((log(y) - alphaL - betaL*log(x))^2))
par2 <- c(exp(alphaL + sdL^2/2), betaL, sdL)
testthat::expect_equal(par1, par2, tolerance = 1e-3)

