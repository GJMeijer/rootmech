# data ####
y0 <- 20
beta <- -0.5
range <- 10
x <- seq(1, 8, l = 101)
y <- stats::runif(length(x), (y0 - 0.5*range)*x^beta, (y0 + 0.5*range)*x^beta)

# test loglikelihood ####
ft <- powerlaw_uniform_fit(x, y)
logL1 <- sum(stats::dunif(
  y,
  min = (ft$multiplier - 0.5*ft$width)*x^ft$exponent - 1e-9,
  max = (ft$multiplier + 0.5*ft$width)*x^ft$exponent + 1e-9,
  log = TRUE
))
logL2 <- powerlaw_uniform_loglikelihood(ft$exponent, x, y)
testthat::expect_equal(logL1, logL2)

# test first derivative of loglikelihood ####
eps <- 1e-9
f0 <- powerlaw_uniform_loglikelihood(beta, x, y)
f1 <- powerlaw_uniform_loglikelihood(beta + eps, x, y)
J_num <- (f1 - f0)/eps
J_ana <- powerlaw_uniform_loglikelihood(beta, x, y, deriv = 1)
testthat::expect_equal(J_num, J_ana, tolerance = 0.001)
