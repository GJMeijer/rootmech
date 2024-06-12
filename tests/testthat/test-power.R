# Power - test derivative of loglikelihood function ####
# generate some data
b <- 2
x <- seq(1, 7, l = 101)^(1/b)

# derivative of loglikelihood
eps <- 1e-6
f0 <- power_loglikelihood(b, x)
f1 <- power_loglikelihood(b + eps, x)
J_num <- (f1 - f0)/eps
J_ana <- power_root(b, x)
testthat::expect_equal(J_num, J_ana, tolerance = 0.001)
