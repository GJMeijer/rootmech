# data ####
shape <- 4
x <- stats::rweibull(50, shape, 2)
w <- stats::runif(length(x))

# check derivative of weibull root solving ####
eps <- 1e-6
v0 <- weibull_root(shape, x, weights = w)
v1 <- weibull_root(shape + eps, x, weights = w)
J_num <- (v1 - v0)/eps
J_ana <- weibull_root_jacobian(shape, x, weights = w)
testthat::expect_equal(J_num, J_ana, tolerance = 0.001)
