#' Power-law fit with scaled uniformly distributed residuals
#'
#' @md
#' @description
#' Power-law fit of (x, y) data. For each value of x, the mean
#' of y is described by
#'
#'   mu = y0*x^beta
#'
#' where y0 is the power law multiplier, and beta the power law exponent.
#'
#' The width of the uniform distribution scales with the mean, i.e.
#'
#'   w = w0*x^beta
#'
#' where w0 is the power-law multiplier for the width of the uniform
#' distribution at x = 1.
#'
#' The optimal fitting parameters y0, beta and w0 are found by maximising
#' the (weighted) loglikelihood.
#'
#' @inheritParams powerlaw_gamma_fit
#' @param method if `method == "chull"`, optional values for the power-law
#'   exponent are acquired from the vertices of a convex hull circumscribing
#'   the log-transformed x,y data. The value minimising the loglikelihood
#'   function is then used. If `method == `uniroot`, the best fitting value
#'   is found using bisection of first derivative of the loglikelihood
#'   function. `"chull"` is more accurate (no approximation) but may be slower
#'   for large datasets.
#' @return a list containing the fields
#'   * `loglikelihood`: the log-likelihood score of the fit
#'   * `multiplier`: fitted power-law multiplier
#'   * `exponent`: fitted power coefficient
#'   * `width`: width of the uniform distribution at x = 1
#' @export
#' @examples
#' y0 <- 20
#' beta <- -0.5
#' range <- 10
#' x <- seq(1, 8, l = 101)
#' y <- stats::runif(length(x), (y0 - 0.5*range)*x^beta, (y0 + 0.5*range)*x^beta)
#'
#' ft <- powerlaw_uniform_fit(x, y)
#' ft
#' xp <- seq(min(x), max(x), l = 251)
#' yp <- ft$multiplier*xp^ft$exponent
#'
#' plot(x, y)
#' lines(xp, yp, col = "red")
#' lines(xp, yp*(1 - 0.5*ft$width/ft$multiplier), col = "blue")
#' lines(xp, yp*(1 + 0.5*ft$width/ft$multiplier), col = "blue")
#'
powerlaw_uniform_fit <- function(
    x,
    y,
    weights = rep(1, length(x)),
    method = "chull",
    start = NULL,
    range = 0.5*c(-1, 1)
) {
  # Convex hull method
  if (method == "chull") {
    # matrix with unique x,y positions
    M <- unique(cbind(x, y, deparse.level = 0))
    # create a convex hull i log-space (clockwise order)
    lx <- log(M[, 1])
    ly <- log(M[, 2])
    i <- grDevices::chull(lx, ly)
    # differences in log(x) and log(y) for points on hull
    dlx <- diff(lx[c(i, i[1])])
    dly <- diff(ly[c(i, i[1])])
    # best fit must be parallel to one of the hull vertices. For each vertex,
    # calculate the intercept of the vertex line with the y-axis, and also the
    # intercepts of lines with the same gradient through all points on the hull.
    # The max difference in gradient is positively related to the likelihood.
    dist <- apply(
      cbind(dlx, dly, deparse.level = 0),
      1,
      function(z) {
        tmp <- ly[i] - lx[i]*z[2]/z[1]
        max(tmp, na.rm = TRUE) - min(tmp, na.rm = TRUE)
      }
    )
    # minimum distance between intercepts = best fit (thinnest band)
    imin <- which.min(dist)
    # calculate exponent
    beta <- dly[imin]/dlx[imin]
  # Bisection method, using uniroot
  } else if (method == "bisection") {
    # intial guess - log-log transform
    if (is.null(start)) {
      ft0 <- stats::lm(log(y) ~ log(x), weights = weights)
      start <- as.numeric(ft0$coefficients[2])
    }
    # find exponent using root solving
    beta <- stats::uniroot(
      powerlaw_uniform_loglikelihood,
      start + range,
      x = x,
      y = y,
      weights = weights,
      deriv = 1,
      extendInt = "downX"
    )$root
  }
  # find min and max of domain (at x = 1)
  lower <- min(y/x^beta, na.rm = TRUE)
  upper <- max(y/x^beta, na.rm = TRUE)
  # loglikelihood
  logL <- powerlaw_uniform_loglikelihood(beta, x, y, weights = weights)
  # return
  list(
    loglikelihood = logL,
    multiplier = 0.5*(lower + upper),
    exponent = beta,
    width = upper - lower
  )
}


#' Calculate power-law normal (uniform) loglikelihood
#'
#' @description
#' Calculate the weighted loglikelihood for a power-law fit with uniform
#' residuals. Can also be used to calculate the first partial
#' derivative with respect to fitting parameters.
#'
#' @inheritParams powerlaw_uniform_fit
#' @param beta power-law exponent
#' @param deriv order of partial derivative requested
#' @return loglikelihood, or its partial derivatives to order `deriv`
#' @examples
#' # generate some data
#' y0 <- 20
#' beta <- -0.5
#' range <- 10
#' x <- seq(1, 8, l = 101)
#' y <- stats::runif(length(x), (y0 - 0.5*range)*x^beta, (y0 + 0.5*range)*x^beta)
#'
#' # fit
#' ft <- powerlaw_uniform_fit(x, y)
#'
#' # check loglikelihood
#' powerlaw_uniform_loglikelihood(ft$exponent, x, y)
#' sum(stats::dunif(
#'   y,
#'   min = (ft$multiplier - 0.5*ft$width)*x^ft$exponent - 1e-9,
#'   max = (ft$multiplier + 0.5*ft$width)*x^ft$exponent + 1e-9,
#'   log = TRUE
#' ))
#'
powerlaw_uniform_loglikelihood <- function(
    beta,
    x,
    y,
    weights = rep(1, length(x)),
    deriv = 0
) {
  if (deriv == 0) {
    lower <- min(y/x^beta, na.rm = TRUE)
    upper <- max(y/x^beta, na.rm = TRUE)
    -log(upper - lower)*sum(weights) - beta*sum(weights*log(x))
  } else if (deriv == 1) {
    il <- which.min(y/x^beta)
    iu <- which.max(y/x^beta)
    lower <- y[il]/x[il]^beta
    upper <- y[iu]/x[iu]^beta
    dlower_dbeta <- -y[il]*log(x[il])/x[il]^beta
    dupper_dbeta <- -y[iu]*log(x[iu])/x[iu]^beta
    (dlower_dbeta - dupper_dbeta)/(upper - lower)*sum(weights) - sum(weights*log(x))
  }
}


#' Calculate Kolmogorov-Smirnov parameters for power-law fit + uniform
#'
#' @description
#' Calculate Kolmogorov-Smirnov parameter for power-law fit with uniformly
#' distributed residuals
#'
#' @md
#' @inheritParams powerlaw_normal_fit
#' @param multiplier,exponent multiplier and exponent for power-law fit
#'   describing the mean
#' @param width with of the distribution, at x = 1
#' @return list with fields
#'   * `ks_distance`: Kolmogorov-Smirnov distance
#' @examples
#' y0 <- 20
#' beta <- -0.5
#' range <- 10
#' x <- seq(1, 8, l = 51)
#' y <- stats::runif(length(x), (y0 - 0.5*range)*x^beta, (y0 + 0.5*range)*x^beta)
#'
#' ft <- powerlaw_uniform_fit(x, y)
#'
#' powerlaw_uniform_ks(x, y, ft$multiplier, ft$exponent, ft$width)
#'
powerlaw_uniform_ks <- function(
    x,
    y,
    multiplier,
    exponent,
    width
) {
  # prediction
  C2 <- sort(stats::punif(
    y,
    (multiplier - 0.5*width)*x^exponent,
    (multiplier + 0.5*width)*x^exponent
  ))
  # real cumulative
  C1_lower <- seq(length(x))/(1 + length(x))
  C1_upper <- C1_lower + 1/(1 + length(x))
  # distance
  ks_distance <- max(abs(c(C1_lower - C2, C1_upper - C2)))
  # return list (in case of future expansion)
  list(
    ks_distance = ks_distance
  )
}

