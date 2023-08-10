# General biomechanics fitting and prediction functions
# 08/08/2023 - G. J. Meijer


#' Wrapper function for fitting biomechanical parameters
#'
#' @description
#' Generate fits for root tensile strain to failure and tensile strength
#' (+ correlation between normalised strain and strength). Uses power-law
#' fitting function `power_weibull_fit()` for power-law fitting and the
#' function `copula_gaussian_fit()` to fit the correlation coefficient
#' for the copula describnig normalised strain and strength relationships.
#'
#' @param dr,epsru,tru vectors with measured diameters, tensile strain to peak
#'   and tensile strength
#' @param dr0 root reference diameter
#' @param weights_power array with weights to use for power-law fitting
#' @param weights_copula array with weights to use for copula fitting
#' @return a dataframe with key parameters for tensile strain to peak fitting
#'   (power law parameters `epsru0` and `betaeps`, Weibull shape parameter
#'   `kappaeps`), tensile strength fitting (power law parameters `tru0` and
#'   `betat`, Weibull shape parameter `kappat`) and copula fitting (
#'   correlation coefficient `rho`)
#' @examples
#' # generate data
#' n <- 101
#' dr <- seq(1, 7, l = n)
#' df <- biomech_random(dr, 0.2, -0.2, 4, 10, -0.5, 6, 0.5)
#'
#' # fit
#' ft <- biomech_fit(df$dr, df$epsru, df$tru)
#'
#' # power-law predictions
#' dr_pred <- seq(min(dr), max(dr), l = 251)
#' epsru_pred <- ft$epsru0*dr_pred^ft$betaeps
#' tru_pred <- ft$tru0*dr_pred^ft$betat
#'
#' # prediction interval
#' df_ell <- bivariate_normal_predictioninterval(rho = ft$rho)
#' df_pred <- bivariate_normal_toweibull(df_ell$x, df_ell$y, ft$kappaeps, ft$kappat)
#'
#' # plot - tensile strain
#' plot(dr, df$epsru)
#' lines(dr_pred, epsru_pred)
#'
#' # plot - tensile strength
#' plot(dr, df$tru)
#' lines(dr_pred, tru_pred)
#'
#' # plot - relationship epsru/tru
#' plot(
#'   df$epsru/(ft$epsru0*dr^ft$betaeps),
#'   df$tru/(ft$tru0*dr^ft$betat)
#' )
#' lines(df_pred$x, df_pred$y)
#'
biomech_fit <- function(
    dr,
    epsru,
    tru,
    dr0 = 1,
    weights_power = dr^2,
    weights_copula = dr^2
) {
  # fit power-law to ultimate strain data
  ft_epsru <- power_weibull_fit(dr/dr0, epsru, weights = weights_power)
  # fit power-law to strength data
  ft_tru <- power_weibull_fit(dr/dr0, tru, weights = weights_power)
  # normalise data
  epsru_rel <- epsru/(ft_epsru$multiplier*(dr/dr0)^ft_epsru$power)
  tru_rel <- tru/(ft_tru$multiplier*(dr/dr0)^ft_tru$power)
  # cumulative probability densities, required for copula fitting
  epsru_rel_cum <- stats::pweibull(
    epsru_rel,
    shape = ft_epsru$shape,
    scale = 1/gamma(1 + 1/ft_epsru$shape)
  )
  tru_rel_cum <- stats::pweibull(
    tru_rel,
    ft_tru$shape,
    scale = 1/gamma(1 + 1/ft_tru$shape)
  )
  # fit gaussian copula
  ft_copula <- copula_gaussian_fit(
    epsru_rel_cum,
    tru_rel_cum,
    weights = weights_copula
  )
  # return dataframe with fitting variable
  data.frame(
    epsru0 = ft_epsru$multiplier,
    betaeps = ft_epsru$power,
    kappaeps = ft_epsru$shape,
    tru0 = ft_tru$multiplier,
    betat = ft_tru$power,
    kappat = ft_tru$shape,
    rho = ft_copula$rho
  )
}


#' Generate random root properties given biomechanical fit
#'
#' @description
#' Generate random biomechanical properties (strain to failure, strength) for
#' known root diameters, using biomechanical fitting parameters
#'
#' @param dr array with known root diameters
#' @param epsru0 diameter-tensile strain power-law multiplier
#' @param betaeps diameter-tensile strain power-law power coefficient
#' @param kappaeps diameter-tensile strain Weibull shape parameter
#' @param tru0 diameter-tensile strength power-law multiplier
#' @param betat diameter-tensile strength power-law power coefficient
#' @param kappat diameter-tensile strength Weibull shape parameter
#' @param rho copula correlation coefficient
#' @param dr0 root reference diameter
#' @return dataframe with root diameters (`dr`), tensile strain to peak
#' (`epsru`) and tensile strength (`tru`)
#' @examples
#' # generate data
#' dr <- seq(2, 9, l = 100)
#' df <- biomech_random(dr, 0.2, -0.3, 4, 10, -0.2, 6, 0.4)
#'
#' # plot
#' plot(df$dr, df$epsru)
#' plot(df$dr, df$tru)
#'
biomech_random <- function(
    dr,
    epsru0,
    betaeps,
    kappaeps,
    tru0,
    betat,
    kappat,
    rho,
    dr0 = 1
) {
  # draw bivariate gaussian data
  r <- bivariate_normal_random(length(dr), mu = c(0, 0), sd = c(1, 1), rho = rho)
  # convert to cumulative density
  c1 <- stats::pnorm(r[, 1], mean = 0, sd = 1)
  c2 <- stats::pnorm(r[, 2], mean = 0, sd = 1)
  # convert back to weibull observations
  x1 <- stats::qweibull(c1, shape = kappaeps, scale = 1/gamma(1 + 1/kappaeps))
  x2 <- stats::qweibull(c2, shape = kappat, scale = 1/gamma(1 + 1/kappat))
  # calculate epsru and tru
  data.frame(
    dr = dr,
    epsru = x1*epsru0*(dr/dr0)^betaeps,
    tru = x2*tru0*(dr/dr0)^betat
  )
}
