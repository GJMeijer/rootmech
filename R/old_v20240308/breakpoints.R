#' Select breakpoints for initial root count guess guess
#'
#' @description
#' Generate a series of optional x-breakpoints on a known interval `xmin` to
#' `xmax`. The domain is split into equal intervals, where `n0` indicates the
#' extra number of optional points on top of `n`. All combinations of `n`
#' breakpoints from the optional set with length `n + n0` are returned as a
#' matrix
#'
#' Specific breakpoints may be fixed using the `fixed` argument.
#'
#' @param n number of breakpoints to select
#' @param xmin,xmax min and max value of domain
#' @param n0 total number of extra optional locations
#' @param fixed array with values of fixed breakpoints. Should contain `NA` for
#'   each breakpoint that is unknown
#' @return matrix with breakpoint options
#' @export
#' @examples
#' # input
#' n <- 3
#' xmin <- 0
#' xmax <- 10
#'
#' # no fixes
#' breakpoint_options(n, xmin, xmax, n0 = 2)
#' # fix in middle
#' breakpoint_options(n, xmin, xmax, fixed = c(NA, 6, NA), n0 = 4)
#' # fix at end
#' breakpoint_options(n, xmin, xmax, fixed = c(NA, NA, 6), n0 = 4)
#'
breakpoint_options <- function(n, xmin, xmax, n0 = 4, fixed = rep(NA, n)) {
  if (all(!is.na(fixed))) {
    cmb <- matrix(fixed, nrow = 1)
  } else {
    # only keep regions in which a fit should be made (must contain at least one NA)
    xbt <- c(xmin, fixed, xmax)
    ina <- which(is.na(xbt))
    t1 <- xbt[ina - 1]
    t2 <- xbt[ina + 1]
    ds <- data.frame(
      xmin = t1[!is.na(t1)],
      xmax = t2[!is.na(t2)],
      nb = which(diff(is.na(xbt)) == -1) - which(diff(is.na(xbt)) == 1)
    )
    # options per segment - fit as evenly as possible
    nt <- with(ds, (xmax - xmin)/(sum(xmax - xmin))*n0)
    ds$ns <- round(ds$nb + nt)
    # options
    Lopts <- lapply(
      1:nrow(ds),
      function(i) seq(ds$xmin[i], ds$xmax[i], l = ds$ns[i] + 2)[2:(ds$ns[i] + 1)]
    )
    # draws per segment
    Lperm <- lapply(
      1:nrow(ds),
      function(i) if (ds$nb[i] == 1) {
        matrix(Lopts[[i]], ncol = 1)
      } else {
        t(utils::combn(Lopts[[i]], ds$nb[i]))
      }
    )
    # all combs of unknown
    co <- Reduce(function(...) merge(..., by = NULL), Lperm)
    # add fixed values back in
    cmb <- matrix(rep(fixed, nrow(co)), nrow = nrow(co), byrow = TRUE)
    cmb[, ina - 1] <- as.matrix(co)
  }
  # return matrix
  cmb
}
